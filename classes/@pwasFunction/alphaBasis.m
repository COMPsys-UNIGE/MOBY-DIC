function alpha = alphaBasis(object,x,nzero)
% alphaBasis   Evaluates the alpha bases functions
%
% ALPHA = eval(OBJ,X,NZERO)
% Evaluates all the alpha bases at the points specified by matrix X. 
% X must be a NPOINTS x NX matrix. ALPHA is a NPOINTS x NV sparse matrix. 
% NX and NV are the domain dimensions and the number of partition
% vertices, respectively. NZERO is the estimated number of non zero
% elements in the sparse matrix ALPHA.
% OBJ is the pwasFunction object.
%
% Mex files are used to speedup the computation.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Matteo Lodi (matteo.lodi@edu.unige.it)
%
% Copyright (C) 2016 University of Genoa, Italy.

% Legal note:
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330,
% Boston, MA  02111-1307  USA

nv = object.getNumberOfVertices();

% Check on x
if size(x,1) ~= object.nx && size(x,2) ~= object.nx
    error(['Matrix x must be of a matrix with ',num2str(object.nx),' rows']);
end

% If the input x is a npoints*nx array it is converted
% in a nx*npoints array
if size(x,2) == object.nx && size(x,1) ~= object.nx
    x = x';
end

% Number of points x
npoints = size(x,2);

% Retrieve rounding tolerance
pars = getMOBYDICpars();
tol = pars.roundtol;

np = object.getNumberOfPartitions();
P = object.getPartition();

usemex = 1;

% Do not use mex files
if ~usemex
       
    rows = zeros(nzero,1);
    cols = zeros(nzero,1);
    vals = zeros(nzero,1);
    
    kk = 1;
    
    for ii = 1:nv
        
        % Indices of the subdivisions corresponding to the vertex of the
        % alpha basis
        part = linear2matrix(ii,np+1);
        
        % Domain of the i-th alpha basis
        xmin = zeros(object.nx,1);
        xmax = zeros(object.nx,1);
        for i = 1:object.nx
            if part(i) == 1
                xmin(i) = P{i}(1);
            else
                xmin(i) = P{i}(part(i)-1);
            end
            if part(i) == np(i)+1
                xmax(i) = P{i}(np(i)+1);
            else
                xmax(i) = P{i}(part(i)+1);
            end
        end
        
        walpha = zeros(nv,1);
        walpha(ii) = 1;
        
        Hd = [eye(object.nx);-eye(object.nx)];
        Kd = [xmax;-xmin];
        
        % Initialize arrays
        int = zeros(object.nx,1);
        idx = zeros(object.nx,1);
        dec = zeros(object.nx,1);
        
        % Partition
        P = object.P;
        
        % Length of subdivisions
        dP = cell(object.nx,1);
        for i = 1:object.nx
            dP{i} = [diff(P{i}) 1];
        end
        
        % Number of subdivisions per dimension
        np = object.np;
        
        % Loop on number of points
        for i = 1:npoints
            
            % Current point
            xx = x(:,i);
            
            % Check if the point is outside the domain
            if any(Hd*xx > Kd+tol)
                
                
            else
                
                % Find bottom-left vertex
                for j = 1:object.nx
                    idx(j) = find(xx(j) >= P{j},1,'last');
                    int(j) = P{j}(idx(j));
                    % Relative position
                    dec(j) = (xx(j)-int(j))/dP{j}(idx(j));
                end
                
                % Sort relative position
                decs = sort(dec,1,'descend');
                
                % Find vertices of simplex containing xx
                a = zeros(object.nx,object.nx+1);
                a(:,1) = 0;
                for j = 1:object.nx
                    a(:,j+1) = dec >= decs(j);
                end
                xs = repmat(idx,1,object.nx+1);
                xs = xs+a;
                
                % Addresses of weights
                addr = findIndex(xs,object.np);
                
                % mu coefficients
                mu = [-diff([1; decs]); decs(end)];
                
                % Indices of non zero mu coefficients
                ind = mu ~= 0;
                
                % Retrieve weights (only corresponding to non zero mu
                % coefficients)
                w = walpha(addr(ind));
                
                % Evaluate function
                val = mu(ind)'*w;
                
                if val ~= 0
                    rows(kk) = i;
                    cols(kk) = ii;
                    vals(kk) = val;
                    kk = kk+1;
                end
                
            end
        end
    end
    
    rows = rows(1:kk-1);
    cols = cols(1:kk-1);
    vals = vals(1:kk-1);
    
    alpha = sparse(rows,cols,vals,npoints,nv);
    
    
    % Use mex files
else
    
    % Number of subdivisions per dimension
    np = object.np;
    % Maximum number of subdivisions
    maxnp = max(np)+1;
    
    % Put all subdivisions in a unique [nx x maxnp] matrix
    P = zeros(object.nx,maxnp);
    for i = 1:object.nx
        P(i,1:np(i)+1) = object.P{i};
    end
    
    % Launch mex function
    % P : n_dim x n_part
    % np : n_dim x 1
    % w: n_w x n_fun
    % x: n_dim x n_points
    
    [rows, cols, vals, nval] = pwasFunction_alphamex(P,np,nzero,x);
    rows = rows(1:nval);
    cols = cols(1:nval);
    vals = vals(1:nval);
    
    alpha = sparse(rows,cols,vals,npoints,nv);
    
end

    function addr = findIndex(x,np)
        coeff = 1;
        addr = 0;
        for k = 1:numel(np)
            addr = addr+coeff*(x(k,:)-1);
            coeff = coeff*(np(k)+1);
        end
        addr = addr+1;
    end

end


