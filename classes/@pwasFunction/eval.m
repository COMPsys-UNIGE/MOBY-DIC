function u = eval(object,x,f)
% eval   Evaluates the PWAS function
%
% Y = eval(OBJ,X)
% Evaluates all NY components of the PWAS function at the points
% specified by matrix X. X must be a NX x NPOINTS matrix. Y is a
% NY x NPOINTS matrix. NX and NY are the domain and codomain dimensions,
% respectively. OBJ is the pwasFunction object.
%
% Y = eval(OBJ,X,IDX)
% Evaluates only the function components defined by IDX. IDX must be a vector
% of indices in the range 1, NY.
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

if isempty(object.w)
    error('No weights have been set. Cannot evaluate function');
end

if ~exist('f','var')
    f = 1:object.ny;
end

if floor(f) ~= f
    error('IDX must be an integer number');
end
if any(f < 1) || any(f > object.ny)
    error(['IDX must be an integer number between 1 and ',num2str(object.ny)]);
end

ny = numel(f);

% Check on x
if size(x,1) ~= object.nx
    error(['Matrix X must be of a matrix with ',num2str(object.nx),' rows']);
end

if any(isnan(x(:))) || any(isinf(x(:)))
    error('X cannot be NaN nor Inf');
end
    

% Number of points x
npoints = size(x,2);

% Retrieve rounding tolerance
pars = getMOBYDICpars();
tol = pars.roundtol;

usemex = 1;

% Do not use mex files
if ~usemex
    
    % Extract domain matrices
    Hd = object.domain.Hd;
    Kd = object.domain.Kd;
    
    % Initialize arrays
    u = zeros(npoints,ny);
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
            u(i,:) = NaN;
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
            addr = findIndex(xs,np);
            
            % mu coefficients
            mu = [-diff([1; decs]); decs(end)];
            
            % Indices of non zero mu coefficients
            ind = mu ~= 0;
            
            % Retrieve weights (only corresponding to non zero mu
            % coefficients)
            w = object.w(addr(ind),f);
            
            % Evaluate function
            u(i,:) = mu(ind)'*w;
            
        end
    end
    
    
    % Use mex files
else
    
    % Retrieve domain boundaries
    Kd = object.domain.Kd;
    xmin = -Kd(object.nx+1:end);
    xmax = Kd(1:object.nx);
    
    % Number of subdivisions per dimension
    np = object.np;
    % Maximum number of subdivisions
    maxnp = max(np)+1;
    
    % Put all subdivisions in a unique [nx x maxnp] matrix
    P = zeros(object.nx,maxnp);
    for i = 1:object.nx
        P(i,1:np(i)+1) = object.P{i};
    end
    
    % Extract weights
    w = object.w(:,f);
    
    % Launch mex function
    % P : n_dim x n_part
    % np : n_dim x 1
    % w: n_w x n_fun
    % x: n_dim x n_points
    u = pwasFunction_evalmex(P,np,w,xmin,xmax,x);
    
    u(u > 1e90) = NaN;
    
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


