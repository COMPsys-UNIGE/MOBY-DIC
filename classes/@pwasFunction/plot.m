function plot(object,varargin)
% plot   Plots the PWAS function
%
% plot(OBJ)
% Plots all the components of the vector PWAS function OBJ.
%
% plot(OBJ,IDX)
% Plots the components of the vector PWAS function specified in the
% vector of indices IDX. IDX must be a vector of indices in the range 1, NY.
%
% plot(OBJ,SEC)
% Plots all the components of the vector PWAS function OBJ in a 
% at most 2-dimensional section defined by SEC. SEC is a structure with
% fields dim and val. Consider a function with 4 inputs x1, x2, x3 and x4.
% In order to plot the function in the section x2 = 1 and x3 = 0 you have
% to set SEC.dim = [2 3] and SEC.val = [1 0].
%
% plot(OBJ,SEC,IDX)
% Plots the components of the vector PWAS function specified in the
% vector of indices IDX in a at most 2-dimensional section defined by SEC.
%
% This functions exploits MPT3 to compute the bounding box of the domain.

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

if nargin == 1
    sec.dim = [];
    sec.val = [];
    idx = [];
elseif nargin == 2
    if isstruct(varargin{1})
        sec = varargin{1};
        idx = [];
    else
        idx = varargin{1};
        sec.dim = [];
        sec.val = [];
    end
elseif nargin == 3
    sec = varargin{1};
    idx = varargin{2};
end

if ~isfield(sec,'dim')
    error('SEC must be a structure with fields dim and val');
end
if ~isfield(sec,'val')
    error('SEC must be a structure with fields dim and val');
end

nsec = numel(sec.dim);
if numel(sec.val) ~= nsec
    error('Wrong structure SEC.');
end

if any(sec.dim < 1 | sec.dim > object.nx)
    error(['SEC.dim must be an integer between 0 and ',num2str(object.nx)]);
end

if isempty(idx)
    idx = 1:object.ny;
end

if floor(idx) ~= idx
    error('IDX must be an integer number');
end
if any(idx < 1) || any(idx > object.ny)
    error(['IDX must be a vector of integer numbers between 1 and ',num2str(object.ny)]);
end

if object.nx > 3
    if isempty(sec)
        error('SEC must be provided!');
    end
end
if object.nx-nsec > 2
    error('Cannot project function');
end

% If no section is provided...
if isempty(sec.dim)
    
    % Weights
    w = object.w;
    
    if isempty(w)
        error('Weights w have not been set. Cannot plot!')
    end
    
    % Extract partition
    P = object.P;
    
    if object.nx == 2
        % Extract vertices
        V = object.getVertices();
        
        % Rotate vertices
        Vr = [V(:,1) -V(:,2)];
        
        % Use MATLAB function delaunay
        tri = delaunay(Vr);
        
        % Extract weights
        w = w(:,idx);
        
        figure
        for i = 1:numel(idx)
            subplot(1,numel(idx),i)
            tmp = w(:,i);
            tmp = reshape(tmp,object.np(:)'+1);
            trisurf(tri,V(:,1),V(:,2),tmp(:),'FaceColor',[0.7 0.7 0.7])
            xlabel(object.xnames{1})
            ylabel(object.xnames{2})
            zlabel(object.ynames{idx(i)})
        end
    else
        % Extract weights
        w = w(:,idx);
        
        figure
        for i = 1:numel(idx)
            subplot(1,numel(idx),i)
            plot(P{1},w(:,i),'k','linewidth',2)
            xlabel(object.xnames{1})
            ylabel(object.ynames{idx(i)})
        end
    end
    
% If a section is provided...
else
    
    % Names
    xnames = object.xnames;
    ynames = object.ynames;
    
    % Retrieve the bounds...
    if isfield(object.domain,'xmin') && isfield(object.domain,'xmax')
        xmin = object.domain.xmin;
        xmax = object.domain.xmax;
        
        % ... otherwise, compute first the bounding box
    else
        
        % Create polyhedron for the domain
        P = Polyhedron(object.domain.Hd,object.domain.Kd);
        
        % Compute bounding box of domain
        B = outerApprox(P);
        
        % Retrieve bounds
        xmin = B.Internal.lb;
        xmax = B.Internal.ub;
    end
    
    % Create grid of points
    if object.nx-nsec == 1
        pvec = linspace(xmin(1),xmax(1),100);
        p = pvec;
    else
        [x, y] = ndgrid(linspace(xmin(1),xmax(1),100),linspace(xmin(2),xmax(2),100));
        p = [x(:) y(:)];
        p = p';
    end
    
    string = '';
    if ~isempty(sec.dim)
        
        npoints = size(p,2);
        
        dims = sec.dim(:);
        vals = sec.val(:);
        vals = repmat(vals,1,npoints);
        
        pext = zeros(object.nx,npoints);
        pext(dims,:) = vals;
        pext(setdiff(1:object.nx,dims),:) = p;
        
        p = pext;
        
        k = 1;
        h = 1;
        for i = 1:object.nx
            if ~any(sec.dim == i)
                xnames{k} = object.xnames{i};
                k = k+1;
            else
                string = [string,object.xnames{i},' = ',num2str(sec.val(h)),'   ']; %#ok<AGROW>
                h = h+1;
            end
        end
        
        
    end
    
    % Evaluate function
    u = object.eval(p,idx);
    
    if object.nx-nsec == 1
        figure
        for i = 1:numel(idx)
            subplot(1,numel(idx),i)
            plot(pvec,u(i,:),'k','linewidth',2);
            xlabel(xnames{1})
            ylabel(ynames{idx(i)})
            title(string)
        end
        
    else
        figure
        for i = 1:numel(idx)
            subplot(1,numel(idx),i)
            mesh(x,y,reshape(u(i,:),100,100));
            xlabel(xnames{1})
            ylabel(xnames{2})
            zlabel(ynames{idx(i)})
            title(string)
        end
        
    end
end