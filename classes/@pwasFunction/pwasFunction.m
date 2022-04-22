classdef pwasFunction < MOBYFunction
    % pwasFunction   Piece-Wise Affine Simplicial function
    %
    % This object represents a piecewise affine function defined over a
    % uniform or non-uniform simplicial domain partition. A PWAS function
    % is defined as follows:
    %
    %              _nv_
    %              \
    %       u(x) = /___  w_j alpha_j(x)
    %              j = 1
    %
    % being nv the number of partition vertices and functions alpha_j
    % defined as:
    %
    %                 /
    %                 | 1, if j = k
    % alpha_j(v_k) = <
    %                 | 0, otherwise
    %                 \
    %
    % where v_k , k = 1, ..., nv are the partition vertices. Once the
    % structure of the partition is fixed, weights w_j uniquely define the
    % shape of the PWAS function. The function domain is hyper rectangular
    % with nx dimensions.
    %
    % OBJ = pwasFunction()
    % Builds an empty pwasFunction object OBJ.
    %
    % OBJ = pwasFunction(DOMAIN,NP)
    % Builds a pwasFunction object OBJ by specifying the hyper-rectangular
    % function domain (DOMAIN) and the number of subdivisions per dimension
    % (NP). DOMAIN is a structure with the following fields:
    % - xmin: array with nx elements specifying the lower domain boundaries
    % - xmax: array with nx elements specifying the upper domain boundaries
    % NP is an array with the same number of components as the domain
    % dimensions which defines the number of subdivisions per dimension.
    % The resulting simplicial partition is uniform. If NP is a scalar, all
    % domain dimensions are subdivided into the same number of segments.
    % For example, if NP = [3 2], the resulting partition is the following
    % (each rectangle is in turn divided into two triangles):
    %     ___________
    %    |   |   |   |
    %    |   |   |   |
    %    |___|___|___|
    %    |   |   |   |
    %    |   |   |   |
    %    |___|___|___|
    %
    % Since the weights are not provided, the shape of the function is not
    % defined.
    %
    % OBJ = pwasFunction(DOMAIN,P)
    % Builds a pwasFunction object OBJ by specifying the hyper-rectangular
    % function domain (DOMAIN) and the subdivisions per dimension (P).
    % DOMAIN is a structure with the following fields:
    % - xmin: array with nx elements specifying the lower domain boundaries
    % - xmax: array with nx elements specifying the upper domain boundaries
    % P is a cell-array with the same number of components as the domain
    % dimensions. The j-th element of P contains the coordinates of the
    % subdivisions along the j-th dimension.
    % For example if P{1} = [1 3 4] and P{2} = [1 2 4], the resulting
    % partition is the following (each rectangle is in turn divided into
    % two triangles):
    %
    %  4  ___________
    %    |       |   |
    %    |       |   |
    %    |       |   |
    %  2 |_______|___|
    %    |       |   |
    %  1 |_______|___|
    %    1       3   4
    %
    % The resulting simplicial partition is in general non uniform. Since
    % the weights are not provided, the shape of the function is not defined.
    %
    % OBJ = pwasFunction(DOMAIN,NP,W)
    % Builds a pwasFunction object OBJ by specifying the hyper-rectangular
    % function domain (DOMAIN), the number of subdivisions per dimension NP
    % and the weights W. W must be a vector with nv (number of vertices)
    % rows and ny (number of outputs) columns. The ordering of the weights
    % is as follows:
    %
    % w_7 ____w_8____ w_9
    %    |       |   |
    %    |       |   |
    %    |       |   |
    % w_4|____w_5|___|w_6
    %    |       |   |
    %    |_______|___|
    % w_1     w_2     w_3
    %
    % Since the weights are provided, the shape of the function is defined.
    %
    % OBJ = pwasFunction(DOMAIN,P,W)
    % Builds a pwasFunction object OBJ by specifying the hyper-rectangular
    % function domain (DOMAIN), the subdivisions per dimension P and the
    % weights W. Since the weights are provided, the shape of the function
    % is defined.
    %
    % pwasFunction methods:
    %   alphaBasis - evaluates the alpha bases functions.
    %   disp - displays some information about the pwasFunction object.
    %   eval - evaluates the PWAS function.
    %   generateC - Generates C files for the circuit implementation of the PWAS function on microcontroller
    %   generateVHDL - Generates VHDL files for the circuit implementation of the PWAS function on FPGA
    %   getNumberOfPartitions - gets the number of subdivisions per dimension
    %   getNumberOfSimplices - gets the number of simplices of the domain partition
    %   getNumberOfVertices - gets the number of partition vertices
    %   getPartition - gets the subdivisions per dimension
    %   getSimplices - Gets the vertices of all simplices of the simplicial partition
    %   getVertices - gets the vertices of the simplices of the domain partition
    %   getWeights - gets the weights w of the PWAS function
    %   isUniform - returns true if the simplicial partition is uniform
    %   plot - plots the PWAS function.
    %   plotPartition - plots the simplicial domain partition.
    %   setWeights - sets the weights w of the PWAS function
    %
    % The pwasFunction object is derived from MOBYFunction and inherits
    % all its methods.
    %
    % See also MOBYFunction, pwagFunction.
    
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
    
    properties (Access = protected)
        
        % Number of vertices
        nv = 0;
        % Number of subdivisions per dimension
        np = 0;
        % Subdivisions per dimension
        P = [];
        % Weights
        w = [];
        % Uniform or non uniform partition
        uniform = [];
        
    end
    
    methods
        
        function object = pwasFunction(varargin)
            if nargin == 0
                
                object.nv = 0;
                object.np = 0;
                object.P = [];
                object.w = [];
                object.uniform = [];
                
            elseif nargin == 2 || nargin == 3
                
                if ~isstruct(varargin{1})
                    error('First argument must be a struct');
                end
                domain = varargin{1};
                if ~isfield(domain,'xmin')
                    error('domain must have fields xmin and xmax');
                end
                if ~isfield(domain,'xmax')
                    error('domain must have fields xmin and xmax');
                end
                
                % Transform into column vectors
                domain.xmin = domain.xmin(:);
                domain.xmax = domain.xmax(:);
                
                % Number of dimensions
                nx = numel(domain.xmin);
                
                if numel(domain.xmax) ~= nx
                    error('domain.xmin and domain.xmax must have the same number of elements');
                end
                
                if any(domain.xmin >= domain.xmax)
                    error('domain.xmin must be lower than domain.xmax');
                end
                
                % Create matrices to express the domain
                Hd = [eye(nx);-eye(nx)];
                Kd = [domain.xmax;-domain.xmin];
                
                if isnumeric(varargin{2})
                    % Number of subdivisions per dimension
                    np = varargin{2};
                    
                    np = np(:);
                    
                    if numel(np) == 1
                        np = repmat(np,nx,1);
                    end
                    
                    if numel(np) ~= nx
                        error(['np must be a scalar or a vector with ',num2str(nx),' elements']);
                    end
                    
                    P = cell(nx,1);
                    for i = 1:nx
                        P{i} = linspace(domain.xmin(i),domain.xmax(i),np(i)+1);
                    end
                    
                    uniform = 1;
                    
                elseif iscell(varargin{2})
                    
                    P = varargin{2};
                    
                    if numel(P) ~= nx
                        error(['P must be a cell array with ',num2str(nx),' elements']);
                    end
                    
                    np = zeros(nx,1);
                    for i = 1:nx
                        if any(P{i} ~= sort(P{i}))
                            error('Elements of P must be sorted in ascending order');
                        end
                        if P{i}(1) ~= domain.xmin(i)
                            error('First element of P must coincide with domain boundaries');
                        end
                        if P{i}(end) ~= domain.xmax(i)
                            error('Last element of P must coincide with domain boundaries');
                        end
                        np(i) = numel(P{i})-1;
                    end
                    
                    unif = 1;
                    for i = 1:nx
                        if any(P{i} ~= linspace(P{i}(1),P{i}(end),numel(P{i})))
                            unif = 0;
                        end
                    end
                    
                    uniform = unif;
                    
                else
                    
                    error('Second argument must be a vector or a cell array');
                    
                end
                
                nv = prod(np+1);
                
                if nargin == 2
                    w = [];
                    ny = 0;
                    
                    xnames = cell(nx,1);
                    for i = 1:nx
                        xnames{i} = ['x_',num2str(i)];
                    end
                    
                    ynames = [];
                    if ~isempty(w)
                        ynames = cell(ny,1);
                        for i = 1:ny
                            ynames{i} = ['y_',num2str(i)];
                        end
                    end
                    
                else
                    w = varargin{3};
                    
                    if size(w,1) ~= nv
                        error(['w must have ',num2str(nv),' rows']);
                    end
                    
                    ny = size(w,2);
                    
                    xnames = cell(nx,1);
                    for i = 1:nx
                        xnames{i} = ['x_',num2str(i)];
                    end
                    
                    ynames = [];
                    if ~isempty(w)
                        ynames = cell(ny,1);
                        for i = 1:ny
                            ynames{i} = ['y_',num2str(i)];
                        end
                    end
                end
                
                
                
                % Fill fields
                object.nx = nx;
                object.ny = ny;
                
                object.domain.Hd = Hd;
                object.domain.Kd = Kd;
                object.xnames = xnames;
                object.ynames = ynames;
                
                object.nv = nv;
                object.np = np;
                object.P = P;
                object.w = w;
                object.uniform = uniform;
                
            else
                error('Wrong input number');
            end
            
        end
    end
    
    methods
        % Get Methods
        np = getNumberOfPartitions(object);
        P = getPartition(object);
        nv = getNumberOfVertices(object);
        w = getWeights(object);
        ns = getNumberOfSimplices(object);
        V = getVertices(object);
        
        % Is Methods
        uniform = isUniform(object);
        
        % Set Methods
        object = setWeights(object,w);
        
        % Other
        plot(object,varargin)
        plotPartition(object)
        y = eval(object,x,f);
        varargout = generateC(object, varargin);
        varargout = generateVHDL(object, varargin);
        
    end
    
end
