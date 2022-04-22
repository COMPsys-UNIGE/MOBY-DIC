classdef nonlinFunction < MOBYFunction
    
    % nonlinFunction   Generic non linear function
    %
    % This object represents a generic non linear vector function f defined
    % over a polytopic domain.
    %
    % OBJ = nonlinFunction()
    % Builds an empty nonLinFunction object OBJ.
    %
    % OBJ = nonlinFunction(NX, NY, FHANDLE, DOMAIN)
    % Builds a nonLinFunction object OBJ by specifying the number of domain
    % dimensions (NX), the number of codomain dimensions (NY), a handle
    % FHANDLE to a user-defined function (it must be defined on a separate
    % m-file belonging to the MATLAB path) and the function domain
    % (DOMAIN). The user-defined function must be in the form:
    %
    % function Y = FUNCTION_NAME(X)
    % ...
    %
    % Input X is a NX x NPOINTS matrix, therefore output Y must be a
    % NY x NPOINTS matrix.
    %
    % An example of function definition, with NX = 3 and NY = 2 is the
    % following:
    %
    % function Y = PLANE_PARABOLA(X)
    %   npts = size(X,2); % Number of points
    %   ny = 2; % Number of codomain dimensions
    %   Y = zeros(ny,npts); % Initialize outputs
    %   Y(1,:) = X(1,:)+X(2,:)+X(3,:); % Compute outputs
    %   Y(2,:) = X(1,:).^2+X(2,:).^2+X(3,:).^2;
    %
    % If the domain is hyper rectangular, DOMAIN can be a
    % structure with fields xmin and xmax, providing the domain boundaries:
    %
    %    xmin <= x <= xmax
    %
    % Otherwise DOMAIN must be a struct with fields Hd and Kd, which define
    % the polytopic domain as follows:
    %
    %    Hd x <= Kd
    %
    % nonlinFunction methods:
    %   disp - displays some information about the object.
    %   eval - evaluates the non linear function.
    %   plot - plots the non linear function.
    %
    % The nonlinFunction object is derived from MOBYFunction and inherits
    % all its methods.
    %
    % See also MOBYFunction.
    
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
    
    % Properties
    
    properties (Access = private)
        fhandle = [];   % Handle to the user-defined function
    end
    
    
    % Methods
    
    methods
        
        % Constructor
        function object = nonlinFunction(varargin)
            
            % Empty object
            if nargin == 0
                object.nx = 0;
                object.ny = 0;
                object.domain.Hd = [];
                object.domain.Kd = [];
                object.xnames = [];
                object.ynames = [];
                object.fhandle = [];
                
            elseif nargin == 4
                nx = varargin{1};
                ny = varargin{2};
                fhandle = varargin{3};
                
                domain = varargin{4};
                
                if isfield(domain,'xmin')
                    if isfield(domain,'Hd') || isfield(domain,'Kd')
                        error('DOMAIN must be a struct with fields "xmin, xmax" or "Hd, Kd"');
                    end
                    if ~isfield(domain,'xmax')
                        error('DOMAIN must be a struct with fields xmin, xmax or Hd, Kd');
                    end
                elseif isfield(domain,'Hd')
                    if isfield(domain,'xmin') || isfield(domain,'xmax')
                        error('DOMAIN must be a struct with fields "xmin, xmax" or "Hd, Kd"');
                    end
                    if ~isfield(domain,'Kd')
                        error('DOMAIN must be a struct with fields xmin, xmax or Hd, Kd');
                    end
                end
                
                if isfield(domain,'xmin')
                    
                    % Transform into column vectors
                    domain.xmin = domain.xmin(:);
                    domain.xmax = domain.xmax(:);
                    
                    if numel(domain.xmax) ~= nx || numel(domain.xmin) ~= nx
                        error(['DOMAIN.xmin and DOMAIN.xmax must have ',num2str(nx),' elements']);
                    end
                    
                    if any(domain.xmin >= domain.xmax)
                        error('DOMAIN.xmin must be lower than DOMAIN.xmax');
                    end
                    
                    % Create matrices to express the domain
                    domain.Hd = [eye(nx);-eye(nx)];
                    domain.Kd = [domain.xmax;-domain.xmin];
                end
                
                Hd = domain.Hd;
                Kd = domain.Kd;
                
                if size(Hd,1) ~= size(Kd,1)
                    error('DOMAIN.Hd and DOMAIN.Kd must have the same number of rows')
                end
                
                % Check domain edges (if provided)
                if ~isempty(domain)
                    if size(Hd,2) ~= nx
                        error(['DOMAIN.Hd must have ',num2str(nx),' columns']);
                    end
                    if size(Kd,2) ~= 1
                        error('DOMAIN.Kd must have 1 column');
                    end
                end
                
                % Check if the domain is bounded
                P = Polyhedron(Hd,Kd);
                isbounded = P.isBounded();
                
                if ~isbounded
                    disp('WARNING: the function domain is unbounded.');
                end
                
                xnames = cell(nx,1);
                ynames = cell(ny,1);
                for i = 1:nx
                    xnames{i} = ['x_',num2str(i)];
                end
                for i = 1:ny
                    ynames{i} = ['y_',num2str(i)];
                end
                
                
                
                % Test function calling it on 10 points
                x = zeros(nx,20);
                try
                    y = fhandle(x);
                catch err
                    error(['Error calling the function handle you provided!',...
                        'Please check the function and the dimension of domain and codomain!',...
                        ' The error is the following: ',err.message])
                end
                if size(y,1) ~= ny || size(y,2) ~= 20
                    error(['The output of the function defined by the handle must be a matrix with size ',num2str(ny),' x NPOINTS']);
                end
                
                % Fill object fields
                
                object.nx = nx;
                
                object.ny = ny;
                
                object.fhandle = fhandle;
                
                object.xnames = xnames;
                object.ynames = ynames;
                
                object.domain = domain;
                
            else
                error('Wrong input arguments.');
            end
            
        end
        
        % Get methods
        
        y = eval(object,x,varargin);
        
        plot(object,f);
        
        disp(object);
        
    end
end


