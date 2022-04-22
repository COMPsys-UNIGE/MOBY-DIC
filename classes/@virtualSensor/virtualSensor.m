classdef virtualSensor
    % virtualSensor   Piece-Wise Affine Simplicial virtual sensor
    %
    % TO DO
    %
    % OBJ = virtualSensor()
    % Builds an empty virtualSensor object OBJ
    %
    % OBJ = virtualSensor(NU,NY,NZ,MU,MY,MZ)
    % Builds a virtual sensor object OBJ by specifying the number of system
    % inputs (NU), the number of measurable system outputs (NY), the number
    % of unmeasurable outputs (NZ), the input time window (MU), the measurable
    % output time window (MY) and the autoregressive time window (MZ).
    %
    % OBJ = virtualSensor(NU,NY,MU,MY,MZ,OPTS)
    % Allows specifying some options. OPTS is a structure with the following
    % fields:
    % - reducedComplexity - if it is set to 1 a reduced complexity virtual
    %                       sensor is obtained (default value: 0)
    % - current - if it is set to 1 the values of u and y at the current time
    %             instant (u_k, y_k) are used to estimate z at the same
    %             instant (z_k); if it is 0, z_k is estimated starting from
    %             u_{k-1}, y_{k-1}. Defalut value: 1.
    %
    % virtualSensor methods:
    %   disp - Displays some information about the virtualSensor object
    %   eval - Evaluates the unmeasurable output z in correspondence of given
    %          system inputs and measurable outputs
    %   getAutoregressiveTimeWindow - Gets the autoregressive time window mz
    %   getFunction - Gets the pwasFunction object defining the virtual sensor
    %   getIdentificationInformation - Gets information about ridge regression
    %                                  process used for sensor identification
    %   getInputNames - gets the names of the system inputs
    %   getInputTimeWindow - Gets the input time window mu
    %   getMeasurableOutputNames - gets the names of the system measurable outputs    
    %   getNumberOfInputs - Gets the number of input variables
    %   getNumberOfMeasurableOutputs - Gets the number of measurable output
    %                                  variables
    %   getNumberOfUnmeasurableOutputs - Gets the number of unmeasurable
    %                                    output variables
    %   getNumberOfPartitions - Gets the number of partitions for each
    %                           dimension of the pwas virtual sensor domain
    %   getOutputTimeWindow - Gets the output time window my
    %   getUnmeasurableOutputNames - gets the names of the system unmeasurable outputs
    %   identify - Identifies the virtual sensor starting from the measured
    %              system inputs and outputs
    %   isCurrent - Indicates if the values of the inputs and measured outputs
    %               at current time instant are used to estimate the unmeasurable output
    %   isIdentified - Returns 1 if the virtual sensor has been identified, 0
    %                  otherwise
    %   generateC - Generates C files for the circuit implementation of
    %               the PWAS function on microcontroller
    %   generateVHDL - Generates VHDL files for the circuit implementation of
    %                  the PWAS function on FPGA
    %   setInputNames - sets the names of the system inputs
    %   setMeasurableOutputNames - sets the names of the system measurable outputs
    %   setUnmeasurableOutputNames - sets the names of the system unmeasurable outputs
    %   validate - Validates the virtual sensor through a test set of data

    % Contributors:
    %
    % Alberto Oliveri (alberto.oliveri@unige.it)
    %
    % Copyright (C) 2015 University of Genoa, Italy.
    
    % Legal note:
    %          This program is free software; you can redistribute it and/or
    %          modify it under the terms of the GNU General Public
    %          License as published by the Free Software Foundation; either
    %          version 2.1 of the License, or (at your option) any later version.
    %
    %          This program is distributed in the hope that it will be useful,
    %          but WITHOUT ANY WARRANTY; without even the implied warranty of
    %          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    %          General Public License for more details.
    %
    %          You should have received a copy of the GNU General Public
    %          License along with this library; if not, write to the
    %          Free Software Foundation, Inc.,
    %          59 Temple Place, Suite 330,
    %          Boston, MA  02111-1307  USA
       
    properties
        nu = 0; % Number of system inputs
        ny = 0; % Number of measurable system outputs
        nz = 0; % Number of unmeasurable system outputs
        mu = 0; % Input time window
        my = 0; % Output time window
        mz = 0; % Autoregressive time window
        np = 0; % Domain partition
        Ts = 0; % Sampling time
        domain = struct('umin',[],'umax',[],'ymin',[],'ymax',[],...
            'zmin',[],'zmax',[]);    % Virtual sensor domain
        fpwas = pwasFunction(); % pwasFunction object
        reducedComplexity = 0;  % Flag indicating if the virtual sensor has 
                                % reduced complexity
        current = 1;    % Flag indicating if u and y at current time step 
                        % are used to estimate z
        info = [];  % Information on virtual sensor identification
        identified = 0; % Flag indicating whether the virtual sensor has 
                        % been identified or not
        unames = {};    % Input names
        ynames = {};    % Measurable output names
        znames = {};    % Unmeasurable output names
    end
    
    methods
        
        % Constructor
        function vs = virtualSensor(varargin)
            
            % vs = virtualSensor()
            if nargin == 0
                vs.nu = 0;
                vs.ny = 0;
                vs.mu = 0;
                vs.my = 0;
                vs.mz = 0;
                vs.np = 0;
                vs.Ts = 0;
                vs.domain = struct('umin',[],'umax',[],'ymin',[],'ymax',[],...
                    'zmin',[],'zmax',[]); 
                vs.fpwas = pwasFunction();
                vs.reducedComplexity = 0;
                vs.current = 1;
                vs.info = [];
                vs.identified = 0;
                vs.unames = {};
                vs.ynames = {};
                vs.znames = {};
                
            elseif nargin == 6 || nargin == 7
                               
                nu = varargin{1};
                ny = varargin{2};
                nz = varargin{3};
                mu = varargin{4};
                my = varargin{5};
                mz = varargin{6};
                
                if ~isa(nu,'double')
                    error('NU must be a non-negative integer value');
                end
                if nu ~= floor(nu) || nu < 0
                    error('NU must be a non-negative integer value');
                end
                if ~isa(ny,'double')
                    error('NY must be a non-negative integer value');
                end
                if ny ~= floor(ny) || ny <= 0
                    error('NY must be a non-negative integer value');
                end
                if ~isa(nz,'double')
                    error('NZ must be a positive integer value');
                end
                if nz ~= floor(nz) || nz <= 0
                    error('NZ must be a positive integer value');
                end
                if nz > 1
                    error('Only NZ = 1 is supported, at present.');
                end
                
                if ~isa(mu,'double')
                    error('MU must be a non-negative integer value');
                end
                if mu ~= floor(mu) || mu < 0
                    error('MU must be a non-negative integer value');
                end
                if numel(mu) == 1
                    mu = repmat(mu,nu,1);
                end
                if ~isa(my,'double')
                    error('MY must be a non-negative integer value');
                end
                if my ~= floor(my) || my < 0
                    error('MY must be a non-negative integer value');
                end
                if numel(my) == 1
                    my = repmat(my,ny,1);
                end
                if ~isa(mz,'double')
                    error('MZ must be a non-negative integer value');
                end
                if mz ~= floor(mz) || mz < 0
                    error('MZ must be a non-negative integer value');
                end
                if numel(mz) ~= 1
                    error('MZ must be a scalar');
                end             
                
                if numel(mu) ~= nu
                    error('MU must be either a scalar or an array with NU elements');
                end
                if numel(my) ~= ny
                    error('MY must be either a scalar or an array with NY elements');
                end
                
                if isempty(mu)
                    mu = 0;
                end
                if isempty(my)
                    my = 0;
                end
                if isempty(mz)
                    mz = 0;
                end
                
                % Fill object properties
                vs.nu = nu;
                vs.ny = ny;
                vs.nz = nz;
                vs.mu = mu;
                vs.my = my;
                vs.mz = mz;
                vs.np = 0;
                vs.Ts = 0;
                vs.domain = struct('umin',[],'umax',[],'ymin',[],'ymax',[],...
                    'zmin',[],'zmax',[]); 
                vs.fpwas = pwasFunction();
                vs.info = [];
                vs.identified = 0;
                
                % Check options
                if nargin == 7
                    options = varargin{7};
                    
                    if ~isfield(options,'reducedComplexity') || isempty(options.reducedComplexity)
                        vs.reducedComplexity = 0;
                    else
                        if options.reducedComplexity == 0 || options.reducedComplexity == 1
                            vs.reducedComplexity = options.reducedComplexity;
                        else
                            error('Field reducedComplexity must be 0 or 1');
                        end
                    end
                    
                    if ~isfield(options,'current') || isempty(options.current)
                        vs.current = 1;
                    else
                        if options.current == 0 || options.current == 1
                            vs.current = options.current;
                        else
                            error('Field current must be 0 or 1');
                        end
                    end
                    
                else
                    
                    vs.reducedComplexity = 0;
                    vs.current = 1;
                    
                end
                
                unames = cell(nu,1);
                for i = 1:nu
                    unames{i} = ['u_',num2str(i)];
                end
                
                ynames = cell(ny,1);
                for i = 1:ny
                    ynames{i} = ['y_',num2str(i)];
                end
                
                znames = cell(nz,1);
                for i = 1:nz
                    znames{i} = ['z_',num2str(i)];
                end
                
                vs.unames = unames;
                vs.ynames = ynames;
                vs.znames = znames;
                
                
            else
                
                error('Wrong input arguments');
            end
            
        end % constructor
        
    end % methods
    
    methods
        
        % Get methods
        nu = getNumberOfInputs(object);
        ny = getNumberOfMeasurableOutputs(object);
        nz = getNumberOfUnmeasurableOutputs(object);
        mu = getInputTimeWindow(object);
        my = getOutputTimeWindow(object);
        mz = getAutoregressiveTimeWindow(object);
        np = getNumberOfPartitions(object);
        fpwas = getFunction(object);
        info = getIdentificationInformation(object);
        
        % Is methods
        identified = isIdentified(object);
              
        % Other methods
        disp(object)
        object = identify(object,u,y,z,varargin);
        zh = eval(object,varargin);
        [zh, err] = validate(object,varargin);
        varargout = generateC(object,varargin);
        varargout = generateVHDL(object,varargin);
        
    end
    
end % classdef


