classdef implicitMPCctrl < controller
    % implicitMPCctrl   Implicit MPC controller
    %
    % This object represents an implicit MPC controller for a LTI
    % system. The controller can be designed either for regulation to a
    % constant reference state (or output) or for tracking an external signal. 
    % The control function is designed to satisfy constraints on
    % the system variables.
    %
    % See the User's Guide for a more detailed explaination of this object.
    %
    % OBJ = implicitMPCctrl()
    % Builds an empty implicitMPCctrl object OBJ.
    %
    % OBJ = implicitMPCctrl(SYS,Ts,CONSTR,OPTS)
    % Builds an implicitMPCctrl object starting from a ltiSys object SYS
    % and a constraints object CONSTR. Ts is the sampling time of the
    % controller. SYS represents the system to regulate and CONSTR the
    % inequality constraints to be fulfilled. OPTS is a structure with the
    % following fields:
    % - N, Nu: prediction and control horizon (default: Nu = N)
    % - P, Q, R, rho: matrices defining the MPC cost function
    %                 (default: P = Q)
    % - norm: norm used in the cost function; available values are 1, 2 or
    %         inf (default: norm = 2)
    % - K, O: define the control function after the control horizon Nu.
    %         u = Kx + O (default: K = 0, O = 0)
    % - ctrlvariable: string indicating if you want to regulate (or track)
    %                 system states or outputs. Possible values are:
    %                 'state', 'output'.
    % - tracking: boolean indicating if the controller is for tracking
    %             (default: false)
    % - trackvariable: indices of the state variables to track (only for
    %                  tracking controllers)
    % - ref, uref: reference state or output (depending on ctrlvariable) 
    %              and reference input (default: ref = 0, uref = 0)
    % - constantInputAfterNu : boolean indicating if the control must be 
    %                          u = Kx + O or u(i) = u(Nu) after Nu
    % - algorithm : optimization algorithm (quadprog or admm)
    % - defaultOutput : starting point for the optimization algorithm
    %                   (default: defaultOutput = 0)
    % - trajectoryTracking : boolean indicating if the reference trajectory
    %                        is known (only for simulation)
    %
    %
    % implicitMPCctrl methods:
    %   disp - displays some information about the implicitMPCctrl object.
    %   eval - evaluates the MPC controller.
    %   generateC - Generates C files for the circuit implementation of the MPC
    %               controller on microcontroller
    %   generateVHDL - Generates VHDL files for the circuit implementation 
    %                  of the MPC controller on FPGA
    %   getInformation - gets the information associated to the implicitMPCctrl object
    %   computeQP - computes the matrices of the QP associated to the
    %               implicit MPC problem
    %
    % The implicitMPCctrl object is derived from controller and inherits
    % all its methods.
    %
    % See also ltiSys, constraints.
    
    % Contributors:
    %
    % Alessandro Ravera (alessandro.ravera@edu.unige.it)
    % Alberto Oliveri (alberto.oliveri@unige.it)
    %
    % Copyright (C) 2021 University of Genoa, Italy.
    
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
        
        sys = [];   % dynSys object representing the system to regulate
        constr = [];    % constraints object representing the constraints
        options = [];   % MPC settings
        
    end    
    
    methods
        
        function object = implicitMPCctrl(varargin)
            
            if nargin == 0
                
                object.sys = [];
                object.constr = [];
                object.options = [];

            % Create an explicit MPC controller starting from the MPC
            % controller generated with MATLAB MPC toolbox
            % TO BE TESTED!            
            elseif (nargin == 1 || nargin == 2) && isa(varargin{1},'mpc')
                    
                mpcobj = varargin{1};
                if exist('varargin{2}')
                    options = varargin{2};
                end

                nu = size(mpcobj.ManipulatedVariables, 2);
                nd = size(mpcobj.DisturbanceVariables, 2);

                A = mpcobj.model.Plant.A;
                B = mpcobj.model.Plant.B(:, 1 : nu);
                C = mpcobj.model.Plant.C;
                D = mpcobj.model.Plant.D(:, 1 : nu);
                Fx = mpcobj.model.Plant.B(:, nu + 1 : nu + nd);
                Fy = mpcobj.model.Plant.D(:, nu + 1 : nu + nd);

                nx = size(A, 2);
                ny = size(C, 1);
                np = 0;

                Ts = mpcobj.Ts;

                object.sys = ltiSys(nx, nu, ny, np, nd, 'dt', Ts);
                object.sys = object.sys.setMatrices('A', A);
                object.sys = object.sys.setMatrices('B', B);
                object.sys = object.sys.setMatrices('C', C);
                object.sys = object.sys.setMatrices('D', D);
                object.sys = object.sys.setMatrices('Fx', Fx);
                object.sys = object.sys.setMatrices('Fy', Fy);

%                     for i = 1 : nx
%                         xnames{i} = mpcobj.StateVariables(i).Name;
%                     end
%                     object.sys = object.sys.setStateNames(xnames);
                for i = 1 : nu
                    unames{i} = mpcobj.ManipulatedVariables(i).Name;
                    umin(i) = mpcobj.ManipulatedVariables(i).Min;
                    umax(i) = mpcobj.ManipulatedVariables(i).Max;
                end
                object.sys = object.sys.setInputNames(unames);
                for i = 1 : ny
                    ynames{i} = mpcobj.OutputVariables(i).Name;
                    ymin(i) = mpcobj.OutputVariables(i).Min;
                    ymax(i) = mpcobj.OutputVariables(i).Max;
                end
                object.sys = object.sys.setOutputNames(ynames);
                for i = 1 : nd
                    dnames{i} = mpcobj.DisturbanceVariables(i).Name;
                end
                object.sys = object.sys.setUnmeasurableInputNames(dnames);

                N = mpcobj.PredictionHorizon;
                Nu = mpcobj.ControlHorizon;

                Q = diag(mpcobj.Weights.OutputVariables);

                if ~isstruct('options')
                    options = [];
                end

                if ~isfield(options,'tracking') || isempty(options.tracking)
                        options.tracking = true;
                end

                if ~isfield(options,'trackvariable') || isempty(options.trackvariable)
                        options.trackvariable = 1:ny;
                end

                if options.tracking == 0
                    R = diag(mpcobj.Weights.ManipulatedVariables);

                    if ~isfield(options,'ref') || isempty(options.ref)
                        options.ref = zeros(ny, 1);
                    end
                    if ~isfield(options,'uref') || isempty(options.uref)
                        options.ref = zeros(nu, 1);
                    end
                else
                    R = diag(mpcobj.Weights.ManipulatedVariablesRate);
                end

                nr = ny;
                constr = constraints(nx, nu, ny, np, nd, nr, N);

                constr = constr.setConstraints('u', umin, umax);
                constr = constr.setConstraints('y', ymin, ymax);

                options.ctrlvariable = 'output';
                options.Q = Q;
                options.P = Q;
                options.R = R;
                options.N = N;
                options.Nu = Nu;

                object = implicitMPCctrl(object.sys, Ts, constr, options);

            elseif nargin == 4
                    
                if ~isa(varargin{1},'ltiSys')
                    error('First argument must be a ''ltiSys'' object');
                end
                    
                if ~isnumeric(varargin{2})
                    error('Second argument must be a double');
                end
                if ~isa(varargin{3},'constraints')
                    error('Third argument must be a ''constraints'' object');
                end
                if ~isstruct(varargin{4})
                    error('Fourth argument must be a struct');
                end

                % Extract inputs
                sys = varargin{1};
                TsCtrl = varargin{2};
                constr = varargin{3};
                options = varargin{4};

                % Check inputs
                if TsCtrl <= 0
                    error('Sampling time must be > than 0');
                end

                if sys.isDiscreteTime() && sys.getSamplingTime ~= TsCtrl
                    error('Sampling time Ts is different from the system sampling time');
                end

                % If the system is continuous-time, discretize it
                if sys.isContinuousTime
                    disp(' ')
                    disp('Discretizing system...')
                    sys = sys.discretize(TsCtrl);
                end

                % Number of states, parameters, inputs and outputs
                nx = sys.getNumberOfStates();
                np = sys.getNumberOfParameters();
                nu = sys.getNumberOfInputs();
                nd = sys.getNumberOfUnmeasurableInputs();
                ny = sys.getNumberOfOutputs();

                if constr.getNumberOfStates ~= nx || ...
                        constr.getNumberOfParameters ~= np || ...
                        constr.getNumberOfInputs ~= nu || ...
                        constr.getNumberOfUnmeasurableInputs ~= nd || ...
                        constr.getNumberOfOutputs ~= ny
                    error('''dynSys'' and ''constraints'' objects are not compatible');
                end

                % System sampling time
                Ts = sys.getSamplingTime();

                % Names
                xnames = cell(nx,1);
                unames = cell(nu,1);
                ynames = cell(ny,1);
                pnames = cell(np,1);
                dnames = cell(nd,1);
                for i = 1:nx
                    xnames{i} = sys.getStateNames{i};
                end
                for i = 1:ny
                    ynames{i} = sys.getOutputNames{i};
                end
                for i = 1:nu
                    unames{i} = sys.getInputNames{i};
                end
                for i = 1:np
                    pnames{i} = sys.getParameterNames{i};
                end
                for i = 1:nd
                    dnames{i} = sys.getUnmeasurableInputNames{i};
                end

                % Set default options
                options = MPCset(sys,options);

                % Check trackvariable
                if options.tracking
                    tracking = true;
                    trackvar = options.trackvariable;
                    xreference = [];
                else
                    tracking = false;
                    trackvar = [];
                    xreference = options.ref;
                end

                ctrlvar = options.ctrlvariable;

                % Number of tracking variables
                ntrack = numel(trackvar);

                % Names of the reference states
                rnames = cell(ntrack,1);
                for i = 1:ntrack
                    rnames{i} = ['ref ',xnames{trackvar(i)}];
                end

                if options.N ~= constr.getTimeHorizon()
                    error('The prediction horizon set in the ''constraints'' object is different from OPTS.N.');
                end
                
                % Update object fields
                object.sys = sys;
                object.constr = constr;
                if options.tracking == false
                    options.trackvariable = [];
                end
                object.options = options;
                
                object.nx = nx;
                object.nu = nu;
                object.ny = ny;
                object.np = np;
                object.nd = nd;
                object.Ts = Ts;
                object.xnames = xnames;
                object.unames = unames;
                if np > 0
                    object.pnames = pnames;
                end
                if ny > 0
                    object.ynames = ynames;
                end
                if nd > 0
                    object.dnames = dnames;
                end
                
                object.tracking = tracking;
                object.trackvar = trackvar;
                object.ctrlvar = ctrlvar;
                object.ref = xreference;
                
                disp(' ')
                disp('Done.')
                disp(' ')
            else
                error('Wrong input arguments');
            end
            
        end
        
        varargout = eval(object,x,varargin);
        
        disp(object);
        
        sys = getSystem(object);
        
        constr = getConstraints(object);
        
        options = getOptions(object);
        
        info = getInformation(object);
       
        
    end
    
end
