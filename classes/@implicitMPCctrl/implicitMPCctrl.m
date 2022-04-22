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
    % - P, Q, R: matrices defining the MPC cost function (default: P = Q)
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
    %
    % The implicitMPCctrl object is derived from controller and inherits
    % all its methods.
    %
    % See also ltiSys, constraints.
    
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
                object.options = [];constantInputAfterNu
                
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