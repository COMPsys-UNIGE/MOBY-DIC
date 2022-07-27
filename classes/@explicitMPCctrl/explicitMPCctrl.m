classdef explicitMPCctrl < controller
    % explicitMPCctrl   Explicit MPC controller
    %
    % This object represents an explicit MPC controller for a LTI
    % system. The controller can be designed either for regulation to a
    % constant reference state (or output) or for tracking an external signal. 
    % The control function u = f(x,p,d,ref) is a piecewise-affine function of
    % the system state (x), parameters (p), unmeasurable inputs (d) and,
    % only for tracking, of the reference signal (ref). The function is
    % defined over a generic polytopic domain partition (pwagFunction
    % object). The control function is designed to satisfy constraints on
    % the system variables.
    %
    % See the User's Guide for a more detailed explaination of this object.
    %
    % OBJ = explicitMPCctrl()
    % Builds an empty explicitMPCctrl object OBJ.
    %
    % OBJ = explicitMPCctrl(MPC_CTRL)
    % Builds a explicitMPCctrl object starting from a MPCController object of
    % toolbox MPT3, representing an implicit MPC controller.
    %
    % OBJ = explicitMPCctrl(EMPC_CTRL)
    % Builds a explicitMPCctrl object starting from a EMPCController object of
    % toolbox MPT3, representing an explicit MPC controller.
    %
    % OBJ = explicitMPCctrl(SYS,Ts,CONSTR,OPTS)
    % Builds a explicitMPCctrl object starting from a ltiSys object SYS
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
    % The MPC controller is designed through Yalmip and MPT3.
    %
    % explicitMPCctrl methods:
    %   computeTree - computes the binary search tree associated to the polytopic partition
    %   disp - displays some information about the explicitMPCctrl object.
    %   eval - evaluates the MPC controller.
    %   generateC - Generates C files for the circuit implementation of the MPC
    %               controller on microcontroller
    %   generateVHDL - Generates VHDL files for the circuit implementation 
    %                  of the MPC controller on FPGA
    %   getConstraints - gets the constraints object representing the constraints
    %   getInformation - gets the information associated to the explicitMPCctrl object
    %   getNumberOfRegions - gets the number of polytopes of the domain partition
    %   getOptions - gets the options structure used to design the MPC controller
    %   getSystem - gets the dynSys object representing the system to control
    %   hasTree - returns true if the binary search tree has been computed
    %   plot - plots the MPC control function.
    %   plotPartition - plots the polytopic domain partition.
    %   plotTree - plots the binary search tree associated to the domain partition.
    %
    % The explicitMPCctrl object is derived from controller and inherits
    % all its methods.
    %
    % See also pwagFunction, ltiSys, pwaSys, constraints.
    
    % Contributors:
    %
    % Alberto Oliveri (alberto.oliveri@unige.it)
    % Mateo Lodi (matteo.lodi@edu.unige.it)
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
        
        sys = [];   % dynSys object representing the system to regulate
        constr = [];    % constraints object representing the constraints
        options = [];   % MPC settings
        
    end
    
    
    methods
        
        function object = explicitMPCctrl(varargin)
            
            if nargin == 0
                
                object.sys = [];
                object.constr = [];
                object.options = [];constantInputAfterNu
                
            elseif nargin == 1 || nargin == 4
                
                % Create an explicit MPC controller starting from the MPC
                % controller generated with MATLAB MPC toolbox
                % TO BE TESTED!
                if isa(varargin{1},'mpc')
                    
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

                    object = explicitMPCctrl(object.sys, Ts, constr, options);
                end

                if isa(varargin{1},'MPCController')
                    
                    disp(' ')
                    disp('Converting MPC controller to explicit form...')
                    varargin{1} = varargin{1}.toExplicit();
                    disp('Done.')
                    disp(' ')
                    
                end
                
                if isa(varargin{1},'EMPCController')
                    
                    % TO DO
                    % Vedere il caso di controllori per l'output

                    % EMPCController object
                    ctrl = varargin{1};
                    
                    % Number of states
                    nx = ctrl.nx;
                    % Number of inputs
                    nu = ctrl.nu;
                    % Number of parameters
                    np = 0;
                    % Number of unmeasurable disturbances
                    nd = 0;constantInputAfterNu
                    % Number of outputs
                    ny = 0;
                    % Sampling time
                    Ts = ctrl.model.Ts;
                    
                    xnames = cell(nx,1);
                    unames = cell(nu,1);
                    pnames = cell(np,1);
                    dnames = cell(nd,1);
                    ynames = cell(nd,1);
                    for i = 1:nx
                        xnames{i} = ['x_',num2str(i)];
                    end
                    for i = 1:nu
                        unames{i} = ['u_',num2str(i)];
                    end
                    for i = 1:np
                        pnames{i} = ['p_',num2str(i)];
                    end
                    for i = 1:nd
                        dnames{i} = ['d_',num2str(i)];
                    end
                    for i = 1:ny
                        dnames{i} = ['y_',num2str(i)];
                    end
                    
                    % Check if the controller is for tracking
                    format = ctrl.xinitFormat.names;
                    if any(ismember(format,'x.reference'))
                        tracking = true;
                        trackvar = 1:nx;
                        rnames = cell(nx,1);
                        for i = 1:nx
                            rnames{i} = ['ref ',xnames{i}];
                        end
                    else
                        tracking = false;
                        trackvar = [];
                        rnames = cell(0,1);
                    end

                    ctrlvar = 'state';
                    
                    % TO DO
                    % Si riesce a recuperare l'informazione su ref?
                    
                    % Extract polytopes
                    polytopes = ctrl.feedback.Set;
                    
                    % Extract functions defined over polytopes
                    functions = polytopes.getFunction('primal');
                    
                    % Number of polytopes
                    nr = numel(polytopes);
                    
                    % Number of functions
                    Ftest = functions(1).F;
                    nfun = size(Ftest,1)/nu;
                    
                    % Create structure "regions"
                    regions = struct('H',cell(nr,1),'K',cell(nr,1),...
                        'F',cell(nr,1),'G',cell(nr,1));
                    
                    % Fill regions
                    for i = 1:nr
                        regions(i).H = polytopes(i).H(:,1:end);
                        regions(i).K = polytopes(i).H(:,end);
                        regions(i).F = functions(i).F;
                        regions(i).G = functions(i).g;
                    end
                    
                    % Create domain
                    HKd = ctrl.feedback.Domain.H;
                    Hd = HKd(:,1:end-1);
                    Kd = HKd(:,end);
                    domain.Hd = Hd;
                    domain.Kd = Kd;
                    
                    % Create pwagFunction object
                    fun = pwagFunction(regions,domain);
                    
                    % Assign function names
                    fun = fun.setInputNames([xnames;pnames;dnames;rnames]);
                    fun = fun.setOutputNames(repmat(unames,nfun,1));
                    
                    object.sys = [];
                    object.constr = [];
                    object.options = [];
                    
                elseif isa(varargin{1},'ltiSys')
                    
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
                    
                    if constr.getNumberOfReferences ~= ntrack
                        error('Number of reference states in ''constraints'' object is not compatible with OPTS.trackvariable');
                    end
                    
                    % Names of the reference states
                    rnames = cell(ntrack,1);
                    for i = 1:ntrack
                        rnames{i} = ['ref ',xnames{trackvar(i)}];
                    end
                    
                    if options.N ~= constr.getTimeHorizon()
                        error('The prediction horizon set in the ''constraints'' object is different from OPTS.N.');
                    end
                                           
                    % Retrieve system matrices
                    A = sys.getMatrices('A');
                    B = sys.getMatrices('B');
                    Ex = sys.getMatrices('Ex');
                    Fx = sys.getMatrices('Fx');
                    Gx = sys.getMatrices('Gx');
                    C = sys.getMatrices('C');
                    D = sys.getMatrices('D');
                    Ey = sys.getMatrices('Ey');
                    Fy = sys.getMatrices('Fy');
                    Gy = sys.getMatrices('Gy');
                    
                    % Retrieve hard constraint matrices
                    [H, K] = constr.getAllConstraints('hard');
                                        
                    % Separate constraint matrix H = [Hx Hu Hp Hd]
                    Hx = cell(options.N+1,1);
                    Hu = cell(options.N+1,1);
                    Hy = cell(options.N+1,1);
                    Hp = cell(options.N+1,1);
                    Hd = cell(options.N+1,1);
                    Hr = cell(options.N+1,1);
                    for i = 1:options.N+1
                        Hx{i} = [];
                        Hu{i} = [];
                        Hy{i} = [];
                        Hp{i} = [];
                        Hd{i} = [];
                        Hr{i} = [];
                    end
                    warn = 0;
                    % TO DO
                    % occhio perche' i vincoli sullo stato possono anche
                    % essere imposti al tempo N+1
                    
                    % ycon indicates if output constraints are imposed
                    ycon = 0;
                    for i = 1:options.N+1
                        if ~isempty(H{i})
                            if i > options.Nu
                                warn = 1;
                            end
                            Hx{i} = H{i}(:,1:nx);
                            Hu{i} = H{i}(:,nx+1:nx+nu);
                            Hp{i} = H{i}(:,nx+nu+1:nx+nu+np);
                            Hd{i} = H{i}(:,nx+nu+np+1:nx+nu+np+nd);
                            Hy{i} = H{i}(:,nx+nu+np+nd+1:nx+nu+np+nd+ny);
                            if any(Hy{i})
                                ycon = 1;
                            end
                        end
                    end

                    if warn
                        warning('Some constraints are imposed after control horizon Nu. This may cause feasibility issues.');
                        disp(' ')
                    end
                    
                    % Retrieve soft constraint matrices
                    [Hs, Ks] = constr.getAllConstraints('soft');
                    
                    % Separate constraint matrix H = [Hx Hu Hp Hd]
                    Hxs = cell(options.N+1,1);
                    Hus = cell(options.N+1,1);
                    Hys = cell(options.N+1,1);
                    Hps = cell(options.N+1,1);
                    Hds = cell(options.N+1,1);
                    Hrs = cell(options.N+1,1);
                    % Number of needed slack variables
                    ns = 0;
                    for i = 1:options.N+1
                        Hxs{i} = [];
                        Hus{i} = [];
                        Hys{i} = [];
                        Hps{i} = [];
                        Hds{i} = [];
                        Hrs{i} = [];
                    end
                    for i = 1:options.N+1
                        if ~isempty(Hs{i})
                            Hxs{i} = Hs{i}(:,1:nx);
                            Hus{i} = Hs{i}(:,nx+1:nx+nu);
                            Hps{i} = Hs{i}(:,nx+nu+1:nx+nu+np);
                            Hds{i} = Hs{i}(:,nx+nu+np+1:nx+nu+np+nd);
                            Hys{i} = Hs{i}(:,nx+nu+np+nd+1:nx+nu+np+nd+ny);
                            if any(Hys{i})
                                ycon = 1;
                            end
                        end
                    end
                    
                    % yref indicates if the controller is for output
                    % regulation (or tracking)
                    if strcmpi(options.ctrlvariable,'output')
                        yref = 1;
                    else
                        yref = 0;
                    end

                    for i = 1:options.N+1
                        if ~isempty(H{i})
                            if yref
                                Hr{i} = zeros(size(H{i},1),ny);
                                Hr{i}(:,trackvar) = H{i}(:,nx+nu+np+nd+ny+1:nx+nu+np+nd+ny+ntrack);
                            else
                                Hr{i} = zeros(size(H{i},1),nx);
                                Hr{i}(:,trackvar) = H{i}(:,nx+nu+np+nd+ny+1:nx+nu+np+nd+ny+ntrack);
                            end
                        end
                        if ~isempty(Hs{i})
                            if yref
                                Hrs{i} = zeros(size(Hs{i},1),ny);
                                Hrs{i}(:,trackvar) = Hs{i}(:,nx+nu+np+nd+ny+1:nx+nu+np+nd+ny+ntrack);
                            else
                                Hrs{i} = zeros(size(Hs{i},1),nx);
                                Hrs{i}(:,trackvar) = Hs{i}(:,nx+nu+np+nd+ny+1:nx+nu+np+nd+ny+ntrack);
                            end
                            ns = ns+size(Hs{i},1);
                        end
                    end
                    
                    ns = 1;
                    
                    disp(' ')
                    disp('Designing explicit MPC controller through MPT3...')
                    
                    % Set up problem in YALMIP
                    x = sdpvar(nx,options.N+1);
                    u = sdpvar(nu,options.N);
                    p = sdpvar(np,1);
                    d = sdpvar(nd,1);
                    slack = sdpvar(ns,1);
                    
                    if ~yref
                        ref = sdpvar(nx,1);
                    else
                        ref = sdpvar(ny,1);
                    end
                    
                    if yref || ycon
                        y = sdpvar(ny,options.N+1);
                    end
                    
                    con = [];
                    obj = 0;
                    idx = 1;
                    
                    % Check if there is some state with constant dynamics
                    M1 = [A B Ex Fx Gx];
                    M2 = [eye(nx) zeros(nx,nu+np+nd+size(Gx,2))];
                    if any(all(M1 == M2,2))
                        warning(['Same states have a constant dynamics. ',...
                            'This may lead to issues in the solution of the MPC ',...
                            'optimization problem. It is suggested to consider ',...
                            'these states as parameters or unmeasurable inputs']);
                    end
                    
                    for i = 1:options.N
                        % prediction model
                        con = [con;   x(:,i+1) == A*x(:,i)+B*u(:,i)+Ex*p+Fx*d+Gx]; %#ok<*AGROW>
                        
                        if yref || ycon
                            % outputs
                            con = [con;   y(:,i) == C*x(:,i)+D*u(:,i)+Ey*p+Fy*d+Gy];
                        end
                        
                        % inequality constraints
                        if yref || ycon
                            if ~isempty(H{i})
                                con = [con;   Hx{i}*x(:,i)+Hu{i}*u(:,i)+...
                                    Hp{i}*p+Hd{i}*d+Hy{i}*y(:,i)+Hr{i}*ref <= K{i}];
                            end
                            if ~isempty(Hs{i})
                                con = [con;   Hxs{i}*x(:,i)+Hus{i}*u(:,i)+...
                                    Hps{i}*p+Hds{i}*d+Hys{i}*y(:,i)+Hrs{i}*ref <=...
                                    Ks{i}+slack(idx:idx+size(Ks{i},1)-1)];
                                idx = idx+size(Ks{i},1);
                            end
                        else
                            if ~isempty(H{i})
                                con = [con;   Hx{i}*x(:,i)+Hu{i}*u(:,i)+...
                                    Hp{i}*p+Hd{i}*d+Hr{i}*ref <= K{i}];
                            end
                            if ~isempty(Hs{i})
                                con = [con;   Hxs{i}*x(:,i)+Hus{i}*u(:,i)+...
                                    Hps{i}*p+Hds{i}*d+Hrs{i}*ref <=...
                                    Ks{i}+slack];%(idx:idx+size(Ks{i},1)-1)];
                                idx = idx+size(Ks{i},1);
                            end
                        end
                        
                        % equality constraints (for i > Nu)
                        if i > options.Nu
                            if options.constantInputAfterNu
                                con = [con; u(:,i) == u(:,i-1)];
                            else
                                con = [con; u(:,i) == options.K*x(:,i)+options.O];
                            end
                        end
                        
                        if options.norm == 2
                            if ~yref
                                obj = obj + (x(:,i)-ref)'*options.Q*(x(:,i)-ref);
                            else
                                obj = obj + (y(:,i)-ref)'*options.Q*(y(:,i)-ref);
                            end
                            obj = obj + (u(:,i)-options.uref)'*options.R*(u(:,i)-options.uref);
                        else
                            if ~yref
                                obj = obj + norm(options.Q*(x(:,i)-ref),options.norm);
                            else
                                obj = obj + norm(options.Q*(y(:,i)-ref),options.norm);
                            end
                            obj = obj + norm(options.R*(u(:,i)-options.uref),options.norm);
                        end
                    end
                    
                    % last output
                    if yref || ycon
                        if any(D)
                            error(['If output constraints are imposed or',...
                                ' the controller is for output regulation'...
                                ' (or tracking), matrix D must be null']);
                        end
                        con = [con;   y(:,options.N+1) == C*x(:,options.N+1)+Ey*p+Fy*d+Gy];
                    end

                    % Terminal cost
                    if options.norm == 2
                        if ~yref
                            obj = obj + (x(:,options.N+1)-ref)'*options.P*(x(:,options.N+1)-ref);
                        else
                            obj = obj + (y(:,options.N+1)-ref)'*options.P*(y(:,options.N+1)-ref);
                        end
                    else
                        if ~yref
                            obj = obj + norm(options.P*(x(:,options.N+1)-ref),options.norm);
                        else
                            obj = obj + norm(options.P*(y(:,options.N+1)-ref),options.norm);
                        end
                    end
                    
                    pars = getMOBYDICpars();
                                   
                    % Put together all non zero elements of matrices P, Q and R
                    vals = [options.Q(:); options.R(:); options.P(:)];
                    vals(vals == 0) = [];
                    
                    ratio = pars.slack/min(vals);
                    
                    if ratio > 1e4
                        warning(['The ratio between global option MOBYDICpars.slack ',...
                            'and the minimum value in matrices P, Q and R is greater ',...
                            'than 10000. This may cause numerical problem in the ',...
                            'solution of the MPC optimization problem']);
                    end
                            
                    
                    % Add penalty on slack variables
                    if ns > 0
                        if options.norm == 2
                            obj = obj + slack'*pars.slack*eye(ns)*slack;
%                             con = [con; slack >= zeros(ns,1)];
                        else
                            con = [con; slack >= zeros(ns,1)];
                            obj = obj + norm(pars.slack*eye(ns)*slack,options.norm);
                        end
                    end
                    
                    % Terminal constraint
                    if yref || ycon
                        if ~isempty(H{options.N+1})
                            con = [con;   Hx{options.N+1}*x(:,options.N+1)+...
                                Hp{options.N+1}*p+Hd{options.N+1}*d+...
                                Hy{options.N+1}*y(:,options.N+1)+Hr{options.N+1}*ref <= K{options.N+1}];
                        end
                    else
                        if ~isempty(H{options.N+1})
                            con = [con;   Hx{options.N+1}*x(:,options.N+1)+...
                                Hp{options.N+1}*p+Hd{options.N+1}*d+...
                                Hr{options.N+1}*ref <= K{options.N+1}];
                        end
                    end
                    
                    % Fix the elements in ref which do not correspond to
                    % trackvar
                    if options.tracking
                        if yref
                            for i = 1:ny
                                idx = i == trackvar;
                                if all(idx == 0)
                                    con = [con; ref(i) == options.ref(i)];
                                end
                            end
                        else
                            for i = 1:nx
                                idx = i == trackvar;
                                if all(idx == 0)
                                    con = [con; ref(i) == options.ref(i)];
                                end
                            end
                        end
                    else
                        con = [con; ref == options.ref];
                    end
                                       
                    % Create MPT3 object Opt
                    if options.tracking
                        problem = Opt(con, obj, [x(:,1);p;d;ref(trackvar)], u(:));
                    else
                        problem = Opt(con, obj, [x(:,1);p;d], u(:));
                    end
                    
                    % Find explicit solution
                    try
                        res = problem.solve;
                    catch err
                        disp(' ')
                        disp('MPT3 failed to solve the optimization problem')
                        disp('and returned the following status:')
                        disp(err.message)
                        error('Failed!')
                    end
                    
                    
                    if ~strcmpi(res.how,'ok')
                        disp(' ')
                        disp('MPT3 failed to solve the optimization problem')
                        disp('and returned the following status:')
                        disp(res.how)
                        error('Failed!')
                    end
                    
                    disp(' ')
                    disp('Creating explicitMPCctrl object...')
                    
                    % Extract PolyUnion object
                    union = res.xopt;
                    
                    if options.tracking
                        if union.Dim ~= nx+np+nd+numel(trackvar)
                            error('PolyUnion object is not compatible with model')
                        end
                    else
                        if union.Dim ~= nx+np+nd
                            error('PolyUnion object is not compatible with model')
                        end
                    end
                    
                    % Extract polytopes
                    polytopes = union.Set;
                    
                    % Extract functions defined over polytopes
                    functions = polytopes.getFunction('primal'); 
                    
                    % Number of polytopes
                    nr = numel(polytopes);
                    
                    % Number of functions
                    Ftest = functions(1).F;
                    nfun = size(Ftest,1)/nu;
                    
                    % Create structure "regions"
                    regions = struct('H',cell(nr,1),'K',cell(nr,1),...
                        'F',cell(nr,1),'G',cell(nr,1));
                    
                    % Fill regions
                    for i = 1:nr
                        regions(i).H = polytopes(i).H(:,1:end-1);
                        regions(i).K = polytopes(i).H(:,end);
                        regions(i).F = functions(i).F;
                        regions(i).G = functions(i).g;
                    end
                    
                    % Set domain
                    HKd = union.Domain.H;
                    Hd = HKd(:,1:end-1);
                    Kd = HKd(:,end);
                    domain.Hd = Hd;
                    domain.Kd = Kd;
                    
                    % Create pwagFunction object
                    fun = pwagFunction(regions,domain);
                    
                    % Assign function names
                    fun = fun.setInputNames([xnames;pnames;dnames;rnames]);
                    fun = fun.setOutputNames(repmat(unames,nfun,1));
                    
                    object.sys = sys;
                    object.constr = constr;
                    if options.tracking == false
                        options.trackvariable = [];
                    end
                    object.options = options;
                    
                end
                
                % Update object fields
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
                
                object.nfun = nfun;
                
                object.fun = fun;
                
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
        
        object = computeTree(object,varargin);
        
        varargout = eval(object,x,varargin);
        
        disp(object);
        
        plot(object,k);
        
        plotPartition(object);
        
        nreg = getNumberOfRegions(object);
        
        answ = hasTree(object);
        
        plotTree(object);
        
        sys = getSystem(object);
        
        constr = getConstraints(object);
        
        options = getOptions(object);
        
        info = getInformation(object);
       
        
    end
    
end
