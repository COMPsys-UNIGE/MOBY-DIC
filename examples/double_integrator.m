%% Double integrator system
% MPC is applied for the regulation of a double integrator system

% Define system matrices
A = [ 0, 1; 0, 0 ];
  
B = [ 0; 1 ];
    
C = [ 1, 0 ];

D = 0;

% Create a continuous-time LTIsys object
nx = 2;    % Number of states
nu = 1;    % Number of inputs
ny = 1;    % Number of outputs
np = 0;    % Number of parameters
nd = 0;    % Number of unmeasurable inputs
nxref = 1; % Number of reference states

% Create continuous-time ltiSys object
ctsys = ltiSys(nx, nu, ny, np, nd, 'ct');

% Set system matrices
ctsys = ctsys.setMatrices('A', A);
ctsys = ctsys.setMatrices('B', B);
ctsys = ctsys.setMatrices('C', C);
ctsys = ctsys.setMatrices('D', D);

% Set name of variables:
% - the first element of the state is the position 
% - the second element of the state is the velocity 
% - the input is the acceleration 
% - the output is the position 
ctsys = ctsys.setStateNames({'s', 'v'});
ctsys = ctsys.setInputNames({'a'});
ctsys = ctsys.setOutputNames({'s'});

%% MPC controller design
% Select MPC type (implicit, explicit, approximate)
mpcType = "implicit";

% Define the sampling time
Ts = 1;

N = 5;  % Prediction horizon
Nu = 3; % Control horizon

% MPC weight matrices
P = [ 1, 0; 0, 0 ];
Q = [ 1, 0; 0, 0 ];
R = 10;

% Constraints
xmin = [ -8, -1 ];
xmax = [ 8, 1 ];
umin = -0.3;
umax = 0.3;

% Create constraints object
constr = constraints(nx, nu, ny, np, nd, nxref, N);

% Set state constraints over the control horizon
constr = constr.setConstraints('x', xmin, xmax, 1:Nu-1);

% Set input constraints at current time instant
constr = constr.setConstraints('u', umin, umax, 0);

% Explicit MPC requires box contraints also on reference
rmin = xmin(1);
rmax = xmax(1);
constr = constr.setConstraints('r', rmin, rmax, 0);

% Controller options
options = struct( ...
             'P', P,          ... % matrix P of the cost function
             'Q', Q,          ... % matrix Q of the cost function
             'R', R,          ... % matrix R of the cost function
             'N', N,          ... % prediction horizon
            'Nu', Nu,         ... % control horizon
             'K', [],         ... % u = Kx+O, after the control horizon
             'O', [],         ... % u = Kx+O, after the control horizon
      'tracking', true,       ... % regulation or tracking problem
 'trackvariable', 1,          ... % name or index of the variable to track (only for tracking)
          'xref', [],         ... % constant reference state (only for regulation)
          'uref', []          ... % constant reference input (only for regulation)
    );

% Design controller
switch mpcType
    
    case "explicit"
        ctrl = explicitMPCctrl(ctsys, Ts, constr, options);
        
    case "approximate"
        ctrl = explicitMPCctrl(ctsys, Ts, constr, options);
        approxOpts.nsamples = 4;
        approxOpts.np = [8 9 8];
        ctrl = ApproxMPCctrl(ctrl,approxOpts);
    
    otherwise % implicit
        ctrl = implicitMPCctrl(ctsys, Ts, constr, options);
end

% % If explicit or approximate MPC with nx+np+nd<=2, control function
% % and polytopic partition can be plotted
% if mpcType ~= "implicit"
%     % Plot control function
%     ctrl.plot;
%     
%     % Plot polytopic partition
%     ctrl.plotPartition;
% end

% Associate controller to LTI system
ctSysReg = ctsys.setController(ctrl);


%% Closed-loop simulation with MPC controller

% Define an initial condition 
x0 = [ 0; 0 ];

% Define a reference position
sref = 5;

% Set a simulation time
T = 30;

% Set constraints for simulation
opts.constraints = constr;

% Closed-loop simulation (continuous time)
ctSysReg.simplot(T, x0, sref, opts);


%% Observer design

% Select the observer type (predictor or filter)
obsType = "predictor";

% Observer sampling time
Tobs = 0.5;

% Noise covariance matrices
Q = 0.1*eye(nx);
R = 0.1*eye(ny);

% Design observer
switch obsType
    
    case "predictor"
        obs = kalmanPredictor(ctsys, Tobs, Q, R);
    
    otherwise % filter
        obs = kalmanFilter(ctsys, Tobs, Q, R);
end

% Associate observer to LTI system
ctSysObs = ctSysReg.setObserver(obs);

%% Closed-loop simulation with MPC controller and observer

% Set initial observer state (optional)
opts.xest0 = x0;

% Closed-loop simulation
ctSysObs.simplot(T, x0, sref, opts);


%% C and VHDL code generation
% The toolbox allows for:
% - generating circuit implementation of controller or observer or controller+observer on microcontroller (C code)
% - generating circuit implementation of controller or observer or controller+observer on FPGA (VHDL code)
% - generating Simulink model for closed-loop simulation of controller or controller+observer
% - generating Simulink model for hardware-in-the-loop simulation of controller or controller+observer (the system runs on MATLAB, the controller/observer in an FPGA)

% Set the signals range for circuit implementation
range.xmin = xmin;
range.xmax = xmax;
range.umin = umin;
range.umax = umax;
range.pmin = [];
range.pmax = [];
range.dmin = [];
range.dmax = [];
range.xrefmin = xmin(1);
range.xrefmax = xmax(1);
range.ymin = xmin(1);
range.ymax = xmax(1);

% Set options for a microcontroller implementation.
% If the MPC controller is implicit, the parameters of Alternating Direction Method of Multipliers (ADMM) algorithm must be set.
if mpcType == "implicit"
    
    ADMMparameters.regPar = 32;   % Regularization parameter
    ADMMparameters.maxIter = 30; % Number of iterations
    
    cOptions = struct( ...
               'range', range,                  ... % signals range
               'ADMMparameters', ADMMparameters ... % parameters for ADMM algorithm
        );

else
    
    cOptions = struct( ...
               'range', range ... % signals range
        );
    
end

% Set options for an FPGA implementation.
% If the MPC controller is implicit, some additional options have to be set:
% - the integer part of the fixed point representation
% - the FPGA device
% - the parameters of ADMM
if mpcType == "implicit"
    
    vhdlOptions = struct( ...
               'architecture', 'fast',            ... % hardware architecture (can be 'fast' or 'small')
               'inputResolution', 12,             ... % resolution of the inputs
               'inputRepresentation', 'unsigned', ... % representation of the inputs
               'coeffResolution', 18,             ... % resolution of the coefficients
               'coeffIntResolution', 6,           ... % integer part of the coefficients (only for implicit MPC)
               'outputResolution', 12,            ... % resolution of the outputs
               'outputRepresentation', 'unsigned',... % representation of the outputs
               'frequency', 1e5,                  ... % FPGA working frequency
               'range', range,                    ... % signals range
               'fpgaBoard', 'xc7z020clg484-1',    ... % target FPGA (only for implicit MPC)
               'ADMMparameters', ADMMparameters   ... % parameters for ADMM algorithm (only for implicit MPC)
        );

else
    
    vhdlOptions = struct( ...
               'architecture', 'fast',            ... % hardware architecture (can be 'fast' or 'small')
               'inputResolution', 12,             ... % resolution of the inputs
               'inputRepresentation', 'unsigned', ... % representation of the inputs
               'coeffResolution', 12,             ... % resolution of the coefficients
               'outputResolution', 12,            ... % resolution of the outputs
               'outputRepresentation', 'unsigned',... % representation of the outputs
               'frequency', 1e5,                  ... % FPGA working frequency
               'range', range                     ... % signals range
        );
    
end

% Create an embeddedSystem object (controller and observer)
es = embeddedSystem(ctSysObs, range);

% Select whether to generate only controller (ctrl)
% or both controller and observer (es)
obj = ctrl;

if isa(obj, "controller")
    sys = ctSysReg;
else % embeddedSystem
    sys = ctSysObs;
end

%% Generate controller/observer C files
obj.generateC(cOptions);

%% Generate controller/observer VHDL files
obj.generateVHDL(vhdlOptions);

%% Generate Simulink model
sys.generateSimulinkModel();

%% Generate Simulink model for hardware-in-the-loop simulation (to run this section you must close MATLAB and open "Model Composer and System Generator")
% Toggle 'simVHDL' for hardware-in-the-loop simulation
simulinkOptions = struct( ...
                'simVHDL', 1 ...
    );

sys.generateSimulinkModel(simulinkOptions, vhdlOptions);
