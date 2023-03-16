%% DC-DC buck converter
% MPC is applied for the regulation of a DC-DC buck converter

% Define system parameters
vin = 6;       % input voltage [V]
L = 47e-6;     % inductance [H]
RL = 400e-3;   % inductor resistance [ohm]
Cap = 1320e-6; % capacitance [F]
vd = 0.9;      % diode forward voltage [V]
f = 30e3;      % switching frequency

% Define system matrices of the linear averaged model
A = [-RL/L, -1/L; 1/Cap, 0];
B = -[(vin + vd)/L; 0];
C = [0, 1];
D = 0;
Ex = [0; -1/Cap];
Gx = [-vd/L; 0] - B;

% Create a continuous-time LTIsys object
nx = 2;    % number of states (inductor current, output voltage)
nu = 1;    % number of inputs (duty cycle)
ny = 1;    % number of outputs (output voltage)
np = 1;    % number of parameters (load current)
nd = 0;    % number of unmeasurable inputs
nxref = 1; % number of reference states (output voltage)

% Create continuous-time ltiSys object
ctsys = ltiSys(nx,nu,ny,np,nd,'ct');

% Set system matrices
ctsys = ctsys.setMatrices('A',A);
ctsys = ctsys.setMatrices('B',B);
ctsys = ctsys.setMatrices('C',C);
ctsys = ctsys.setMatrices('D',D);
ctsys = ctsys.setMatrices('Ex',Ex);
ctsys = ctsys.setMatrices('Gx',Gx);

% Set name of variables:
% - the first element of the state is the inductor current 
% - the second element of the state is the output voltage 
% - the input is the duty cycle
% - the output is the output voltage
% - the tunable parameter is the load current
ctsys = ctsys.setStateNames({'i','vout'});
ctsys = ctsys.setInputNames({'d'});
ctsys = ctsys.setOutputNames({'vout'});
ctsys = ctsys.setParameterNames({'iout'});

%% MPC controller design

% Select MPC type (implicit, explicit, approximate)
mpcType = "explicit";

% Define the sampling time
Ts = 1/f;

% Prediction horizon
N = 5;
% Control horizon
Nu = 3;
% Constraints horizon
Nc = 1;

% Constraints
imin = 0;              % minimum inductor current
imax = 3;              % maximum inductor current
voutmin = -0.1;        % minimum output voltage
voutmax = 3.75;        % maximum output voltage
dmin = 0.2;            % minimum duty cycle
dmax = 0.8;            % maximum duty cycle
ioutmin = 0;           % minimum load current
ioutmax = 4.5;         % maximum load current
xmin = [imin voutmin]; % state bounds
xmax = [imax voutmax]; % state bounds
umin = dmin;           % input bounds
umax = dmax;           % input bounds
rmin = voutmin;        % reference bounds
rmax = voutmax;        % reference bounds
pmin = ioutmin;        % parameter bounds
pmax = ioutmax;        % parameter bounds

% Create constraints object
constr = constraints(nx, nu, ny, np, nd, nxref, N);

% Set hard input constraints
constr = constr.setConstraints('u', umin, umax, 0);

% Set hard state constraints
H = [[1, -umax/(2*f*L), 0], zeros(1,np+ny+nxref)];
K = xmax(1) - (umax/(2*f*L))*vin;
constr = constr.setConstraints(H, K, 1);

% Explicit MPC requires box constraints on every variable
constr = constr.setConstraints('r', rmin, rmax, 0);
constr = constr.setConstraints('p', pmin, pmax, 0);
constr = constr.setConstraints('x1', xmin(1), xmax(1), 0);
constr = constr.setConstraints('x2', xmin(2), xmax(2), 0);

% MPC weight matrices
P = [ 0, 0; 0, 2 ];
Q = [ 0, 0; 0, 2 ];
R = 1;

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
       'trackvariable', 2,          ... % name or index of the variable to track (only for tracking)
                'xref', [],         ... % constant reference state (only for regulation)
                'uref', [],         ... % constant reference input (only for regulation)
'constantInputAfterNu', true        ... % u(i) = u(Nu), i=Nu+1,...,N
    );

% Design controller
switch mpcType    
    case "explicit"
        ctrl = explicitMPCctrl(ctsys, Ts, constr, options);
        
    case "approximate"
        ctrl = explicitMPCctrl(ctsys, Ts, constr, options);
        approxOpts.nsamples = 4;
        approxOpts.np = [9 9 1 9];
        ctrl = ApproxMPCctrl(ctrl,approxOpts);
    
    otherwise % implicit
        ctrl = implicitMPCctrl(ctsys, Ts, constr, options);
end

% Associate controller to LTI system
ctSysReg = ctsys.setController(ctrl);

%% Closed-loop simulation with MPC controller

% Define an initial condition
x0 = [0 0]';

% Load current
iout = 0.8;

% Reference output voltage
vref = 1.8;

% Simulation time
T = 10e-3;

% Set constraints for simulation
opts.constraints = constr;

% Closed-loop simulation (continuous time)
ctSysReg.simplot(T, x0, iout, [], vref, opts);

%% Observer design

% Select the observer type (predictor or filter)
obsType = "predictor";

% Observer sampling time
Tobs = Ts;

% Noise covariance matrices
Q = 0.1*eye(nx+nd);
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
ctSysObs.simplot(T, x0, iout, [], vref, opts);

%% C and VHDL code generation
% The toolbox allows for:
% - generating circuit implementation of controller or observer or 
%   controller+observer on microcontroller (C code)
% - generating circuit implementation of controller or observer or 
%   controller+observer on FPGA (VHDL code)
% - generating Simulink model for closed-loop simulation of controller or
%   controller+observer
% - generating Simulink model for hardware-in-the-loop simulation of 
%   controller or controller+observer (the system runs on MATLAB, the 
%   controller/observer in an FPGA)

% Set the signals range for circuit implementation
range.xmin = xmin;
range.xmax = xmax;
range.umin = umin;
range.umax = umax;
range.pmin = pmin;
range.pmax = pmax;
range.dmin = [];
range.dmax = [];
range.xrefmin = rmin;
range.xrefmax = rmax;
range.ymin = xmin(2);
range.ymax = xmax(2);

% Set options for a microcontroller implementation.
% If the MPC controller is implicit, the parameters of Alternating 
% Direction Method of Multipliers (ADMM) algorithm must be set.
if mpcType == "implicit"    
    ADMMparameters.regPar = 2;   % Regularization parameter
    ADMMparameters.maxIter = 40; % Number of iterations
    
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
               'coeffIntResolution', 7,           ... % integer part of the coefficients (only for implicit MPC)
               'outputResolution', 12,            ... % resolution of the outputs
               'outputRepresentation', 'unsigned',... % representation of the outputs
               'frequency', 5e7,                  ... % FPGA working frequency
               'useADC', 1,                       ... % manage ADC scalings
               'useDAC', 1,                       ... % manage DAC scalings
               'range', range,                    ... % signals range
               'fpgaBoard', 'xc7z020clg484-1',    ... % target FPGA (only for implicit MPC)
               'ADMMparameters', ADMMparameters,  ... % parameters for ADMM algorithm (only for implicit MPC)
               'defaultOutput', 0.5               ... % output initial value
        );
else
    vhdlOptions = struct( ...
               'architecture', 'fast',            ... % hardware architecture (can be 'fast' or 'small')
               'inputResolution', 12,             ... % resolution of the inputs
               'inputRepresentation', 'unsigned', ... % representation of the inputs
               'coeffResolution', 12,             ... % resolution of the coefficients
               'outputResolution', 12,            ... % resolution of the outputs
               'outputRepresentation', 'unsigned',... % representation of the outputs
               'frequency', 5e7,                  ... % FPGA working frequency
               'useADC', 1,                       ... % manage ADC scalings
               'useDAC', 1,                       ... % manage DAC scalings
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

%% Generate Simulink model for hardware-in-the-loop simulation
% (to run this section you must close MATLAB and open Model Composer and System Generator)
% Toggle 'simVHDL' for hardware-in-the-loop simulation
simulinkOptions = struct( ...
                'simVHDL', 1 ...
    );

sys.generateSimulinkModel(simulinkOptions, vhdlOptions);