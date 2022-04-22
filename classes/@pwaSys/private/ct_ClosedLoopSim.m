function signals = ct_ClosedLoopSim(object,x0,p,d,xref,time,options)
% ct_ClosedLoopSim    Simulates the continuous-time closed-loop 
%                     dynamical system
%
% This is a private method called by method ltiSys/sim.
%
% See also: ltiSys/sim.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
%
% Copyright (C) 2015 University of Genoa, Italy.

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


ODEsolver = options.ODEsolver;
ODEoptions = options.ODEoptions; %#ok<NASGU>

% Current initial condition
x0cur = x0; %#ok<NASGU>

% Controller
ctrl = object.controller;

if ctrl.isTracking && isempty(xref)
    error('xref must be provided for tracking controllers');
end

% Sampling time of the controller
Tsctrl = ctrl.getSamplingTime();

% Time instants in which control function is applied
Tctrl = unique([0:Tsctrl:time time]);

% Number of time instants where control function is applied
Nctrl = numel(Tctrl)-1;

% Current initial condition
x0cur = x0;

% Initialize times, states and inputs
T = [];
X = [];
U = [];
Y = [];

h = waitbar(0,'Simulating...');

% Simulation loop
for i = 1:Nctrl
    
    % Start and stom time for simulation
    tstart = Tctrl(i); %#ok<NASGU>
    tend = Tctrl(i+1); %#ok<NASGU>
    
    % Control action
    if ctrl.isTracking
        u = ctrl.eval(x0cur,p,d,xref);
    else
        u = ctrl.eval(x0cur,p,d); 
    end
    
    if isnan(u)
        warning('MPC controller returned NaN. State is possibly outside the feasibility domain');
        break;
    end
    
    % String representing the call to the ODE solver
    command = [ODEsolver,'(@sys,[tstart tend],x0cur,ODEoptions);'];
    
    % Call ODE solver
    [Tcur,Xcur] = eval(command);
    
    % Number of points in which solution is provided
    NT = numel(Tcur);
    
    % Replicate input in order to have it for all time instants
    Ucur = repmat(u,NT,1);
    
    % Compute system outputs
    Ycur = object.computeOutput(Xcur,Ucur,p,d);
    
    % Initial conditions at next iteration
    x0cur = Xcur(end,:);
    
    % Update times, states and inputs
    T = [T; Tcur];
    X = [X; Xcur];
    Y = [Y; Ycur];
    U = [U; Ucur];
    
    waitbar(i/Nctrl)
    
end

close(h)

signals.time = T;
signals.state = X;
if object.nu > 0
    signals.input = U;
end
signals.output = Y;


    function dxdt = sys(t,x) %#ok<INUSL,DEFNU>
        
        % Compute derivative
        dxdt = object.evaluate(x,u,p,d);
    end

end







