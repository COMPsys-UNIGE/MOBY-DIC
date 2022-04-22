function signals = dt_ClosedLoopSim(object,x0,p,d,xref,time,options)
% dt_ClosedLoopSim    Simulates the discrete-time closed-loop 
%                     dynamical system
%
% This is a private method called by method ltiSys/sim.
%
% See also: ltiSys/sim.

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


% Controller
ctrl = object.controller;

if ctrl.isTracking && isempty(xref)
    error('xref must be provided for tracking controllers');
end

% Sampling time of the controller
Tsctrl = ctrl.getSamplingTime();

% Sampling time of the system
Tssys = object.getSamplingTime();

% Time instants in which the controller is applied
Tctrl = unique([0:Tsctrl:time time]);

% Time instants in which the system is integrated, within the controller
% sampling time
Tsys = 0:Tssys:Tsctrl;

% Number of time instants where controller is applied
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
    
    % Input corresponding to initial condition
    if ctrl.isTracking
        u = ctrl.eval(x0cur,p,d(:)',xref);
    else
        u = ctrl.eval(x0cur,p,d(:)');
    end
    
    if isnan(u)
        warning('MPC controller returned NaN. State is possibly outside the feasibility domain');
        break;
    end
       
    % Times
    Tcur = Tctrl(i)+Tsys;
    Tcur = Tcur';
    
    % Remove times over maximum time
    Tcur(Tcur > time) = [];
    
    % Number of iterations
    Nit = numel(Tcur);
    
    % Initialize vector containing state evolution
    Xcur = zeros(Nit,object.nx);
    
    % Current initial condition
    Xcur(1,:) = x0cur;
    
    % Iterate system
    for k = 2:Nit
        Xcur(k,:) = object.evaluate(Xcur(k-1,:),u,p,d);
    end
    
    x0cur = Xcur(end,:);
    
    % Number of points in which solution is provided
    NT = numel(Tcur);
    
    % Replicate input in order to have it for all time instants
    Ucur = repmat(u,NT,1);
       
    % Compute system outputs
    Ycur = object.computeOutput(Xcur,Ucur,p,d);
    
    % Update times, states and inputs
    X = [X; Xcur];
    Y = [Y; Ycur];
    U = [U; Ucur];
    T = [T; Tcur];
    
    waitbar(i/Nctrl)
    
end

close(h)

signals.time = T;
signals.state = X;
if object.nu > 0
    signals.input = U;
end
signals.output = Y;

end







