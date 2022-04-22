function signals = dt_ClosedLoopObsSim(object,x0,p,d,xref,time,options)
% dt_ClosedLoopObsSim    Simulates the discrete-time closed-loop
%                        dynamical system with observer
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


if ~isfield(options,'xest0')
    options.xest0 = zeros(object.nx,1);
end
if ~isfield(options,'dest0')
    options.dest0 = zeros(object.nd,1);
end

if numel(options.xest0) ~= object.nx
    error(['options.xest0 must be a vector with ',num2str(object.nx),' elements']);
end
if numel(options.dest0) ~= object.nd
    error(['options.dest0 must be a vector with ',num2str(object.nd),' elements']);
end

x0est = [options.xest0(:);options.dest0(:)];

% Controller
ctrl = object.controller;

if ctrl.isTracking && isempty(xref)
    error('xref must be provided for tracking controllers');
end

% Observer
obs = object.observer;

if isa(obs,'kalmanFilter')
    if any(object.D ~= 0)
        error(['Matrix D is non-zero and a Kalman fIlter is used.',...
            'This leads to algebraic loops. Consider using a Kalman predictor.'])
    end
end

% Sampling time of the observer
Tsobs = obs.getSamplingTime();

% Sampling time of the controller
Tsctrl = ctrl.getSamplingTime();

% Sampling time of the system
Tssys = object.getSamplingTime();

% Time instants in which the controller is applied
Tctrl = unique([0:Tsctrl:time time]);

% Time instants in which the observer is applied, within the controller
% sampling time
Tobs = 0:Tsobs:Tsctrl;

% Time instants in which the system is integrated, within the observer
% sampling time
Tsys = 0:Tssys:Tsobs;

% Number of time instants where controller is applied
Nctrl = numel(Tctrl)-1;

% Number of time instants where observer is applied
Nobs = numel(Tobs)-1;

% Current initial condition
x0cur = x0;

% Initialize control at 0
u = zeros(1,object.nu);

% Actual oputput (notice that the control is 0 but the output should not
% depend on u, otherwise an algebraic loop will occur)
y0 = object.computeOutput(x0,u,p,d);
y0cur = y0;

xest = x0est;

% Initialize times, states and inputs
T = [];
X = [];
U = [];
Y = [];
Xest = [];
Yest = [];

% Current estimated state and output (notice that the the first time control is 0 but the
% estimated state should not depend on u, otherwise an algebraic loop will occur)
[xestupd, yestupd] = obs.computeOutput(xest,u,p,y0cur);

xctrl = xestupd(1:object.nx);
dctrl = xestupd(object.nx+1:end);

% Control action
if ctrl.isTracking
    u = ctrl.eval(xctrl(:)',p,dctrl(:)',xref);
else
    u = ctrl.eval(xctrl(:)',p,dctrl(:)');
end

if isnan(u)
    warning('MPC controller returned NaN. State is possibly outside the feasibility domain');
end

h = waitbar(0,'Simulating...');

% Simulation loop
for i = 1:Nctrl
    
    for j = 1:Nobs
        
        % Predict state at following time instant based on current
        % measurements
        xestpred = obs.predict(xest,u,p,y0cur);
        
        % Times
        Tcur = Tctrl(i)+Tobs(j)+Tsys;
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
        
        % Number of points in which solution is provided
        NT = numel(Tcur);
        
        % Replicate input in order to have it for all time instants
        Ucur = repmat(u,NT,1);
        
        % Replicate predicted state in order to have it for all time instants
        Xestcur = repmat(xestupd',NT,1);
        Yestcur = repmat(yestupd',NT,1);
        
        % Compute system outputs
        Ycur = object.computeOutput(Xcur,Ucur,p,d);
        
        % Initial conditions at next iteration
        % Initial conditions at next iteration
        x0cur = Xcur(end,:);
        y0cur = Ycur(end,:);
        
        % Update times, states, inputs, predicted states and predicted
        % outputs
        T = [T; Tcur];
        X = [X; Xcur];
        Y = [Y; Ycur];
        U = [U; Ucur];
        Xest = [Xest; Xestcur];
        Yest = [Yest; Yestcur];
        
        xest = xestpred;
        
        % Current estimated state and output (notice that the the first time control is 0 but the estimated state should not depend on u, otherwise an algebraic loop will occur)
        [xestupd, yestupd] = obs.computeOutput(xest,u,p,y0cur);
        
    end
    
    xctrl = xestupd(1:object.nx);
    dctrl = xestupd(object.nx+1:end);
    
    % Control action
    if ctrl.isTracking
        u = ctrl.eval(xctrl(:)',p,dctrl(:)',xref);
    else
        u = ctrl.eval(xctrl(:)',p,dctrl(:)');
    end
    
    if isnan(u)
        warning('MPC controller returned NaN. State is possibly outside the feasibility domain');
        break;
    end
    
    waitbar(i/Nctrl)
    
end

close(h);

signals.time = T;
signals.state = X;
if object.nu > 0
    signals.input = U;
end
signals.output = Y;
signals.est_state = Xest;
signals.est_output = Yest;

end







