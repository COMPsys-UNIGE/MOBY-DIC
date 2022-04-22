function signals = dt_OpenLoopObsSim(object,x0,p,d,time,options)
% dt_OpenLoopObsSim    Simulates the discrete-time open-loop 
%                      dynamical system with observer
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

% Observer
obs = object.observer;

% Sampling time of the observer
Tsobs = obs.getSamplingTime();

% Sampling time of the system
Tssys = object.getSamplingTime();

% Time instants in which the observer is applied
Tobs = unique([0:Tsobs:time time]);

% Time instants in which the system is integrated, within the controller
% sampling time
Tsys = 0:Tssys:Tsobs;

% Number of time instants where observer is applied
Nobs = numel(Tobs)-1;

% Current initial condition
x0cur = x0;

% Control action (equal to 0 since this is a open-loop
% simulation)
if object.nu == 0
    u = [];
else
    u = zeros(1,object.nu);
end

% Actual oputput
y0 = object.computeOutput(x0,u,p,d);
y0cur = y0;

% Initialize times, states and inputs
T = [];
X = [];
U = [];
Y = [];
Xest = [];
Yest = [];

xest = x0est;

h = waitbar(0,'Simulating...');

% Simulation loop
for i = 1:Nobs
    
    % Current estimated state and output
    [xestupd, yestupd] = obs.computeOutput(xest,u,p,y0cur);
    
    % Predict state at following time instant based on current
    % measurements
    xestpred = obs.predict(xest,u,p,y0cur);  
       
    % Times
    Tcur = Tobs(i)+Tsys;
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
    
    % Replicate predicted state in order to have it for all time instants
    Xestcur = repmat(xestupd',NT,1);
    Yestcur = repmat(yestupd',NT,1);  
    
    % Replicate input in order to have it for all time instants
    Ucur = repmat(u,NT,1);
    
    % Compute system outputs
    Ycur = object.computeOutput(Xcur,Ucur,p,d);
    
    % Initial conditions at next iteration
    x0cur = Xcur(end,:);
    y0cur = Ycur(end,:);
       
    % Update times, states and inputs
    X = [X; Xcur];
    Y = [Y; Ycur];
    U = [U; Ucur];
    T = [T; Tcur];
    Xest = [Xest; Xestcur];
    Yest = [Yest; Yestcur];
    
    xest = xestpred;
    
    waitbar(i/Nobs)

end

close(h)

signals.time = T;
signals.state = X;
if object.nu > 0
    signals.input = U;
end
signals.output = Y;
signals.est_state = Xest;
signals.est_output = Yest;

end







