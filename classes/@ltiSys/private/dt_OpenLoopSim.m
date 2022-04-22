function signals = dt_OpenLoopSim(object,x0,p,d,time,options)
% dt_OpenLoopSim    Simulates the discrete-time open-loop dynamical system
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



% Number of iterations
Nit = floor(time/object.Ts);

% Rounded stop time
time = object.Ts*Nit;

% Start and stop time for simulation
tstart = 0;
tend = time;

% Create vector T
T = tstart:object.Ts:tend;
T = T';

% Number of iterations
nit = numel(T);

% Create vector X
X = zeros(nit,object.nx);

% Set vector U
U = zeros(nit,object.nu);

% initial condition
X(1,:) = x0;

h = waitbar(0,'Simulating...');

for k = 2:nit
    X(k,:) = object.evaluate(X(k-1,:),U(k,:),p,d);
    waitbar(k/nit);
end
close(h)

% Compute system outputs
Y = object.computeOutput(X,U,p,d);

signals.time = T;
signals.state = X;
if object.nu > 0
    signals.input = U;
end
signals.output = Y;

end







