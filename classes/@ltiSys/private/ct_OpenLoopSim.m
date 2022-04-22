function signals = ct_OpenLoopSim(object,x0,p,d,time,options) %#ok<INUSL>
% ct_OpenLoopSim    Simulates the continuous-time open-loop dynamical system
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


ODEsolver = options.ODEsolver;
ODEoptions = options.ODEoptions; %#ok<NASGU>

% Start and stom time for simulation
tstart = 0; %#ok<NASGU>
tend = time; %#ok<NASGU>

% Control action (equal to 0 since this is a open-loop
% simulation)
if object.nu == 0
    u = [];
else
    u = zeros(1,object.nu);
end

% String representing the call to the ODE solver
command = [ODEsolver,'(@sys,[tstart tend],x0,ODEoptions);'];

% Call ODE solver
[T,X] = eval(command);

% Set vector U
U = zeros(numel(T),object.nu);

% Compute system outputs
Y = object.computeOutput(X,U,p,d);

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







