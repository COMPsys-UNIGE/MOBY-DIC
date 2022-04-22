function dtSys = discretize(object,Ts)
% deltau   Converts the system from continuous to discrete-time 
%
% dtSys = discretize(OBJ,Ts)
% Converts the system from continuous to discrete-time, with zero-order-hold 
% method. Ts is the sampling time of the resulting discrete-time system.
% The conversion is performed by exploiting 'c2d' function of MATLAB
% Control System Toolbox. OBJ is the ltiSys object.

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


% Check if system is already discrete-time
if strcmp(object.mode,'dt')
    error('System is already discrete-time!')
end

% Check sampling time
if Ts <= 0
    error('Sampling time must be > 0')
end

% Extract continuous-time matrices
Ac = object.A;
Bc = object.B;
Cc = object.C;
Dc = object.D;
Exc = object.Ex;
Eyc = object.Ey;
Fxc = object.Fx;
Fyc = object.Fy;
Gxc = object.Gx;
Gyc = object.Gy;

% Put together matrices B, Ex, Fx and Gx
Bcext = [Bc Exc Fxc Gxc];
Dcext = [Dc Eyc Fyc Gyc];

% Create a continuous-time ss object
ctsys = ss(Ac, Bcext, Cc, Dcext);

% Convert ss object into discrete-time
dtsys = c2d(ctsys,Ts);

% Extract discrete-time matrices from ss object
[Ad, Bdext, Cd, Ddext] = ssdata(dtsys);
% Ad=eye(size(Ac,1))+Ts*Ac;
% Bdext=Ts*Bcext;
% Cd=Cc;
% Ddext=Dcext;

% Separate matrices B, W and f
Bd = Bdext(:,1:object.nu);
Exd = Bdext(:,object.nu+1:object.nu+object.np);
Fxd = Bdext(:,object.nu+object.np+1:object.nu+object.np+object.nd);
Gxd = Bdext(:,end);
% Take only first nu columns from D
Dd = Ddext(:,1:object.nu);
Eyd = Ddext(:,object.nu+1:object.nu+object.np);
Fyd = Ddext(:,object.nu+object.np+1:object.nu+object.np+object.nd);
Gyd = Ddext(:,end);

% Create discrete-time ltiSys object
dtSys = ltiSys(object.nx,object.nu,object.ny,object.np,object.nd,'dt',Ts);
dtSys = dtSys.setMatrices('A',Ad);
dtSys = dtSys.setMatrices('B',Bd);
dtSys = dtSys.setMatrices('C',Cd);
dtSys = dtSys.setMatrices('D',Dd);
dtSys = dtSys.setMatrices('Ex',Exd);
dtSys = dtSys.setMatrices('Ey',Eyd);
dtSys = dtSys.setMatrices('Fx',Fxd);
dtSys = dtSys.setMatrices('Fy',Fyd);
dtSys = dtSys.setMatrices('Gx',Gxd);
dtSys = dtSys.setMatrices('Gy',Gyd);

dtSys = dtSys.setStateNames(object.getStateNames());
dtSys = dtSys.setInputNames(object.getInputNames());
dtSys = dtSys.setOutputNames(object.getOutputNames());
dtSys = dtSys.setParameterNames(object.getParameterNames());
dtSys = dtSys.setUnmeasurableInputNames(object.getUnmeasurableInputNames());

% TO DO
% controllare se il controllore va bene
if ~isempty(object.controller)
    dtSys = dtSys.setController(object.controller);
end
if ~isempty(object.observer)
    dtSys = dtSys.setObserver(object.observer);
end




