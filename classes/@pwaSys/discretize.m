function dtSys = discretize(object,Ts)
% discretize  - Converts the system from continuous to discrete-time 
% dtSys = discretize(OBJ,Ts)
% Converts the system from continuous to discrete-time, with zero-order-hold 
% method. Ts is the sampling time of the resulting discrete-time system.
% The conversion is performed by exploiting 'c2d' function of MATLAB
% Control System Toolbox. OBJ is the ltiSys object.

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


% Check if system is already discrete-time
if strcmp(object.mode,'dt')
    error('System is already discrete-time!')
end

% Check sampling time
if Ts <= 0
    error('Sampling time must be > 0')
end

nDyn = object.nDyn;

A = object.A;
B = object.B;
C = object.C;
D = object.D;
Ex = object.Ex;
Ey = object.Ey;
Fx = object.Fx;
Fy = object.Fy;
Gx = object.Gx;
Gy= object.Gy;

H = object.H;
K = object.K;

% Create PWA system
dtSys = pwaSys(object.nx,object.nu,object.ny,object.np,object.nd,'dt',Ts);


for i=1:nDyn
% Extract continuous-time matrices for each dynamic
Ac = A{i};
Bc = B{i};
Cc = C{i};
Dc = D{i};
Exc = Ex{i};
Eyc = Ey{i};
Fxc = Fx{i};
Fyc = Fy{i};
Gxc = Gx{i};
Gyc = Gy{i};

HH = H{i};
KK = K{i};

% Put together matrices B, Ex, Fx and Gx
Bcext = [Bc Exc Fxc Gxc];
Dcext = [Dc Eyc Fyc Gyc];

% Create a continuous-time ss object
ctsys = ss(Ac, Bcext, Cc, Dcext);

% Convert ss object into discrete-time
dtsys = c2d(ctsys,Ts);

% Extract discrete-time matrices from ss object
[Ad, Bdext, Cd, Ddext] = ssdata(dtsys);

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

dtSys = dtSys.addDynamics(HH,KK);

% Create discrete-time ltiSys object
dtSys = dtSys.setMatrices(i,'A',Ad);
dtSys = dtSys.setMatrices(i,'B',Bd);
dtSys = dtSys.setMatrices(i,'C',Cd);
dtSys = dtSys.setMatrices(i,'D',Dd);
dtSys = dtSys.setMatrices(i,'Ex',Exd);
dtSys = dtSys.setMatrices(i,'Ey',Eyd);
dtSys = dtSys.setMatrices(i,'Fx',Fxd);
dtSys = dtSys.setMatrices(i,'Fy',Fyd);
dtSys = dtSys.setMatrices(i,'Gx',Gxd);
dtSys = dtSys.setMatrices(i,'Gy',Gyd);

end
dtSys = dtSys.setStateNames(object.getStateNames());
dtSys = dtSys.setInputNames(object.getInputNames());
dtSys = dtSys.setOutputNames(object.getOutputNames());
dtSys = dtSys.setParameterNames(object.getParameterNames());
dtSys = dtSys.setUnmeasurableInputNames(object.getUnmeasurableInputNames());



