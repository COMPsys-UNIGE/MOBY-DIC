function duSys = deltau(object)
% deltau   Converts the system to delta-u formulation
%
% duSys = deltau(OBJ)
% Converts the system to delta-u formulation. The original input becomes a
% system state and the input is replaced by its derivative (for
% continuous-time systems) or increment (for discrete-time systems).

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



% Extract matrices
A = object.A;
B = object.B;
C = object.C;
D = object.D;
Ex = object.Ex;
Ey = object.Ey;
Fx = object.Fx;
Fy = object.Fy;
Gx = object.Gx;
Gy = object.Gy;

H = object.H;
K = object.K;

% Extract dimensions
nx = object.nx;
nu = object.nu;
ny = object.ny;
np = object.np;
nd = object.nd;

ndu = nu;

nDyn = object.nDyn;

% Extract names
xnames = object.xnames;
unames = object.unames;
ynames = object.ynames;
pnames = object.pnames;
dnames = object.dnames;

nxnew = nx+nu;
nunew = ndu;
npnew = np;
ndnew = nd;
nynew = ny;

% Create PWA system object
if object.isContinuousTime
    duSys = pwaSys(nxnew,nunew,nynew,npnew,ndnew,'ct');
else
    duSys = pwaSys(nxnew,nunew,nynew,npnew,ndnew,'dt',object.Ts);
end


for j=1:nDyn
    
    if object.isContinuousTime
        Anew = [A{j} B{j}; zeros(nu,nx+nu)];
    else
        Anew = [A{j} B{j}; zeros(nu,nx) eye(nu)];
    end
    
    Bnew = [zeros(nx,ndu); eye(nu,ndu)];
    Cnew = [C{j} D{j}];
    Dnew = zeros(ny,ndu);
    Exnew = [Ex{j}; zeros(nu,np)];
    Eynew = Ey{j};
    Fxnew = [Fx{j}; zeros(nu,nd)];
    Fynew = Fy{j};
    Gxnew = [Gx{j}; zeros(nu,1)];
    Gynew = Gy{j};
    
    Hnew = [H{j} zeros(size(H{j},1),nu)];
    Knew = K{j};
    
    xnamesnew = [xnames unames];
    unamesnew = unames;
    for i = 1:ndu
        unamesnew{i} = ['delta-',unames{i}];
    end
    ynamesnew = ynames;
    pnamesnew = pnames;
    dnamesnew = dnames;
    
    duSys = duSys.addDynamics(Hnew,Knew);
    
    duSys = duSys.setMatrices(j,'A',Anew);
    duSys = duSys.setMatrices(j,'B',Bnew);
    duSys = duSys.setMatrices(j,'C',Cnew);
    duSys = duSys.setMatrices(j,'D',Dnew);
    duSys = duSys.setMatrices(j,'Ex',Exnew);
    duSys = duSys.setMatrices(j,'Ey',Eynew);
    duSys = duSys.setMatrices(j,'Fx',Fxnew);
    duSys = duSys.setMatrices(j,'Fy',Fynew);
    duSys = duSys.setMatrices(j,'Gx',Gxnew);
    duSys = duSys.setMatrices(j,'Gy',Gynew);
    
end


duSys = duSys.setStateNames(xnamesnew);
duSys = duSys.setInputNames(unamesnew);
duSys = duSys.setOutputNames(ynamesnew);
duSys = duSys.setParameterNames(pnamesnew);
duSys = duSys.setUnmeasurableInputNames(dnamesnew);



