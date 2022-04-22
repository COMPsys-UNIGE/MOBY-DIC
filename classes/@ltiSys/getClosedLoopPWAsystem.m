function CL_PWAsystem = getClosedLoopPWAsystem(object)
% getClosedLoopPWAsystem    Gets the closed-loop PWA system
%
% CL_PWASYS = getClosedLoopPWAsystem(OBJ)
% Gets a pwaSystem object defining the closed loop system comprising the 
% LTI system (plant) and the explicit PWA MPC controller. 

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

controller = object.getController;

if object.isContinuousTime
    Ts = object.getController.getSamplingTime;
        warning(['Closed loop need discrete time system! The system will be discretize with sample time ',num2str(Ts),'!']);
    dtSys = object.discretize(Ts);
else
    Ts = object.getSamplingTime;
    dtSys = object;
end

if isempty(controller)
    error 'No close loop controllere set';
end;

pwagFun = controller.getFunction;

nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nu = object.getNumberOfInputs;
ny = object.getNumberOfOutputs;
nxref = numel(controller.getTrackingVariable);
if isa(pwagFun,'pwagFunction')
nRegion = pwagFun.getNumberOfRegions;
regions = pwagFun.getRegions;
elseif isa(pwagFun,'pwasFunction')
    nRegion = pwagFun.getNumberOfSimplices;
    S = pwagFun.getSimplices;
    for i=1:nRegion
        P = Polyhedron(S{i});
        y = pwagFun.eval(S{i});
        x = [S{i}, ones(nx+np+nd+nxref+1,1)];
        FG = linsolve(x,y);
        F = FG(1:end-1)';
        G = FG(end);
        regions(i).F = F;
        regions(i).G = G;
        regions(i).H = P.A;
        regions(i).K = P.b;
    end
else
    error('Unknown control function');
end

CL_PWAsystem = pwaSys(nx,nu,ny,np,nd,'dt',Ts);

[A, B, C, D, Ex, Ey, Fx, Fy, Gx, Gy] = dtSys.getMatrices();

for i=1:nRegion
    reg = regions(i);
    
    FFx = reg.F(1:nu,1:nx);
    if np > 0
        FFp = reg.F(1:nu,nx+1:nx+np);
    else
        FFp = zeros(nu,0);
    end
    if nd > 0
        FFd = reg.F(1:nd,nx+np+1:end);
    else
        FFd = zeros(nu,0);
    end
    
    A_cl = A+B*FFx;
    B_cl = B;
    Ex_cl = Ex+B*FFp;
    Fx_cl = Fx+B*FFd;
    Gx_cl = Gx+B*reg.G(1:nu);
    
    C_cl = C+D*FFx;
    D_cl = D;
    Ey_cl = Ey+D*FFp;
    Fy_cl = Fy+D*FFd;
    Gy_cl = Gy+D*reg.G(1:nu);
    
    
    Hx_cl = reg.H(:,1:nx);
    Hp_cl = reg.H(:,nx+1:nx+np);
    Hd_cl = reg.H(:,nx+np+1:end);
    K_cl = reg.K;
    
    CL_PWAsystem = CL_PWAsystem.addDynamics([Hx_cl Hp_cl Hd_cl],K_cl);
    % Set system matrices
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'A',A_cl);
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'B',B_cl);
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'C',C_cl);
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'D',D_cl);
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'Ex',Ex_cl);
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'Fx',Fx_cl);
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'Gx',Gx_cl);
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'Ey',Ey_cl);
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'Fy',Fy_cl);
    CL_PWAsystem = CL_PWAsystem.setMatrices(i,'Gy',Gy_cl);
    
    
end

end



