function range = getRangeFromController(object)
%GETRANGEFROMCONTROLLER
% Get range of signals involved in exact or approximate explicit MPC controller
% 
% This is a private method.

ctrl = object.getController;

np = object.np;
nd = object.nd;
nxref = numel(object.controller.getTrackingVariable);

% get controller constraints and fill range struct
% state constraints
nfun = object.getController.getFunction.getCodomainDimensions;
fun = object.getController.getFunction;
nDim = fun.getDomainDimensions;
nx = object.nx;
nu = object.getController.getNumberOfInputs;

idx = 1:nu;

[Hd Kd] = fun.getDomain;
P = Polyhedron(Hd,Kd);

if ~P.isBounded
    error(['PWAG function domain must be bounded in order'...
            ' to generate circuit file']);
end
Dp = P.outerApprox();
DVert = Dp.V;

minV = min(DVert);
maxV = max(DVert);

range.xmin = minV(1:nx);
range.xmax = maxV(1:nx);

range.pmin = minV(nx+1:nx+np);
range.pmax = maxV(nx+1:nx+np);

range.dmin = minV(nx+np+1:nd+nx+np);
range.dmax = maxV(nx+np+1:nd+nx+np);

range.xrefmin = minV(nx+np+nd+1:nd+nx+np+nxref);
range.xrefmax = maxV(nx+np+nd+1:nd+nx+np+nxref);


% get codomain limits

VV = fun.getVertices;
V = [];
if isa(fun,'pwagFunction')
    for i=1:numel(VV)
        V = [V; VV{i}];
    end
elseif isa(fun,'pwasFunction')
    V = VV;
elseif isa(fun,'discontinuousPwasFunction')
    V = VV;
else
    error(['Cannot find codomain of ',class(fun),' function']);
end

evalV = fun.eval(V');

umin = min(evalV,[],2);
umax = max(evalV,[],2);

range.umin = umin(idx);
range.umax = umax(idx);
    
end