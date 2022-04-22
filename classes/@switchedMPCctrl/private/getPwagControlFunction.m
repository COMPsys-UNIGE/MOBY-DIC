function pwag = getPwagControlFunction(object)

ndyn = object.getNumberOfDynamics;

ctrl = object.controllers(1).MPCcontroller;
[Hd Kd] = ctrl.getDomain;
ny = ctrl.getFunction.getCodomainDimensions;
reg = [];
for i=1:ndyn
    
    H = [object.controllers(i).H, zeros(1,numel(ctrl.getTrackingVariable))];
    K = object.controllers(i).K;
       
    pwagFun = object.controllers(i).MPCcontroller.getFunction;
   
    for j=1:pwagFun.getNumberOfRegions
currentReg = pwagFun.getRegions(j);
        Hp = [Hd;H;currentReg.H];
        Kp = [Kd;K;currentReg.K];
        pp = Polyhedron(Hp,Kp);
        
        if ~pp.isEmptySet
            currentReg.H = pp.A;
            currentReg.K = pp.b;
            reg = [reg; currentReg];
        end
    
    end
    
end

D.Hd = Hd;
D.Kd = Kd;

pwag = pwagFunction(reg,D);

end