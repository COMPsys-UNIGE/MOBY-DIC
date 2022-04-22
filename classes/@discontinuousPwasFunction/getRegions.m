function regions = getRegions(object)
%GETREGIONS Summary of this function goes here
%   Detailed explanation goes here

% TO DO


S = object.getSimplices;
ind = 1;
for j=1:object.nDyn
    Hd = [object.domain.Hd; object.functions(j).H];
    Kd = [object.domain.Kd; object.functions(j).K];
    
    for i=1:numel(S)
        Ptmp = Polyhedron(S{i});
        P = Polyhedron([Hd;Ptmp.A],[Kd;Ptmp.b]);
        if ~P.isEmptySet 
        y = object.functions(j).pwasFunction.eval(S{i});
        x = [S{i}, ones(object.getDomainDimensions+1,1)];
        FG = linsolve(x,y);
        F = FG(1:end-1)';
        G = FG(end);
        regions(ind).F = F;
        regions(ind).G = G;
        regions(ind).H = P.A;
        regions(ind).K = P.b;
        ind = ind+1;
        end
    end
    
end

end

