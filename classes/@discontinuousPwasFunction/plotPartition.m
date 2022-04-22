function plotPartition(object)
% plotPartition   Plots the polytopic domain partition
%
% plotPartition(OBJ)
% Plots the polytopic domain partition of the discontinue PWAS function OBJ.
%
% This function exploits MPT3 to perform the plot.

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

% TO DO

func = object.functions;

figure
[Hd , Kd] = object.getDomain;



P = [];
for i=1:object.nDyn
    fi = func(i);
    
    Pd = Polyhedron([fi.H;Hd],[fi.K;Kd]);
    
    S = fi.pwasFunction.getSimplices;
    for j=1:numel(S)
    Ps = Polyhedron(S{j});
    Pint = Ps.intersect(Pd);
    P  = [P;Pint];
    end
end
P.plot();

if P(1).Dim == 2
    
    xlabel(object.xnames{1})
    ylabel(object.xnames{2})
    
elseif P(1).Dim == 3
    
    xlabel(object.xnames{1})
    ylabel(object.xnames{2})
    zlabel(object.xnames{3})
    
end
% regions = object.getRegions();
% 
% for i = 1:object.nr
%     H = regions(i).H;
%     K = regions(i).K;
%     P(i) = Polyhedron(H,K);
% end
% 
% figure
% P.plot();
% 
% if P(1).Dim == 2
%     
%     xlabel(object.xnames{1})
%     ylabel(object.xnames{2})
%     
% elseif P(1).Dim == 3
%     
%     xlabel(object.xnames{1})
%     ylabel(object.xnames{2})
%     zlabel(object.xnames{3})
%     
% end

