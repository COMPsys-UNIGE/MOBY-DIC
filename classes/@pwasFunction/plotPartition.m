function plotPartition(object)
% plotPartition   Plots the simplicial domain partition
%
% plotPartition(OBJ)
% Plots the simplicial domain partition of the PWAS function OBJ.
%
% This functions exploits MATLAB functions delaunay and triplot.

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


if object.nx ~= 2
    error('Only two-dimensional partitions can be plotted')
end

% Extract vertices
V = object.getVertices();

% Rotate vertices
Vr = [V(:,1) -V(:,2)];

% Use MATLAB function delaunay
tri = delaunay(Vr);

% Plot partition
figure
triplot(tri,V(:,1),V(:,2),'color',[0 0 0],'linewidth',2)

xlabel(object.xnames{1});
ylabel(object.xnames{2});