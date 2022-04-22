function disp(object)
% disp   Displays some information about the pwagFunction object
%
% disp(OBJ)
% OBJ is the pwagFunction object.

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


if object.nx == 0
    disp('Empty pwagFunction object');
    disp(' ')
else
    disp('pwagFunction object');
    disp(' ')
    disp([' - Domain dimensions: ',num2str(object.nx)]);
    disp([' - Codomain dimensions: ',num2str(object.ny)]);
    disp(' ')
    disp([' - Number of regions: ',num2str(object.nr)]);
    disp([' - Number of edges: ',num2str(size(object.getEdges,1))]);
    disp([' - Number of linear functions: ',num2str(size(object.getFunctions,1))]);
    xstr = [];
    for i = 1:object.nx-1
        xstr = [xstr object.xnames{i} ', ']; %#ok<*AGROW>
    end
    xstr = [xstr object.xnames{end}];
    ystr = [];
    for i = 1:object.ny-1
        ystr = [ystr object.ynames{i} ', '];
    end
    ystr = [ystr object.ynames{end}];
    disp(' ')
    disp([' - Input names: ',xstr])
    disp([' - Output names: ',ystr])
    disp(' ')
    
    if ~isempty(object.tree)
        disp(' - Tree computed:');
        disp(['     number of nodes: ',num2str(object.tree.numberOfNodes)]);
        disp(['     minimum depth: ',num2str(object.tree.minDepth)]);
        disp(['     mean depth: ',num2str(object.tree.meanDepth)]);
        disp(['     maximum depth: ',num2str(object.tree.maxDepth)]);
        disp(['     balance: ',num2str(object.tree.balance)]);
    end
end
end