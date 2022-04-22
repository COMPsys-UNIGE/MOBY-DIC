function info_out = treeStepPlot(nodes,node,info_in,ascissas,depthmax)
% treeStepPlot   Explores the binary search tree in order to retrieve information for the
% plotting of the tree

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


if node.leaf == 0
    name = info_in.prevName;
    followNames = [2*name+1 2*name+2];
    info_in.names(info_in.index+1) = name;
    info_in.points(name+1,:) = [ascissas{node.depth+1}(name-2^node.depth+2) 1-node.depth/depthmax];
    info_in.index = info_in.index+1;
    info_in.leaves(name+1) = node.leaf;
    for i = 1:2
        info_in.prevName = followNames(i);
        % if node.children(i) ~= 0
        info_in = treeStepPlot(nodes,nodes(node.children(i)),info_in,ascissas,depthmax);
	% end
    end
else
    name = info_in.prevName;
    info_in.names(info_in.index+1) = name;
    info_in.points(name+1,:) = [ascissas{node.depth+1}(name-2^node.depth+2) 1-node.depth/depthmax];
    info_in.index = info_in.index+1;
    info_in.leaves(name+1) = node.leaf;
end
info_out = info_in;