function tree = getTree(object)
% getTree   Gets the binary search tree associated to the domain partition
%
% TREE = getTree(OBJ)
% OBJ is the pwagFunction object. TREE is a structure with the following
% fields:
%  - numberOfNodes: number of nodes of the tree
%  - minDepth: minimum depth of the tree
%  - meanDepth: mean depth of the tree
%  - maxDepth: maximum depth of the tree
%  - nodes: array of structures corresponding to all nodes of the tree
%  - balance: estimation of the tree balance (computed as the variance of
%             the node depths)
%
% See also: pwagFunction/computeTree

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

tree = object.tree;