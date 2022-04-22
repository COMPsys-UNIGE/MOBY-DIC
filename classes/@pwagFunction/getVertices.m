function V = getVertices(object,varargin)
% getVertices  Gets the vertices of the polytopes of the domain partition
%
% V = getVertices(OBJ)
% Gets all the vertices of the domain partition. V is a cell array with as 
% many elements as the number of polytopes, whose i-th element is a 
% N_VERTICES x NX matrix containing the N_VERTICES vertices of the i-th 
% region. OBJ is the pwagFunction object.
%
% V = getVertices(OBJ,REG)
% Gets the vertices of the regions whose indices are specified in array REG.
% V is a cell array with as many elements as the number of indexes in REG.
%
% Vertices are extracted through MPT3.

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



% Check on input number
if nargin < 1 || nargin > 2
    error('Wrong input number')
end

% If reg is not provided, all regions are considered
if nargin == 1
    reg = 1:object.nr;
else
    reg = varargin{1};
end

% Consider array reg as column vector
reg = reg(:);

% Number of regions where to compute vertices
nr = numel(reg);

% initialize cell array containing vertices
V = cell(nr,1);

% Loop on regions
for i = 1:nr
    
    [H, K] = object.getEdges(reg(i));
    
    % Find vertices through MPT3
    P = Polyhedron(H,K);
    V{i} = P.V;
end