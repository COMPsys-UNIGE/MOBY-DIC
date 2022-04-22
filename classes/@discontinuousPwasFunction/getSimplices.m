function S = getSimplices(object)

% GET SIMPLICES
% Gets the vertices of all simplices of the simplicial partition
%
% S = getSimplices(OBJ)
% S is a cell-array with as many elements as the number of simplices. The
% i-th element of S is a matrix with dimensions NX+1 X NX, whose rows
% correspond to the simplex vertices. OBJ is the pwasFunction object.
%
% ACKNOWLEDGEMENTS
%
% Contributors:
% 
% * Alberto Oliveri (alberto.oliveri@unige.it)
% 
% Copyright is with:
% 
% * Copyright (C) 2011 University of Genoa, Italy.

% -------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% -------------------------------------------------------------------------

% Number of dimensions
ndim = object.getDomainDimensions();

% Get partition
P = object.getPartition();

% Number of subdivisions per dimension
Nx = object.getNumberOfPartitions();

% Find distance between points subdividing each dimension
dPt = cell(1,ndim);
for i = 1:ndim
    dPt{i} = [diff(P{i})];
end

% Number of hyper rectangles of the simplicial partition
nr = prod(Nx);

% Convert cell array into a matrix. Each row of the matrix contains the
% distances between the current vertex and the vertices close to it
[dPt{1:ndim}] = ndgrid(dPt{:});
dP = zeros(nr,ndim);
for i = 1:ndim
    dP(:,i) = dPt{i}(:);
end

% Remove last vertex from each component of P
Pt = P;
for i = 1:ndim
    Pt{i} = P{i}(1:end-1);
end

% All vertices of simplicial partition excluded the ones at the boundary of
% the domain
[Pt{1:ndim}] = ndgrid(Pt{:});
Vs = zeros(nr,ndim);
for i = 1:ndim
    Vs(:,i) = Pt{i}(:);
end

% All permutations of vector 1:ndim. They are used to locate all simplices
% within a hyper rectangle
p = perms(1:ndim);

% Number of simplices in a hyper-rectangle
np = size(p,1);

% Total number of simplices
NS = prod(Nx)*factorial(ndim);

% Init cell array containing all simplices
S = cell(1,NS);

Vtemp = zeros(ndim+1,ndim);
t = 1;
% Loop on all hyper rectangles
for i = 1:nr
    % Loop on all simplices in a given hyper rectangle
    for j = 1:np
        % Loop on all simplex vertices
        Vtemp(1,:) = Vs(i,:);
        for k = 1:ndim
            Vtemp(k+1,:) = Vtemp(k,:);
            Vtemp(k+1,p(j,k)) = Vtemp(k,p(j,k))+dP(i,p(j,k));
        end
        S{t} = Vtemp;
        t = t+1;
    end 
end
