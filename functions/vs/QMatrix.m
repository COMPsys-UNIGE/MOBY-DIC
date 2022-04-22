%% Q MATRIX
% Computes the Q matrix for order 1 Tikhonov regularization
%
% SYNTAX
%
% Q = QMatrix(D,P)
%
% D is a matrix specifying the domain in the form:
%  _                                 _
% |  x_min(1) x_min(2) ... x_min(nx)  |
% |_ x_max(1) x_max(2) ... x_max(nx) _|
%
% P can be an array containing the number of subdivisions per dimensions
% (in case of uniform partition) or a cell array whose i-th element contains
% the i-th component of the vertices of the simplicial partition
% (for non-uniform partition).
%
% ACKNOWLEDGEMENTS
%
% Contributors:
%
% * Tomaso Poggi (tpoggi@essbilbao.org)
%
% Copyright is with:
%
% * Copyright (C) 2010 University of Genoa, Italy.


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

function Q = QMatrix(D,P)

% Number of dimensions
n = numel(P);

% Number of subdivisions per dimension
np = zeros(1,n);

if iscell(P)
    % Number of vertices of the simplicial partition
    Nbs = 1;
    for i = 1:numel(P)
        Nbs = Nbs * numel(P{i});
        np(i) = numel(P{i})-1;
    end
else
    % Number of vertices of the simplicial partition
    Nbs = prod(P+1);
    np = P;
    % Convert P into a cell array defining the vertices of the partition
    P = cell(n);
    for i = 1:n
        P{i} = linspace(D(1,i),D(2,i),np(i)+1);
    end
end

% Generate vertices of a fictitious simplicial partition in a normalized
% domain whose i-th dimension ranges from 0 to P(i)
v = cell(n,1);
for i = 1:n
    v{i} = 0:np(i);
end
if n > 1
    [v{1:n}] = ndgrid(v{:});
end
Jv = zeros(Nbs,n);
for i = 1:n
    Jv(:,i) = v{i}(:);
end

d = cell(n,1);
for i = 1:n
    d{i} = [diff(P{i}) NaN];
end
if n > 1
    [d{1:n}] = ndgrid(d{:});
end
Jd = zeros(Nbs,n);
for i = 1:n
    Jd(:,i) = d{i}(:);
end

% Initialize Q matrix
Q = spalloc(Nbs,Nbs,Nbs*(n+1));

for k = 1:Nbs
    
    % It is the reciprocal of the square of the distance between the 2
    % adjacent vertices
    entry = 1./Jd(k,:).^2;
    
    % Loop on the number of dimensions
    for i = 1:n
        % Initialize vector tmp: [0 0 ... 0]
        tmp = zeros(1,n);
        
        % Put a "1" in i-th position: [0 0 ... 0 1 0 ... 0]
        tmp(i) = 1;
        
        % Add k-th vertex to tmp
        % What you obtain at the end is that tmp is the vertex adjacent to
        % the current vertex Jv, along the i-th dimension.
        tmp = Jv(k,:)+tmp;
        
        % This part of code is performed only if the adjacent vertex is
        % inside the domain
        if tmp(i) <= np(i)
            
            % Find the index of the adjacent vertex
            dij=findIndex(tmp(:),np); % it works because tmp is 1xn
            
            % Fill Q matrix
            Q(k,k)=Q(k,k)+entry(i);
            Q(dij,dij)=Q(dij,dij)+entry(i);
            Q(k,dij)=Q(k,dij)-entry(i);
            Q(dij,k)=Q(dij,k)-entry(i);
        end
    end
end

    function addr = findIndex(x,np)
        coeff = 1;
        addr = 0;
        for kk = 1:numel(np)
            addr = addr+coeff*(x(kk,:));
            coeff = coeff*(np(kk)+1);
        end
        addr = addr+1;
    end
    
end