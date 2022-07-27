function dyn = findDynamics(object,x,p,d)
% findDynamics   Finds the index of the system dynamics
%
% DYN = findDynamics(OBJ,X,P,D) 
% Finds the index of the dynamics related to state X, 
% parameters P and unmeasurable inputs D. Based on these values,
% all inequalities H [x p d]' <= K are checked until one of them is
% satisfied. The index DYN of such inequality is returned. If no inequality
% is satisfied, DYN = 0. OBJ is the pwaSys object.
%
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

if size(x,1) ~= object.nx
    if size(x,2) == object.nx
        x = x';
    else
    error(['X must be a matrix with  ',num2str(object.nx),' columns'])
    end
end

% Number of points
npts = size(x,2);

if ~exist('p','var') && object.np > 0
    error('Vector of parameters P must be provided')
end

if ~exist('d','var') && object.nd > 0
    error('Vector of unmeasurable inputs D must be provided')
end

if ~exist('p','var')
    p = [];
end

if ~exist('d','var')
    d = [];
end

if size(p,1) ~= object.np
    if size(p,2) == object.np
        p = p';
    else
    error(['P must be a matrix with  ',num2str(object.np),' columns'])
    end
end
if size(d,1) ~= object.nd
    if size(d,2) == object.nd
        d = d';
    else
    error(['D must be a matrix with  ',num2str(object.nd),' columns'])
    end
end

if numel(d) == 0
    d = zeros(0,npts);
end
if numel(p) == 0
    p = zeros(0,npts);
end

if size(p,2) == 1
    p = repmat(p,size(p,2),npts);
end

if size(d,2) == 1
    d = repmat(d,size(d,2),npts);
end


dyn = zeros(npts,1);

H = object.H;
K = object.K;

for n =1:npts
    for i=1:object.nDyn
        if H{i}*[x(:,n); p(:,n); d(:,n)] <= K{i}
            break
        end
    end
    
    dyn(n) = i;
end

end