%% findDynamics Find the index of the dynamic related to state value x
%
% Y = findDynamics(OBJ,X)
% where X is the system state and OBJ is the discontinue pwas object.
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


function dyn = findDynamics(object,x)

if size(x,1) ~= object.nx
    if size(x,2) == object.nx
        x = x';
    else
        error(['X must be a matrix with  ',num2str(object.nx),' columns'])
    end
end

% Number of points
npts = size(x,2);
dyn = zeros(npts,1);

for n =1:npts
    for i=1:object.nDyn
        if object.functions(i).H*[x(:,n)] <= object.functions(i).K
            break
        end
    end
    dyn(n) = i;
end

end
