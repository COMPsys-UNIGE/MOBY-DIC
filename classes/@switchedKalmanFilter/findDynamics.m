function dyn = findDynamics(object,x,p)
% findDynamics   Find the index of the dynamicd related to state value x, 
% parameters value p and unmeasurable inputs value d
%
% Y = computeOutput(OBJ,X,P,D) 
% where X, U, P and D are the system state, input, parameters and 
% unmeasurable inputs, respectively. OBJ is the ltiSys object.
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

if size(x,1) ~= object.nx+object.nd
    if size(x,2) == object.nx+object.nd
        x = x';
    else
    error(['X must be a matrix with  ',num2str(object.nx+object.nd),' columns'])
    end
end

% Number of points
npts = size(x,2);

if ~exist('p','var') && object.np > 0
    error('Vector of parameters P must be provided')
end



if ~exist('p','var')
    p = [];
end


if size(p,1) ~= object.np
    if size(p,2) == object.np
        p = p';
    else
    error(['P must be a matrix with  ',num2str(object.np),' columns'])
    end
end


if numel(p) == 0
    p = zeros(0,npts);
end

if size(p,2) == 1
    p = repmat(p,size(p,1),npts);
end


dyn = zeros(npts,1);

for n =1:npts
    for i=1:object.nDyn
        if object.filters(i).H*[x(1:object.nx,n); p(:,n); x(object.nx+1:end,n)] <= object.filters(i).K
            break
        end
    end
    
    dyn(n) = i;
end

end