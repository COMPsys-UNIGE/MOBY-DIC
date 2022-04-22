%% eval   Evaluates the discontinuous PWAS function

% Y = eval(OBJ,X)
% Evaluates all NY components of the discontinuous PWAS function at the points
% specified by matrix X. X must be a NPOINTS x NX matrix. Y is a
% NPOINTS x NY matrix. NX and NY are the domain and codomain dimensions,
% respectively. OBJ is the discontinuous pwasFunction object.
%
% Y = eval(OBJ,X,IDX)
% Evaluates only the function components defined by IDX. IDX must be a vector
% of indices in the range 1, NY.
%
% Mex files are used to speedup the computation.

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


function u = eval(object,x)

if size(x,1) ~= object.getDomainDimensions
    if size(x,2) == object.getDomainDimensions
        x = x';
    else
    error(['X must be a matrix with  ',num2str(object.getDomainDimensions),' rows'])
    end
end

% Number of points
npts = size(x,2);


dyn = object.findDynamics(x);
for i=1:npts
u = object.functions(dyn(i)).pwasFunction.eval(x(:,i));
end
end
