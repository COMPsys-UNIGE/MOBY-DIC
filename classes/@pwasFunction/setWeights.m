function object = setWeights(object,w)
% setWeights  Sets the weights w of the PWAS function
%
% OBJ = setWeights(OBJ,W)
% Sets the weights W defining the PWAS function in the form:
%              _nv_
%              \
%       u(x) = /___  w_j alpha_j(x)
%              j = 1
%
% W is a NV x NY matrix, being NV the number of vertices and NY the
% number of codomain dimensions. OBJ is the pwasFunction object.
% The ordering of the weights is as follows:
%
% w_7 ____w_8____ w_9
%    |       |   |
%    |       |   |
%    |       |   |
% w_4|____w_5|___|w_6
%    |       |   |
%    |_______|___|
% w_1     w_2     w_3

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

if size(w,1) ~= object.nv
    error(['Matrix w must have ',num2str(object.nv),' rows'])
end
object.w = w;

ny = size(w,2);

object.ny = ny;

ynames = cell(ny,1);
for i = 1:ny
    ynames{i} = ['y_',num2str(i)];
end
object.ynames = ynames;

end