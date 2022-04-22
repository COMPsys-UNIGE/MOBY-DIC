function object = addDynamics(object,H,K)
% addDynamics   Adds a dynamics to the pwaSys object
%
% object = addDynamics(object,H,K)
% Adds to the system the dynamics defined by matrices H and K, as follows:
%    _   _
%   |  x  |
% H |  p  | <= K
%   |_ d _| 
% 
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

object.nDyn = object.nDyn+1;

object.H{end+1} = H;
object.K{end+1} = K;

object.A{end+1} = zeros(object.nx);
object.C{end+1} = zeros(object.ny,object.nx);
object.B{end+1} = zeros(object.nx,object.nu);
object.D{end+1} = zeros(object.ny,1);
object.Ex{end+1} = zeros(object.nx,object.np);
object.Ey{end+1} = zeros(object.ny,object.np);
object.Fx{end+1} = zeros(object.nx,object.nd);
object.Fy{end+1} = zeros(object.ny,object.nd);
object.Gx{end+1} = zeros(object.nx,1);
object.Gy{end+1} = zeros(object.ny,1);

end