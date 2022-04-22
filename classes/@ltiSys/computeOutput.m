function y = computeOutput(object,x,u,p,d)
% computeOutput   Computes the output of the system
%
% Y = computeOutput(OBJ,X,U,P,D)
% Returns the system output defined as:
%
%   Y = C X + D U + E_y P + F_y D + G_y
% 
% where X, U, P and D are the system state, input, parameters and 
% unmeasurable inputs, respectively. OBJ is the ltiSys object.
%
% Y = computeOutput(OBJ,X,U,P)
% Computes the system output when the system has no unmeasurable inputs.
%
% Y = computeOutput(OBJ,X,U)
% Computes the system output when the system has no parameters and 
% unmeasurable inputs.

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

if size(x,2) ~= object.nx
    error(['X must be a matrix with  ',num2str(object.nx),' columns'])
end
if size(u,2) ~= object.nu
    error(['U must be a matrix with  ',num2str(object.nu),' columns'])
end

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

if size(p,2) ~= object.np
    error(['P must be a matrix with  ',num2str(object.np),' columns'])
end
if size(d,2) ~= object.nd
    error(['D must be a matrix with  ',num2str(object.nd),' columns'])
end

% Number of points
npts = size(x,1);

if size(u,1) == 1
    u = repmat(u,npts,1);
end
if size(p,1) == 1
    p = repmat(p,npts,1);
end
if size(d,1) == 1
    d = repmat(d,npts,1);
end
if size(u,1) ~= npts
    error(['U must be a matrix with either 1 or ',num2str(npts),' rows'])
end
if size(p,1) ~= npts
    error(['P must be a matrix with either 1 or ',num2str(npts),' rows'])
end
if size(d,1) ~= npts
    error(['D must be a matrix with either 1 or ',num2str(npts),' rows'])
end

y = x*object.C'+u*object.D'+p*object.Ey'+d*object.Fy'+repmat(object.Gy',npts,1);