function xnext = evaluate(object,x,u,p,d)
% evaluate    Computes the derivative of the state (if the system is 
%             continuous-time) or the state at next time instant (if the 
%             system is discrete-time)
%
% DXDT = evaluate(OBJ,X,U,P,D)
% If the system is continuous-time, the time derivative dxdt at current
% time instant is computed as:
%
%   DXDT = A X + B U + E_x P + F_x D + G_x
%
% where X, U, P and D are the system state, input, parameters and 
% unmeasurable inputs at current time instant, respectively. 
% OBJ is the ltiSys object.
%
% XNEXT = evaluate(OBJ,X,U,P,D)
% If the system is discrete-time, the system state at next time instant 
% is computed as:
%
%   XNEXT = A X + B U + E_x P + F_x D + G_x
%
% where X, U, P and D are the system state, input, parameters and 
% unmeasurable inputs at current time instant, respectively. 
% OBJ is the ltiSys object.
%
% If the system has no parameters and/or unmeasurable inputs, the following
% compact syntaxes can be used:
% DXDT = evaluate(OBJ,X,U,P)
% DXDT = evaluate(OBJ,X,U)
% XNEXT = evaluate(OBJ,X,U,P)
% XNEXT = evaluate(OBJ,X,U)

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

if numel(x) ~= object.nx
    error(['X must be a vector with  ',num2str(object.nx),' elements'])
end
if numel(u) ~= object.nu
    error(['U must be a vector with  ',num2str(object.nu),' elements'])
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

if numel(p) ~= object.np
    error(['P must be a vector with  ',num2str(object.np),' elements'])
end
if numel(d) ~= object.nd
    error(['D must be a vector with  ',num2str(object.nd),' elements'])
end

xnext = object.A*x(:)+object.B*u(:)+object.Ex*p(:)+object.Fx*d(:)+object.Gx;