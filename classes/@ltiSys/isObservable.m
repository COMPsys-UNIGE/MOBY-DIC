function answ = isObservable(object)
% isObservable    Checks if the system is observable
%
% The LTI system is in this form (consider only continuous-time case,
% analog considerations hold also for discrete-time systems):
%   _
%  |
%  |  dx
%  | ---- = A x(t) + B u(t) + E_x p(t) + F_x d(t) + G_x
% <   dt
%  |
%  | y(t) = C x(t) + D u(t) + E_y p(t) + F_y d(t) + G_y
%  |_
%
% The input d(t) is not measurable, therefore it must be estimated
% and is therefore considered as a system state. Conversely, p(t) is a
% measurable parameter and is then considered as an input.
%
% Call    _    _              _    _
%        | x(t) |            | u(t) |
% z(t) = |      | and v(t) = |      |
%        |_d(t)_|            |_p(t)_|
%
% the system can be reformulated as follows:
%   _
%  |
%  |  dz    
%  | ---- = A' z(t) + B' v(t) + G_x'
% <   dt    
%  |
%  | y(t) = C' z(t) + D' v(t) + G_y'
%  |_
%
% where,
%        _      _         _      _          _   _
%       | A  F_x |       | B  E_x |        | G_x |
%  A' = |        |, B' = |        |, G_x'= |     |
%       |_0   0 _|       |_0   0 _|        |_ 0 _|
%
%  C' = [ C  F_y ], D' = [ D  E_y ], G_y' = G_y
%
% The system is observable if the following matrix has full rank (n is the 
% dimension of matrix A'):
%  _             _ 
% | C'            |
% | C' A'         |
% | C' (A')^2     |
% |    ...        |
% |_C' (A')^(n-1)_|
%
% The syntax of the method is the following:
%
% ANSW = isObservable(OBJ)
% OBJ is the ltiSys object.

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

A = object.A;
Fx = object.Fx;
C = object.C;
Fy = object.Fy;

nx = object.nx;
nd = object.nd;

An = [A Fx;zeros(nd,nx) zeros(nd,nd)];
Cn = [C Fy];

O = obsv(An,Cn);

if rank(O) == nx
    answ = true;
else
    answ = false;
end
