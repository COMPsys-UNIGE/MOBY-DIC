function plot(object,k)
% disp   Plots the MPC control function
%
% plot(OBJ)
% Plots the switched MPC control function at current time step. OBJ is the 
% switched MPCctrl object.
%
% plot(OBJ,IDX)
% Plots the switched MPC control function at time step k.
%
% This function exploits MPT3 to perform the plots.

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

if ~exist('k','var')
    k = 1:object.nu;
end

if floor(k) ~= k
    error('k must be an integer number');
end
if any(k < 1) || any(k > object.nu*object.nfun)
    error(['k must be an integer number between 1 and ',num2str(object.nfun)]);
end

object.fun.plot(k);