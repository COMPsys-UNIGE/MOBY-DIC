function [xcur, ycur] = computeOutput(object,x,u,p,y)
% computeOutput   Computes the Kalman filter output
%
% [XO, YO] = computeOutput(OBJ,X,U,P,Y)
% Computes the estimated updated state (XO) and output (YO) at current 
% time instant based on the estimated state (X), input (U) and parameter 
% (P) and on the measured output (Y) at current time instant. 
% OBJ is the kalmanFilter object. 

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


    if numel(x) ~= object.nx+object.nd
        error(['xcur must be a vector with ',num2str(object.nx+object.nd),' elements']);
    end
    if numel(u) ~= object.nu
        error(['u must be a vector with ',num2str(object.nu),' elements']);
    end
    if numel(p) ~= object.np
        error(['p must be a vector with ',num2str(object.np),' elements']);
    end

    ii = object.findDynamics(x,p);
    
    [xcur, ycur] = object.filters(ii).kFilter.computeOutput(x,u,p,y);

end