function ndim = getNumberOfDimensions(object)
% getNumberOfDimensions   Gets the number of dimensions of the control function 
%
% NDIM = getNumberOfDimensions(OBJ)
% If the controller is for regulation, the number of dimensions 
% (i.e., the number of inputs to the control function) is equal to the sum
% of the number of states, parameters and unmeasurable inputs. If the
% controller is for tracking, the number of variables to track must be
% added. OBJ is the controller object.

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


if object.tracking
    ndim = object.nx+object.np+object.nd+numel(object.trackvar);
else
    ndim = object.nx+object.np+object.nd;
end
    
end