function ctrlvar = getControlVariable(object)
% getControlVariable   Indicates if the controller regulates (or tracks)
%                      system states or outputs
%
% CTRLVAR = getControlVariable(OBJ)
% String CTRLVAR can be 'state' or 'output' and indicates if the controller 
% is designed for the regulation (or tracking) of system states or outputs. 
% OBJ is the controller object.
%
% See also controller/isTracking. 

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



ctrlvar = object.ctrlvar;

end