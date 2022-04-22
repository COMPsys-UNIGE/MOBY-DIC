function object = setInputNames(object,inputNames)
% setInputNames   Gives a mnemonical name to the system inputs
%
% OBJ = setInputNames(OBJ,INPUTNAMES)
% INPUTNAMES must be a cell array whose elements are strings representing
% the name of the corresponding system input. OBJ is the dynSys object.
%
% See also dynSys/setOutputNames, dynSys/setParameterNames, 
% dynSys/setStateNames, dynSys/setUnmeasurableInputNames.
 
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

if numel(inputNames) ~= object.nu
    error('Number of names is different from number of inputs!');
end

object.unames = inputNames(:)';