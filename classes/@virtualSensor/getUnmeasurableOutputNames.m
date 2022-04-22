function names = getUnmeasurableOutputNames(object)
% getUnmeasurableOutputNames   Gets the mnemonical name of the system 
%                              unmeasurable outputs
%
% NAMES = getUnmeasurableOutputNames(OBJ)
% NAMES is a cell array whose elements are strings representing
% the name of the corresponding system unmeasurable output. OBJ is the 
% virtualSensor object.
 
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

names = object.znames;