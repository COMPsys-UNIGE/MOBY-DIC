function varargout = getGain(object)
% getGain   Gets the gains of the Kalman filter
%
% K = getGain(OBJ)
% Gets the gain K of the Kalman filter. OBJ is the kalmanFilter object.
%
% [Kx, Kd] = getGain(OBJ)
% Gets gain K split into its two components:
%      _  _
%     | Kx |
% K = |    |
%     |_Kd_|
%
% where Kx is related to the system states and Kd to the unmeasurable
% inputs.

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

if nargout == 0 || nargout == 1
    
    varargout{1} = object.K;
    
elseif nargout == 2
    
    K = object.K;
    
    Kx = K(1:object.nx,:);
    Kd = K(object.nx+1:end,:);
    
    varargout{1} = Kx;
    varargout{2} = Kd;
    
else
    error('Wrong number of outputs');
end

end