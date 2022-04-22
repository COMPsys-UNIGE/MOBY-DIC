%% GET IDENTIFICATION INFORMATION
% Gets information about ridge regression process used for sensor
% identification
%
% SYNTAX
%
% info = getIdentificationInformation(object)
%
% info is a structure with the following fields:
%
% * wnorm: it is the 2 norm of the weights vector w (for 0 order Tikhonov
%           regularization) or the 2 norm of L w (for 1 order Tikhonov
%           regularization)
% * residual: it is the 2 norm of the residual $G w \textrm{--} d$ ( i.e. $f(X) \textrm{--} Y$ )
% * test_error: it is the 2 norm of the residual computed in the test set
%                (if provided), i.e. $f(Xt) \textrm{--} Yt$.
%
% see function ridgeRegression for further details.
%
% ACKNOWLEDGEMENTS
%
% Contributors:
% 
% * Alberto Oliveri (alberto.oliveri@unige.it)
% 
% Copyright is with:
% 
% * Copyright (C) 2011 University of Genoa, Italy.

% -------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% -------------------------------------------------------------------------

function info = getIdentificationInformation(object)
    % getIdentificationInformation  Gets information about ridge regression 
    %                               process used for sensor identification
    %
    % INFO = getIdentificationInformation(OBJ)
    % INFO is a structure with the following fields:
    %  - wnorm: it is the 2 norm of the weights vector w (for 0 order Tikhonov
    %           regularization) or the 2 norm of L w (for 1 order Tikhonov
    %           regularization)
    %  - residual: it is the 2 norm of the residual $G w \textrm{--} d$ ( i.e. $f(X) \textrm{--} Y$ )
    %  - test_error: it is the 2 norm of the residual computed in the test set
    %                (if provided), i.e. $f(Xt) \textrm{--} Yt$.
 
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
    
info = object.info;

end