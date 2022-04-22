function varargout = getMatrices(object,varargin)
% getMatrices    Gets the matrices defining the LTI system
%
% [A, B, C, D, Ex, Ey, Fx, Fy, Gx, Gy, H, K] = getMatrices(OBJ)
% Gets matrices A, B, C, D, Ex, Ey, Fx, Fy, Gx and Gy (see ltiSys for
% an explanation of their meaning). OBJ is the ltiSys object.
%
% M = getMatrices(OBJ,MATR)
% Gets only the matrix M specified in label MATR (allowed values are 'A', 
% 'B', 'C', 'D', 'Ex', 'Ey', 'Fx', 'Fy', 'Gx' or 'Gy').

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


if nargin == 1
    varargout{1} = object.A;
    varargout{2} = object.B;
    varargout{3} = object.C;
    varargout{4} = object.D;
    varargout{5} = object.Ex;
    varargout{6} = object.Ey;
    varargout{7} = object.Fx;
    varargout{8} = object.Fy;
    varargout{9} = object.Gx;
    varargout{10} = object.Gy;
    varargout{11} = object.H;
    varargout{12} = object.K;
elseif nargin == 2
    if strcmpi(varargin{1},'A')
        varargout{1} = object.A;
    elseif strcmpi(varargin{1},'B')
        varargout{1} = object.B;
    elseif strcmpi(varargin{1},'C')
        varargout{1} = object.C;
    elseif strcmpi(varargin{1},'D')
        varargout{1} = object.D;
    elseif strcmpi(varargin{1},'Ex')
        varargout{1} = object.Ex;
    elseif strcmpi(varargin{1},'Ey')
        varargout{1} = object.Ey;
    elseif strcmpi(varargin{1},'Fx')
        varargout{1} = object.Fx;
    elseif strcmpi(varargin{1},'Fy')
        varargout{1} = object.Fy;
    elseif strcmpi(varargin{1},'Gx')
        varargout{1} = object.Gx;
    elseif strcmpi(varargin{1},'Gy')
        varargout{1} = object.Gy;
        elseif strcmpi(varargin{1},'H')
        varargout{1} = object.H;
        elseif strcmpi(varargin{1},'K')
        varargout{1} = object.K;
    else
        error(['Input argument matr must be either ''A'', ''B'', ''C''',...
        '''D'', ''Ex'', ''Ey'', ''Fx'', ''Fy'', ''Gx'', ''Gy'', ''H'' or ''K''']);
    end
else
    error('Wrong number of inputs');
end

end