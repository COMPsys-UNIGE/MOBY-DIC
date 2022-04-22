function varargout = getMatrices(object,varargin)
% getMatrices   Gets the matrices defining the Kalman filter
%
% [A B C D Gx Gy H K] = getMatrices(OBJ)
% Gets matrices A, B, C, D, Gx and Gy defining the Kalman filter
% discrete-time dynamical system. OBJ is the kalmanFilter object.
%
% M = getMatrices(OBJ,MATR)
% Gets only the matrix M specified in label MATR (allowed values are 'A',
% 'B', 'C', 'D', 'Gx', 'Gy', 'H' or 'K').

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


nDyn = object.nDyn;

A = cell(nDyn,1);
B = cell(nDyn,1);
C = cell(nDyn,1);
D = cell(nDyn,1);
Gx = cell(nDyn,1);
Gy = cell(nDyn,1);
H = cell(nDyn,1);
K = cell(nDyn,1);

for i=1:nDyn
    [A{i}, B{i}, C{i}, D{i}, Gx{i}, Gy{i}] = object.filters(i).kFilter.getMatrices;
    H{i} = object.filters(i).H;
    K{i} = object.filters(i).K;
end

if nargin == 1
    varargout{1} = A;
    varargout{2} = B;
    varargout{3} = C;
    varargout{4} = D;
    varargout{5} = Gx;
    varargout{6} = Gy;
    varargout{7} = H;
    varargout{8} = K;
elseif nargin == 2
    if strcmpi(varargin{1},'A')
        varargout{1} = A;
    elseif strcmpi(varargin{1},'B')
        varargout{1} = B;
    elseif strcmpi(varargin{1},'C')
        varargout{1} = C;
    elseif strcmpi(varargin{1},'D')
        varargout{1} = D;
    elseif strcmpi(varargin{1},'Gx')
        varargout{1} = Gx;
    elseif strcmpi(varargin{1},'Gy')
        varargout{1} = Gy;
    elseif strcmpi(varargin{1},'H')
        varargout{1} = H;
    elseif strcmpi(varargin{1},'K')
        varargout{1} = K;
    else
        error(['Input argument matr must be either ''A'', ''B'', ''C''',...
            '''D'', ''Gx'' or ''Gy''']);
    end
else
    error('Wrong number of inputs');
end

end