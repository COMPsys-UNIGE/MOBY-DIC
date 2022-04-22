% eval   Evaluates the approximate MPC control function that have to be
%        applied  for the specific dynamic
%
% U = eval(OBJ,X)
% Evaluates the control function at the system states specified by matrix
% X. X must be a NPOINTS x NX matrix. U is a NPOINTS x NU matrix.
% NX and NU are the number of system states and inputs, respectively.
% OBJ is the MPCctrl object.
%
% U = eval(OBJ,X,XREF)
% If the controller is designed for tracking, the reference state XREF must
% also be provided. XREF must be a matrix with as many columns as the
% number of variables to track (defined in OPTS.trackvar) and as many rows
% as the number of rows of X. If XREF has only one row, the same value is
% used for all points.
%
% U = eval(OBJ,X,P,D)
% Evaluates the control function at the system states specified by matrix
% X, with the parameters specified by matrix P and the unmeasurable inputs
% specified by matrix D. X must be a NPOINTS x NX matrix. P must be either
% a 1 x NP or a NPOINTS x NP matrix, being NP the number of system
% parameters. If P has one row, the same value is replicated for all
% points. P must be either a 1 x ND or a NPOINTS x ND matrix, being ND the
% number of unmeasurable inputs. If D has one row, the same value is
% replicated for all points.
%
% U = eval(OBJ,X,P,D,XREF)
% If the controller is designed for tracking, the reference state XREF must
% also be provided.
%
% If the "usemex" option is enabled in MOBY-DIC settings (see
% MOBYDICsettings), mex files are used to speedup the computation.

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
% -------------------------------------------------------------------------

function varargout = eval(object,x,varargin)

if nargin == 2
    p = [];
    d = [];
    xref = [];
elseif nargin == 3
    p = [];
    d = [];
    xref = varargin{1};
elseif nargin == 4
    p = varargin{1};
    d = varargin{2};
    xref = [];
elseif nargin == 5
    p = varargin{1};
    d = varargin{2};
    xref = varargin{3};
else
    error('Wrong input arguments');
end

nxref = numel(object.getTrackingVariable);

if size(x,1) ~= object.nx
    if size(x,2) == object.nx
        x = x';
    else
    error(['X must be a matrix with  ',num2str(object.nx),' rows'])
    end
end

if size(p,1) ~= object.np
    if size(p,2) == object.np
        p = p';
    else
    error(['P must be a matrix with  ',num2str(object.np),' rows'])
    end
end

if size(d,1) ~= object.nd
    if size(d,2) == object.nd
        d = d';
    else
    error(['D must be a matrix with  ',num2str(object.nd),' rows'])
    end
end

if size(xref,1) ~= nxref
    if size(xref,2) == nxref
        xref = xref';
    else
    error(['XREF must be a matrix with  ',num2str(nxref),' rows'])
    end
end

npts = size(x,2);

if numel(d) == 0
    d = zeros(0,npts);
end

if numel(p) == 0
    p = zeros(0,npts);
end

if numel(xref) == 0
    xref = zeros(0,npts);
end

if size(p,2) == 1
    p = repmat(p,size(p,1),npts);
end

if size(d,2) == 1
    d = repmat(d,size(d,1),npts);
end

if size(xref,2) == 1
    d = repmat(xref,size(xref,1),npts);
end


y = object.getFunction.eval([x;p;d;xref]);

if nargout < 2
    varargout{1} = y(1:object.nu,:);
else
    varargout{1} = y(1:object.nu,:);
    varargout{2} = y(object.nu:end,:);
end

end

