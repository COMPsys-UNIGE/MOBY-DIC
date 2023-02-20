function varargout = eval(object,x,varargin)
% eval   Evaluates the MPC control function
%
% U = eval(OBJ,X)
% Evaluates the control function at the system states specified by matrix 
% X. X must be a NX x NPOINTS matrix. U is a NU x NPOINTS matrix.
% NX and NU are the number of system states and inputs, respectively.
% OBJ is the explicitMPCctrl object.
%
% U = eval(OBJ,X,REF)
% If the controller is designed for tracking, the reference state or output 
% REF must also be provided. X must be a NX x NPOINTS matrix. REF must be a 
% matrix with as many rows as the number of variables to track (defined in 
% OPTS.trackvar in the class constructor) and NPOINTS columns. If REF has 
% only one column, the same value is used for all points. U is a 
% NU x NPOINTS matrix. NX and NU are the number of system states and inputs, 
% respectively.
%
% U = eval(OBJ,X,P,D)
% Evaluates the control function at the system states specified by matrix 
% X, with the parameters specified by matrix P and the unmeasurable inputs 
% specified by matrix D. X must be a NX x NPOINTS matrix. P must be either
% a NP x 1 or a NP x NPOINTS matrix. If P has one column, the same value is 
% replicated for all points. D must be either a ND x 1 or a ND x NPOINTS
% matrix. If D has one column, the same value is replicated for all points. 
% U is a NU x NPOINTS matrix. NX, NU, NP and ND are the number of system 
% states, inputs, parameters and unmeasurable inputs, respectively. 
%
% U = eval(OBJ,X,P,D,REF)
% If the controller is designed for tracking, the reference state or output
% REF must also be provided.
%
% Mex files are used to speedup the computation.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Mateo Lodi (matteo.lodi@edu.unige.it)
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

if nargin == 2
    p = [];
    d = [];
    ref = [];
elseif nargin == 3
    p = [];
    d = [];
    ref = varargin{1};
elseif nargin == 4
    p = varargin{1};
    d = varargin{2};
    ref = [];
elseif nargin == 5
    p = varargin{1};
    d = varargin{2};
    ref = varargin{3};
else
    error('Wrong input arguments');
end

if object.tracking && isempty(ref)
    error('ref must be provided for tracking controllers');
end

% Number of variables
nx = object.nx;
np = object.np;
nd = object.nd;
nu = object.nu;
nfun = object.nfun;

% Check on x
if size(x,1) ~= object.nx
    error(['Matrix x must be of a matrix with ',num2str(nx),' rows']);
end

% Check on p
if size(p,1) ~= object.np
    error(['Matrix p must be of a matrix with ',num2str(np),' rows']);
end

% Check on d
if size(d,1) ~= object.nd
    error(['Matrix d must be of a matrix with ',num2str(nd),' rows']);
end

% Check on ref
if ~isempty(ref)
    if size(ref,1) ~= numel(object.trackvar)
        error(['ref must be a matrix with ',num2str(numel(object.trackvar)),' rows']);
    end
end

% Number of points x
npoints = size(x,2);

% If only one parameter or unmeasurable input is provided, it is replicated
% for each x
if size(p,2) == 1
    p = repmat(p,1,npoints);
end
if size(d,2) == 1
    d = repmat(d,1,npoints);
end
if object.tracking
    if size(ref,2) == 1 || size(ref,2) ~= npoints
        ref = repmat(ref(:,1),1,npoints);
    end
end

% Number of points
nptp = size(p,2);
nptd = size(d,2);
nptref = size(ref,2);

if ~isempty(p)
    if nptp ~= npoints
        error(['p must be a matrix with either 1 or ',num2str(npoints),' columns'])
    end
end
if ~isempty(d)
    if nptd ~= npoints
        error(['d must be a matrix with either 1 or ',num2str(npoints),' columns'])
    end
end
if ~isempty(ref)
    if nptref ~= npoints
        error(['ref must be a matrix with either 1 or ',num2str(npoints),' columns'])
    end
end

% Extended state
x = [x; p; d; ref];

if nargout < 2
    u = object.fun.eval(x,1:nu);
    varargout{1} = u;
else
    u = object.fun.eval(x);
    varargout{1} = u(:,1:nu);
    varargout{2} = cell(nfun,1);
    for i = 1:nfun
        varargout{2}{i} = u(:,(1:nu)+(i-1)*nu);
    end
end
end