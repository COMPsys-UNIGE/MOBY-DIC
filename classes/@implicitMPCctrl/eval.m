function u = eval(object,x,varargin)
% eval   Evaluates the MPC control function
%
% U = eval(OBJ,x)
% Evaluates the control function at the current system state x.
% x must be a vector of legnth NX. U is a vector of length NU.
% NX and NU are the number of system states and inputs, respectively.
% OBJ is the MPCctrl object.
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
ny = object.ny;
N = object.options.N;
Nu = object.options.Nu;

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

% Check the algorithm for solving QP (quadprog or admm)
if isfield(object.options, 'algorithm')
    algorithm = object.options.algorithm;
    if algorithm ~= "ADMM"
        warning('%s is not a supported algorithm!', algorithm);
    end
else
    algorithm = "";
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
    if size(ref,2) == 1
        ref = repmat(ref,1,npoints);
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

% Define matrices of QP
if object.options.tracking
    persistent uOld;
    if isempty(uOld)
        if ~isfield(object.options, 'defaultOutput')
            uOld = zeros(nu,1);
        else
            if isempty(object.options.defaultOutput)
                uOld = zeros(nu,1);
            else
                uOld = object.options.defaultOutput;
            end
        end
    end
    
    if algorithm == "ADMM"
        [QP1, QP2, ineqL, ineqR, eqL, eqR] = computeQP(object, x, p, d, ref, uOld, "ADMM");
    else
        [QP1, QP2, ineqL, ineqR, eqL, eqR] = computeQP(object, x, p, d, ref, uOld);
    end
else
    if algorithm == "ADMM"
        [QP1, QP2, ineqL, ineqR, eqL, eqR] = computeQP(object, x, p, d, ref, "ADMM");
    else
        [QP1, QP2, ineqL, ineqR, eqL, eqR] = computeQP(object, x, p, d, ref);
    end
end

% Solve QP with admm
if algorithm == "ADMM"
    
    if isfield(object.options, 'algoParameters')
        regPar = object.options.algoParameters.regPar;
        admmIter = object.options.algoParameters.maxIter;
    else
        regPar = 2;
        admmIter = 40;
    end

    if object.isTracking
        x = [x; uOld];
    end
    
    [M1, M2, M3, M4, v1, v2, v3] = object.computeADMMmatrices(QP1, QP2, ineqL, ineqR, eqL, eqR, regPar);
    
    z = object.admm(M1, M2, M3, M4, v1, v2, v3, regPar, admmIter);
    
% Solve QP with quadprog
else
    evalc('z = quadprog(QP1, QP2, ineqL, ineqR, eqL, eqR);');
end
        
try
    u = z(1:nu);
catch
    error('The problem is infeasible, try to consider soft constraints on state variables.');
end

if object.options.tracking
    u = u + uOld;
    uOld = u;
end

end