function optout = MPCset(sys,options)
% MPCset   Checks and sets the options for MPCctrl object
%
% OPTOUT = MPCset(SYS,OPTS)
% SYS is the ltiSys or pwaSys object to control, OPTS is the structure
% defining the options. See the documentation of MPCctrl for a list of all
% structure fields.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Tomaso Poggi (tpoggi@essbilbao.org)
% Alessandro Ravera (alessandro.ravera@edu.unige.it)
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


nx = sys.getNumberOfStates();
nu = sys.getNumberOfInputs();
ny = sys.getNumberOfOutputs();
np = sys.getNumberOfParameters();
nd = sys.getNumberOfUnmeasurableInputs();
% Names
xnames = cell(nx,1);
for i = 1:nx
    xnames{i} = sys.getStateNames{i};
end

if ~isfield(options,'ctrlvariable')
    optout.ctrlvariable = 'state';
elseif isempty(options.ctrlvariable)
    optout.ctrlvariable = 'state';
else
    if ~strcmpi(options.ctrlvariable,'state') &&...
            ~strcmpi(options.ctrlvariable,'output')
        error('ctrlvariable must be either ''state'' or ''output''');
    end
    optout.ctrlvariable = options.ctrlvariable;
end

% Set default values
if ~isfield(options,'Q')
    error('Matrix Q must be provided!')
elseif isempty(options.Q)
    error('Matrix Q must be provided!')
else
    optout.Q = options.Q;
end

if strcmpi(optout.ctrlvariable,'state')
    if any(size(optout.Q) ~= nx)
        error(['Q must be a ',num2str(nx),'x',num2str(nx),' matrix']);
    end
else
    if any(size(optout.Q) ~= ny)
        error(['Q must be a ',num2str(ny),'x',num2str(ny),' matrix']);
    end
end

if ~isfield(options,'R')
    error('Matrix R must be provided!')
elseif isempty(options.R)
    error('Matrix R must be provided!')
else
    optout.R = options.R;
end

if any(size(optout.R) ~= nu)
    error(['R must be a ',num2str(nu),'x',num2str(nu),' matrix']);
end

if strcmpi(optout.ctrlvariable,'state')
    if ~isfield(options,'ref')
        optout.ref = zeros(nx,1);
    elseif isempty(options.ref)
        optout.ref = zeros(nx,1);
    else
        optout.ref = options.ref;
        optout.ref = optout.ref(:);
    end
    if numel(optout.ref) ~= nx
        error(['ref must be a vector with ',num2str(nx),' elements']);
    end
else
    if ~isfield(options,'ref')
        optout.ref = zeros(ny,1);
    elseif isempty(options.ref)
        optout.ref = zeros(ny,1);
    else
        optout.ref = options.ref;
        optout.ref = optout.ref(:);
    end
    if numel(optout.ref) ~= ny
        error(['ref must be a vector with ',num2str(ny),' elements']);
    end
end

if ~isfield(options,'uref')
    optout.uref = zeros(nu,1);
elseif isempty(options.uref)
    optout.uref = zeros(nu,1);
else
    optout.uref = options.uref;
    optout.uref = optout.uref(:);
end

if numel(optout.uref) ~= nu
    error(['uref must be a vector with ',num2str(nu),' elements']);
end

if ~isfield(options,'norm')
    optout.norm = 2;
elseif isempty(options.norm)
    optout.norm = 2;
else
    optout.norm = options.norm;
end

A = sys.getMatrices('A');
B = sys.getMatrices('B');
C = sys.getMatrices('C');

if ~isfield(options,'N')
    error('Prediction horizon N must be provided!')
elseif isempty(options.N)
    error('Prediction horizon N must be provided!')
else
    optout.N = options.N;
end

if ~isfield(options,'Nu')
    optout.Nu = optout.N;
elseif isempty(options.Nu)
    optout.Nu = optout.N;
else
    optout.Nu = options.Nu;
end

if optout.Nu < 1 || optout.Nu > optout.N
    error('Nu cannot be lower than 1 or greater than N')
end

if ~isfield(options,'tracking')
    optout.tracking = false;
elseif isempty(options.tracking)
    optout.tracking = false;
else
    optout.tracking = boolean(options.tracking);
end

if ~isfield(options,'trackvariable')
    if optout.tracking
        error('Tracking variable trackvariable must be provided')
    else
        optout.trackvariable = [];
    end
elseif isempty(options.trackvariable)
    if optout.tracking
        error('Tracking variable trackvariable must be provided')
    else
        optout.trackvariable = [];
    end
else
    optout.trackvariable = options.trackvariable;
end

if optout.tracking
    if isnumeric(options.trackvariable)
        trackvar = options.trackvariable;
        if any(floor(trackvar) ~= trackvar)
            error('trackvariable must be an integer array')
        end
        if any(trackvar < 1 | trackvar > nx)
            error('The state variables you want to track do not exist');
        end
    else
        error('trackvariable must be either an integer or a string')
    end
    if numel(trackvar) > nx
        error(['trackvariable must be a vector with atmost ', num2str(nx), 'elements']);
    end
    optout.trackvariable = trackvar;
end

if ~isfield(options,'constantInputAfterNu')
    optout.constantInputAfterNu = false;
else
    optout.constantInputAfterNu = options.constantInputAfterNu;
end

if ~isfield(options,'P')
    optout.P = optout.Q;
elseif isempty(options.P)
    optout.P = optout.Q;
else
    optout.P = options.P;
end

if ~isfield(options,'K')
    optout.K = zeros(nu,nx);
elseif isempty(options.K)
    optout.K = zeros(nu,nx);
else
    optout.K = options.K;
end

if ~isfield(options,'O')
    optout.O = zeros(nu,1);
elseif isempty(options.O)
    optout.O = zeros(nu,1);
else
    optout.O = options.O;
end

if size(optout.K) ~= [nu nx]
    error(['K must be a ',num2str(nu),'x',num2str(nx),' matrix']);
end

if size(optout.O) ~= [nu 1]
    error(['O must be a ',num2str(nu),'x1 matrix']);
end

if strcmpi(optout.ctrlvariable,'state')
    if any(size(optout.P) ~= nx)
        error(['P must be a ',num2str(nx),'x',num2str(nx),' matrix']);
    end
else
    if any(size(optout.P) ~= ny)
        error(['P must be a ',num2str(ny),'x',num2str(ny),' matrix']);
    end
end

if isfield(options,'rho')
    optout.rho = options.rho;
end

if isfield(options,'defaultOutput')
    optout.defaultOutput = options.defaultOutput;
end

if isfield(options,'algorithm')
    optout.algorithm = options.algorithm;
    if isfield(options,'algoParameters')
        optout.algoParameters = options.algoParameters;
    end
end

if ~isfield(options,'trajectoryTracking')
    optout.trajectoryTracking = false;
else
    optout.trajectoryTracking = options.trajectoryTracking;
end

end
