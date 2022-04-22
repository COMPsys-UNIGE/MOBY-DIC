function optout = embeddedSystemVHDLset(fun,opts)
% MPCset   Checks and sets the options for the generation of VHDL code for
%          the PWAGFunction object
%
% OPTOUT = pwagVHDLset(FUN,OPTS)
% FUN is a PWAGFunction object, OPTS is the structure defining the options.
% See the documentation of PWAGFunction.generateVHDL for a list of all
% structure fields.

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

ctrl = fun.getController;

if ~isfield(opts,'architecture') || isempty(opts.architecture)
    optout.architecture = 'fast';
else
    if ~strcmpi(opts.architecture,'small') && ~strcmpi(opts.architecture,'fast')
        error('OPTS.architecture must be eithrer ''small'' or ''fast''');
    end
    optout.architecture = opts.architecture;
end

if isa(ctrl, 'explicitMPCctrl')
    nx = fun.nx;
    nd = fun.nd;
    
    % Number of domain dimensions
    nin = fun.ny+fun.np+numel(fun.getController.getTrackingVariable);
    nout = fun.nu;
else
    nx = fun.getNumberOfStates;
    np = fun.getNumberOfParameters;
    nu = fun.getNumberOfInputs;
    nd = fun.getNumberOfUnmeasurableInputs;
    ny = fun.getNumberOfOutputs;
    nxref = numel(ctrl.getTrackingVariable);
    nin = np + ny + nxref;
    nout = nu;
end

if ~isfield(opts,'inputResolution') || isempty(opts.inputResolution)
    optout.inputResolution = 12;
else
    if ~isnumeric(opts.inputResolution)
        error('OPTS.inputResolution must be an integer');
    end
    if opts.inputResolution ~= floor(opts.inputResolution)
        error('OPTS.inputResolution must be an integer');
    end
    if opts.inputResolution < 1
        error('OPTS.inputResolution must be greater than zero');
    end
    optout.inputResolution = opts.inputResolution;
end

if ~isfield(opts,'inputRepresentation') || isempty(opts.inputRepresentation)
    optout.inputRepresentation = 'signed';
else
    if ~strcmpi(opts.inputRepresentation,'unsigned') &&...
            ~strcmpi(opts.inputRepresentation,'signed')
        error('OPTS.inputRepresentation must be either ''unsigned'' or ''signed''');
    end
    optout.inputRepresentation = opts.inputRepresentation;
end

% Number of bits for integer and decimal part
nbit = optout.inputResolution;
nint = nbit;
ndec = nbit-nint;

if ~isfield(opts,'defaultOutput') || isempty(opts.defaultOutput)
    optout.defaultOutput = zeros(1,nout);
else
    if numel(opts.defaultOutput) ~= nout
        error(['defaultOutput must be a vector with ',num2str(nout),' elements']);
    else
         optout.defaultOutput =  opts.defaultOutput; 
    end
end

% Check if number of bits is sufficient
if any(ndec < 0)
    error('Number of bits is too low. Increase OPTS.inputResolution')
end

if ~isfield(opts,'initialState') || isempty(opts.initialState)
    optout.initialState = zeros(nx+nd,1);
else
    optout.initialState = opts.initialState;
end


if ~isfield(opts,'inputRange') || isempty(opts.inputRange)
    if strcmpi(optout.inputRepresentation,'unsigned')
        optout.inputRange.min = 0;
        optout.inputRange.max = 2^optout.inputResolution-1;
    else
        optout.inputRange.min = -2^(optout.inputResolution-1);
        optout.inputRange.max = 2^(optout.inputResolution-1)-1;
    end
else
    if ~isstruct(opts.inputRange)
        error('OPTS.inputRange must be a structure with fields ''min'' and ''max''');
    end
    if ~isfield(opts.inputRange,'min') || ~isfield(opts.inputRange,'max') 
        error('OPTS.inputRange must be a structure with fields ''min'' and ''max''');
    end
    if numel(opts.inputRange.min) ~= numel(opts.inputRange.max)
        error('OPTS.inputRange.min and OPTS.inputRange.max must have the same number of elements');
    end
    if any(opts.inputRange.min ~= floor(opts.inputRange.min)) ||...
            any(opts.inputRange.max ~= floor(opts.inputRange.max))
        error('OPTS.inputRange.min and OPTS.inputRange.max must be integer numbers');
    end
    if strcmpi(optout.inputRepresentation,'unsigned')
        if any(opts.inputRange.min < 0) || any(opts.inputRange.max < 0)
            error(['The input must be in the range [',num2str(0),' ',num2str(2^optout.inputResolution-1),']']);
        end
        if any(opts.inputRange.min > 2^optout.inputResolution-1) || any(opts.inputRange.max > 2^optout.inputResolution-1)
            error(['The input must be in the range [',num2str(0),' ',num2str(2^optout.inputResolution-1),']']);
        end
    else
        if any(opts.inputRange.min < -2^(optout.inputResolution-1)) || any(opts.inputRange.max < -2^(optout.inputResolution-1))
            error(['The input must be in the range [',num2str(-2^(optout.inputResolution-1)),' ',num2str(2^(optout.inputResolution-1)-1),']']);
        end
        if any(opts.inputRange.min > 2^(optout.inputResolution-1)-1) || any(opts.inputRange.max > 2^(optout.inputResolution-1)-1)
            error(['The input must be in the range [',num2str(-2^(optout.inputResolution-1)),' ',num2str(2^(optout.inputResolution-1)-1),']']);
        end
    end
    optout.inputRange = opts.inputRange;
end

if numel(optout.inputRange.min) == 1
    optout.inputRange.min = repmat(optout.inputRange.min,nin,1);
    optout.inputRange.max = repmat(optout.inputRange.max,nin,1);
end

if ~isfield(opts,'coeffResolution') || isempty(opts.coeffResolution)
    optout.coeffResolution = optout.inputResolution;
else
    if ~isnumeric(opts.coeffResolution)
        error('OPTS.coeffResolution must be an integer');
    end
    if opts.coeffResolution ~= floor(opts.coeffResolution)
        error('OPTS.coeffResolution must be an integer');
    end
    if opts.coeffResolution < 1
        error('OPTS.coeffResolution must be greater than zero');
    end
    optout.coeffResolution = opts.coeffResolution;
end

if ~isfield(opts,'outputResolution') || isempty(opts.outputResolution)
    optout.outputResolution = optout.coeffResolution;
else
    if ~isnumeric(opts.outputResolution)
        error('OPTS.outputResolution must be an integer');
    end
    if opts.outputResolution ~= floor(opts.outputResolution)
        error('OPTS.outputResolution must be an integer');
    end
    if opts.outputResolution < 1
        error('OPTS.outputResolution must be greater than zero');
    end

    if opts.outputResolution > optout.coeffResolution
        error('OPTS.outputResolution must be lower than or equal to OPTS.coeffResolution');
    end
    if opts.inputResolution > optout.coeffResolution
        error('OPTS.inputResolution must be lower than or equal to OPTS.coeffResolution');
    end
    optout.outputResolution = opts.outputResolution;
end

if ~isfield(opts,'outputRepresentation') || isempty(opts.outputRepresentation)
    optout.outputRepresentation = 'signed';
else
    if ~strcmpi(opts.outputRepresentation,'unsigned') &&...
            ~strcmpi(opts.outputRepresentation,'signed')
        error('OPTS.outputRepresentation must be either ''unsigned'' or ''signed''');
    end
    optout.outputRepresentation = opts.outputRepresentation;
end

if ~isfield(opts,'outputRange') || isempty(opts.outputRange)
    if strcmpi(optout.outputRepresentation,'unsigned')
        optout.outputRange.min = 0;
        optout.outputRange.max = 2^optout.outputResolution-1;
    else
        optout.outputRange.min = -2^(optout.outputResolution-1);
        optout.outputRange.max = 2^(optout.outputResolution-1)-1;
    end
else
    if ~isstruct(opts.outputRange)
        error('OPTS.outputRange must be a structure with fields ''min'' and ''max''');
    end
    if ~isfield(opts.outputRange,'min') || ~isfield(opts.outputRange,'max') 
        error('OPTS.outputRange must be a structure with fields ''min'' and ''max''');
    end
    if numel(opts.outputRange.min) ~= numel(opts.outputRange.max)
        error('OPTS.outputRange.min and OPTS.outputRange.max must have the same number of elements');
    end
    if any(opts.outputRange.min ~= floor(opts.outputRange.min)) ||...
            any(opts.outputRange.max ~= floor(opts.outputRange.max))
        error('OPTS.outputRange.min and OPTS.outputRange.max must be integer numbers');
    end
    if strcmpi(optout.outputRepresentation,'unsigned')
        if any(opts.outputRange.min < 0) || any(opts.outputRange.max < 0)
            error(['The output must be in the range [',num2str(0),' ',num2str(2^optout.outputResolution-1),']']);
        end
        if any(opts.outputRange.min > 2^optout.outputResolution-1) || any(opts.outputRange.max > 2^optout.outputResolution-1)
            error(['The output must be in the range [',num2str(0),' ',num2str(2^optout.outputResolution-1),']']);
        end
    else
        if any(opts.outputRange.min < -2^(optout.outputResolution-1)) || any(opts.outputRange.max < -2^(optout.outputResolution-1))
            error(['The output must be in the range [',num2str(-2^(optout.outputResolution-1)),' ',num2str(2^(optout.outputResolution-1)-1),']']);
        end
        if any(opts.outputRange.min > 2^(optout.outputResolution-1)-1) || any(opts.outputRange.max > 2^(optout.outputResolution-1)-1)
            error(['The output must be in the range [',num2str(-2^(optout.outputResolution-1)),' ',num2str(2^(optout.outputResolution-1)-1),']']);
        end
    end
    optout.outputRange = opts.outputRange;
end

if numel(optout.outputRange.min) == 1
    optout.outputRange.min = repmat(optout.outputRange.min,nout,1);
    optout.outputRange.max = repmat(optout.outputRange.max,nout,1);
end

if ~isfield(opts,'frequency') || isempty(opts.frequency)
    optout.frequency = 50000000;
else
    if ~isnumeric(opts.frequency)
        error('OPTS.frequency must be a number');
    end
    if opts.frequency < 0
        error('OPTS.frequency must be greater than zero');
    end    
    optout.frequency = opts.frequency;
end

if ~isfield(opts,'generateInterface') || isempty(opts.generateInterface)
    optout.generateInterface = true;
else
    optout.generateInterface = opts.generateInterface;
end

% Set or create folder in which to store the VHDL files
if isfield(opts,'folder')
    folder = opts.folder;
else
    folder = '';
end
if ~isempty(folder)
    if strcmp(folder(end),'\')
        folder(end) = '/';
    elseif ~strcmp(folder(end),'/')
        folder = [folder,'/'];
    end
else
    created = 0;
    i = 1;
    while ~created
        testfolder = [pwd,'/FPGA_circuit_',num2str(i),'/'];
        if ~exist(testfolder,'dir')
            folder = testfolder;
            created = 1;
        else
            i = i+1;
        end
    end
end
optout.folder = folder;

if ~isfield(opts,'useADC') || isempty(opts.useADC)
    optout.useADC = 1;
else
    if (opts.useADC ~= 0) && (opts.useADC ~= 1)
        error('OPTS.useADC must be either 0 or 1');
    end
    optout.useADC = opts.useADC;
end

if ~isfield(opts,'useDAC') || isempty(opts.useDAC)
    optout.useDAC = 1;
else
    if (opts.useDAC ~= 0) && (opts.useDAC ~= 1)
        error('OPTS.useDAC must be either 0 or 1');
    end
    optout.useDAC = opts.useDAC;
end

if isa(ctrl, 'implicitMPCctrl')
    if ~isfield(opts,'fpgaBoard') || isempty(opts.fpgaBoard)
        error('OPTS.fpgaBoard is requested for RTL synthesis')
    else
        optout.fpgaBoard = opts.fpgaBoard;
    end
    
    if ~isfield(opts,'coeffIntResolution') || isempty(opts.coeffIntResolution)
        error('OPTS.coeffIntResolution must be an integer');
    else
        if ~isnumeric(opts.coeffIntResolution)
            error('OPTS.coeffIntResolution must be an integer');
        end
        if opts.coeffIntResolution ~= floor(opts.coeffIntResolution)
            error('OPTS.coeffIntResolution must be an integer');
        end
        if opts.coeffIntResolution < 1
            error('OPTS.coeffIntResolution must be greater than zero');
        end
        if opts.coeffIntResolution > opts.coeffIntResolution
            error('OPTS.coeffIntResolution must be smaller than OPTS.coeffResolution');
        end
        optout.coeffIntResolution = opts.coeffIntResolution;
    end
    
    if ~isfield(opts,'ADMMparameters') || isempty(opts.ADMMparameters)
        optout.ADMMparameters.maxIter = 40;
        optout.ADMMparameters.regPar = 2;
    else
        if ~isstruct(opts.ADMMparameters)
            error('OPTS.ADMMparameters must be a structure with fields ''maxIter'' and ''regPar''');
        end
        if ~isfield(opts.ADMMparameters,'maxIter')
            optout.ADMMparameters.maxIter = 40;
        end
        if ~isfield(opts.ADMMparameters,'regPar')
            optout.ADMMparameters.regPar = 2;
        end
        optout.ADMMparameters = opts.ADMMparameters;
    end
    
    if ~isfield(opts,'range') || isempty(opts.range)
        error('Range field of options structure must be provided');
    else
        r = opts.range;
        if (~ctrl.isTracking)
            if ~isfield(r,'xmin')||~isfield(r,'xmax')||...
                    ~isfield(r,'umin')||~isfield(r,'umax')||...
                    ~isfield(r,'pmin')||~isfield(r,'pmax')||...
                    ~isfield(r,'dmin')||~isfield(r,'dmax')||...
                    ~isfield(r,'ymin')||~isfield(r,'ymax')
                error(['opts.range must be a struct with fields xmin,xmax,umin,umax,'...
                    'pmin,pmax,dmin,dmax,ymin,ymax']);
            else
                if numel(r.xmin) == 1
                    optout.range.xmin = repmat(r.xmin,1,nx);
                    optout.range.xmax = repmat(r.xmax,1,nx);
                else
                    optout.range.xmin = r.xmin;
                    optout.range.xmax = r.xmax;
                end
                if numel(r.umin) == 1
                    optout.range.umin = repmat(r.umin,1,nu);
                    optout.range.umax = repmat(r.umax,1,nu);
                else
                    optout.range.umin = r.umin;
                    optout.range.umax = r.umax;
                end
                if numel(r.dmin) == 1
                    optout.range.dmin = repmat(r.dmin,1,nd);
                    optout.range.dmax = repmat(r.dmax,1,nd);
                else
                    optout.range.dmin = r.dmin;
                    optout.range.dmax = r.dmax;
                end
                if numel(r.pmin) == 1
                    optout.range.pmin = repmat(r.pmin,1,np);
                    optout.range.pmax = repmat(r.pmax,1,np);
                else
                    optout.range.pmin = r.pmin;
                    optout.range.pmax = r.pmax;
                end
                if numel(r.ymin) == 1
                    optout.range.ymin = repmat(r.ymin,1,ny);
                    optout.range.ymax = repmat(r.ymax,1,ny);
                else
                    optout.range.ymin = r.ymin;
                    optout.range.ymax = r.ymax;
                end
            end
        else
            if ~isfield(r,'xmin')||~isfield(r,'xmax')||...
                    ~isfield(r,'umin')||~isfield(r,'umax')||...
                    ~isfield(r,'pmin')||~isfield(r,'pmax')||...
                    ~isfield(r,'dmin')||~isfield(r,'dmax')||...
                    ~isfield(r,'ymin')||~isfield(r,'ymax')||...
                    ~isfield(r,'xrefmin')||~isfield(r,'xrefmax')
                error(['opts.range must be a struct with fields xmin,xmax,umin,umax,'...
                    'pmin,pmax,dmin,dmax,ymin,ymax,xremax,xrefmin']);
            else
                if numel(r.xmin) == 1
                    optout.range.xmin = repmat(r.xmin,1,nx);
                    optout.range.xmax = repmat(r.xmax,1,nx);
                else
                    optout.range.xmin = r.xmin;
                    optout.range.xmax = r.xmax;
                end
                if numel(r.umin) == 1
                    optout.range.umin = repmat(r.umin,1,nu);
                    optout.range.umax = repmat(r.umax,1,nu);
                else
                    optout.range.umin = r.umin;
                    optout.range.umax = r.umax;
                end
                if numel(r.dmin) == 1
                    optout.range.dmin = repmat(r.dmin,1,nd);
                    optout.range.dmax = repmat(r.dmax,1,nd);
                else
                    optout.range.dmin = r.dmin;
                    optout.range.dmax = r.dmax;
                end
                if numel(r.pmin) == 1
                    optout.range.pmin = repmat(r.pmin,1,np);
                    optout.range.pmax = repmat(r.pmax,1,np);
                else
                    optout.range.pmin = r.pmin;
                    optout.range.pmax = r.pmax;
                end
                if numel(r.ymin) == 1
                    optout.range.ymin = repmat(r.ymin,1,ny);
                    optout.range.ymax = repmat(r.ymax,1,ny);
                else
                    optout.range.ymin = r.ymin;
                    optout.range.ymax = r.ymax;
                end
                if numel(r.xrefmin) == 1
                    optout.range.xrefmin = repmat(r.xrefmin,1,nxref);
                    optout.range.xrefmax = repmat(r.xrefmax,1,nxref);
                else
                    optout.range.xrefmin = r.xrefmin;
                    optout.range.xrefmax = r.xrefmax;
                end
            end
        end
    end
end