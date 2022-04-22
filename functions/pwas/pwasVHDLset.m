function optout = pwasVHDLset(fun,opts,idx)
% MPCset   Checks and sets the options for the generation of VHDL code for
%          the PWASFunction object
%
% OPTOUT = pwasVHDLset(FUN,OPTS)
% FUN is a PWASFunction object, OPTS is the structure defining the options.
% See the documentation of PWASFunction.generateVHDL for a list of all
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

if ~isfield(opts,'architecture') || isempty(opts.architecture)
    optout.architecture = 'fast';
else
    if ~strcmpi(opts.architecture,'small') && ~strcmpi(opts.architecture,'fast')
        error('OPTS.architecture must be eithrer ''small'' or ''fast''');
    end
    optout.architecture = opts.architecture;
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

% Number of domain dimensions
ndim = fun.getDomainDimensions();
nfun = numel(idx);
% Number of subdivisions per dimension
np = fun.getNumberOfPartitions();

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
    optout.inputRange.min = repmat(optout.inputRange.min,ndim,1);
    optout.inputRange.max = repmat(optout.inputRange.max,ndim,1);
end

if numel(optout.inputRange.min) ~= ndim || numel(optout.inputRange.max) ~= ndim
    error(['OPTS.inputRange.min and .max must have 1 or ',num2str(ndim),' elements']);
end


% Number of actual bits within the input range
nbiteff = ceil(log2(abs(optout.inputRange.max-optout.inputRange.min)));

% Number of bits for integer and decimal part
nint = ceil(log2(np(:)+1));
ndec = nbiteff(:)-nint(:);

% Check if number of bits is sufficient
if any(ndec) < 0
    error('Number of bits is too low. Increase OPTS.inputResolution or OPTS.inputRange')
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
    optout.outputRange.min = repmat(optout.outputRange.min,nfun,1);
    optout.outputRange.max = repmat(optout.outputRange.max,nfun,1);
end

if ~isfield(opts,'frequency') || isempty(opts.frequency)
    optout.frequency = 50000000;
else
    if ~isnumeric(opts.frequency)
        error('OPTS.frequency must be a number');
    end
    if opts.frequency < 1
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