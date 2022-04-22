function optout = MPCctrlCset(fun,opts)
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

% if ~isfield(opts,'architecture') || isempty(opts.architecture)
%     optout.architecture = 'fast';
% else
%     if ~strcmpi(opts.architecture,'small') && ~strcmpi(opts.architecture,'fast')
%         error('OPTS.architecture must be eithrer ''small'' or ''fast''');
%     end
%     optout.architecture = opts.architecture;
% end

if ~isfield(opts,'generate_main')
    opt.generate_main = true;
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

if ~isfield(opts,'inputRange') || isempty(opts.inputRange)
        optout.inputRange.min = -2^(optout.inputResolution-1);
        optout.inputRange.max = 2^(optout.inputResolution-1)-1;
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
   
        if any(opts.inputRange.min < -2^(optout.inputResolution-1)) || any(opts.inputRange.max < -2^(optout.inputResolution-1))
            error(['The input must be in the range [',num2str(-2^(optout.inputResolution-1)),' ',num2str(2^(optout.inputResolution-1)-1),']']);
        end
        if any(opts.inputRange.min > 2^(optout.inputResolution-1)-1) || any(opts.inputRange.max > 2^(optout.inputResolution-1)-1)
            error(['The input must be in the range [',num2str(-2^(optout.inputResolution-1)),' ',num2str(2^(optout.inputResolution-1)-1),']']);
        end

    optout.inputRange = opts.inputRange;
end

if numel(optout.inputRange.min) ~= 1
    error('inputRange.min and inputRange.max must be scalars');
end

if ~isfield(opts,'outputResolution') || isempty(opts.outputResolution)
    optout.outputResolution = optout.inputResolution;
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
    optout.outputResolution = opts.outputResolution;
end

if ~isfield(opts,'outputRange') || isempty(opts.outputRange)  
        optout.outputRange.min = -2^(optout.outputResolution-1);
        optout.outputRange.max = 2^(optout.outputResolution-1)-1;
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
    
        if any(opts.outputRange.min < -2^(optout.outputResolution-1)) || any(opts.outputRange.max < -2^(optout.outputResolution-1))
            error(['The output must be in the range [',num2str(-2^(optout.outputResolution-1)),' ',num2str(2^(optout.outputResolution-1)-1),']']);
        end
        if any(opts.outputRange.min > 2^(optout.outputResolution-1)-1) || any(opts.outputRange.max > 2^(optout.outputResolution-1)-1)
            error(['The output must be in the range [',num2str(-2^(optout.outputResolution-1)),' ',num2str(2^(optout.outputResolution-1)-1),']']);
        end
    optout.outputRange = opts.outputRange;
end

if ~isfield(opts,'initialCondition') || isempty(opts.initialCondition)
    optout.initialCondition = 0;
else
    optout.initialCondition = opts.initialCondition;
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
        testfolder = [pwd,'/C_circuit_',num2str(i),'/'];
        if ~exist(testfolder,'dir')
            folder = testfolder;
            created = 1;
        else
            i = i+1;
        end
    end
end
optout.folder = folder;

if ~isfield(opts,'generate_main')
    optout.generate_main = 'true';
else
    optout.generate_main = opts.generate_main;
end