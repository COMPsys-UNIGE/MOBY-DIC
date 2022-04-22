function optout = PWASset(fun,options)
% MPCset   Checks and sets the options for PWAS approximation
%
% OPTOUT = PWASset(FUN,OPTS)
% FUN is a MOBYFunction object, OPTS is the structure defining the options.
% See the documentation of MOBYFunction.approximate for a list of all
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


if ~isstruct(options)
    error('OPTS must be a structure');
end

% Domain dimension
nx = fun.getDomainDimensions();

% Set default values
if ~isfield(options,'P') && ~isfield(options,'np')
    error('Either P or np must be provided!')
end

if isfield(options,'P') && isfield(options,'np')
    if isempty(options.P) && isempty(options.np)
        error('Either P or np must be provided!')
    end
end

if ~isfield(options,'nsamples')
    optout.nsamples = 5;
elseif isempty(options.nsamples)
    optout.nsamples = 5;
else
    if ~isnumeric(options.nsamples)
        error('OPTS.nsamples must be a scalar number');
    end
    if numel(options.nsamples) > 1
        error('OPTS.nsamples must be a scalar number');
    end
    if options.nsamples < 1
        error('OPTS.nsamples must be greater than 0');
    end
    optout.nsamples = options.nsamples;
end



if isfield(options,'P')
    if ~iscell(options.P)
        error('OPTS.P must be a cell array');
    end
    if isfield(options,'np')
        error('OPTS.P and OPTS.NP cannot be set together');
    end
    
    if numel(options.P) ~= nx
        error(['OPTS.P must be a cell array with ',num2str(nx),' elements']);
    end
    optout.P = options.P;
    optout.np = [];
else
    if numel(options.np) == 1
        options.np = repmat(options.np,nx,1);
    end
    
    if numel(options.np) ~= nx
        error(['OPTS.np must be an array with ',num2str(nx),' elements']);
    end
    optout.np = options.np;
    optout.P = [];
end

if ~isfield(options,'domain')
    optout.domain.xmin = [];
    optout.domain.xmax = [];
elseif isempty(options.domain)
    optout.domain.xmin = [];
    optout.domain.xmax = [];
else
    if ~isstruct(options.domain)
        error('OPTS.domain must be a structure with fields xmin and xmax');
    end
    if ~isfield(options.domain,'xmin') || ~isfield(options.domain,'xmax')
        error('OPTS.domain must be a structure with fields xmin and xmax');
    end
    
    % Check domain
    if numel(options.domain.xmin) ~= nx
        error(['OPTS.domain.xmin must be a vector with ',num2str(nx),' elements']);
    end
    if numel(options.domain.xmax) ~= nx
        error(['OPTS.domain.xmax must be a vector with ',num2str(nx),' elements']);
    end
    if any(options.domain.xmax < options.domain.xmin)
        error('OPTS.domain.xmax must be greater than OPTS.domain.xmin');
    end
    
    optout.domain.xmin = options.domain.xmin(:);
    optout.domain.xmax = options.domain.xmax(:);
end




