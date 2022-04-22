function optout = VSset(vs,options)
% VSset   Checks and sets the options for virtual sensor identification
%
% OPTOUT = VSset(OPTS)
% OPTS is the structure defining the options.
% See the documentation of virtualSensor.identify for a list of all
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

if ~isfield(options,'lambda')
    optout.lambda = 1e-6;
elseif isempty(options.lambda)
    optout.lambda = 1e-6;
else
    optout.lambda = options.lambda;
end

if ~isfield(options,'gamma')
    optout.gamma = 1.1;
elseif isempty(options.gamma)
    optout.gamma = 1.1;
else
    optout.gamma = options.gamma;
end

if ~isfield(options,'order')
    optout.order = 0;
elseif isempty(options.order)
    optout.order = 0;
else
    optout.order = options.order;
end

if ~isfield(options,'np')
    error('OPTS.np must be provided');
elseif isempty(options.order)
    error('OPTS.np must be provided');
else
    optout.np = options.np;
end

if ~isfield(options,'Ts')
    error('OPTS.Ts must be provided');
elseif isempty(options.order)
    error('OPTS.Ts must be provided');
else
    optout.Ts = options.Ts;
end

nu = vs.getNumberOfInputs();
ny = vs.getNumberOfMeasurableOutputs();
mz = vs.getAutoregressiveTimeWindow();

if ~isfield(options,'domain')
    optout.domain = [];
elseif isempty(options.domain)
    optout.domain = [];
else
    if ~isstruct(options.domain)
        error('OPTS.domain must be a struct');
    end
    if nu > 0
        if ~isfield(options.domain,'umin') || ~isfield(options.domain,'umax')
            error('OPTS.domain must be a struct with fields umin and umax');
        end
        umin = options.domain.umin;
        umax = options.domain.umax;
        if numel(umin) == 1
            umin = repmat(umin,nu,1);
        end
        if numel(umax) == 1
            umax = repmat(umax,nu,1);
        end
        if numel(umin) ~= nu || numel(umax) ~= nu
            error(['OPTS.domain.umin and OPTS.domain.umax must be either scalars or arrays with ',num2str(nu),' elements']);
        end
        optout.domain.umin = umin;
        optout.domain.umax = umax;
    end
    
    if ny > 0
        if ~isfield(options.domain,'ymin') || ~isfield(options.domain,'ymax')
            error('OPTS.domain must be a struct with fields ymin and ymax');
        end
        ymin = options.domain.ymin;
        ymax = options.domain.ymax;
        if numel(ymin) == 1
            ymin = repmat(ymin,ny,1);
        end
        if numel(ymax) == 1
            ymax = repmat(ymax,ny,1);
        end
        if numel(ymin) ~= ny || numel(ymax) ~= ny
            error(['OPTS.domain.ymin and OPTS.domain.ymax must be either scalars or arrays with ',num2str(ny),' elements']);
        end
        optout.domain.ymin = ymin;
        optout.domain.ymax = ymax;
    end
    
    if mz > 0
        if ~isfield(options.domain,'zmin') || ~isfield(options.domain,'zmax')
            error('OPTS.domain must be a struct with fields zmin and zmax');
        end
        zmin = options.domain.zmin;
        zmax = options.domain.zmax;
        if numel(zmin) == 1
            zmin = repmat(zmin,nz,1);
        end
        if numel(zmax) == 1
            zmax = repmat(zmax,nz,1);
        end
        if numel(zmin) ~= nz || numel(zmax) ~= nz
            error(['OPTS.domain.zmin and OPTS.domain.zmax must be either scalars or arrays with ',num2str(nz),' elements']);
        end
        optout.domain.zmin = zmin;
        optout.domain.zmax = zmax;
    end
    
end



