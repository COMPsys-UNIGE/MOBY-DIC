%% VALIDATE
% Validates the virtual sensor through a test set of data
%
% SYNTAX
%
% [zh err] = validate(object,yt,zt,[options])
%
% yt and zt are matrices (or cell arrays of matrices) defining the test set. 
% options is a structure with the following fields:
%
% * plot: if it is 1, a plot of z, zt and the relative and absolute error 
%          is performed. Defalut value: 1.
% * initialCondition: initial condition for the estimate of z. This
%                      condition is used until all the data at needed time
%                      instants are available. This field can be an array,
%                      indicating a different initial condition for all the
%                      elements of the cell arrays zt, or a scalar, if the
%                      same initial condition must be used for all elements
%                      of the cell arrays. Default initial condition is 0.
%
% zh is the estimated output and err is a struct with fields MAX and RMSE
% containing the maximum and root mean square error respectively.
%
% [zh err] = validate(object,ut,yt,zt,[options])
%
% ACKNOWLEDGEMENTS
%
% Contributors:
% 
% * Tomaso Poggi (tpoggi@essbilbao.org)
% 
% Copyright is with:
% 
% * Copyright (C) 2010 University of Genoa, Italy.

% -------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% -------------------------------------------------------------------------

function [zh err] = validate(object,varargin)

% Check on inputs
if nargin == 3
    u = [];
    y = varargin{1};
    z = varargin{2};
    options.plot = 1;
elseif nargin == 4
    if isstruct(varargin{3})
        u = [];
        y = varargin{1};
        z = varargin{2};
        options = varargin{3};
    else
        u = varargin{1};
        y = varargin{2};
        z = varargin{3};
        options.plot = 1;
    end
elseif nargin == 5
    u = varargin{1};
    y = varargin{2};
    z = varargin{3};
    options = varargin{4};
else
    error('Wrong input arguments.');
end

if ~object.isIdentified()
    error('Virtual sensor is not identified.');
end

if ~isfield(options,'plot')
    options.plot = 1;
end

if ~iscell(y)
    yy = cell(1,1);
    yy{1} = y;
    y = yy;
    
    uu = cell(1,1);
    uu{1} = u;
    u = uu;
    
    zz = cell(1,1);
    zz{1} = z;
    z = zz;
    
else
    
    if numel(y) ~= numel(u) && ~isempty(u)
        error('Cell arrays u and y must have the same number of elements');
    end
    
    if isempty(u)
        uu = cell(size(y));
        for i = 1:numel(uu)
            uu{i} = [];
        end
        u = uu;
    end
    
    if numel(y) ~= numel(z)
        error('Cell arrays y and z must have the same number of elements');
    end
end

if ~isfield(options,'initialCondition')
    options.initialCondition = zeros(1,numel(y));
end

if ~isnumeric(options.initialCondition)
    error('Initial condition must be an array of doubles');
end
if isscalar(options.initialCondition)
    options.initialCondition = repmat(options.initialCondition,1,numel(y));
else
    if numel(options.initialCondition) ~= numel(y)
        error(['Initial conditions must be an array of ', num2str(numel(y)),' elements']);
    end
end

z0 = options.initialCondition;

zh = cell(1,numel(y));

for ii = 1:numel(y)
    
    if ~isempty(u{ii})
        if size(u{ii},2) ~= object.nu
            error(['Input u must be a matrix with ',num2str(object.nu),' columns']);
        end
    end
    if size(y{ii},2) ~= object.ny
        error(['Input y must be a matrix with ',num2str(object.ny),' columns']);
    end
    
    % Number of samples of system measurable output
    T = size(y{ii},1);
    
    if ~isempty(u{ii})
        if size(u{ii},1) ~= T
            error('Inputs u and y must have the same number of rows');
        end
    end
    
    if size(z{ii},1) ~= T
        error('Inputs y and z must have the same number of rows');
    end
    
    % Estimate unmeasurable output
    zh{ii} = object.eval(u{ii},y{ii},z0(ii));
    
    if ~object.current
        zh{ii}(end) = [];
        z{ii}(1) = [];
    end
    
    % Compute errors. The first T/3 elements are not considered because, for
    % these data the data at previous time instants are not available and so
    % they are set to 0. For this reason the error could be altered.
    tau = round(T/3);
    err(ii).RMSE = sqrt(mean((zh{ii}(tau:end)-z{ii}(tau:end)).^2));
    err(ii).MAX = max(zh{ii}(tau:end)-z{ii}(tau:end));
    
    if options.plot
        % plot estimated output vs real output
        figure
        subplot(3,1,1)
        hold on
        plot(z{ii})
        plot(zh{ii},'r')
        xlabel('t')
        ylabel('z')
        legend('test set','estimation')
        % plot absolute error
        subplot(3,1,2)
        hold on
        plot(abs(z{ii}-zh{ii}))
        xlabel('t')
        ylabel('absolute error')
        % plot relative error
        subplot(3,1,3)
        hold on
        plot(abs(z{ii}-zh{ii})./max(abs(z{ii}))*100)
        xlabel('t')
        ylabel('% relative error')
    end
    
end
