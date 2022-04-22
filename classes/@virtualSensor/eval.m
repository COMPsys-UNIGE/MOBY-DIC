%% EVAL
% Evaluates the unmeasurable output z in correspondence of given system
% inputs and measurable outputs
%
% SYNTAX
%
% zh = eval(object,y,[z0])
%
% Only the measures of system outputs y are provided (if the system does
% not have any input). z0 (optional) is the initial condition for the
% estimate z. If it is not provided it is set to 0.
%
% zh = eval(object,u,y,[z0])
%
% The measures of system inputs u are also provided.
% 
% ACKNOWLEDGEMENTS
%
% Contributors:
% 
% * Alberto Oliveri (alberto.oliveri@unige.it)
% 
% Copyright is with:
% 
% * Copyright (C) 2011 University of Genoa, Italy.

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

function z = eval(object,varargin)

% Check on inputs
if nargin == 2
    y = varargin{1};
    u = [];
    z0 = 0;
elseif nargin == 3
    if isscalar(varargin{2})
        y = varargin{1};
        z0 = varargin{2};
        u = [];
    else
        u = varargin{1};
        y = varargin{2};
        z0 = 0;
    end
elseif nargin == 4
    u = varargin{1};
    y = varargin{2};
    z0 = varargin{3};
else
    error('Wrong input arguments.');
end

if size(u,2) ~= object.nu
    error(['Input u must be a matrix with ',num2str(object.nu),' columns']);
end

if size(y,2) ~= object.ny
    error(['Inputs y must be a matrix with ',num2str(object.ny),' columns']);
end

% Number of samples of system measurable output
T = size(y,1);

if ~isempty(u)
    if size(u,1) ~= T
        error('Inputs u and y must have the same number of rows');
    end
end

if ~object.isIdentified()
    error('Virtual sensor is not identified.');
end

% Time windows
mu = object.mu;
my = object.my;
mz = object.mz;

% Maximum time window
l = max([max(mz)+1,max(mu),max(my)]);

% Number of dimensions of the space
Ndim = mz+sum(my)+sum(mu);

% pwas object
fpwas = object.fpwas;


if ~isempty(u)
    
    u = [nan(l-1,object.nu); u];
    
end

y = [nan(l-1,object.ny); y];

% This is necessary for reduced complexity virtual sensors
% The i-th element of cell array indices will contain the indices of the
% columns of X corresponding to time instant k-i+1 (k-i if current time
% instant is not considered)
indices = cell(1,l);

% Initialize indices with empty elements
for i = 1:l
    indices{i} = [];
end

% X will contain all points u and y at the proper time instants used to
% estimate z
X = zeros(size(y,1)-l+1,Ndim-mz);

% Indices of the starting column
kstart = 1;

% Loop on the number of system inputs
for j = 1:object.nu
    
    % Index used to indicize cell array indices
    h = 1;
    
    % Loop on the number of past time instants to consider for the
    % current input
    for k = kstart:kstart+mu(j)-1
        
        kk = k-kstart+1;
        % Fill matrix X with input data in the current time instant
        % This is u_j at instant t-kk+1
        X(:,k) = u(l-kk+1:end-kk+1,j);
        
        % Update indices
        indices{h} = [indices{h} k];
        h = h+1;
        
    end
    
    if isempty(k)
        k = kstart-1;
    end
    
    % Update index of the starting column
    kstart = k+1;
    
end

% Loop on the number of system outputs
for j = 1:object.ny
    
    % Index used to indicize cell array indices
    h = 1;
    
    % Loop on the number of past time instants to consider for the
    % current output
    for k = kstart:kstart+my(j)-1
        
        kk = k-kstart+1;
        % Fill matrix X with output data in the current time instant
        % This is y_j at instant t-kk+1
        X(:,k) = y(l-kk+1:end-kk+1,j);
        
        % Update indices
        indices{h} = [indices{h} k];
        h = h+1;
        
    end
    
    if isempty(k)
        k = kstart-1;
    end
    
    % Update index of the starting column
    kstart = k+1;
    
end

ndata = size(X,1);

% If a reduced complexity virtual sensor is used, you need to split the
% dataset into cell arrays, one for each time instant. The indices of the
% columns of X corresponding to a given time instant are contained in
% 'indices'.
if object.reducedComplexity
    
    % XX will contain the dataset
    XX = cell(1,l);
    
    for i = 1:l
        
        XX{i} = X(:,indices{i});
        
    end
    
    X = XX;
    clear XX
    
end

% Buffer containing the values of the estimated output for all needed time
% instants
zbuffer = z0*ones(1,object.mz);

% z will contain the estimated output
z = zeros(ndata,1);

inan = 0;

% Estimate z
for i = 1:ndata
    
    if ~object.reducedComplexity
        
        % Estimate z at next time step
        if any(isnan(X(i,:)))
            z(i) = z0;
        else
            z(i) = fpwas.eval([X(i,:) zbuffer]);
        end
        
    else
        
        zacc = 0;
        
        for j = 1:numel(fpwas)
            
            % Take the value at previous step from zbuffer
            if object.current
                if j-1 <= numel(zbuffer) && j ~= 1
                    zz = zbuffer(j-1);
                else
                    zz = [];
                end
            else
                if j <= numel(zbuffer)
                    zz = zbuffer(j);
                else
                    zz = [];
                end
            end
            
            if any(isnan(X{j}(i,:)))
                zacc = z0;
                break
            end
            
            % Estimate z at next time step
            zacc = zacc+fpwas(j).eval([X{j}(i,:) zz]');
            
            if any(isnan(zacc))
                inan = 1;
                warning('Virtual sensor returned NaN. Validation did not complete.');
                break;
            end
            
        end
        
        if inan
            break;
        end
        
        % Accumulate z
        z(i) = zacc;
        
    end
    
    if inan
        break;
    end
    
    % Shift left zbuffer
    zbuffer = [z(i) zbuffer];
    zbuffer(end) = [];
       
end

end
