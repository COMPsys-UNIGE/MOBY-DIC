function signals = sim(object,T,x0,varargin)
% sim    Simulates the dynamical system
%
% If the system has not a controller, an open-loop simulation is performed,
% otherwise the control action of the controller is also applied.
% If the system has an observer, the observed state variable is fed into
% the controller and both real and observed variables are monitored.
%
% SIGNALS = sim(OBJ,T,X0)
% Simulates the LTI system OBJ, starting from initial condition X0 for a 
% time T. Output SIGNALS is a structure with the following fields:
% - time: time instants in which variables are computed
% - state: system states
% - input: system inputs provided by the controller associated to the 
%          system (if present)
% - output: system outputs
% - est_state: system states estimated by the observer associated to the
%              system (if present)
% - est_output: system outputs estimated by the observer associated to the
%               system  (if present)
% - unmeasurable_input: system unmeasurable inputs (if present)
% - parameter: system parameters (if present)
% - ref: reference state or output where the system is regulated (if
%        present)
%
% SIGNALS = sim(OBJ,T,X0,OPTS)
% If the system is continuous-time, a structure OPTS can be provided with
% the following fields:
% - ODEsolver: string indicating the solver used to integrate the
%              differential equations. It can be either 'ode45','ode23',
%              'ode113','ode15s','ode23s','ode23t' or 'ode23tb'. See the
%              documentation of ode45 to get help (doc ode45).
%              Default: 'ode45'.
% - ODEoptions: options to provide to the ODE solver. Type doc odeset to
%               get help. Default: odeset.
%
% SIGNALS = sim(OBJ,T,X0,REF)
% SIGNALS = sim(OBJ,T,X0,REF,OPTS)
% If a tracking controller is associated to the system, the reference state
% or output REF must be provided.
%
% SIGNALS = sim(OBJ,T,X0,P,D)
% SIGNALS = sim(OBJ,T,X0,P,D,OPTS)
% SIGNALS = sim(OBJ,T,X0,P,D,REF)
% SIGNALS = sim(OBJ,T,X0,P,D,REF,OPTS)
% If the system has parameters and/or unmeasurable inputs the must be 
% provided through vectors P and D.
%
% See also: ltiSys/simplot.

% Contributors:
%
% Alberto Oliveri (alberto.oliveri@unige.it)
% Matteo Lodi (matteo.lodi@edu.unige.it)
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


if nargin == 3
    p = [];
    d = [];
    ref = [];
elseif nargin == 4
    if isstruct(varargin{1})
        p = [];
        d = [];
        ref = [];
        options = varargin{1};
    else
        p = [];
        d = [];
        ref = varargin{1};
    end
elseif nargin == 5
    if isstruct(varargin{2})
        p = [];
        d = [];
        ref = varargin{1};
        options = varargin{2};
    else
        p = varargin{1};
        d = varargin{2};
        ref = [];
    end
elseif nargin == 6
    if isstruct(varargin{3})
        p = varargin{1};
        d = varargin{2};
        ref = [];
        options = varargin{3};
    else
        p = varargin{1};
        d = varargin{2};
        ref = varargin{3};
    end
elseif nargin == 7
    p = varargin{1};
    d = varargin{2};
    ref = varargin{3};
    options = varargin{4};
else
    error('Wrong input arguments');
end

% Check options structure
if ~exist('options','var')
    options.ODEoptions = odeset;
    options.ODEsolver = 'ode45';
end

if ~isstruct(options)
    error('options must be a structure');
end

if object.isContinuousTime
    if ~isfield(options,'ODEoptions')
        options.ODEoptions = odeset;
    end
    if isempty(options.ODEoptions)
        options.ODEoptions = odeset;
    end
    if ~isstruct(options.ODEoptions)
        error('ODEoptions must be a structure');
    end
    
    if ~isfield(options,'ODEsolver')
        options.ODEsolver = 'ode45';
    end
    if isempty(options.ODEsolver)
        options.ODEsolver = 'ode45';
    end
    availsolvers = {'ode45','ode23','ode113','ode15s','ode23s','ode23t','ode23tb'};
    if ~any(ismember(availsolvers,options.ODEsolver))
        error('The specified ODEsolver is not available. Type ''doc ode45'' to get help');
    end
end

% Check initial condition
x0 = x0(:);

if numel(x0) ~= object.nx
    error(['X0 must be a vector with ',num2str(object.nx),' elements']);
end

% Check parameters
if object.np == 0
    if ~isempty(p)
        disp('WARNING: system has no parameters. Input ''P'' is ignored');
    end
    p = zeros(object.np,1);
else
    
    p = p(:);
    
    if numel(p) ~= object.np
        error(['p must be a vector with ',num2str(object.np),' elements']);
    end
end

% Check unmeasurable inputs
if object.nd == 0
    if ~isempty(d)
        disp('WARNING: system has no unmeasurable inputs. Input ''D'' is ignored');
    end
    d = zeros(object.nd,1);
else
    
    d = d(:);
    
    if numel(d) ~= object.nd
        error(['d must be a vector with ',num2str(object.nd),' elements']);
    end
end

disp('Simulating...')
disp(' ')

% Simulation of continuous-time system
if object.isContinuousTime
    
    % Check time
    if T <= 0
        error('T must be a positive number');
    end
    
    % Open-loop system
    if ~object.hasController()
        
        if ~object.hasObserver()
            
            % Without observer
            signals = ct_OpenLoopSim(object,x0,p,d,T,options);
            
        else
            
            % With observer
            signals = ct_OpenLoopObsSim(object,x0,p,d,T,options);
        end
        
        % Closed-loop system
    else
        
        % Without observer
        if ~object.hasObserver()
            
            signals = ct_ClosedLoopSim(object,x0,p,d,ref,T,options);
            
            % With observer
        else
            
            signals = ct_ClosedLoopObsSim(object,x0,p,d,ref,T,options);
            
        end
        
    end
    
    % Simulation of discrete-time system
else
    
    % Check time
    if T <= 0
        error('T must be a positive number');
    end
    
    % Open-loop system
    if ~object.hasController()
        
        if ~object.hasObserver()
            
            % Without observer
            signals = dt_OpenLoopSim(object,x0,p,d,T,options);
            
        else
            
            % With observer
            signals = dt_OpenLoopObsSim(object,x0,p,d,T,options);
        end
        
        % Closed-loop system
    else
        
        % Without observer
        if ~object.hasObserver()
            
            signals = dt_ClosedLoopSim(object,x0,p,d,ref,T,options);
            
            % With observer
        else
            
            signals = dt_ClosedLoopObsSim(object,x0,p,d,ref,T,options);
            
        end
        
    end
    
end

% Add parameter, unmeasurable input and reference state to signals
signals.unmeasurable_input = d;
signals.parameter = p;
signals.ref = ref;

disp('Done.')
disp(' ')


end







