classdef ltiSys < dynSys
    % ltiSys   Linear time-invariant dynamical system object
    %
    % The ltiSys object can rapresent either an affine continuous-time
    % system in the form:
    %   _
    %  |
    %  |  dx
    %  | ---- = A x(t) + B u(t) + E_x p(t) + F_x d(t) + G_x
    % <   dt
    %  |
    %  | y(t) = C x(t) + D u(t) + E_y p(t) + F_y d(t) + G_y
    %  |_
    %
    % or a discrete-time affine system in the form:
    %   _
    %  |
    %  | x(k+1) = A x(k) + B u(k) + E_x p(k) + F_x d(k) + G_x
    % <
    %  |   y(k) = C x(k) + D u(k) + E_y p(k) + F_y d(k) + G_y
    %  |_
    %
    % where k denotes the discrete-time instant, with sampling time Ts.
    %
    % The quantities have the following meaning:
    %
    % - x: states (dimension nx)
    % - u: inputs (dimension nu)
    % - y: outputs (dimension ny)
    % - p: parameters (dimension np)
    % - d: unmeasurable inputs (dimension nd)
    %
    % A, B, C, D, E_x, E_y, F_x, F_y, G_x, G_y are time-invariant matrices 
    % of appropriate dimension.
    %
    % OBJ = ltiSys()
    % Builds an empty ltiSys object OBJ.
    %
    % OBJ = ltiSys(NX,NU,NY,'ct')
    % Builds a continuous-time ltiSys object by specifying the number of 
    % states (NX), inputs (NU) and outputs (NY). The number of parameters 
    % and unmeasurable inputs is assumed to be 0. String 'ct' indicates that
    % the system is continuous-time.
    %
    % OBJ = ltiSys(NX,NU,NY,'dt',Ts)
    % Builds a discrete-time ltiSys object by specifying the number of 
    % states (NX), inputs (NU) and outputs (NY). The number of parameters 
    % and unmeasurable inputs is assumed to be 0. String 'dt' indicates that
    % the system is discrete-time and Ts is the sampling time.
    %
    % OBJ = ltiSys(NX,NU,NY,NP,ND,'ct')
    % Builds a continuous-time ltiSys object by specifying the number of 
    % states (NX), inputs (NU), outputs (NY), parameters (NP) and 
    % unmeasurable inputs (ND). 
    %
    % OBJ = ltiSys(NX,NU,NY,NP,ND,'dt',Ts)
    % Builds a discrete-time ltiSys object by specifying the number of 
    % states (NX), inputs (NU), outputs (NY), parameters (NP) and 
    % unmeasurable inputs (ND). 
    %
    % ltiSys methods:
    %   computeOutput - computes the output of the system
    %   deltau - converts the system to delta-u formulation
    %   discretize - converts the system from continuous to discrete-time
    %   disp - displays some information about the ltiSys object
    %   evaluate - Computes the derivative of the state (if the system is continuous-time)
    %              or the state at next time instant (if the system is discrete-time)
    %   generateSimulinkModel - generates a Simulink model for the system simulation
    %   getClosedLoopPWAsystem - gets the matrices defining the LTI system controlled in closed   
    %   getMatrices - gets the A, B, C, D, W and f matrices defining the LTI system
    %   isObservable - checks if the system is observable
    %   setMatrices - sets the A, B, C, D, W and f matrices defining the LTI system
    %   sim - simulates the dynamical system.
    %   simplot - simulates the dynamical system and plots time evolution of
    %             states and inputs.
    %
    % The ltiSys object is derived from dynSys and inherits all its methods.
    %
    % See also dynSys, pwaSys
    
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
    
    
    
    
    % Properties
    
    properties (Access = private)
        A = [];
        B = [];
        C = [];
        D = [];
        Ex = [];
        Ey = [];
        Fx = [];
        Fy = [];
        Gx = [];
        Gy = [];
    end
    
    % Methods
    
    methods
        
        % Constructor
        function sys = ltiSys(varargin)
            
            if nargin == 0
                sys.A = [];
                sys.B = [];
                sys.C = [];
                sys.D = [];
                sys.Ex = [];
                sys.Ey = [];
                sys.Fx = [];
                sys.Fy = [];
                sys.Gx = [];
                sys.Gy = [];
                
            elseif nargin == 4
                sys.nx = varargin{1};
                sys.nu = varargin{2};
                sys.ny = varargin{3};
                
                if sys.nx <= 0
                    error('Number of states must be > 0')
                end
                if sys.nu < 0
                    error('Number of inputs must be >= 0')
                end
                if sys.ny <= 0
                    error('Number of outputs must be > 0')
                end
                
                mode = varargin{4};
                if ~strcmpi(mode,'ct')
                    error('Sampling time Ts must be provided!')
                end
                sys.mode = lower(mode);
                
                sys.A = zeros(sys.nx);
                sys.C = zeros(sys.ny,sys.nx);
                sys.B = zeros(sys.nx,sys.nu);
                sys.D = zeros(sys.ny,1);
                sys.Ex = zeros(sys.nx,sys.np);
                sys.Ey = zeros(sys.ny,sys.np);
                sys.Fx = zeros(sys.nx,sys.nd);
                sys.Fy = zeros(sys.ny,sys.nd);
                sys.Gx = zeros(sys.nx,1);
                sys.Gy = zeros(sys.ny,1);
                
                sys.Ts = 0;
                for i = 1:sys.nx
                    sys.xnames{i} = ['x_',num2str(i)];
                end
                for i = 1:sys.nu
                    sys.unames{i} = ['u_',num2str(i)];
                end
                for i = 1:sys.ny
                    sys.ynames{i} = ['y_',num2str(i)];
                end
                for i = 1:sys.np
                    sys.pnames{i} = ['p_',num2str(i)];
                end
                for i = 1:sys.nd
                    sys.dnames{i} = ['d_',num2str(i)];
                end
                
            elseif nargin == 5
                sys.nx = varargin{1};
                sys.nu = varargin{2};
                sys.ny = varargin{3};
                
                if sys.nx <= 0
                    error('Number of states must be > 0')
                end
                if sys.nu < 0
                    error('Number of inputs must be >= 0')
                end
                if sys.ny <= 0
                    error('Number of outputs must be > 0')
                end
                
                mode = varargin{4};
                if ~strcmpi(mode,'dt')
                    error('Sampling time Ts is not needed for continuous-time systems!')
                end
                sys.mode = lower(mode);
                
                sys.Ts = varargin{5};
                
                if sys.Ts <= 0
                    error('Sampling time must be > 0')
                end
                
                sys.A = zeros(sys.nx);
                sys.C = zeros(sys.ny,sys.nx);
                sys.B = zeros(sys.nx,sys.nu);
                sys.D = zeros(sys.ny,1);
                sys.Ex = zeros(sys.nx,sys.np);
                sys.Ey = zeros(sys.ny,sys.np);
                sys.Fx = zeros(sys.nx,sys.nd);
                sys.Fy = zeros(sys.ny,sys.nd);
                sys.Gx = zeros(sys.nx,1);
                sys.Gy = zeros(sys.ny,1);
                
                for i = 1:sys.nx
                    sys.xnames{i} = ['x_',num2str(i)];
                end
                for i = 1:sys.nu
                    sys.unames{i} = ['u_',num2str(i)];
                end
                for i = 1:sys.ny
                    sys.ynames{i} = ['y_',num2str(i)];
                end
                for i = 1:sys.np
                    sys.pnames{i} = ['p_',num2str(i)];
                end
                for i = 1:sys.nd
                    sys.dnames{i} = ['d_',num2str(i)];
                end
                
            elseif nargin == 6
                sys.nx = varargin{1};
                sys.nu = varargin{2};
                sys.ny = varargin{3};
                sys.np = varargin{4};
                sys.nd = varargin{5};
                
                if sys.nx <= 0
                    error('Number of states must be > 0')
                end
                if sys.nu < 0
                    error('Number of inputs must be >= 0')
                end
                if sys.ny <= 0
                    error('Number of outputs must be > 0')
                end
                if sys.np < 0
                    error('Number of parameters must be >= 0')
                end
                if sys.nd < 0
                    error('Number of unmeasurable inputs must be >= 0')
                end
                
                mode = varargin{6};
                if ~strcmpi(mode,'ct')
                    error('Sampling time Ts must be provided!')
                end
                sys.mode = lower(mode);
                
                sys.Ts = 0;
                
                sys.A = zeros(sys.nx);
                sys.C = zeros(sys.ny,sys.nx);
                sys.B = zeros(sys.nx,sys.nu);
                sys.D = zeros(sys.ny,1);
                sys.Ex = zeros(sys.nx,sys.np);
                sys.Ey = zeros(sys.ny,sys.np);
                sys.Fx = zeros(sys.nx,sys.nd);
                sys.Fy = zeros(sys.ny,sys.nd);
                sys.Gx = zeros(sys.nx,1);
                sys.Gy = zeros(sys.ny,1);
                
                for i = 1:sys.nx
                    sys.xnames{i} = ['x_',num2str(i)];
                end
                for i = 1:sys.nu
                    sys.unames{i} = ['u_',num2str(i)];
                end
                for i = 1:sys.ny
                    sys.ynames{i} = ['y_',num2str(i)];
                end
                for i = 1:sys.np
                    sys.pnames{i} = ['p_',num2str(i)];
                end
                for i = 1:sys.nd
                    sys.dnames{i} = ['d_',num2str(i)];
                end
                
            elseif nargin == 7
                sys.nx = varargin{1};
                sys.nu = varargin{2};
                sys.ny = varargin{3};
                sys.np = varargin{4};
                sys.nd = varargin{5};
                
                if sys.nx <= 0
                    error('Number of states must be > 0')
                end
                if sys.nu < 0
                    error('Number of inputs must be >= 0')
                end
                if sys.ny <= 0
                    error('Number of outputs must be > 0')
                end
                if sys.np < 0
                    error('Number of parameters must be >= 0')
                end
                if sys.nd < 0
                    error('Number of unmeasurable inputs must be >= 0')
                end
                
                mode = varargin{6};
                if ~strcmpi(mode,'dt')
                    error('Sampling time Ts is not needed for continuous-time systems!')
                end
                sys.mode = lower(mode);
                
                sys.Ts = varargin{7};
                
                if sys.Ts <= 0
                    error('Sampling time must be > 0')
                end
                
                sys.A = zeros(sys.nx);
                sys.C = zeros(sys.ny,sys.nx);
                sys.B = zeros(sys.nx,sys.nu);
                sys.D = zeros(sys.ny,1);
                sys.Ex = zeros(sys.nx,sys.np);
                sys.Ey = zeros(sys.ny,sys.np);
                sys.Fx = zeros(sys.nx,sys.nd);
                sys.Fy = zeros(sys.ny,sys.nd);
                sys.Gx = zeros(sys.nx,1);
                sys.Gy = zeros(sys.ny,1);
                
                for i = 1:sys.nx
                    sys.xnames{i} = ['x_',num2str(i)];
                end
                for i = 1:sys.nu
                    sys.unames{i} = ['u_',num2str(i)];
                end
                for i = 1:sys.ny
                    sys.ynames{i} = ['y_',num2str(i)];
                end
                for i = 1:sys.np
                    sys.pnames{i} = ['p_',num2str(i)];
                end
                for i = 1:sys.nd
                    sys.dnames{i} = ['d_',num2str(i)];
                end
                
            else
                error('Wrong input arguments for ltiSys object constructor');
            end
        end
        
        % Set methods
        
        object = setMatrices(object,varargin);
        
        % Get methods
        
        varargout = getMatrices(object,varargin);
        
        % Is methods
        answ = isObservable(object);
        
        % Other methods
        
        signals = sim(object,mode,x0,p,d,N,ctrl);
        
        signals = simplot(object,mode,x0,p,d,N,ctrl);
        
        duSys = deltau(object);
        
        y = computeOutput(object,x,u,p,d);
        
        xnext = evaluate(object,x,u,p,d);
        
        dynSys = discretize(object,Ts);
        
        disp(object);
        
        
    end
    
    methods (Access = private)
        
        signals = ct_OpenLoopSim(object,x0,p,d,time,options);
        signals = dt_OpenLoopSim(object,x0,p,d,time,options);
        signals = ct_OpenLoopObsSim(object,x0,p,d,time,options);
        signals = dt_OpenLoopObsSim(object,x0,p,d,time,options);
        signals = ct_ClosedLoopSim(object,x0,p,d,xref,time,options);
        signals = dt_ClosedLoopSim(object,x0,p,d,xref,time,options);
        signals = ct_ClosedLoopObsSim(object,x0,p,d,xref,time,options);
        signals = dt_ClosedLoopObsSim(object,x0,p,d,xref,time,options);
        
        ct_simulinkOpenLoopSim(object,options);
        dt_simulinkOpenLoopSim(object,options);
        ct_simulinkOpenLoopObsSim(object,options);
        dt_simulinkOpenLoopObsSim(object,options);
        ct_simulinkClosedLoopSim(object,options);
        dt_simulinkClosedLoopSim(object,options);
        ct_simulinkClosedLoopObsSim(object,options);
        dt_simulinkClosedLoopObsSim(object,options);
        
    end
    
end


