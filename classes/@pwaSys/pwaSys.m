classdef pwaSys < dynSys
    % Piecewise affine discrete-time dynamical system object.
    %
    % The PWA object can rapresent either a picewise affine continuous-time
    % pwa system in the form:
    %   _
    %  |
    %  |  dx
    %  | ---- = A^{(i)} x(t) + B^{(i)} u(t) + E_x^{(i)} p(t) + F_x^{(i)} d(t) + G_x^{(i)}
    % <   dt
    %  |                                                                                  ,if H^{(i)} [x(t) p(t) d(t)]' <= K^{(i)}
    %  | y(t) = C^{(i)} x(t) + D^{(i)} u(t) + E_y^{(i)} p(t) + F_y^{(i)} d(t) + G_y^{(i)}
    %  |_
    %
    % or a discrete-time affine system in the form:
    %   _
    %  |
    %  | x(k+1) = A^{(i)} x(k) + B^{(i)} u(k) + E_x^{(i)} p(k) + F_x^{(i)} d(k) + G_x^{(i)}
    % <
    %  |                                                                                  ,if H^{(i)} x(t) <= K^{(i)}
    %  |   y(k) = C^{(i)} x(k) + D^{(i)} u(k) + E_y^{(i)} p(k) + F_y^{(i)} d(k) + G_y^{(i)}
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
    % A^{(i)}, B^{(i)}, C^{(i)}, D^{(i)}, E_x^{(i)}, E_y^{(i)}, F_x^{(i)},
    % F_y^{(i)}, G_x^{(i)}, G_y^{(i)} are time-invariant matrices
    % of appropriate dimension.
    %
    % OBJ = pwaSys()
    % Builds an empty ltiSys object OBJ.
    %
    % OBJ = pwaSys(NX,NU,NY,'ct')
    % Builds a continuous-time ltiSys object by specifying the number of
    % states (NX), inputs (NU) and outputs (NY). The number of parameters
    % and unmeasurable inputs is assumed to be 0. String 'ct' indicates that
    % the system is continuous-time.
    %
    % OBJ = pwaSys(NX,NU,NY,'dt',Ts)
    % Builds a discrete-time ltiSys object by specifying the number of
    % states (NX), inputs (NU) and outputs (NY). The number of parameters
    % and unmeasurable inputs is assumed to be 0. String 'dt' indicates that
    % the system is discrete-time and Ts is the sampling time.
    %
    % OBJ = pwaSys(NX,NU,NY,NP,ND,'ct')
    % Builds a continuous-time ltiSys object by specifying the number of
    % states (NX), inputs (NU), outputs (NY), parameters (NP) and
    % unmeasurable inputs (ND).
    %
    % OBJ = pwaSys(NX,NU,NY,NP,ND,'dt',Ts)
    % Builds a discrete-time ltiSys object by specifying the number of
    % states (NX), inputs (NU), outputs (NY), parameters (NP) and
    % unmeasurable inputs (ND).
    %
    % pwaSys methods:
    %   computeOutput - computes the output of the system
    %   deltau - converts the system to delta-u formulation
    %   discretize - converts the system from continuous to discrete-time
    %   disp - displays some information about the ltiSys object
    %   evaluate - Computes the derivative of the state (if the system is continuous-time)
    %              or the state at next time instant (if the system is discrete-time)
    %   generateSimulinkModel - generates a Simulink model for the system simulation
    %   getMatrices - gets the A, B, C, D, W and f matrices defining the LTI system
    %   isObservable - checks if the system is observable
    %   setMatrices - sets the A, B, C, D, W and f matrices defining the LTI system
    %   sim - simulates the dynamical system.
    %   simplot - simulates the dynamical system and plots time evolution of
    %             states and inputs.
    %
    % The pwaSys object is derived from dynSys and inherits all its methods.
    %
    % See also dynSys, ltiSys
    
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
    
    
    
    % Properties
    
    properties (Access = private)
        
        nDyn = 0;
        
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
        
        H = [];
        K = [];
    end
    
    % Methods
    
    methods
        
        % Constructor
        function sys = pwaSys(varargin)
            
            if nargin == 0
                sys.nDyn = 0;
                
                sys.A = cell(0);
                sys.B = cell(0);
                sys.C = cell(0);
                sys.D = cell(0);
                sys.Ex = cell(0);
                sys.Ey = cell(0);
                sys.Fx = cell(0);
                sys.Fy = cell(0);
                sys.Gx = cell(0);
                sys.Gy = cell(0);
                
            elseif nargin == 4
                sys.nx = varargin{1};
                sys.nu = varargin{2};
                sys.ny = varargin{3};
                sys.nDyn = 0;
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
                
                
                
                sys.A = cell(0);
                sys.C = cell(0);
                sys.B = cell(0);
                sys.D = cell(0);
                sys.Ex = cell(0);
                sys.Ey = cell(0);
                sys.Fx = cell(0);
                sys.Fy = cell(0);
                sys.Gx = cell(0);
                sys.Gy = cell(0);
                
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
                sys.nDyn = 0;
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
                
                sys.A = cell(0);
                sys.C = cell(0);
                sys.B = cell(0);
                sys.D = cell(0);
                sys.Ex = cell(0);
                sys.Ey = cell(0);
                sys.Fx = cell(0);
                sys.Fy = cell(0);
                sys.Gx = cell(0);
                sys.Gy = cell(0);
                
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
                sys.nDyn = 0;
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
                
                sys.A = cell(0);
                sys.C = cell(0);
                sys.B = cell(0);
                sys.D = cell(0);
                sys.Ex = cell(0);
                sys.Ey = cell(0);
                sys.Fx = cell(0);
                sys.Fy = cell(0);
                sys.Gx = cell(0);
                sys.Gy = cell(0);
                
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
                sys.nDyn = 0;
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
                
                sys.A = cell(0);
                sys.C = cell(0);
                sys.B = cell(0);
                sys.D = cell(0);
                sys.Ex = cell(0);
                sys.Ey = cell(0);
                sys.Fx = cell(0);
                sys.Fy = cell(0);
                sys.Gx = cell(0);
                sys.Gy = cell(0);
                
                sys.H = cell(0);
                sys.K = cell(0);
                
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
        
        object = addDynamics(object,H,K);
        
        % Set methods
        
        object = setMatrices(object,dyn,varargin);
        
        % Get methods
        
        nDyn = getNumberOgDynamics(object);
        
        varargout = getMatrices(object,varargin);
               
        % Other methods
        
        signals = sim(object,mode,x0,p,d,N,ctrl);
        
        signals = simplot(object,mode,x0,p,d,N,ctrl);
        
        duSys = deltau(object); 
        
        y = computeOutput(object,x,u,p,d);
        
        xnext = evaluate(object,x,u,p,d);
        
        dynSys = discretize(object,Ts);
        
        dyn = findDynamics(object,x,p,d);
        
        disp(object);
        
        plotDynamicPartition(object);
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


