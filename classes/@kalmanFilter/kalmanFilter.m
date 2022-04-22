classdef kalmanFilter < observer
    % kalmanFilter   Kalman filter for the estimation of the states 
    %                of a LTI system
    %
    % This object represents a Kalman filter for the estimation of
    % unmeasurable states of a LTI system (ltiSys object). The knowledge of
    % all system states is necessary for the design of the MPC controller.
    % The Kalman filter is a dynamical system in the form:
    %
    % The Kalman filter is a discrete-time dynamical system, whose state
    % evolution follows the evolution of the real system state.
    % Refer to MOBY-DIC Toolbox User's Guide for a detailed explanation of
    % the Kalman filter.
    %
    % OBJ = kalmanFilter()
    % Builds an empty kalmanFilter object OBJ.
    %
    % OBJ = kalmanFilter(SYS,Ts,K)
    % Builds a kalmanFilter object OBJ by specifying the ltiSys object
    % (SYS) whose states need to be estimated, the observer sampling time
    % Ts and the Kalman gain K. The gain K weights the difference between
    % the measured system output and the estimated one. K should be the
    % solution of the Riccati equation.
    %
    % OBJ = kalmanFilter(SYS,Ts,Q,R)
    % Builds a kalmanFilter object OBJ by specifying the ltiSys object
    % (SYS) whose states need to be estimated, the observer sampling time
    % Ts and the covariance matrices Q and R of the model and measurement
    % noises, respectively. If the elements in Q are higher than the 
    % elements in R, the measurements are trusted more than the model and
    % viceversa. Q must be a square matrix with dimension NX+ND, being NX
    % the number of system states and ND the number of unmeasurable inputs.
    % R must be a matrix of dimension NY, being NY the number of system
    % outputs. The Kalman gain K is computed through MATLAB function
    % 'kalman'.
    %
    % kalmanFilter methods:
    %   computeOutput - Computes the Kalman filter output
    %   disp - displays information about kalmanFilter object
    %   getCovarianceMatrix - gets the asymptotic covariance matrix
    %   getGain - gets the gain of the Kalman filter
    %   getMatrices - gets the matrices defining the Kalman filter
    %   predict - predicts the system state at next time instant based on current
    %             measurements
    %   update - updates the system state based on current outputs
    %
    % The kalmanFilter object is derived from observer and inherits
    % all its methods.
    %
    % See also kalmanPredictor.
    
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
        A = [];
        B = [];
        C = [];
        D = [];
        Gx = [];
        Gy = [];
        
        K = [];
        Q = [];
        R = [];
        
        P = [];
        
    end
    
    methods
        
        function obs = kalmanFilter(varargin)
            
            if nargin == 0
                
                obs.A = [];
                obs.B = [];
                obs.C = [];
                obs.D = [];
                obs.Gx = [];
                obs.Gy = [];
                
                obs.K = [];
                obs.Q = [];
                obs.R = [];
                obs.P = [];
                
            elseif nargin == 3
                
                sys = varargin{1};
                obs.Ts = varargin{2};
                obs.K = varargin{3};
                
                % Check inputs
                if ~isa(sys,'ltiSys')
                    error('First argument must be LTIsys object');
                end
                
                if obs.Ts <= 0
                    error('Sampling time must be > than 0');
                end
                
                if LTIsys.isDiscreteTime() && LTIsys.getSamplingTime ~= obs.Ts
                    error('Sampling time Ts is different from the system sampling time');
                end
                
                % Set number of variables
                obs.nx = sys.getNumberOfStates();
                obs.nu = sys.getNumberOfInputs();
                obs.ny = sys.getNumberOfOutputs();
                obs.np = sys.getNumberOfParameters();
                obs.nd = sys.getNumberOfUnmeasurableInputs();
                
                % Check gain dimensions
                if any(size(obs.K) ~= [ obs.nx+obs.nd obs.ny ])
                    error(['Gain K must be a matrix with ',num2str(obs.nx+obs.nd),' rows and ',num2str(obs.ny),' columns']);
                end
                
                % If the system is continuous-time, discretize it
                if sys.isContinuousTime 
                    disp(' ')
                    disp('Discretizing system...')
                    sys = sys.discretize(obs.Ts);
                end
                
                % Extract system matrices
                [Asys, Bsys, Csys, Dsys, Exsys, Eysys, Fxsys, Fysys, Gxsys, Gysys] = sys.getMatrices();
                
                % Put together state and unmeasurable input this way:
                %
                % X = [ x 
                %       d ];
                %
                % and input and parameters this way:
                %
                % U = [ u 
                %       p ];
                % Put together state and unmeasurable input                
                Ap = [Asys Fxsys; zeros(obs.nd,obs.nx) eye(obs.nd)];
                Bp = [Bsys Exsys; zeros(obs.nd,obs.nu) zeros(obs.nd,obs.np)];
                Gxp = [Gxsys; zeros(obs.nd,1)];
                Cp = [Csys Fysys];
                Dp = [Dsys Eysys];
                Gyp = Gysys;
                
               
                % Set observer matrices
                obs.A = Ap-Ap*obs.K*Cp;
                obs.B = [Bp-Ap*obs.K*Dp Ap*obs.K];
                obs.Gx = Gxp-Ap*obs.K*Gyp;
                obs.C = [Cp-Cp*obs.K*Cp;eye(obs.nx+obs.nd)-obs.K*Cp];
                obs.D = [[Dp-Cp*obs.K*Dp Cp*obs.K];-obs.K*Dp obs.K];
                obs.Gy = [Gyp-Cp*obs.K*Gyp;-obs.K*Gyp];
                
                % Set names
                obs.xnames = sys.getStateNames();
                obs.unames = sys.getInputNames();
                obs.ynames = sys.getOutputNames();
                obs.pnames = sys.getParameterNames();
                obs.dnames = sys.getUnmeasurableInputNames();
                
            elseif nargin == 4
                
                sys = varargin{1};
                obs.Ts = varargin{2};
                obs.Q = varargin{3};
                obs.R = varargin{4};
                
                % Check inputs
                if ~isa(sys,'ltiSys')
                    error('First argument must be a LTIsys object');
                end
                
                if obs.Ts <= 0
                    error('Sampling time must be > than 0');
                end
                
                if sys.isDiscreteTime() && sys.getSamplingTime ~= obs.Ts
                    error('Sampling time Ts is different from the system sampling time');
                end
                
                % Set number of variables
                obs.nx = sys.getNumberOfStates();
                obs.nu = sys.getNumberOfInputs();
                obs.ny = sys.getNumberOfOutputs();
                obs.np = sys.getNumberOfParameters();
                obs.nd = sys.getNumberOfUnmeasurableInputs();
                
                % Check gain dimensions
                if any(size(obs.Q) ~= [ obs.nx+obs.nd obs.nx+obs.nd ])
                    error(['Matrix Q must be a matrix with ',num2str(obs.nx+obs.nd),' rows and ',num2str(obs.nx+obs.nd),' columns']);
                end
                if any(size(obs.R) ~= [ obs.ny obs.ny ])
                    error(['Matrix R must be a matrix with ',num2str(obs.ny),' rows and ',num2str(obs.ny),' columns']);
                end
                
                % If the system is continuous-time, discretize it
                if sys.isContinuousTime 
                    disp(' ')
                    disp('Discretizing system...')
                    sys = sys.discretize(obs.Ts);
                end
                
                % Extract system matrices
                [Asys, Bsys, Csys, Dsys, Exsys, Eysys, Fxsys, Fysys, Gxsys, Gysys] = sys.getMatrices();
                
                % Put together state and unmeasurable input this way:
                %
                % X = [ x 
                %       d ];
                %
                % and input and parameters this way:
                %
                % U = [ u 
                %       p ];
                % Put together state and unmeasurable input
                Ap = [Asys Fxsys; zeros(obs.nd,obs.nx) eye(obs.nd)];
                Bp = [Bsys Exsys; zeros(obs.nd,obs.nu) zeros(obs.nd,obs.np)];
                Gxp = [Gxsys; zeros(obs.nd,1)];
                Cp = [Csys Fysys];
                Dp = [Dsys Eysys];
                Gyp = Gysys;
                
                % Matrix multiplying the model disturbance
                W = eye(obs.nx+obs.nd);
                
                % Matrix multiplying the model disturbance in the output
                % function
                V = zeros(obs.ny,obs.nx+obs.nd);
                
                % Create ss object for kalman function
                dt_sys = ss(Ap,[Bp Gxp W],Cp,[Dp Gyp V],obs.Ts);
                
                % Design Kalman filter
                disp(' ')
                disp('Designing Kalman filter...')
                [kSys,Ktmp,obs.P] = kalman(dt_sys,obs.Q,obs.R,'current');
                disp(' ')
                disp('Done.')
                disp(' ')
                               
                % Kalman gain, computed by hand because Ktmp (computed by
                % kalman function) is premultiplied by A
                obs.K = (obs.P*Cp')/(Cp*obs.P*Cp'+obs.R);
                
                % Set observer matrices
                obs.A = Ap-Ap*obs.K*Cp;
                obs.B = [Bp-Ap*obs.K*Dp Ap*obs.K];
                obs.Gx = Gxp-Ap*obs.K*Gyp;
                obs.C = [Cp-Cp*obs.K*Cp;eye(obs.nx+obs.nd)-obs.K*Cp];
                obs.D = [[Dp-Cp*obs.K*Dp Cp*obs.K];-obs.K*Dp obs.K];
                obs.Gy = [Gyp-Cp*obs.K*Gyp;-obs.K*Gyp];
                
                % Set names
                obs.xnames = sys.getStateNames();
                obs.unames = sys.getInputNames();
                obs.ynames = sys.getOutputNames();
                obs.pnames = sys.getParameterNames();
                obs.dnames = sys.getUnmeasurableInputNames();
                
            else
                
                error('Wrong input arguments.');
            end
            
        end
        
        % Get methods
        
        P = getCovarianceMatrix(object);
        
        varargout = getGain(object);
        
        varargout = getMatrices(object,varargin)
        
        % Other methods
        
        xpred = predict(object,xcur,u,p,y)
        
        [xcur, ycur] = computeOutput(object,x,u,p,y)
        
        disp(object);
        
        varargout = generateVHDL(object,circuit_parameters,options);
        
        
    end
end