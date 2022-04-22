classdef switchedKalmanPredictor < observer
    % switchedKalmanPredictor switched  Kalman predictor for the estimation
    %                   of the states of a PWA system
    %
    % This object represents a switched Kalman predictor for the estimation of
    % unmeasurable states of a pwa system (pwaSys object). The knowledge of
    % all system states is necessary for the design of the MPC controller.
    % The Kalman predictor is a dynamical system in the form:
    %
    % The Kalman predictor is a discrete-time dynamical system, whose state
    % evolution follows the evolution of the real system state.
    % Refer to MOBY-DIC Toolbox User's Guide for a detailed explanation of
    % the Kalman predictor.
    %
    % OBJ = switchedKalmanPredictor()
    % Builds an empty kalmanPredictor object OBJ.
    %
    % OBJ = switchedKalmanPredictor(SYS,Ts,K)
    % Builds a kalmanPredictor object OBJ by specifying the ltiSys object
    % (SYS) whose states need to be estimated, the observer sampling time
    % Ts and the Kalman gain K. The gain K weights the difference between
    % the measured system output and the estimated one. K should be the
    % solution of the Riccati equation. K could be a matrices or a cell
    % array with nDyn x 1 elements; if it is a matrix the value of K is the
    % same for each dynamic, otherwise the i-th element of the cell array
    % define the k matrix for the i-th dynamic.
    %
    % OBJ = switchedKalmanPredictor(SYS,Ts,Q,R)
    % Builds a kalmanPredictor object OBJ by specifying the ltiSys object
    % (SYS) whose states need to be estimated, the observer sampling time
    % Ts and the covariance matrices Q and R of the model and measurement
    % noises, respectively. Q and R could be a matrices or a cell array with
    % nDyn x 1 elements; if they are matrices the value of Q and R are the
    % same for each dynamic, otherwise the i-th element of the cell array
    % define the covariance matrices for the i-th dynamic.
    % If the elements in Q{i} are higher than the
    % elements in R{i}, the measurements are trusted more than the model and
    % viceversa. Q{i} must be a square matrix with dimension NX+ND, being NX
    % the number of system states and ND the number of unmeasurable inputs.
    % R{i} must be a matrix of dimension NY, being NY the number of system
    % outputs. The Kalman gain K is computed through MATLAB function
    % 'kalman'.
    %
    % switchedKalmanPredictor methods:
    %   computeOutput - Computes the Kalman predictor output
    %   disp - displays information about kalmanPredictor object
    %   getCovarianceMatrix - gets the asymptotic covariance matrix
    %   getGain - gets the gain of the Kalman predictor
    %   getMatrices - gets the matrices defining the Kalman predictor
    %   predict - predicts the system state at next time instant based on current
    %             measurements
    %
    % The kalmanPredictor object is derived from observer and inherits
    % all its methods.
    %
    % See also switchedKalmanFilter.
    
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
        filters = [];
        nDyn = 0;        
    end
    
    methods
        
        function obs = switchedKalmanPredictor(varargin)
            
            if nargin == 0
                
                obs.filters = [];
                obs.nDyn = 0;
            elseif nargin == 3
                
                pwasys = varargin{1};
                obs.Ts = varargin{2};
                kk = varargin{3};
                
                % Check inputs
                if ~isa(pwasys,'pwaSys')
                    error('First argument must be pwaSys object');
                end
                
                if obs.Ts <= 0
                    error('Sampling time must be > than 0');
                end
                
                if pwasys.isDiscreteTime() && pwasys.getSamplingTime ~= obs.Ts
                    error('Sampling time Ts is different from the system sampling time');
                end
                
                % Set number of variables
                obs.nx = pwasys.getNumberOfStates();
                obs.nu = pwasys.getNumberOfInputs();
                obs.ny = pwasys.getNumberOfOutputs();
                obs.np = pwasys.getNumberOfParameters();
                obs.nd = pwasys.getNumberOfUnmeasurableInputs();
                
                nDyn = pwasys.getNumberOfDynamics;
                
                obs.nDyn = nDyn;
                
                % Check if k is cell array or a matrix
                if iscell(kk)
                    if ~(numel(kk) == nDyn)
                        error('K must be a nDyn x 1 cell array or a matrix');
                    end
                else
                    tmp = kk;
                    kk = cell(nDyn,1);
                    for i=1:nDyn
                        kk{i} = tmp;
                    end
                end
                
                % If the system is continuous-time, discretize it
                if pwasys.isContinuousTime
                    disp(' ')
                    disp('Discretizing system...')
                    pwasys = pwasys.discretize(obs.Ts);
                end
                
                
                % Extract system matrices
                [Asys, Bsys, Csys, Dsys, Exsys, Eysys, Fxsys, Fysys, Gxsys, Gysys, H, K] = pwasys.getMatrices();
                
                for i=1:nDyn
                    
                    % Check gain dimensions
                    if any(size(kk{i}) ~= [ obs.nx+obs.nd obs.ny ])
                        error(['Gain K{i} must be a matrix with ',num2str(obs.nx+obs.nd),' rows and ',num2str(obs.ny),' columns']);
                    end
                    
                    % Create continuous-time ltiSys object
                    linsys = ltiSys(obs.nx,obs.nu,obs.ny,obs.np,obs.nd,'dt',obs.Ts);
                    
                    % Set system matrices
                    linsys = linsys.setMatrices('A',Asys{i});
                    linsys = linsys.setMatrices('B',Bsys{i});
                    linsys = linsys.setMatrices('C',Csys{i});
                    linsys = linsys.setMatrices('D',Dsys{i});
                    linsys = linsys.setMatrices('Ex',Exsys{i});
                    linsys = linsys.setMatrices('Ey',Eysys{i});
                    linsys = linsys.setMatrices('Fx',Fxsys{i});
                    linsys = linsys.setMatrices('Fy',Fysys{i});
                    linsys = linsys.setMatrices('Gx',Gxsys{i});
                    linsys = linsys.setMatrices('Gy',Gysys{i});
                    
                    % Set variable names
                    linsys = linsys.setStateNames(pwasys.getStateNames());
                    linsys = linsys.setInputNames(pwasys.getInputNames());
                    linsys = linsys.setParameterNames(pwasys.getParameterNames());
                    linsys = linsys.setUnmeasurableInputNames(pwasys.getUnmeasurableInputNames());
                    linsys = linsys.setOutputNames(pwasys.getOutputNames());
                    
                    obs.filters(i).H = H{i};
                    obs.filters(i).K = K{i};
                    obs.filters(i).kPredictor = kalmanPredictor(linsys,obs.Ts,kk{i});
                    
                    obs.filters(i) = fil_tmp;
                end
                % Set names
                obs.xnames = pwasys.getStateNames();
                obs.unames = pwasys.getInputNames();
                obs.ynames = pwasys.getOutputNames();
                obs.pnames = pwasys.getParameterNames();
                obs.dnames = pwasys.getUnmeasurableInputNames();
                
            elseif nargin == 4
                
                pwasys = varargin{1};
                obs.Ts = varargin{2};
                Q = varargin{3};
                R = varargin{4};
                
                % Check inputs
                if ~isa(pwasys,'pwaSys')
                    error('First argument must be a pwaSys object');
                end
                
                if obs.Ts <= 0
                    error('Sampling time must be > than 0');
                end
                
                if pwasys.isDiscreteTime() && pwasys.getSamplingTime ~= obs.Ts
                    error('Sampling time Ts is different from the system sampling time');
                end
                
                % Set number of variables
                obs.nx = pwasys.getNumberOfStates();
                obs.nu = pwasys.getNumberOfInputs();
                obs.ny = pwasys.getNumberOfOutputs();
                obs.np = pwasys.getNumberOfParameters();
                obs.nd = pwasys.getNumberOfUnmeasurableInputs();
                
                % If the system is continuous-time, discretize it
                if pwasys.isContinuousTime
                    disp(' ')
                    disp('Discretizing system...')
                    pwasys = pwasys.discretize(obs.Ts);
                end
                
                nDyn = pwasys.getNumberOfDynamics;
                obs.nDyn = nDyn;
                % Check if k is cell array or a matrix
                if iscell(Q) && iscell(R)
                    if ~(numel(Q) == nDyn) || ~(numel(R) == nDyn)
                        error('Q and R must be a nDyn x 1 cell array');
                    end
                elseif ~iscell(Q) && ~iscell(R)
                    tmp1 = Q;
                    tmp2 = R;
                    
                    Q = cell(nDyn,1);
                    R = cell(nDyn,1);
                    for i=1:nDyn
                        Q{i} = tmp1;
                        R{i} = tmp2;
                    end
                else
                    error('Q and R must be cell array nDyn x 1 or matrices');
                end
                
                % Extract system matrices
                [Asys, Bsys, Csys, Dsys, Exsys, Eysys, Fxsys, Fysys, Gxsys, Gysys, H ,K] = pwasys.getMatrices();
                
                
                
                
                for i=1:nDyn
                    
                    % Check gain dimensions
                    if any(size(Q{i}) ~= [ obs.nx+obs.nd obs.nx+obs.nd ])
                        error(['Matrix Q{i} must be a matrix with ',num2str(obs.nx+obs.nd),' rows and ',num2str(obs.nx+obs.nd),' columns']);
                    end
                    if any(size(R{i}) ~= [ obs.ny obs.ny ])
                        error(['Matrix R must be a matrix with ',num2str(obs.ny),' rows and ',num2str(obs.ny),' columns']);
                    end
                    
                    
                    
                    % Create continuous-time ltiSys object
                    linsys = ltiSys(obs.nx,obs.nu,obs.ny,obs.np,obs.nd,'dt',obs.Ts);
                    
                    % Set system matrices
                    linsys = linsys.setMatrices('A',Asys{i});
                    linsys = linsys.setMatrices('B',Bsys{i});
                    linsys = linsys.setMatrices('C',Csys{i});
                    linsys = linsys.setMatrices('D',Dsys{i});
                    linsys = linsys.setMatrices('Ex',Exsys{i});
                    linsys = linsys.setMatrices('Ey',Eysys{i});
                    linsys = linsys.setMatrices('Fx',Fxsys{i});
                    linsys = linsys.setMatrices('Fy',Fysys{i});
                    linsys = linsys.setMatrices('Gx',Gxsys{i});
                    linsys = linsys.setMatrices('Gy',Gysys{i});
                    
                    % Set variable names
                    linsys = linsys.setStateNames(pwasys.getStateNames());
                    linsys = linsys.setInputNames(pwasys.getInputNames());
                    linsys = linsys.setParameterNames(pwasys.getParameterNames());
                    linsys = linsys.setUnmeasurableInputNames(pwasys.getUnmeasurableInputNames());
                    linsys = linsys.setOutputNames(pwasys.getOutputNames());
                    
                    obs.filters(i).H = H{i};
                    obs.filters(i).K = K{i};
                    obs.filters(i).kPredictor = kalmanPredictor(linsys,obs.Ts,Q{i},R{i});
                    
                end
                % Set names
                obs.xnames = pwasys.getStateNames();
                obs.unames = pwasys.getInputNames();
                obs.ynames = pwasys.getOutputNames();
                obs.pnames = pwasys.getParameterNames();
                obs.dnames = pwasys.getUnmeasurableInputNames();
                
            else
                
                error('Wrong input arguments.');
            end
            
        end
        
        % Get methods
        
        P = getCovarianceMatrix(object);% TO DO
        
        varargout = getGain(object);% TO DO
        
        varargout = getMatrices(object,varargin)% TO DO
        
        % Other methods
        
        [xcur, ycur] = computeOutput(object,x,u,p,y)
        
        xpred = predict(object,xcur,u,p,y)
        
        disp(object); % TO DO
        
        dyn = findDynamics(object,x,p)
        
        
    end
end