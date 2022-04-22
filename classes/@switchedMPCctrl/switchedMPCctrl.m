classdef switchedMPCctrl < controller
    
    % switchedMPCctrl   Explicit switchedMPC controller
    %
    % This object represents an explicit switched MPC controller for a PWA
    % system. The controller can be designed either for regulation to a
    % costant reference state or for tracking an external signal.  The
    % control function u = f(x,p,d,ref) is a piecewise-affine function of
    % the system state (x), parameters (p), unmeasurable inputs (d) and,
    % only for tracking (true/false) , of the reference signal (ref). The
    % function is defined over a generic polytopic domain partition
    % (pwagFunction object).
    % The control function is designed to satisfy constraints on
    % the system variables.
    %
    % See the User's Guide for a more detailed explaination of this object.
    %
    % object = switchedMPCctrl()
    % Builds an empty switchedMPCctrl object.
    %
    % object = switchedMPCctrl(SYS,Ts,CONSTR,OPTS)
    % Builds a switchedMPCctrl object starting from a pwaSys object SYS
    % and a constraints object CONSTR. Ts is the sampling time of the
    % controller. SYS represents the system to regulate and CONSTR the
    % inequality constraints to be fulfilled. OPTS is a structure with the
    % following fields:
    % - N, Nu: prediction and control horizon (default: Nu = N)
    % - P, Q, R: matrices defining the MPC cost function (default: P = Q)
    % - norm: norm used in the cost function; available values are 1, 2 or
    %         inf (default: norm = 2)
    % - K, O: define the control function after the control horizon Nu.
    %         u = Kx + O (default: K = 0, O = 0)
    % - rho: penalty for the slack variables of soft constraints (default:
    %        rho = 1e3)
    % - tracking: boolean indicating if the controller is for tracking
    %             (default: false)
    % - trackvariable: indices of the state variables to track (only for
    %                  tracking controllers)
    % - ref, uref: reference state and reference input (default: xref = 0,
    %               uref = 0)
    % - refbounds: bounds for the reference state. The first column
    %              contains the lower bounds, the second the upper bounds.
    %              refbounds must have as many rows as the number of
    %              variables to track (only for tracking controllers)
    % - constantInputAfterNu : boolean indicating if the control must be
    %                          u = Kx + O or u(i) = u(Nu) after Nu
    %
    %
    % The switchedMPC controller is designed through Yalmip and MPT3. (?)
    %
    % switchedMPCctrl methods:
    %   disp - displays some information about the pwagFunction object
    %   getNumberOfDynamics - gets the number of dynamics of the PWA-sys
    %   eval - evaluates the MPC controller
    %
    % The switchedMPCctrl object is derived from controller and inherits
    % all its methods.
    %
    % See also pwagFunction, ltiSys, pwaSys, constraints.
    
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
    
    % TO DO
% gestire controllori per output

    properties (Access = private)
        
        sys = [];
        constr = [];
        options = [];
        controllers = [];
        nDyn = 0;
        
    end
    
    methods
        
        function object = switchedMPCctrl(varargin)
            
            if nargin == 0
                object.controllers = [];
                object.nDyn = 0;
                
            elseif nargin == 4
                
                % Extract inputs
                sys = varargin{1};
                object.Ts = varargin{2};
                constr = varargin{3};
                options = varargin{4};
                
                object.sys = sys;
                object.constr = constr;
                object.options = options;
                
                % Check inputs
                if ~isa(sys,'pwaSys')
                    error('First argument must be pwaSys object');
                end
                
                if object.Ts <= 0
                    error('Sampling time must be > than 0');
                end
                
                if sys.isDiscreteTime() && sys.getSamplingTime ~= object.Ts
                    error('Sampling time Ts is different from the system sampling time');
                end
                
                object.nDyn = sys.getNumberOfDynamics;
                
                nconstrs = numel (constr);
                noptions = numel (options);
                
                if nconstrs == 1
                    constr = repmat (constr,object.nDyn,1);
                elseif nconstrs ~= object.nDyn
                    error ('The number of constrains must be the same of number of dynamics')
                end
                
                if noptions == 1
                    options = repmat (options,object.nDyn,1);
                elseif noptions ~= object.nDyn
                    error ('The number of options must be the same of number of dynamics')
                elseif noptions == object.nDyn
                    truefalse = options(1).tracking;
                    for j = 1:object.nDyn
                        if options(j).tracking ~= truefalse
                            error ('The option tracking must be true/false for all the nDyn dynamics')
                        end
                    end
                end
                
                % If the system is continuous-time, discretize it
                if sys.isContinuousTime
                    disp(' ')
                    disp('Discretizing system...')
                    sys = sys.discretize(object.Ts);
                end
                
                % Set number of variables
                object.nx = sys.getNumberOfStates();
                object.nu = sys.getNumberOfInputs();
                object.np = sys.getNumberOfParameters();
                object.nd = sys.getNumberOfUnmeasurableInputs();
                
                object.nDyn = sys.getNumberOfDynamics();
                
                % Names
                xnames = cell(object.nx,1);
                unames = cell(object.nu,1);
                pnames = cell(object.np,1);
                dnames = cell(object.nd,1);
                for i = 1:object.nx
                    xnames{i} = sys.getStateNames{i};
                end
                for i = 1:object.nu
                    unames{i} = sys.getInputNames{i};
                end
                for i = 1:object.np
                    pnames{i} = sys.getParameterNames{i};
                end
                for i = 1:object.nd
                    dnames{i} = sys.getUnmeasurableInputNames{i};
                end
                
                object.xnames = xnames;
                object.unames = unames;
                if object.np > 0
                    object.pnames = pnames;
                end
                if object.nd > 0
                    object.dnames = dnames;
                end
                
                %Extract system matrices
                [Asys, Bsys, Csys, Dsys, Exsys, Eysys, Fxsys, Fysys, Gxsys, Gysys, H, K] = sys.getMatrices();
                
                for i = 1:object.nDyn
                    % Create continuous-time nDyn-ltiSys object
                    linsys = ltiSys(object.nx,object.nu,sys.getNumberOfOutputs,object.np,object.nd,'dt',object.Ts);
                    
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
                    
                    object.controllers(i).H = H{i};
                    object.controllers(i).K = K{i};
                    
                    object.controllers(i).MPCcontroller = MPCctrl(linsys,object.Ts,constr(i),options(i));
                    
                end
                
                % Check trackvariableIt is also necessary to call function interruptTimerRoutine()
                    if options(1).tracking
                        object.tracking = true;
                        object.trackvar = options(1).trackvariable;
                    else
                        object.tracking = false;
                        object.trackvar = [];
                    end
                
                object.fun = object.getPwagControlFunction();
                object.nfun = object.fun.getCodomainDimensions;
            else
                error ('Wrong input arguments!');
            end
            
            
            
        end
        
        varargout = eval(object,x,varargin);
        
        generateC(object,varargin);
        
        disp(object);
        
        dyn = findDynamics(object,x,p,d);
        
        nDyn = getNumberOfDynamics(object);
        
        inf = getInformation(object);
        
        plotPartition(object);
        
        plot(object);

        
    end
    
    
end


