classdef switchedApproxMPCctrl < controller
    
    % switched_ApproxMPCctrl  Explicit approximate switched_MPC controller
    %
    % This object represents an explicit approximate switched MPC controller
    % for a PWA system.
    % The controller can be designed either for regulation to a
    % costant reference state or for tracking an external signal. The
    % control function u = f(x,p,d,ref) is a piecewise-affine Simplicial
    % function of the system state (x), parameters (p), unmeasurable inputs
    % (d) and, only for tracking (true/false),of the reference signal (ref).
    % The function is defined over a uniform or non-uniform simplicial domain
    % partition (pwasFunction object).
    % The control function is designed to satisfy constraints on
    % the system variables.
    %
    % See the User's Guide for a more detailed explaination of this object.
    %
    % object = switchedApproxMPCctrl()
    % Builds an empty switched_ApproxMPCctrl object .
    %
    % object = switchedApproxMPCctrl(switched_MPCctrl,OPTS)
    % Builds a switchedApproxMPCctrl object starting from a switchedMPCctrl
    % object and approximation options OPTS.
    %
    % The switchedApproxMPCctrl controller is designed through Yalmip and MPT3.
    %
    % switched_ApproxMPCctrl methods:
    %   disp - displays some information about the pwagFunction object
    %   getNumberOfDynamics - gets the number of dynamics of the PWA-sys
    %   eval - evaluates the MPC controller
    %   getInformation - get the information of the pwa sys starting from
    %   the MPc ctrl
    %   findDynamics - find the index of the dynamic
    %
    % The switched_Approx MPCctrl object is derived from controller and
    % inherits all its methods.
    %
    % See also pwasFunction, switched_MPCctrl, pwaSys, options.
    
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
    
    
    properties (Access = private)
        
        sys = [];       % dynSys object representing the system to regulate
        constr = [];    % constraints object representing the constraints
        Approxoptions = [];   % MPC settings
        controllers = [];
        nDyn = 0;
        
    end
    
    methods
        
        function object = switchedApproxMPCctrl(varargin)
            
            if nargin == 0  % object=switched_MPCctrl(switched_MPCctrl,options)
                % where the options are the approx options
                object.sys = [];
                object.constr = [];
                object.options = [];
                object.controllers = [];
                object.nDyn = 0;
                
            elseif nargin == 2
                
                % Extract inputs
                ctrl = varargin{1};
                options = varargin{2};
                
                object.nDyn = getNumberOfDynamics(ctrl);
                TsCtrl = ctrl.getSamplingTime();
                
                % Checkinputs
                if ~isa(varargin{1},'switchedMPCctrl')
                    error('First argument must be a ''switchedMPCctrl'' object');
                end
                
                numopt = numel (options);
                if numopt == 1
                    if ~isstruct(varargin{2})
                        error('Second argument must be a struct');
                    end
                    options = repmat (options,object.nDyn,1);
                elseif numopt ~= object.nDyn
                    error ('The number of options must be the same of number of dynamics');
                end
                
                
                if isfield(options(1),'np')
                np = options(1).np;
                for i=2:object.nDyn
                    if isfield(options(i),'np')
                        if np ~= options(i).np
                             error('All dynamics  must have same partition');
                        end
                    else
                        error('All dynamics  must have same partition');
                    end
                end
                end
                
                if isfield(options(1),'P')
                P = options(1).P;
                for i=2:object.nDyn
                    if isfield(options(i),'P')
                        if P ~= options(i).P
                             error('All dynamics  must have same partition');
                        end
                    else
                        error('All dynamics  must have same partition');
                    end
                end
                end
                
                
                [Hd Kd] = ctrl.getDomain;
                P = Polyhedron(Hd,Kd);
                if ~P.isBounded
                    error('Controller function must be bounded!');
                end
                Dp = P.outerApprox();
                DVert = Dp.V;
                
                minV = min(DVert);
                maxV = max(DVert);
                
                for i=1:object.nDyn
                    options(i).xmin = minV;
                    options(i).xmax = maxV;
                end
                
                
                % Names
                object.nx = ctrl.getNumberOfStates();
                object.nu = ctrl.getNumberOfInputs();
                object.np = ctrl.getNumberOfParameters();
                object.nd = ctrl.getNumberOfUnmeasurableInputs();
                object.Approxoptions = options;
                object.Ts = TsCtrl;
                
                object.tracking = ctrl.isTracking;
                object.trackvar = ctrl.getTrackingVariable;
                
                xnames = cell(object.nx,1);
                unames = cell(object.nu,1);
                pnames = cell(object.np,1);
                dnames = cell(object.nd,1);
                
                for i = 1:object.nx
                    xnames{i} = ctrl.getStateNames{i};
                end
                for i = 1:object.nu
                    unames{i} = ctrl.getInputNames{i};
                end
                for i = 1:object.np
                    pnames{i} = ctrl.getParameterNames{i};
                end
                for i = 1:object.nd
                    dnames{i} = ctrl.getUnmeasurableInputNames{i};
                end
                
                object.xnames = xnames;
                object.unames = unames;
                if object.np > 0
                    object.pnames = pnames;
                end
                if object.nd > 0
                    object.dnames = dnames;
                end
                
                disp(' ')
                disp('Approximating explicit switchedMPC controller...')
                
                
                
                Hmat = ctrl.getInformation.sys.getMatrices('H');
                Kmat = ctrl.getInformation.sys.getMatrices('K');
                
                HforFunc = cell(1,object.nDyn);
                KforFunc = cell(1,object.nDyn);
                
                MPCctrlarray = ctrl.getInformation.controllers;
                f = cell (1,object.nDyn);
                
                for j = 1:object.nDyn
                    
                    info = getInformation (ctrl);
                    object.sys = info;
                    
                    object.controllers(j).H = Hmat{j};
                    object.controllers(j).K = Kmat{j};
                    
                    object.controllers(j).ApproxMPCcontroller = ApproxMPCctrl(MPCctrlarray(j).MPCcontroller,options(j));
                    
                    f{j} = object.controllers(j).ApproxMPCcontroller.getFunction;
                    
                    HforFunc{j} = [Hmat{j}, zeros(1,numel(ctrl.getTrackingVariable))];
                    KforFunc{j} = Kmat{j};
                end
                
                disp(' ')
                disp('Done.')
                disp(' ')
                
                object.nfun = 1;
                object.fun = discontinuousPwasFunction(f,HforFunc,KforFunc);
                
            else
                error('Wrong input arguments!');
            end
            
        end
        
        varargout = eval(object,x,varargin);
        
        varargout = generateC(object,varargin);
        
        disp(object);
        
        nDyn = getNumberOfDynamics(object);
        
        inf = getInformation(object);
        
        dyn = findDynamics(object,x,p,d);
        
        plotPartition(object);
        
        plot(object);
    end
    
end




