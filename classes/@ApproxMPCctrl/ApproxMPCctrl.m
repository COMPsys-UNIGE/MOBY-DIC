classdef ApproxMPCctrl < controller
    % ApproxMPCctrl   Approximate explicit MPC controller
    %
    % This object represents an approximate explicit MPC controller for a
    % LTI system. The controller can be designed either for regulation
    % to a costant reference state or for tracking an external signal. The
    % control function u = f(x,p,d,ref) is a piecewise-affine function of
    % the system state (x), parameters (p), unmeasurable inputs (d) and,
    % only for tracking, of the reference signal (ref). The function is
    % defined over a uniform or non-uniform simplicial domain partition
    % (pwasFunction object). The control function is designed to satisfy
    % constraints on the system variables. Some constraints could be
    % relaxed due to the approximation of the exact control function.
    %
    % See the User's Guide for a more detailed explaination of this object.
    %
    % OBJ = ApproxMPCctrl()
    % Builds an empty ApproxMPCctrl object OBJ.
    %
    % OBJ = ApproxMPCctrl(MPCcontroller,OPTS)
    % Builds an ApproxMPCctrl object starting from an MPCctrl object.
    % OPTS is a structure with the following fields:
    % - nsamples: the integrals needed to set up the L2 norm minimization are
    %             computed numerically on a grid of points. The PWAS function
    %             domain is divided into hyper-rectangles which are in turn
    %             split into simplices. In each hyper-rectangle a regular grid
    %             of points is computed where to evaluate the integrals. The
    %             grid of points is computed by taking nsamples points along
    %             each side of the hyper-rectangle. The greater nsamples, the
    %             more accurate approximation, at the cost of a greater
    %             computation time. By increasing the number of dimensions and
    %             the number of simplices, the computation effor grows very
    %             fast. Low values of nsamples should therefore be used.
    % - np: this is the number of subdivisions per dimensions, which generates
    %       a uniform simplicial partition. np can be an array with as many
    %       elements as the number of domain dimensions, or a scalar. In this
    %       last case the same number of subdivisions is applied to all domain
    %       dimensions. If a non uniform partition is desired, np must be empty
    %       (or not set) and OPTS.P must be provided.
    % - P: this defines any non uniform simplicial partition. P is a cell-array
    %      with the same number of components as the domain dimensions. The
    %      j-th element of P contains the coordinates of the subdivisions along
    %      the j-th dimension.
    %      For example if P{1} = [1 3 4] and P{2} = [1 2 4], the resulting
    %      partition is the following (each rectangle is in turn divided into
    %      two triangles):
    %      4  ___________
    %        |       |   |
    %        |       |   |
    %        |       |   |
    %      2 |_______|___|
    %        |       |   |
    %      1 |_______|___|
    %        1       3   4
    %    If a uniform partition is desired, P must be empty (or not set) and
    %    OPTS.np must be provided.
    % - xmin, xmax: arrays defining the hyper-rectangular domain the
    %               approximate function is defined over. If they are not
    %               provided, the hyper-rectangular domain is automatically
    %               computed by taking the bounding box of the exact function
    %               domain.
    %
    % ApproxMPCctrl methods:
    %   disp - displays some information about the ApproxMPCctrl object.
    %   getNumberOfRegions - gets the number of simplices of the domain partition
    %   eval - evaluates the approximate MPC controller.
    %   plot - plots the approximate MPC control function.
    %   plotPartition - plots the simplicial domain partition.
    %
    % The ApproxMPCctrl object is derived from controller and inherits
    % all its methods.
    %
    % See also pwasFunction, ltiSys, pwaSys, constraints, MPCctrl.
    %
    % ApproxMPCctrl methods:
    %   disp - displays some information about the ApproxMPCctrl object.
    %   getNumberOfRegions - gets the number of simplices of the domain partition
    %   eval - evaluates the approximate MPC controller.
    %   plot - plots the approximate MPC control function.
    %   plotPartition - plots the simplicial domain partition.
    %
    % The ApproxMPCctrl object is derived from controller and inherits
    % all its methods.
    %
    % See also pwasFunction, ltiSys, pwaSys, constraints, MPCctrl.
    
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
        
        sys = [];   % dynSys object representing the system to regulate
        constr = [];    % constraints object representing the constraints
        options = [];   % MPC settings
        
    end
    
    
    methods
        
        function object = ApproxMPCctrl(varargin)
            
            if nargin == 0
                
                object.sys = [];
                object.constr = [];
                object.options = [];
                
            elseif nargin == 2
                
                if ~isa(varargin{1},'explicitMPCctrl')
                    error('First argument must be a ''explicitMPCctrl'' object');
                end
                if ~isstruct(varargin{2})
                    error('Second argument must be a struct');
                end
                
                ctrl = varargin{1};
                opts = varargin{2};
                
                fun = ctrl.fun;
                
                % Set default options
                opts = PWASset(fun,opts);
                
                % Hyper-rectangular domain
                xmin = opts.domain.xmin;
                xmax = opts.domain.xmax;
                
                % Retrieve function domain
                [Hd, Kd] = ctrl.getDomain;
                
                % Compute bounding box
                P = Polyhedron(Hd,Kd);
                B = P.outerApprox;
                
                % PWAG function domain
                xminc = B.Internal.lb;
                xmaxc = B.Internal.ub;
                
                pars = getMOBYDICpars();
                tol = pars.roundtol;
                
                % Round to prevent numerical problems
                xminc = tol*round(xminc/tol);
                xmaxc = tol*round(xmaxc/tol);
                
                % If the domain is not provided set the bounding box of the
                % controller domain
                if isempty(xmin)
                    xmin = xminc;
                    xmax = xmaxc;
                end
                
                % Check if the approximate function domain is equal or smaller
                % with the respect to the controller domain
                if any(xmin > xminc) || any(xmax < xmaxc)
                    warning(['The domain of the approximate controller is ',...
                        'smaller than the domain of the exact controller. ',...
                        'This may cause the trajectories to go outside the domain. ',...
                        'Consider enlarging the domain of the approximate controller.']);
                end
                
                % Create domain structure
                domain.xmin = xmin;
                domain.xmax = xmax;
                
                % Create pwasFunction template
                if ~isempty(opts.P)
                    P = opts.P;
                    fappr = pwasFunction(domain,P);
                else
                    np = opts.np;
                    fappr = pwasFunction(domain,np);
                end
                % Number of vertices
                nv = fappr.getNumberOfVertices();
                
                % Codomain dimensions
                ny = ctrl.getNumberOfInputs;
                
                disp(' ')
                disp('Computing cost function matrices...')
                disp(' ')
                % Compute matrices for the L2 approximation
                
                [H, f, alpha, fopt] = computeMatrices(fappr,ctrl,opts);
                
                % Extract system and constraints
                info = ctrl.getInformation;
                sys = info.sys;
                constr = info.constr;
                
                disp('Computing constraint matrices...')
                disp(' ')
                [Ain, Bin, Aeq, Beq, nsigma, constrred] = constraintMatrices(fappr,ctrl,sys,constr);
                
                disp('Solving QP problem...')
                disp(' ')
                
                % Original matrices
                Horig = H;
                forig = f;
                
                if isempty(Ain) && isempty(Aeq)
                    
                    % Perform least squares optimization
                    w = H\f;
                    
                    % Reshape weights to the correct shape
                    w = reshape(w,nv,ny);
                    
                else
                    
                    pars = getMOBYDICpars();
                                        
                    MOBYpar = getMOBYDICpars();
                    
                    % Add part for Tikhonov regularization
                    H = H+pars.tikhonov*speye(size(H));
                    
                    % Add part managing the slack variables
                    H = blkdiag(H,pars.slack*speye(nsigma));
                    f = [f;zeros(nsigma,1)];
                    
                    if strcmp(MOBYpar.qpsolver,'mpt')
                        % Create Opt object (MPT3 object defining an
                        % optimization problem)
                        problem = Opt('H',H,'f',-f,'A',Ain,'b',Bin,'Ae',Aeq,'be',Beq);
                        % Solve optimization problem through MPT3 interface
                        sol = problem.solve;
                        
                        % Optimal solution
                        xopt = sol.xopt;
                        % Status
                        exitflag = sol.exitflag;
                        how = sol.how;
                    elseif strcmp(MOBYpar.qpsolver,'quadprog')
                        [sol, fval, outFlag] = quadprog(H,-f,Ain,Bin,Aeq,Beq);
                        xopt = sol;
                        exitflag = outFlag;
                        how = ['quadrprog error ',num2str(outFlag)];
                    else
                        error(['Solver ',MOBYpar.qpsolver,' is not allowed'])
                    end
                    if exitflag ~= 1
                        error(['The QP solver returned the following status:', how]);
                    end
                    
                    % Extract weights
                    w = xopt(1:ny*nv);
                    
                    % Extract slack variables
                    sigma = xopt(ny*nv+1:end);
                    
                    if ~isempty(sigma)
                        disp(['Maximum constraints violation: ',num2str(max(sigma))]);
                    end
                    
                    % Reshape weights
                    w = reshape(w,nv,ny);
                    
                    % Optimal value
                    % TO DO
                    % Usare un valore piu significativo
                    err = alpha*w-fopt;
                    rmse = sqrt(sum(err.^2)/size(err,1));
                    
                    disp(['Estimated RMSE: ',num2str(rmse)]);
                    
                    
                    
                end
                
                disp('Done.')
                disp(' ')
                
                % Assign weights to the approximate function
                fappr = fappr.setWeights(w);
                
                % Set number of variables
                nx = ctrl.getNumberOfStates();
                nu = ctrl.getNumberOfInputs();
                np = ctrl.getNumberOfParameters();
                nd = ctrl.getNumberOfUnmeasurableInputs();
                object.nx = nx;
                object.nu = nu;
                object.np = np;
                object.nd = nd;
                
                % Set names
                xnames = ctrl.getStateNames();
                unames = ctrl.getInputNames();
                pnames = ctrl.getParameterNames();
                dnames = ctrl.getUnmeasurableInputNames();
                object.xnames = xnames;
                object.unames = unames;
                object.pnames = pnames;
                object.dnames = dnames;
                object.fun = fappr;
                object.nfun = 1;
                
                object.tracking = ctrl.isTracking();
                object.trackvar = info.options.trackvariable;
                object.Ts = ctrl.getSamplingTime();
                
                object.sys = info.sys;
                object.constr = constrred;
                
                options = info.options;
                options.apporoximation = opts;
                object.options = options;
                
                disp(' ')
                disp('Done.')
                disp(' ')
                
            else
                error('Wrong input arguments!');
            end
            
        end
        
        varargout = eval(object,x,varargin);
        
        disp(object);
        
        plot(object,k);
        
        plotPartition(object);
        
        nreg = getNumberOfRegions(object);
        
        info = getInformation(object);
        
    end
    
end