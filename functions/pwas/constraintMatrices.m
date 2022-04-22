function [Ain, Bin, Aeq, Beq, nsigma, constrred] = constraintMatrices(fappr,ctrl,sys,constr)
% constraintMatrices   Computes the matrices defining the optimization constraints
%
% When approximating an MPC control function with a PWAS function, the same
% constraints imposed to the exact function are also enforced on the
% approximate function. These constraints translate into linear
% inequalities in the form
%
%                           Ain w <= Bin    (1)
%
% being w the weights defining the PWAS function. Soft constraints are
% relaxed through slack variables sigma. Additional equality constraints
% can be imposed to impose that the approximate function coincides with the
% exact one in a given point or region. These constraints translate into
% linear equalities in the form:
%
%                           Aeq w = Beq     (2)
%
% These constraints are imposed in the optimization problem described in
% the documentation of function computeMatrices.m.
%
% Not all the constraints can be applied to the approximate function but
% only the constraints at current step and the constraints at the following 
% step. The constraints at current step k which do not involve the input u(k)
% are discarded as well as the constraints at following step k+1 which do 
% involve the input u(k+1).
%
% [AIN, BIN, AEQ, BEQ, NSIGMA, CONSTR_RED] = ...
%        constraintMatrices(FPWAS,SYS,CONSTR)
% FPWAS is a pwasFunction object where the domain and the simplicial
% partition must be defined. SYS and CONSTR are a ltiSys object and a
% constraints object describing the system to be controlled and the
% constraints to be fulfilled. AIN, BIN, AEQ, BEQ are the constraints
% matrices in (1) and (2). NSIGMA is the number of slack variables and
% CONSTR_RED is a constraints object containing only the constraints which
% are actually applied to the PWAS function.

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
    
if isempty(constr)
    Ain = [];
    Bin = [];
    nsigma = 0;
    constrred = [];
else
    
    % Vertices of the simplicial partition
%     V = fappr.getVertices();
    
    % Evaluate optimal function at the vertices to retrieve the domain
    nx = ctrl.getNumberOfStates();
    np = ctrl.getNumberOfParameters();
    nd = ctrl.getNumberOfUnmeasurableInputs();
    nu = ctrl.getNumberOfInputs();
    ny = ctrl.getNumberOfOutputs();
    nxref = numel(ctrl.getTrackingVariable);
    
    % Control function domain
    [Hd, Kd] = ctrl.getDomain();
    Pd = Polyhedron(Hd,Kd);
    Hd = Hd';
    Kd = Kd';
    
    % Get vertices of each simplex
    S = fappr.getSimplices();
    V = [];
    Vin = [];
    
    count = 0;
    
    % Loop on all simplices
    for i = 1:numel(S)
        
        % Check if the vertices of the simplex are inside the controller
        % domain
        inside = all(S{i}*Hd <= repmat(Kd,nx+np+nd+nxref+1,1),2);
        if all(inside)
            % If all vertices are inside, add the vertices to the list
            Vin = [Vin; S{i}];
            count = count+1;
        else
            % If some vertices are outside, compute the intersection
            % between the simplex and the controller domain
            Ps = Polyhedron(S{i});
            Pint = Pd.intersect(Ps);
            Vint = Pint.V;
            
            % Add the vertices of the intersection to the list
            V = [V; Vint];
        end
    end
    
    V = [V;Vin];
    
    % Flag indicating if all vertices are inside the domain
    if count == numel(S)
        allinside = 1;
    else
        allinside = 0;
    end
    
    pars = getMOBYDICpars;
    tol = pars.roundtol;
    
    % Round vertices to prevent numerical problems
    V = tol*round(V/tol);
    
    % Remove repeated vertices
    V = unique(V,'rows');
    V = V';
    
    % Number of vertices of the simplicial partition
    nvtot = fappr.getNumberOfVertices();
    
    % Number of vertices of the simplicial partition inside the domain
    nv = size(V,2);
    
    % Constraints time horizon
    N = constr.getTimeHorizon();
    
    % Extract all hard constraints at time instant k and k+1
    [H{1}, K{1}] = constr.getAllConstraints('hard',0);
    if N >= 1
        [H{2}, K{2}] = constr.getAllConstraints('hard',1);
    else
        H{2} = [];
        K{2} = [];
    end
    
    % Extract all soft constraints at time instant k and k+1
    [Hs{1}, Ks{1}] = constr.getAllConstraints('soft',0);
    if N >= 1
        [Hs{2}, Ks{2}] = constr.getAllConstraints('soft',1);
    else
        Hs{2} = [];
        Ks{2} = [];
    end
    
    % Number of variables
    nr = constr.getNumberOfReferences();
    
    % Remove constraints at current instant not involving u
    % (they are meaningless)
    Hu = H{1}(:,nx+1:nx+nu);
    idx = any(Hu,2);
    H{1} = H{1}(idx,:);
    K{1} = K{1}(idx);
    
    Hu = Hs{1}(:,nx+1:nx+nu);
    idx = any(Hu,2);
    Hs{1} = Hs{1}(idx,:);
    Ks{1} = Ks{1}(idx);
    
    % Remove constraints at next instant involving u
    % (they cannot be imposed)
    Hu = H{2}(:,nx+1:nx+nu);
    idx = any(Hu,2);
    H{2} = H{2}(~idx,:);
    K{2} = K{2}(~idx);
    
    Hu = Hs{2}(:,nx+1:nx+nu);
    idx = any(Hu,2);
    Hs{2} = Hs{2}(~idx,:);
    Ks{2} = Ks{2}(~idx);
    
    % Create reduced constraints object, with the only constraints which
    % are imposed
    constrred = constraints(nx,nu,ny,np,nd,nr,1);
    constrred = constrred.setConstraints(H{1},K{1},0);
    constrred = constrred.setConstraints(H{2},K{2},1);
    constrred = constrred.setSoftConstraints(Hs{1},Ks{1},0);
    constrred = constrred.setSoftConstraints(Hs{2},Ks{2},1);
    
    % Constraint matrices at the current time instant
    Hu0 = H{1}(:,nx+1:nx+nu);
    Hx0 = H{1}(:,1:nx);
    Hp0 = H{1}(:,nx+nu+1:nx+nu+np);
    Hd0 = H{1}(:,nx+nu+np+1:nx+nu+np+nd);
    Hr0 = H{1}(:,nx+nu+np+nd+1:nx+nu+np+nd+nr);
    K0 = K{1};
    
    Hu0s = Hs{1}(:,nx+1:nx+nu);
    Hx0s = Hs{1}(:,1:nx);
    Hp0s = Hs{1}(:,nx+nu+1:nx+nu+np);
    Hd0s = Hs{1}(:,nx+nu+np+1:nx+nu+np+nd);
    Hr0s = Hs{1}(:,nx+nu+np+nd+1:nx+nu+np+nd+nr);
    K0s = Ks{1};
    
    % Constraint matrices at next time instant
    Hu1tmp = H{2}(:,nx+1:nx+nu);
    Hx1tmp = H{2}(:,1:nx);
    Hp1tmp = H{2}(:,nx+nu+1:nx+nu+np);
    Hd1tmp = H{2}(:,nx+nu+np+1:nx+nu+np+nd);
    Hr1tmp = H{2}(:,nx+nu+np+nd+1:nx+nu+np+nd+nr);
    K1tmp = K{2};
    
    Hu1stmp = Hs{2}(:,nx+1:nx+nu);
    Hx1stmp = Hs{2}(:,1:nx);
    Hp1stmp = Hs{2}(:,nx+nu+1:nx+nu+np);
    Hd1stmp = Hs{2}(:,nx+nu+np+1:nx+nu+np+nd);
    Hr1stmp = Hs{2}(:,nx+nu+np+nd+1:nx+nu+np+nd+nr);
    K1stmp = Ks{2};
    
    if any(Hx1tmp)
        if isempty(sys)
            error('Property ''sys'' is not set in the MPCctrl object. COnstraints at next time instant cannt be imposed.');
        end
        if sys.isContinuousTime()
            error('System must be in discrete-time form');
        end
        % Extract system matrices
        A = sys.getMatrices('A');
        B = sys.getMatrices('B');
        Ex = sys.getMatrices('Ex');
        Fx = sys.getMatrices('Fx');
        Gx = sys.getMatrices('Gx');
        
        Hu1 = Hx1tmp*B;
        Hx1 = Hx1tmp*A;
        Hp1 = Hx1tmp*Ex+Hp1tmp;
        Hd1 = Hx1tmp*Fx+Hd1tmp;
        Hr1 = Hr1tmp;
        K1 = K1tmp - Hx1tmp*Gx;
        
    else
        
        Hu1 = [];
        Hx1 = [];
        Hp1 = [];
        Hd1 = [];
        Hr1 = [];
        K1 = [];
        
    end
    
    if any(Hx1stmp)
        if isempty(sys)
            error('Property ''sys'' is not set in the MPCctrl object. Constraints at next time instant cannt be imposed.');
        end
        if sys.isContinuousTime()
            error('System must be in discrete-time form');
        end
        % Extract system matrices
        A = sys.getMatrices('A');
        B = sys.getMatrices('B');
        Ex = sys.getMatrices('Ex');
        Fx = sys.getMatrices('Fx');
        Gx = sys.getMatrices('Gx');
        
        Hu1s = Hx1stmp*B;
        Hx1s = Hx1stmp*A;
        Hp1s = Hx1stmp*Ex+Hp1stmp;
        Hd1s = Hx1stmp*Fx+Hd1stmp;
        Hr1s = Hr1stmp;
        K1s = K1stmp - Hx1stmp*Gx;
        
    else
        
        Hu1s = [];
        Hx1s = [];
        Hp1s = [];
        Hd1s = [];
        Hr1s = [];
        K1s = [];
    end
    
    % Put matrices together
    Hu = [Hu0;Hu1];
    Hx = [Hx0;Hx1];
    Hp = [Hp0;Hp1];
    Hd = [Hd0;Hd1];
    Hr = [Hr0;Hr1];
    K = [K0;K1];
    
    Hus = [Hu0s;Hu1s];
    Hxs = [Hx0s;Hx1s];
    Hps = [Hp0s;Hp1s];
    Hds = [Hd0s;Hd1s];
    Hrs = [Hr0s;Hr1s];
    Ks = [K0s;K1s];
    
    Hv = [Hx Hp Hd Hr];
    Hvs = [Hxs Hps Hds Hrs];
    
    % Separate constraints involving only inputs from constraints involving
    % also states
    idx = any(Hv,2);
    
    Hvu = Hv(~idx,:);
    Hvx = Hv(idx,:);
    Huu = Hu(~idx,:);
    Hux = Hu(idx,:);
    Ku = K(~idx);
    Kx = K(idx);
    
    idx = any(Hvs,2);
    
    Hvus = Hvs(~idx,:);
    Hvxs = Hvs(idx,:);
    Huus = Hus(~idx,:);
    Huxs = Hus(idx,:);
    Kus = Ks(~idx);
    Kxs = Ks(idx);
    
    % Create matrix alpha on all simplicial partition vertices (the input 
    % constraints are imposed on all vertices in order to saturate the 
    % approximate function, thus preventing the weights to assume huge
    % values
    alpha = speye(nvtot);
    
    % Create matrix alpha on the simplicial partition vertices inside
        % the domain and on the vertices in the intersection between the
        % simplices and the controller domain (the state constraints are 
        % imposed only on these vertices  
        nzero = nv*nvtot;   % TO DO this is an upper bound
        alpha = fappr.alphaBasis(V,nzero);
    
    
    % Create matrices Ain and Bin for constraints involving only u
    Ainu = [];
    Ainus = [];
    Binu = [];
    Binus = [];
    if ~isempty(Huu)
        Ainu = kron(Huu,alpha);
        Ku = repmat(Ku,1,nv);
        Ku = Ku';
        Ku = Ku(:);
        Binu = sparse(Ku); 
    end
    if ~isempty(Huus)
        Ainus = kron(Huus,alpha);
        Kus = repmat(Kus,1,nv);
        Kus = Kus';
        Kus = Kus(:);
        Binus = sparse(Kus);
    end
     
    if ~allinside
        % Create matrix alpha on the simplicial partition vertices inside
        % the domain and on the vertices in the intersection between the
        % simplices and the controller domain (the state constraints are 
        % imposed only on these vertices  
        nzero = nv*nvtot;   % TO DO this is an upper bound
        alpha = fappr.alphaBasis(V,nzero);
    end
    
    % Create matrices Ain and Bin for constraints involving also x
    Ainx = [];
    Ainxs = [];
    Binx = [];
    Binxs = [];

    if ~isempty(Hux)
        tmp = Hvx*V;
        tmp = tmp';
        tmp = tmp(:);
        Kx = repmat(Kx,1,nv);
        Kx = Kx';
        Kx = Kx(:);
        Ainx = kron(Hux,alpha);
        Binx = Kx-tmp;
    end
    if ~isempty(Huxs)
        tmp = Hvxs*V;
        tmp = tmp';
        tmp = tmp(:);
        Kxs = repmat(Kxs,1,nv);
        Kxs = Kxs';
        Kxs = Kxs(:);
        Ainxs = kron(Huxs,alpha);
        Binxs = Kxs-tmp;
    end
%     Ainxs = [];
%     Binxs = [];
    % TO DO
    % Controllare ancora vincoli e variabile di slack (1 o tante)
    
    % Put matrices together
    Ain = [Ainu;Ainx];
    Ains = [Ainus;Ainxs];
    Bin = [Binu;Binx];
    Bins = [Binus;Binxs];
        
    % Add single slack variable to prevent numerical problems
    Ain = [Ain [zeros(size(Ainu,1),size(Ainx,1));-speye(size(Ainx,1))]];
    Ains = [Ains zeros(size(Ains,1),size(Ainx,1))];
    
    % Number of slack variables
    nsigma = size(Ains,1);    
    
    % Add portion to matrix a to manage slack variables
    Ain = [Ain zeros(size(Ain,1),nsigma)];
    Ains = [Ains -speye(nsigma)];
    
    % Also include the single slack variable used to prevent numerical problems
    nsigma = nsigma+size(Ainx,1);
    
    % Put matrices together
    Ain = [Ain; Ains];
    Bin = [Bin; Bins];
    
end

% Number of variables
np = ctrl.getNumberOfParameters();
nd = ctrl.getNumberOfUnmeasurableInputs();
ndim = ctrl.getNumberOfDimensions();

tracking = ctrl.isTracking();

if np > 0 || nd > 0 || tracking
    Aeq = [];
    Beq = [];
else
    
    % COmpute alpha basis in reference state
    xref = ctrl.getReference();
    nzero = 4^ndim;
    alphaxref = fappr.alphaBasis(xref,nzero);
    Aeq = kron(eye(nu),alphaxref);
    Beq = ctrl.eval(xref);
    Beq = Beq(:);
    
    Aeq = [Aeq zeros(size(Aeq,1),nsigma)];
end

