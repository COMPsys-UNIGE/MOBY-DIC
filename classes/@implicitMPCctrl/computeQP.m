function [QP1, QP2, ineqL, ineqR, eqL, eqR] = computeQP(object, x, p, d, ref, varargin)
% computeMPC
% Computes the matrices of the quadratic programming (QP) problem required
% by MPC every time step.
%
% Computes matrices of the constrained quadratic programming (QP) problem
% requested by implicit MPC at the current time. The QP problem is:
%      min_z [(1/2)*z'*QP1*z + z'*QP2]
%       s.t. ineqL*z <= ineqR
%            eqL*z = eqR
%
% OBJ is the implicitMPCctrl object, x, p and d are current system states,
% parameters and unmeasurable inputs, respectively.

% varargin can comprehend the system inputs at the previous time step
% (uOld), in case of tracking controller, and algo="ADMM" in case the QP of
% the implicit MPC controller is solved with ADMM (default is MATLAB
% function quadprog).

if nargin == 5
    uOld = 0;
    algo = "";
elseif nargin == 6
    if isstring(varargin{1})
        uOld = 0;
        algo = varargin{1};
    else
        uOld = varargin{1};
        algo = "";
    end
else
    uOld = varargin{1};
    algo = varargin{2};
end

% Prediction horizon
N = object.options.N;

% Control horizon
Nu = object.options.Nu;

% Number of variables
nx = object.nx;
np = object.np;
nd = object.nd;
nu = object.nu;
ny = object.ny;
nref = length(object.getInformation.options.trackvariable);

% After the control horizon the control action for the prediction is
% u = Kx+O, if K and O are defined, otherwise u = uOld
if isfield(object.options, 'K')
    Ku = object.options.K;
else
    Ku = zeros(nu,nx);
end
if isfield(object.options, 'O')
    Ou = object.options.O;
else
    Ou = zeros(nu,1);
end

sys = object.sys;

% LTI system matrices
[A, B, C, D, Ex, Ey, Fx, Fy, Gx, Gy] = sys.getMatrices;

% MPC weights
Q = object.options.Q;
P = object.options.P;
R = object.options.R;

% If the controller track a reference state, consider augmented state
% [x, uOld] as state, and the variation of input deltaU as input
if object.options.tracking
    sys = sys.deltau();

    % Number of variables
    nx = sys.getNumberOfStates();
    nu = sys.getNumberOfInputs();
    ny = sys.getNumberOfOutputs();

    % Augment state
    x = [x; uOld];

    % Extend matrices
    Ku = [Ku, zeros(nu,nu)];
    Q = [Q, zeros(nx-nu,nu); zeros(nu, nx)];
    P = [P, zeros(nx-nu,nu); zeros(nu, nx)];

    % LTI system matrices
    [A, B, C, D, Ex, Ey, Fx, Fy, Gx, Gy] = sys.getMatrices;
end

% Reference state or output
if strcmpi(object.ctrlvar,'state')
    if ~object.options.tracking
        xref = object.ref;
    else
        xref = zeros(nx,size(ref,2));
        xref(object.options.trackvariable,:) = ref;
    end
else
    yref = zeros(ny,1);
    if ~object.options.tracking
        yref = object.ref;
    else
        yref(object.options.trackvariable) = ref;
    end
end

% Constraints matrices
[H, ~] = object.constr.getAllConstraints;

% Get indices of states with soft constraints (not supported by ADMM)
% [Hs,~] = object.constr.getAllConstraints('soft');
% softIdx = [];
% for i=1:size(Hs{1},2)
%     if any(Hs{1}(:,i)~=0)
%         softIdx = [softIdx, i];
%     end
% end
% nSoft = length(softIdx);

[Hs,~] = object.constr.getAllConstraints('soft');
softIdx = [];
for i=1:size(H{1},1)
%     if any(Hs{1}(i,:)~=0)
    if ismember(H{1}(i, :), Hs{1}, 'rows')
        softIdx = [softIdx, i];
    end
end
%
% for i=1:size(Hs{1},1)
%     for j=1:size(Hs{1},1)
%     if Hs{1}(i,:) == H{1}(i,:)
%         softIdx = [softIdx, i];
%     end
% end
nSoft = length(softIdx);

%%
[H, K] = object.constr.getAllConstraints;
for i = 1:N+1
    H{i} = H{i}(:, 1:end-np-nd-ny-nref);
    Hs{i} = Hs{i}(:, 1:end-np-nd-ny-nref);
end

ineqR = [];
% ineqR_old = [];
% [~, ineqR] = object.constr.getInputConstraints(0); 
% for i = 2:N+1
%     ineqR = [ineqR; K{i}];
% end
for i = 1:N
    
    try
        if i < Nu+1
            firstRow = length(ineqR) + 1;
            [~, ineqR_tmp] = object.constr.getInputConstraints(i-1); 
            ineqR = [ineqR; ineqR_tmp];
            lastRow = length(ineqR);
            
            
            if object.options.tracking
                try
            %             ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, i) = H{i}(:, nx-nu+1:nx);
                    [ineqL(firstRow:lastRow, (i-1)*nu + 1), ~] = object.constr.getInputConstraints(i-1);
                    [ineqL(firstRow:lastRow, Nu*nu + i*nx - nu + 1 : Nu*nu + i*nx), ~] = object.constr.getInputConstraints(i-1);
                catch
                    ineqL(firstRow:lastRow, (i-1)*nu + 1) = 0;
                    ineqL(firstRow:lastRow, Nu*nu + i*nx - nu + 1 : Nu*nu + i*nx) = 0;
                end
            else
                try
            %             ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, i) = H{i}(:, nx+1:nx+nu);
                    [ineqL(firstRow:lastRow:length(K{i})*(i-1)+2*nu, i), ~] = object.constr.getInputConstraints(i-1);
                catch
                    ineqL(firstRow:lastRow:length(K{i})*(i-1)+2*nu, i) = 0;
                end
            end
        end
    catch
    end
        
    try
        firstRow = length(ineqR) + 1;
        [~, ineqR_tmp] = object.constr.getStateConstraints(i-1); 
        ineqR = [ineqR; ineqR_tmp];
        lastRow = length(ineqR);

        if object.options.tracking
            try
                ineqL(firstRow:lastRow, Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx - nu) = H{i}(any((H{i}(:, 1:nx-nu)), 2), 1:nx-nu);
    %             ineqL(length(cell2mat(K(1:i-1)))+1:length(cell2mat(K(1:i))), Nu*nu + i*nx -nu + 1 : Nu*nu + i*nx) = H{i}(:, nx-nu+1:nx);
            catch
%                 ineqL(firstRow+1:lastRow, Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx) = 0;
            end
        else
            try
                if any(H{i}(:, 1:nx))
                    ineqL(firstRow:lastRow, Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx) = H{i}(any((H{i}(:, 1:nx)), 2), 1:nx);
    %             else
    %                ineqL(length(cell2mat(K(1:i)))+1:length(cell2mat(K(1:i+1))), Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx) = 0;
                end
            catch
    %             ineqL(length(cell2mat(K(1:i)))+1:length(cell2mat(K(1:i+1))), Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx) = 0;
            end
        end
    catch
    end
end

% ineqL = zeros(numel(ineqR), Nu*nu + N*nx);
ineqL = [ineqL, zeros(size(ineqL,1), Nu*nu + (N+1)*nx - size(ineqL,2))];

% for i = 1:Nu
%     if object.options.tracking
%         try
% %             ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, i) = H{i}(:, nx-nu+1:nx);
%             [ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, i), ~] = object.constr.getInputConstraints(i-1);
%             [ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, Nu*nu + i*nx - nu + 1 : Nu*nu + i*nx), ~] = object.constr.getInputConstraints(i-1);
%         catch
%             ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, i) = 0;
%             ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, Nu*nu + i*nx - nu + 1 : Nu*nu + i*nx) = 0;
%         end
%     else
%         try
% %             ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, i) = H{i}(:, nx+1:nx+nu);
%             [ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, i), ~] = object.constr.getInputConstraints(i-1);
%         catch
%             ineqL(length(K{i})*(i-1)+1:length(K{i})*(i-1)+2*nu, i) = 0;
%         end
%     end
% end

% for i = 1:N-1
% 
%     if object.options.tracking
%         try
%             ineqL(length(cell2mat(K(1:i)))+1:length(cell2mat(K(1:i+1))), Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx - nu) = H{i+1}(:, 1:nx-nu);
% %             ineqL(length(cell2mat(K(1:i-1)))+1:length(cell2mat(K(1:i))), Nu*nu + i*nx -nu + 1 : Nu*nu + i*nx) = H{i}(:, nx-nu+1:nx);
%         catch
%             ineqL(length(cell2mat(K(1:i-1)))+1:length(cell2mat(K(1:i))), Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx) = 0;
%         end
%     else
%         try
%             if any(H{i+1}(:, 1:nx))
%                 ineqL(length(cell2mat(K(1:i)))+1:length(cell2mat(K(1:i+1))), Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx) = H{i+1}(:, 1:nx);
% %             else
% %                ineqL(length(cell2mat(K(1:i)))+1:length(cell2mat(K(1:i+1))), Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx) = 0;
%             end
%         catch
% %             ineqL(length(cell2mat(K(1:i)))+1:length(cell2mat(K(1:i+1))), Nu*nu + (i-1)*nx + 1 : Nu*nu + i*nx) = 0;
%         end
%     end
% end

% ----- TO DO: ADD SOFT CONSTRAINTS FOR ADMM -----
% Reformulate the constraints for QP in case of soft constraints
if object.options.tracking
    uref = zeros(nu,1);

%     if nSoft == 0
%         xmin = [-K{1}(1:nx-nu); -K{1}(2*(nx-nu)+1:2*(nx-nu)+nu)];
%         xmax = [K{1}(nx-nu+1:2*(nx-nu)); K{1}(2*(nx-nu)+nu+1:2*(nx-nu)+2*nu)];
%     else
%         xmin = [];
%         xmax = [];
%         for i=1:2:2*(nx-nu)
%             xmin = [xmin; -K{1}(i)];
%             xmax = [xmax; K{1}(i+1)];
%         end
%         xmin = [xmin; -K{1}(2*(nx-nu)+1:2*(nx-nu)+nu)];
%         xmax = [xmax; K{1}(2*(nx-nu)+nu+1:2*(nx-nu)+2*nu)];
%     end
%     umin = -K{1}(2*(nx-nu)+1:2*(nx-nu)+nu) - K{1}(2*(nx-nu)+nu+1:2*(nx-nu)+2*nu);
%     umax = K{1}(2*(nx-nu)+nu+1:2*(nx-nu)+2*nu) + K{1}(2*(nx-nu)+1:2*(nx-nu)+nu);
%
%     % Create constraints object
%     constr = constraints(nx,nu,np,nd,ny,N);
%     % Set hard state constraints
%     constr = constr.setConstraints('x',xmin,xmax);
%     % Set hard input constraints
%     constr = constr.setConstraints('u',umin,umax);

else
    
    uref = object.options.uref;
%     constr = object.constr;
%     if nSoft ~= 0
%         xmin = [];
%         xmax = [];
%         for i=1:2:2*nx
%             xmin = [xmin; -K{1}(i)];
%             xmax = [xmax; K{1}(i+1)];
%         end
%         umin = -K{1}(2*nx+1:2*nx+nu);
%         umax = K{1}(2*nx+nu+1:2*nx+2*nu);
%
%         % Create constraints object
%         constr = constraints(nx,nu,np,nd,ny,N);
%         % Set hard state constraints
%         constr = constr.setConstraints('x',xmin,xmax);
%         % Set hard input constraints
%         constr = constr.setConstraints('u',umin,umax);
%     end

end

% Compute QP matrices and constraints
if strcmpi(object.ctrlvar,'state')

    % QP matrices
    R1 = zeros(Nu*nu);
    for i=1:nu:Nu*nu
        R1(i:i+nu-1,i:i+nu-1) = R;
    end
    Q1 = zeros((N+1)*nx);
    for i=1:nx:(N+1)*nx
        Q1(i:i+nx-1,i:i+nx-1) = Q;
    end
    Q1(end-nx+1:end,end-nx+1:end) = P;
    S1 = zeros((N+1)*nx, Nu*nu);
    QP1 = 2.*[R1, S1.'; S1, Q1];

    QP2 = zeros(Nu*nu+(N+1)*nx, 1);
    for i=1:nu:Nu*nu
        QP2(i:i+nu-1) = -2.*R*uref;
    end
    j=1;
    for i=Nu*nu+1:nx:Nu*nu+(N+1)*nx-nx
        QP2(i:i+nx-1) = -2.*Q*xref(:,j);
        j=j+1;
    end
    QP2(end-nx+1:end) = -2.*P*xref(:,j);

    % QP equality constraints
    eqConB = zeros((N+1)*nx, Nu*nu);
    for i=1:nu:Nu*nu
        if i==1
            j = nx+1;
        end
        eqConB(j:j+nx-1,i:i+nu-1) = -B;
        j = j + nx;
    end
    if ~object.isTracking
        for i=Nu*nx+1:nx:(N+1)*nx
            eqConB(i:i+nx-1,Nu*nu:Nu*nu+nu-1) = -B;
        end
    end

    eqConA = zeros((N+1)*nx);
    for i=1:nx:(Nu-1)*nx
        eqConA(i:i+nx-1,i:i+nx-1) = eye(nx);
        j = i+nx;
        eqConA(j:j+nx-1,i:i+nx-1) = -A;
    end
    if object.isTracking && object.options.constantInputAfterNu
        for i=(Nu-1)*nx+1:nx:N*nx
            eqConA(i:i+nx-1,i:i+nx-1) = eye(nx);
            j = i+nx;
            eqConA(j:j+nx-1,i:i+nx-1) = -A+[zeros(nx-nu,nx-nu),A(1:nx-nu,end); zeros(nu,nx)]-B*Ku;
        end
    else
        for i=(Nu-1)*nx+1:nx:N*nx
            eqConA(i:i+nx-1,i:i+nx-1) = eye(nx);
            j = i+nx;
            eqConA(j:j+nx-1,i:i+nx-1) = -A-B*Ku;
        end
    end
    eqConA(end-nx+1:end,end-nx+1:end) = eye(nx);

    eqL = [eqConB eqConA];

    eqR = zeros((N+1)*nx, 1);
%     eqR(1:nx) = A*x + Ex*p + Fx*d + Gx;
    eqR(1:nx) = x;
    for i=nx+1:nx:Nu*nx
        eqR(i:i+nx-1) = Ex*p + Fx*d + Gx;
    end
    for i=Nu*nx+1:nx:(N+1)*nx
        eqR(i:i+nx-1) = Ex*p + Fx*d + Gx + B*Ou;
    end

else
% ----- TO DO: OUTPUT REGULATION AND TRACKING IS NOT YET SUPPORTED -----
%     if ycon == 0
%         % QP matrices
%         R1 = zeros(N*nu);
%         for i=1:nu:(N-1)*nu
%             R1(i:i+nu-1,i:i+nu-1) = R + (D'*Q*D);
%         end
%         R1(end-nu+1:end,end-nu+1:end) = R + (D'*Q*D);
%         Q1 = zeros((N+1)*nx);
%         for i=1:nx:N*nx+1
%             Q1(i:i+nx-1,i:i+nx-1) = C'*Q*C;
%         end
%         Q1(end-nx+1:end,end-nx+1:end) = C'*P*C;
%         S1 = zeros((N+1)*nx, N*nu);
%         j = 1;
%         for i=1:nx:N*nx-1
%             S1(i:i+nx-1,j:j+nu-1) = C'*Q*D;
%             j = j + nu;
%         end
%         S1(end-nx+1:end,end-nu+1:end) = C'*P*D;
%         QP1 = [R1, S1.'; S1, Q1];
%
%         QP2 = zeros(N*(nu+nx)+nx, 1);
%         for i=1:nu:(N-1)*nu
%             QP2(i:i+nu-1) = -D'*Q*yref-R*uref;
%         end
%         QP2((N-1)*nu+1:N*nu) = -D'*P*yref-R*uref;
%         for i=N*nu+1:nx:N*(nu+nx)
%             QP2(i:i+nx-1) = -C'*Q*yref;
%         end
%         QP2(end-nx+1:end) = -C'*P*yref;
%
%         % QP inequality constraints
%         ineqL = zeros(2*(N*(nu+nx)+nx), N*(nu+nx)+nx);
%         ineqR = zeros(2*(N*(nu+nx)+nx), 1);
%
%         k = 0;
%         for i=1:N
%             for j=1:nu
%                 ineqL(2*k+j,k+j) = Hu{i}(2*nx+j,j);
%                 ineqL(2*k+j+nu,k+j) = Hu{i}(2*nx+j+nu,j);
%                 ineqR(2*k+j) = K{i}(2*nx+j,j);
%                 ineqR(2*k+j+nu) = K{i}(2*nx+j+nu,j);
%             end
%             k = k + nu;
%         end
%         for i=1:N+1
%             for j=1:nx
%                 ineqL(2*k+j,k+j) = Hx{i}(j,j);
%                 ineqL(2*k+j+nx,k+j) = Hx{i}(j+nx,j);
%                 ineqR(2*k+j) = K{i}(j);
%                 ineqR(2*k+j+nx) = K{i}(j+nx);
%             end
%             k = k + nx;
%         end
%
%         % QP equality constraints
%         eqConB = zeros((N+1)*nx, N*nu);
%         for i=1:nu:N*nu
%             if i==1
%                 j=i+nx;
%             end
%             eqConB(j:j+nx-1,i:i+nu-1) = B;
%             j = j + nx;
%         end
%
%         eqConA = zeros((N+1)*nx);
%         for i=1:nx:N*nx
%             eqConA(i:i+nx-1,i:i+nx-1) = -eye(nx);
%             j = i+nx;
%             eqConA(j:j+nx-1,i:i+nx-1) = A;
%         end
%         eqConA(end-nx+1:end,end-nx+1:end) = -eye(nx);
%
%         eqL = [eqConB eqConA];
%
%         eqR = zeros((N+1)*nx, 1);
%         eqR(1:nx) = -x;
%         for i=nx+1:nx:(N+1)*nx
%             eqR(i:i+nx-1) = -Gx;
%         end
%
%     else
%         % QP matrices
%         R1 = zeros(N*nu);
%         for i=1:nu:N*nu
%             R1(i:i+nu-1,i:i+nu-1) = R;
%         end
%         Q1 = zeros((N+1)*ny);
%         for i=1:ny:N*ny
%             Q1(i:i+ny-1,i:i+ny-1) = Q;
%         end
%         Q1(end-ny+1:end,end-ny+1:end) = P;
%         QP1 = [R1, zeros(N*nu, (N+1)*(nx+ny)); zeros((N+1)*nx, N*(nu+nx+ny)+nx+ny); zeros((N+1)*ny, N*(nu+nx)+nx), Q1];
%
%         QP2 = zeros(N*(nu+nx+ny)+nx+ny, 1);
%         for i=1:nu:N*nu
%             QP2(i:i+nu-1) = -R*uref;
%         end
%         for i=N*(nu+nx)+nx+1:ny:N*(nu+nx+ny)+nx+ny
%             QP2(i:i+ny-1) = -Q*yref;
%         end
%
%         ineqL = zeros(2*(N*(nu+nx+ny)+nx+ny), N*(nu+nx+ny)+nx+ny);
%         ineqR = zeros(2*(N*(nu+nx+ny)+nx+ny), 1);
%
%         % QP inequality constraints
%         k = 0;
%         for i=1:N
%             for j=1:nu
%                 ineqL(2*k+j,k+j) = Hu{i}(2*nx+j,j);
%                 ineqL(2*k+j+nu,k+j) = Hu{i}(2*nx+j+nu,j);
%                 ineqR(2*k+j) = K{i}(2*nx+j,j);
%                 ineqR(2*k+j+nu) = K{i}(2*nx+j+nu,j);
%             end
%             k = k + nu;
%         end
%         for i=1:N+1
%             for j=1:nx
%                 ineqL(2*k+j,k+j) = Hx{i}(j,j);
%                 ineqL(2*k+j+nx,k+j) = Hx{i}(j+nx,j);
%                 ineqR(2*k+j) = K{i}(j);
%                 ineqR(2*k+j+nx) = K{i}(j+nx);
%             end
%             k = k + nx;
%         end
%         for i=1:N
%             for j=1:ny
%                 ineqL(2*k+j,k+j) = Hy{i}(2*(nx+nu)+j,j);
%                 ineqL(2*k+j+ny,k+j) = Hy{i}(2*(nx+nu)+j+ny,j);
%                 ineqR(2*k+j) = K{i}(2*(nx+nu)+j);
%                 ineqR(2*k+j+ny) = K{i}(2*(nx+nu)+j+ny);
%             end
%             k = k + ny;
%         end
%         ineqL(2*k+1:2*k+ny,k+1:k+ny) = Hy{N+1}(2*nx+1:2*nx+ny,1:ny);
%         ineqL(2*k+ny+1:2*k+2*ny,k+1:k+ny) = Hy{N+1}(2*nx+ny+1:2*nx+2*ny,1:ny);
%         ineqR(2*k+1:2*k+ny) = K{N+1}(2*nx+1:2*nx+ny);
%         ineqR(2*k+ny+1:2*k+2*ny) = K{N+1}(2*nx+ny+1:2*nx+2*ny);
%
%         % QP equality constraints
%         eqConB = zeros((N+1)*nx, N*nu);
%         for i=1:nu:N*nu
%             if i==1
%                 j=i+nx;
%             end
%             eqConB(j:j+nx-1,i:i+nu-1) = B;
%             j = j + nx;
%         end
%
%         eqConA = zeros((N+1)*nx);
%         for i=1:nx:N*nx
%             eqConA(i:i+nx-1,i:i+nx-1) = -eye(nx);
%             j = i+nx;
%             eqConA(j:j+nx-1,i:i+nx-1) = A;
%         end
%         eqConA(end-nx+1:end,end-nx+1:end) = -eye(nx);
%
%         eqConD = zeros((N+1)*ny, N*nu);
%         j = 1;
%         for i=1:nu:N*nu
%             eqConD(j:j+ny-1,i:i+nu-1) = D;
%             j = j + ny;
%         end
%
%         eqConC = zeros((N+1)*ny, (N+1)*nx);
%         j = 1;
%         for i=1:nx:(N+1)*nx
%             eqConC(j:j+ny-1,i:i+nx-1) = C;
%             j = j + ny;
%         end
%
%         eqL = [eqConB eqConA zeros((N+1)*nx, (N+1)*ny); eqConD eqConC -eye((N+1)*ny)];
%
%         eqR = zeros((N+1)*(nx+ny), 1);
%         eqR(1:nx) = -x;
%         for i=nx+1:nx:(N+1)*nx
%             eqR(i:i+nx-1) = -Gx;
%         end
%         for i=(N+1)*nx+1:ny:(N+1)*(nx+ny)
%             eqR(i:i+ny-1) = -Gy;
%         end
%     end
end

% ----- TO DO: ADD SOFT CONSTRAINTS ALSO FOR ADMM -----
% In case of soft constraints, extend the state of QP
if nSoft ~= 0 && algo ~= "ADMM"
    rho = object.options.rho;
    Qsoft = zeros((N+1)*nSoft);
    for i=1:nSoft:(N+1)*nSoft
        Qsoft(i:i+nSoft-1,i:i+nSoft-1) = rho;
    end
    QP1 = [QP1, zeros(size(QP1,1),(N+1)*nSoft); zeros((N+1)*nSoft,size(QP1,2)), Qsoft];

    QP2 = [QP2; zeros((N+1)*nSoft,1)];

    eqL = [eqL, zeros(size(eqL,1),(N+1)*nSoft)];

    ineqLlength = size(ineqL, 1);
    ineqL = [ineqL, zeros(size(ineqL,1), (N+1)*nSoft)];
    for i = 1:(N+1)
        ineqLsoft = zeros(size(H{i}, 1), nSoft);
        for j = 1:size(H{i}, 1)
            if ismember(H{i}(j, :), Hs{i}, 'rows')
                ineqLsoft(j, j) = 1;
            end
        end
        ineqL((i-1)*size(H{i}, 1) + 1 : (i-1)*size(H{i}, 1) + size(H{i}, 1), ineqLlength + (i-1)*size(Hs{i}, 1) + 1 : ineqLlength + (i-1)*size(Hs{i}, 1) + nSoft) = ineqLsoft;
    end

end

end
