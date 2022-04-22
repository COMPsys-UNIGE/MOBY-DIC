%% IDENTIFY
% Identifies the virtual sensor starting from the measured system inputs
% and outputs
%
% SYNTAX
%
% object = identify(object,u,y,z,np)
%
% u is a [npoints x nu] matrix whose i-th row corresponds to the i-th
% measurement of the nu inputs of the system. The same is for y (for the
% measurable outputs) and z (for the unmeasurable outputs).
% u, y and z can also be cell arrays whose i-th element correspond to a
% dataset obtained iterating the system starting from different initial
% conditions.
% np is an array containing the number of subdivisions per dimension of the
% simplicial partition. The domain of the pwas function is estimated from
% data.
% ut, yt and zt (optional) are matrices or cell arrays just as u, y
% and z, but they contain a test dataset which can be used to estimate
% the best value of the regularization parameter lambda. If they are not
% provided the optimal value for lambda is computed with a Generalized
% Cross Validation technique.
%
% object = identify(object,u,y,z,opts)
%
% The domain is provided from outside the function. D is a matrix in the form:
%  _                                  _
% |  x_min(1) x_min(2) ... x_min(nx)   |
% |_  x_max(1) x_max(2) ... x_max(nx) _|
%
% object = identify(object,u,y,z,np,[ut],[yt],[zt],D,options)
%
% options is a structure with the following fields:
%
% * order: Tikhonov regularization order. It can be 0 or 1, default 0.
% * lambda: regularisation weight. It must be a scalar > 0. If it is not
%            provided it is estimated inside the function.
% * nsplits: number of splits of the dataset used to solve the regularised
%             least square problem with an iterative approach (in order to
%             save memory). In practice, instead of solving $|| H w - F ||^2$,
%             you can solve || H_1 w - F_1 ||^2 + || H_2 w - F_2 ||^2 + ... +
%             || H_n w - F_n ||^2, in which H = [H_1 H_2 ... H_n]' and F = [F_1
%             F_2 ... F_n]'. nsplits corresponds to mz.
%             A low value of nsplits makes ridgeRegression faster but it
%             can result in memroy occupation problems.
%             Default value, 1.
% * constraints: only if reducedComplexity = 0.
%                 array of structures defining optional constraints on the
%                 pwas function you want to obtain. Each element of the
%                 array defines different constraints. Tre structures have
%                 the following fields:
%
%                  - type: it can be either 'bounds' or 'equality'. If it
%                          is 'bounds' constraints in the form
%                          lbound <= f(x) <= ubound are imposed. lbound and
%                          ubound can be provided through fields lbound and
%                          ubound. If they are not provided they are
%                          imposed as the minimum and maximum values of the
%                          data Y and Yt contained in the dataset.
%                          If type is equality, equality constraints in the
%                          form f(x) = k can be imposed for x lying on a
%                          hyper-plane parallel to the domain components.
%
%                  - lbound: lower bound for constraints of type 'bounds'
%
%                  - ubound: upper bound for constraints of type 'bounds'
%
%                  - variables: needed for equality constraints. It is an
%                               array of strings indicating for which set of
%                               points you want to impose the constraints. For
%                               example if you want to impose a constraint
%                               on f, for x1 = 3 and x3 = 2, variables must
%                               be ['x1 = 3', 'x3 = 2'].
%
%                  - value: needed for equality constraints. It is the
%                           value you want the function to assume in
%                           correspondence of the points defined in field
%                           variables. To impose the constraint f(x) = 0,
%                           for x2 = 1 and x4 = 5 you gave to set
%                           variables = ['x2 = 1', 'x4 = 5'] and value = 0.
% * gamma: domain expansion parameter. If the domain D is not provided from
%           outside the function, it is automatically extrapolated from data
%           X. The tight domain (Dt) contains exactly data X. Such domain is
%           grown of a factor gamma such as D = gamma Dt. Default value
%           1.1.
% * solver: string specifying the solver you want to use to solve the QP
%            problem (this is necessary only for constrained ridge
%            regression, otherwise the solution is analytical). Possible
%            choices are 'quadprog' (default), 'cvx', 'cplex', 'yalmip'
%            'clp'.
% * verbose: if verbose is set to 1, messages are displayed indicating the
%            status of the ridge regression process.
%
% ACKNOWLEDGEMENTS
%
% Contributors:
%
% * Tomaso Poggi (tpoggi@essbilbao.org)
%
% Copyright is with:
%
% * Copyright (C) 2010 University of Genoa, Italy.

% -------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
%
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the
%          Free Software Foundation, Inc.,
%          59 Temple Place, Suite 330,
%          Boston, MA  02111-1307  USA
%
% -------------------------------------------------------------------------

function object = identify(object,u,y,z,opts)

% Check input arguments
if nargin ~= 5
    error('Wrong input arguments');
end

options = VSset(object,opts);

% Check on input datasets

if iscell(z)
    % Number of datasets
    n_x0 = numel(z);
    
    % If u is not provided, create an empty cell array
    if isempty(u)
        for i = 1:n_x0
            u{i} = [];
        end
    end
    
    % Check dataset consistency
    if ~iscell(y) || ~iscell(u)
        error('Dataset not consistent. u, y and z must be all cell arrays.');
    end
    
    if numel(y) ~= n_x0 || numel(u) ~= n_x0
        error('Dataset not consistent. u, y and z must be cell arrays with the same number of elements.');
    end
    
    for i = 1:n_x0
        % Check if all elements of cell arrays have the same number of
        % columns
        if size(u{i},2) ~= object.nu || size(y{i},2) ~= object.ny || size(z{i},2) ~= object.nz
            error('All components of u, y and z cell arrays must have the same number of columns');
        end
        
        % Check if corresponding elements of cell arrays have the same
        % number of rows
        if isempty(u{i})
            if size(y{i},1) ~= size(z{i},1)
                error('Dataset not consistent. y and z must have the same number of rows.');
            end
        else
            if size(y{i},1) ~= size(z{i},1) || size(y{i},1) ~= size(u{i},1)
                error('Dataset not consistent. u, y and z must have the same number of rows.');
            end
        end
    end
    
    
else
    
    % Check dataset consistency
    if isempty(u)
        if size(y,1) ~= size(z,1)
            error('Dataset not consistent. y and z must have the same number of rows.');
        end
    else
        if size(y,1) ~= size(z,1) || size(y,1) ~= size(u,1)
            error('Dataset not consistent. u, y and z must have the same number of rows.');
        end
    end
    
    if size(z,2) ~= 1
        error('z must have one column');
    end
    
    % Put dataset in a cell array
    uu{1} = u;
    yy{1} = y;
    zz{1} = z;
    u = uu;
    y = yy;
    z = zz;
    
    clear uu yy zz
    
    n_x0 = 1;
end

% Time windows
mu = object.mu;
my = object.my;
mz = object.mz;
np = options.np;

% Maximum time window
l = max([max(mz)+1,max(mu),max(my)]);

% Number of dimensions of the space
Ndim = mz+sum(my)+sum(mu);

if Ndim == 0
    error('Number of dimensions is 0. Increase MU, MY or MZ.');
end

% Reply np for all dimensions of the domain
if numel(np) == 1
    np = repmat(np,1,Ndim);
end

% Matrices who will contain the dataset
X = cell(n_x0,1);
F = cell(n_x0,1);
D = cell(n_x0,1);

% This is necessary for reduced complexity virtual sensors
% The i-th element of cell array indices will contain the indices of the
% columns of X corresponding to time instant k-i+1 (k-i if current time
% instant is not considered)
indices = cell(1,l);

% Initialize indices with empty elements
for i = 1:l
    indices{i} = [];
end

umin = zeros(n_x0,object.nu);
umax = zeros(n_x0,object.nu);
ymin = zeros(n_x0,object.ny);
ymax = zeros(n_x0,object.ny);
zmin = zeros(n_x0,object.nz);
zmax = zeros(n_x0,object.nz);

% Loop on the number of datasets
for i = 1:n_x0
    
    % Initialize i-th element of cell array
    X{i} = zeros(size(y{i},1)-l+1,Ndim);
    F{i} = zeros(size(y{i},1)-l+1,1);
    D{i} = zeros(2,Ndim);
    
    % Indices of the starting column
    kstart = 1;
    
    % Loop on the number of system inputs
    for j = 1:object.nu
        
        % Index used to indicize cell array indices
        h = 1;
        
        % Loop on the number of past time instants to consider for the
        % current input
        for k = kstart:kstart+mu(j)-1
            
            kk = k-kstart+1;
            % Fill matrix X with input data in the current time instant
            % This is u_j at instant t-kk+1
            X{i}(:,k) = u{i}(l-kk+1:end-kk+1,j);
            
            % Update indices
            if i == 1
                indices{h} = [indices{h} k];
                h = h+1;
            end
            
        end
        
        umin(i,j) = min(min(X{i}(:,kstart:kstart+mu(j)-1)));
        umax(i,j) = max(max(X{i}(:,kstart:kstart+mu(j)-1)));
        
        if isempty(options.domain)
            D{i}(:,kstart:kstart+mu(j)-1) = ...
                [min(X{i}(:,kstart:kstart+mu(j)-1));max(X{i}(:,kstart:kstart+mu(j)-1))];
        else
            D{i}(:,kstart:kstart+mu(j)-1) = ...
                [repmat(options.domain.umin(j),1,mu(j));repmat(options.domain.umax(j),1,mu(j))];
        end
        
        if isempty(k)
            k = kstart-1;
        end
        
        % Update index of the starting column
        kstart = k+1;
        
    end
    
    % Loop on the number of system outputs
    for j = 1:object.ny
        
        % Index used to indicize cell array indices
        h = 1;
        
        % Loop on the number of past time instants to consider for the
        % current output
        for k = kstart:kstart+my(j)-1
            
            kk = k-kstart+1;
            % Fill matrix X with output data in the current time instant
            % This is y_j at instant t-kk+1
            X{i}(:,k) = y{i}(l-kk+1:end-kk+1,j);
            
            % Update indices
            if i == 1
                indices{h} = [indices{h} k];
                h = h+1;
            end
            
        end
        
        ymin(i,j) = min(min(X{i}(:,kstart:kstart+my(j)-1)));
        ymax(i,j) = max(max(X{i}(:,kstart:kstart+my(j)-1)));
        
        if isempty(options.domain)
            D{i}(:,kstart:kstart+my(j)-1) = ...
                [min(X{i}(:,kstart:kstart+my(j)-1));max(X{i}(:,kstart:kstart+my(j)-1))];
        else
            D{i}(:,kstart:kstart+my(j)-1) = ...
                [repmat(options.domain.ymin(j),1,my(j));repmat(options.domain.ymax(j),1,my(j))];
        end
        
        if isempty(k)
            k = kstart-1;
        end
        
        % Update index of the starting column
        kstart = k+1;
        
    end
    
    % Index used to indicize cell array indices, it varies according to
    % 'current'
    if object.current
        h = 2;
    else
        h = 1;
    end
    
    % Loop on the number of past time instants to consider for the
    % unmeasured output
    for k = kstart:kstart+mz-1
        
        kk = k-kstart+1;
        % Fill matrix X with unmeasurable output data in the previous time instant
        % This is z_j at instant t-kk
        X{i}(:,k) = z{i}(l-kk:end-kk);
        
        % Update indices
        if i == 1
            indices{h} = [indices{h} k];
            h = h+1;
        end
        
    end   
    
    % Fill matrix F with the unmeasurable output
    F{i} = z{i}(l:end);
    
    % If you do not want to use the data of u and y at the current instant,
    % you have to slightly arrange the matrices of the dataset
    if ~object.current
        
        X{i}(1:end-1,end-mz+1:end) = X{i}(2:end,end-mz+1:end);
        X{i}(end,:) = [];
        
        F{i}(1:end-1) = F{i}(2:end);
        F{i}(end) = [];
        
    end
    
    if mz > 0
        zmin(i) = min(min(X{i}(:,kstart:kstart+mz-1)));
        zmax(i) = max(max(X{i}(:,kstart:kstart+mz-1)));
        
        if isempty(options.domain)
            D{i}(:,kstart:kstart+mz-1) = ...
                [min(X{i}(:,kstart:kstart+mz-1));max(X{i}(:,kstart:kstart+mz-1))];
        else
            D{i}(:,kstart:kstart+mz-1) = ...
                [repmat(options.domain.zmin,1,mz);repmat(options.domain.zmax,1,mz)];
        end
        
    else
        zmin(i) = min(min(F{i}));
        zmax(i) = max(max(F{i}));
    end
    
end

if isempty(options.domain)
    
    umin = min(umin);
    umax = max(umax);
    ymin = min(ymin);
    ymax = max(ymax);
    zmin = min(zmin);
    zmax = max(zmax);
    
    umean = (umin+umax)/2;
    ymean = (ymin+ymax)/2;
    zmean = (zmin+zmax)/2;
    
    if object.mu > 0
        object.domain.umin = umean+options.gamma*(umin-umean);
        object.domain.umax = umean+options.gamma*(umax-umean);
    end
    if object.my > 0
        object.domain.ymin = ymean+options.gamma*(ymin-ymean);
        object.domain.ymax = ymean+options.gamma*(ymax-ymean);
    end
    object.domain.zmin = zmean+options.gamma*(zmin-zmean);
    object.domain.zmax = zmean+options.gamma*(zmax-zmean);
    
else
    
    object.domain = options.domain;
    
end

% Put all the datasets together in a unique matrix
X = cell2mat(X(:));
F = cell2mat(F(:));
D = cell2mat(D(:));
D = [min(D);max(D)];

if isempty(options.domain)
    Dm = mean(D);
    D = [Dm+options.gamma*(D(1,:)-Dm);Dm+options.gamma*(D(2,:)-Dm)];
end

% Do the same procedure for the test set, if provided

% If a reduced complexity virtual sensor is used, the partition must be
% specified as a cell array
if object.reducedComplexity
    
    % PP will contain the partition
    PP = cell(1,l);
    
    for i = 1:l
        
        PP{i} = np(indices{i});
        
    end
    
    np = PP;
    clear PP
    
end

% Check if the number of free paramerters is higher than the number of
% equations
if object.reducedComplexity
    % Number of free parameters
    npars = 0;
    for i = 1:l
        npars = npars + prod(np{i}+1);
    end
else
    npars = prod(np+1);
end

if (npars >= length(X))
    display(['Number of partitons/number of data = ', num2str(npars/length(X))])
    disp('WARNING: too many partitions (basis functions) compared to number of data')
end


% pwas regression: Regularized Least Squares

if object.reducedComplexity
    
    % If a reduced complexity virtual sensor is used, you need to split the
    % dataset into cell arrays, one for each time instant. The indices of the
    % columns of X corresponding to a given time instant are contained in
    % 'indices'.
    
    % XX will contain the dataset
    XX = cell(1,l);
    DD = cell(1,l);
    
    for i = 1:l
        
        XX{i} = X(:,indices{i});
        DD{i} = D(:,indices{i});
        
    end
    
    X = XX;
    D = DD;
    clear XX DD
    
    if isempty(X{l})
        X(l) = [];
    end
    
    if isempty(D{l})
        D(l) = [];
    end
    
    if isempty(np{l})
        np(l) = [];
    end
    
    options.np = np;
    options.D = D;
    
    [object.fpwas, object.info] = hierarchicalRidgeRegression(X,F,options);
    
else
    
    options.D = D;
    [object.fpwas, object.info] = ridgeRegression(X,F,options);
    
end

object.identified = 1;
object.Ts = options.Ts;

return;
