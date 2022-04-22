%% HIERARCHICAL RIDGE REGRESSION
% Performs a hierarchical PWAS constrained ridge regression on a set of
% data
%
% SYNTAX
%
% [fpwas info] = hierarchicalRidgeRegression(X,Y,P)
%
% X must be a cell array of matrices and Y a cell array of arrays.
% The number of cell arrays will coincide with the number of pwas function
% to be added. P defines the simplicial partition you want the fpwas
% functions to be defined on. If P is a scalar, each dimension of the
% domain of all pwas functions is subdivided into P intervals. If it is a
% cell array, you specify individually the number of subdivisions per
% dimension for each pwas function. Each element of cell array P can be an
% array or a cell array itself. The domain of the pwas function is
% automatically extrapolated from input data X and Xt.
% Xt and Yt (optional) represent a test dataset used to estimate the
% optimal value for the Tikhonov parameter lambda. lambda is chosen in
% order to minimize || f(Xt) - Yt ||^2 being f the pwas function obtained
% starting from the training dataset (X,Y). If Xt and Yt are not provided,
% lambda is chosen with a GCV approach. Xt must be a [ndatatest x ndim]
% matrix and Yt a [ndatatest x 1] array.
%
% fpwas is an array of pwas objects defining the pwas functions obtained
% after the regression which must be added to fit data Y.
%
% info is a struct with the following fields:
%
% * wnorm: it is the 2 norm of the weights vector w (for 0 order Tikhonov
%           regularization) or the 2 norm of L w (for 1 order Tikhonov
%           regularization)
% * residual: it is the 2 norm of the residual G w - d ( i.e. f(X) - Y )
% * test_error: it is the 2 norm of the residual computed in the test set
%                (if provided), i.e. f(Xt) - Yt
%
% *[fpwas info] = hierarchicalRidgeRegression(X,Y,P,D)*
%
% As above, but the domain of the pwas functions is passed from outside the
% function. D is a cell array of matrices in the form:
%  _                                 _
% |  x_min(1) x_min(2) ... x_min(nx)  |
% |_ x_max(1) x_max(2) ... x_max(nx) _|
%
% Each element of the cell array is related to a pwas function.
%
% [fpwas info] = hierarchicalRidgeRegression(X,Y,P,options)
%
% options is a structure with the following fields:
%
% * order: Tikhonov regularization order. It can be 0 or 1, default 0.
% * lambda: regularization weight. It must be a scalar > 0. If it is not
%            provided it is estimated inside the function.
% * nsplits: number of splits of the dataset used to solve the regularised
%             least square problem with an iterative approach (in order to
%             save memory). In practice, instead of solving || G w - d ||^2,
%             you can solve || G_1 w - d_1 ||^2 + || G_2 w - d_2 ||^2 + ... +
%             || G_n w - d_n ||^2, in which G = [G_1 G_2 ... G_n]' and d = [d_1
%             d_2 ... d_n]. nsplits corresponds to n.
%             A low value of nsplits makes ridgeRegression faster but it
%             can result in memroy occupation problems.
%             Default value, 1.
%
% * gamma: domain expansion parameter. If the domain D is not provided from
%           outside the function, it is automatically extrapolated from data
%           X (and Xt). The tight domain (Dt) contains exactly data X.
%           Such domain is increased of a factor gamma such as D = gamma Dt.
%           Default value 1.1.
% * solver: string specifying the solver you want to use to solve the QP
%            problem (this is necessary only for constrained ridge
%            regression, otherwise the solution is analytical). Possible
%            choices are 'quadprog' (default), 'cvx', 'cplex', 'yalmip'
%            'clp'.
% * verbose: if verbose is set to 1, messages are displayed indicating the
%            status of the ridge regression process.
%
% *[fpwas info] = hierarchicalRidgeRegression(X,Y,P,D,options)*
%
% All fields explained above are specified.
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

function [fpwas, info] = hierarchicalRidgeRegression(X,Y,opts)

% Check input arguments
if nargin ~= 3
    error('Wrong input arguments');
end
   
% Extract parameters
P = opts.np;
gamma = opts.gamma;
lambda = opts.lambda;
order = opts.order;

% Check on dataset consistency
if iscell(X)
    npwas = numel(X);
    for i = 1:npwas
        if (size(X{i},1) ~= size(Y,1))
            error('The number of data X must be equal to the number of samples Y')
        end
    end
else
    error('X must be a cell array');
end

% Number of dimensions of the domain
nd = zeros(1,npwas);

for i = 1:npwas
    nd(i) = size(X{i},2);
end

% Check partition P
if iscell(P)
    if numel(P) ~= npwas
        error(['P must be a scalar or a cell array with ',num2str(npwas),' elements']);
    end
    for i = 1:npwas
        if ~iscell(P{i}) && numel(P{i}) == 1
            P{i} = repmat(P{i},1,nd(i));
        end
        if ~iscell(P{i}) && numel(P{i}) ~= nd(i)
            error('Elements of P must be cell arrays, scalars or arrays with as many elements as the number of dimensions');
        end
    end
else
    if numel(P) == 1
        tmp = P;
        P = cell(1,npwas);
        for i = 1:npwas
            P{i} = repmat(tmp,1,nd(i));
        end
    else
        error(['P must be a scalar or a cell array with ',num2str(npwas),' elements']);
    end
end

% TO DO
% Custom domain
D = cell(npwas,1);

% % Check if the provided domain is compatible with the dataset
% if ~isempty(D)
%     if ~iscell(D)
%         error(['D must be a cell array']);
%     end
%     
%     if numel(D) ~= npwas
%         error(['D must be a cell array with ',num2str(npwas),' elements']);
%     end
%     
%     for i = 1:npwas
%         if ~isempty(D{i})
%             if size(D{i},1) ~= 2
%                 error('Domain D must have 2 rows');
%             end
%             
%             if size(D{i},2) ~= nd(i)
%                 error(['Domain D must have ',num2str(nd(i)),' columns']);
%             end
%         end
%     end
% else
%     D = cell(1,npwas);
% end


% pwas function domain identification, if necessary

for i = 1:npwas
    if isempty(D{i})
        Xmin = min(X{i});
        Xmax = max(X{i});
            
        % domain that contains the data exactly
        D{i} = [Xmin;Xmax];
        
        % the domain is enlarged by a factor of gamma with respect to its center
        % Dnew = gamma*D;
        Dm = mean(D{i});
        D{i} = [Dm+gamma*(D{i}(1,:)-Dm);Dm+gamma*(D{i}(2,:)-Dm)];
    end
    if any(D{i}(1,:) > min(X{i})) || any(D{i}(2,:) < max(X{i}))
        error('The domain is too small! Some data in the dataset are outside.');
    end
    
    % Check on non uniform partition
    if iscell(P{i})
        for j=1:numel(P{i})
            if P{i}{j}(1) ~= D{i}(1,j) || P{i}{j}(end) ~= D{i}(2,j)
                error(['The first and last vertices of the non-uniform partition',...
                    ' must coincide with the domain boundaries']);
            end
        end
    end
    
end

% Create pwas objects containing only the definition of the domain and of
% the partition
for i = 1:npwas
    domain.xmin = D{i}(1,:);
    domain.xmax = D{i}(2,:);
    fpwas(i) = pwasFunction(domain,P{i});
    % Number of vertices of the simplicial partition
    Nb(i) = fpwas(i).getNumberOfVertices();
end

disp('Computing optimization matrices...');

Q = [];
for i = 1:npwas
    % Create 0 or 1 order Tikhonov regularization matrix
    if order == 1
        Q = blkdiag(Q,QMatrix(D{i},P{i}));
    else
        Q = blkdiag(Q,speye(Nb(i)));
    end        
end

% Compute matrices
d = sparse(Y);

% Compute matrices for all other data chunks
    
G = [];
for i = 1:npwas
    Gtmp = fpwas(i).alphaBasis(X{i},Nb(i)*size(X{i},1));
    G = [G Gtmp];
end

%     notrem=unique([notrem find(sum(alpha)~=0)]);
H = G'*G;
y = G'*d;
    
disp('Optimization matrices computed.');

% Add Tikhonov regularization matrix
H = H+lambda*Q;

disp('Optimization started...');

% if isempty(constraints)
% If there are no constraints, solve least squares regularized problem

w = full(H\y);

for i = 1:npwas
    if i == 1
        idx = 1;
    else
        idx = idx+Nb(i-1);
    end
    fpwas(i) = fpwas(i).setWeights(w(idx:idx+Nb(i)-1));
end

% Fill structure info
if nargout > 1
    info.wnorm = full(sqrt(w'*Q*w));
    info.residual = norm(G*w-d);
end

disp('Optimization terminated.');
disp('Ridge regression performed succesfully.');

return;

