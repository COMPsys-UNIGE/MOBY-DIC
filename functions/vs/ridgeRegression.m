%% RIDGE REGRESSION
% Performs a PWAS constrained ridge regression on a set of data.
%
% SYNTAX
%
% [fpwas info] = ridgeRegression(X,Y,P,[Xt],[Yt])
%
% X must be a [ndata x ndim] matrix and Y a [ndata x 1] array. P
% defines the simplicial partition you want the fpwas function to be
% defined on. If P is a scalar, each dimension of the domain is subdivided
% into P intervals. If it is a vector, you specify individually the number
% of subdivisions per dimension. P can also be a cell array whose i-th
% element contains the i-th component of the vertices of the simplicial
% partition (for non-uniform partition). The domain of the pwas function is
% automatically extrapolated from input data X and Xt.
% Xt and Yt (optional) represent a test dataset used to estimate the
% optimal value for the Tikhonov parameter lambda. lambda is chosen in
% order to minimize || f(Xt) - Yt ||^2 being f the pwas function obtained
% starting from the training dataset (X,Y). If Xt and Yt are not provided,
% lambda is chosen with a GCV approach. Xt must be a [ndatatest x ndim]
% matrix and Yt a [ndatatest x 1] array.
%
% fpwas is a pwas object defining the pwas function obtained after the
% regression.
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
% *[fpwas info] = ridgeRegression(X,Y,P,[Xt],[Yt],D)*
%
% As above, but the domain of the pwas function is passed from outside the
% function. D is a matrix in the form:
%  _                                 _
% |  x_min(1) x_min(2) ... x_min(nx)  |
% |_ x_max(1) x_max(2) ... x_max(nx) _|
%
% [fpwas info] = ridgeRegression(X,Y,P,[Xt],[Yt],options)
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
% * constraints: array of structures defining optional constraints on the
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
% *[fpwas info] = ridgeRegression(X,Y,P,[Xt],[Yt],D,options)*
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

function [fpwas, info] = ridgeRegression(X,Y,opts)

% Check input arguments
if nargin ~= 3
    error('Wrong input arguments');
end

% Extract parameters
P = opts.np;
D = opts.D;
lambda = opts.lambda;
order = opts.order;

% Check on dataset consistency

if (size(X,1) ~= size(Y,1))
    error('The number of data X must be equal to the number of samples Y')
end


% Number of dimensions of the domain
nd = size(X,2);

% If P is a scalar, transform it into a vector of nd equal elements
if numel(P) == 1
    P = repmat(P,1,nd);
end

% Check if the partition is compatible with the domain
if numel(P) ~= nd
    error(['Partition P must have ',num2str(nd),' elements']);
end

if any(D(1,:) > min(X)) || any(D(2,:) < max(X))
    error('The domain is too small! Some data in the dataset are outside.');
end

% Check on non uniform partition
if iscell(P)
    for i=1:numel(P)
        if P{i}(1) ~= D(1,i) || P{i}(end) ~= D(2,i)
            error(['The first and last vertices of the non-uniform partition',...
                ' must coincide with the domain boundaries']);
        end
    end
end

% Create pwas object containing only the definition of the domain and of
% the partition
domain.xmin = D(1,:);
domain.xmax = D(2,:);
fpwas = pwasFunction(domain,P);

% Number of vertices of the simplicial partition
Nb = fpwas.getNumberOfVertices();

% Set constraints
disp('Computing optimization matrices...');

% Create 0 or 1 order Tikhonov regularization matrix
if order == 1
    Q = QMatrix(D,P);
else
    Q = speye(Nb);
end

% Compute matrices
d = sparse(Y);

% Compute matrices for all other data chunks
G = fpwas.alphaBasis(X,Nb*size(X,1));
    
%     notrem=unique([notrem find(sum(alpha)~=0)]);
H = G'*G;
y = G'*d;

disp('Optimization matrices computed.');

% Add Tikhonov regularization matrix
HH = H+lambda*Q;

disp('Optimization started...');


% If there are no constraints, solve least squares regularized problem
w = full(HH\y);

% Fill structure info
if nargout > 1
    info.wnorm = full(sqrt(w'*Q*w));
    info.residual = norm(G*w-d);
end

disp('Optimization terminated.');
disp('Ridge regression performed succesfully.');

% Setting weights to pwas object
%     fpwas = fpwas.setWeights(w);

fpwas = fpwas.setWeights(w);

return;

