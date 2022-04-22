function zOpt = admm(object, M1, M2, M3, M4, v1, v2, v3, p, maxIter)
%ADMM
% Given a Qp defined by the matrices M1, M2, M3, M4, v1, v2, v3 and
% given the current system states x, this function runs the alternating
% direction method of multipliers (ADMM) algorithm.
% zOpt is the optimal primal variable.
% p is a regularization parameter and maxIter is the number of iterations
% to perform the algorithm.
% 
% This is a private method.

% Nuber of state elements
nx = object.getInformation.sys.getNumberOfStates;

% Number of inputs
nu = object.getInformation.sys.getNumberOfInputs;

% Prediction horizion
N = object.getInformation.options.N;

% Control horizion
Nu = object.getInformation.options.Nu;

% warm start
persistent z;
persistent y;
persistent lambda;
if isempty(y)
    z = zeros(size(v1));
    y = zeros(size(v3));
    lambda = zeros(size(v2));
else
    if object.getInformation.options.tracking
        z = [z(nu+1:Nu*nu); zeros(nu,1); z(Nu*nu+nx+nu+1:end); zeros(nx+nu,1)];
    else
        z = [z(nu+1:Nu*nu); zeros(nu,1); z(Nu*nu+nx+1:end); zeros(nx,1)];
    end
end

maxV = [];

for i=1:maxIter
    
    % primal variable update        
    z = v1 + M2*([y;zeros(numel(lambda)-numel(y), 1)] + lambda);
    
    % dual variable update
    y = M3*z - lambda(1:numel(v3)) + v3;
    for j=1:numel(y)
        if y(j) < 0
            y(j) = 0;
        end
    end
    
    % update Lagrange multiplier
    lambda = lambda + M4*z + [y;zeros(numel(lambda)-numel(y), 1)] - v2;
    
%     maxV = [maxV; max(abs([y;z;lambda;M3*z;M4*z]))];
    
end

% maxInt = fix(max([maxV; abs(z); abs(y); abs(lambda); abs(v1); abs(v2); abs(v3)]))

zOpt = z;

end

