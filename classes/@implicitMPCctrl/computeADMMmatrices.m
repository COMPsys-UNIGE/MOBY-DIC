function [M1, M2, M3, M4, v1, v2, v3] = computeADMMmatrices(object, QP1, QP2, ineqL, ineqR, eqL, eqR, regPar)

% computeADMMmatrices
% Computes the matrices for the ADMM algorithm, given the matrices of the
% quadratic programming (QP) problem required by MPC every time step (QP1,
% QP2, ineqL, eqL, eqR) and the regularization paramater of ADMM (regPar).
%
% Given a quadratic programming (QP) problem:
%     min_x [(1/2)*x'*QP1*x + x'*QP2]
%      s.t. ineqL*x <= ineqR
%           eqL*x = eqR,
% this function computes the matrices required by ADMM.
%
% OBJ is the implicitMPCctrl object, regPar is a regularization parameter
% for ADMM algorithm.

% Number of variables
nx = object.getNumberOfStates;
np = object.getNumberOfParameters;
nd = object.getNumberOfUnmeasurableInputs;
nu = object.getNumberOfInputs;
N = object.options.N;
Nu = object.options.Nu;

A = [ineqL; eqL];
B = [eye(size(ineqL, 1)); zeros(size(eqL, 1), size(ineqL, 1))];
c = [ineqR; eqR];

% ----- TO DO: ADD SOFT CONSTRAINTS FOR ADMM -----
% % Get indices of states with soft constraints
% [Hs,~] = object.constr.getAllConstraints('soft');
% softIdx = [];
% for i=1:size(Hs{1},2)
%     if any(Hs{1}(:,i)~=0)
%         softIdx = [softIdx, i];
%     end
% end
% nSoft = length(softIdx);

KKT = -inv(QP1 + (regPar*(A'*A)));
M1 = KKT;
M2 = regPar*M1*(A');
M3 = -ineqL;
M4 = A;

v1 = M1*QP2 - M2*c;
v2 = c;
v3 = ineqR;

end

