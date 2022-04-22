function MOBYDICpars = getMOBYDICpars()

% Tolerance used for equalities. Two numbers which differ for less than
% this tolerance are considered to be equal.
MOBYDICpars.roundtol = 1e-6;

% Tolerance used for polyhedra. A polyhedron is considered to be empty if
% its Chebyshev radius is below this threshold.
MOBYDICpars.emptytol = 1e-6;

% Solver used to solve qp problems in the pwas approximation. Possible
% choices are 'quadprog' or 'mpt'
MOBYDICpars.qpsolver = 'quadprog';

% Coefficient used for the Tikhonov regularization in the PWAS
% approximation. A too high value compromises the approximation accuracy, a
% too low value eliminates the smoothing effects of the Tikhonov
% regularization.
MOBYDICpars.tikhonov = 1e-6;

% Coefficient used too weight the slack variables used to soften the soft
% constraints. A too low value softens the constraints too much, a too high
% value makes the constraints hard.
MOBYDICpars.slack = 1e4;