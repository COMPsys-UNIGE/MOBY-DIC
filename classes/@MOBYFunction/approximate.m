function fappr = approximate(object,opts)

% approximate   approximates the function with a PWAS function
%
% FAPPR = approximate(OBJ,OPTS)
% FAPPR is a pwasFunction object representing a PWA function defined over a
% simplicial partition which approximates the generic function represented
% by MOBYFunction object OBJ. The approximation is performed by minimizing
% the L2 norm of the error between the exact and approximate function.
% OPTS is a structure with the following
% fields:
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
%             Default: nsamples = 5.
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
%      7  ___________
%        |       |   |
%        |       |   |
%        |       |   |
%      5 |_______|___|
%        |       |   |
%      2 |_______|___|
%        1       3   4
%    If a uniform partition is desired, P must be empty (or not set) and
%    OPTS.np must be provided.
% - domain: structure with fields xmin and xmax, defining the 
%           hyper-rectangular domain the approximate function is defined 
%           over. If it is not provided, the hyper-rectangular domain is 
%           automatically computed by taking the bounding box of the exact 
%           function domain. The domain shown in the above figure 
%           corresponds to domain.xmin = [1 2] and domain.xmax = [5 7].

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

% Set default options
opts = PWASset(object,opts);

% Hyper-rectangular domain
xmin = opts.domain.xmin;
xmax = opts.domain.xmax;

% If the domain is not provided...
if isempty(xmin)
    % Retrieve function domain
    Hd = object.domain.Hd;
    Kd = object.domain.Kd;
    
    % Compute bounding box
    P = Polyhedron(Hd,Kd);
    B = P.outerApprox;
    
    % PWAS function domain
    xmin = B.Internal.lb;
    xmax = B.Internal.ub;
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
ny = object.ny;

disp(' ');
disp('Approximating function...');
disp(' ');

% Compute matrices for the L2 approximation
[H, f, alpha, fopt] = computeMatrices(fappr,object,opts);

% Perform least squares optimization
w = H\f;

% Reshape weights to the correct shape
w = reshape(w,nv,ny);

err = alpha*w-fopt;
rmse = sqrt(sum(err.^2)/size(err,1));

% Assign weights to the approximate function
fappr = fappr.setWeights(w);

disp('Done.');
disp(' ');
disp(['Estimated RMSE: ',num2str(rmse)]);
disp(' ');
disp(' ')

% Set names
fappr = fappr.setInputNames(object.xnames);
fappr = fappr.setOutputNames(object.ynames);