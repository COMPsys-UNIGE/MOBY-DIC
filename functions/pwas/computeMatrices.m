function [H, f, alpha, fval] = computeMatrices(fpwas,fun,options)
% computeMatrices   Computes the matrices defining the quadratic cost function
%
% A PWAS function g(x) is defined as
%
%              _nv_
%              \
%       g(x) = /___  w_j alpha_j(x)     (1)
%              j = -1
%
% In order to approximate a generic non linear function u(x) with a PWAS
% function g(x), it is necessary to solve the following optimization
% problem:
%                /
%           min  | ||u(x)-g(x)||^2 dx        (2)
%           g(x) |
%                / D
%
% By substituting (1) in (2) it is possible to recast (2) as a quadratic
% function of weights vector w, as follows:
%
%           min  w'Hw-2f'w               (3)
%            w
%
% Terms H(ij) are defined as the integral of alpha_i x alpha_j, terms f(i)
% are defined as the integral of alpha_i x u.
%
% [H, F] = computeMatrices(FPWAS,FUN,OPTS)
% FPWAS is a pwasFunction object where the domain and the simplicial
% partition must be defined. FUN is either a MOBYFunction object (defining
% the nonlinear function u(x) to be approximated) or a controller object
% (defining the control function u(x) to be approximated. OPTS is a structure with
% the following fields:
% - nsamples: number of samples per dimension used to numerically compute
%             the integrals. The higher nsamples the higher the
%             approximation accuracy at the cost of a higher computation
%             effort.
%
% H and F are the matrices defining the quadratic cost function.

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

nsamples = options.nsamples;

% Extract information from the pwasFunction object
nv = fpwas.getNumberOfVertices();
P = fpwas.getPartition();
np = fpwas.getNumberOfPartitions();

if isa(fun,'MOBYFunction')
    nx = fun.getDomainDimensions();
    ny = fun.getCodomainDimensions();
elseif isa(fun,'controller')
    nx = fun.getNumberOfDimensions();
    ny = fun.getNumberOfInputs();
end

% Generate samples
samples = cell(nx,1);
steps = cell(nx,1);
nsamplestot = 1;
for i = 1:nx
    samples{i} = zeros(1,np(i)*(nsamples-1));
    steps{i} = zeros(1,np(i)*(nsamples-1));
    offset = 1;
    for j = 1:np(i)
        tmp = linspace(P{i}(j),P{i}(j+1),nsamples);
        tmp(end) = [];
        step = tmp(2)-tmp(1);
        samples{i}(offset:offset+nsamples-2) = tmp;
        steps{i}(offset:offset+nsamples-2) = step;
        offset = offset+nsamples-1;
    end
    nsamplestot = nsamplestot*numel(samples{i});
end

% Create samples in domain
[x{1:nx}] = ndgrid(samples{:});
[ddx{1:nx}] = ndgrid(steps{:});
X = zeros(nsamplestot,nx);
dx = zeros(nsamplestot,nx);
for k = 1:nx
    X(:,k) = x{k}(:);
    dx(:,k) = ddx{k}(:);
end
dx = prod(dx,2);

% Evaluate function to approximate
if isa(fun,'MOBYFunction')
    fval = fun.eval(X');
elseif isa(fun,'controller')
    NX = fun.getNumberOfStates();
    NP = fun.getNumberOfParameters();
    ND = fun.getNumberOfUnmeasurableInputs();
    
    XX = X(:,1:NX);
    PP = X(:,NX+1:NX+NP);
    DD = X(:,NX+NP+1:NX+NP+ND);
    XREF = X(:,NX+NP+ND+1:end);
    
    fval = fun.eval(XX',PP',DD',XREF');
end
fval = fval';

ninan = all(~isnan(fval),2);
fval = fval(ninan,:);
X = X(ninan,:);
dx = dx(ninan);

nzero = nv*(2*nsamples)^nx;
alpha = fpwas.alphaBasis(X,nzero);

allzero = ~any(alpha,1);

% Initialize arrays used to create the sparse matrix
rows = zeros(nv*3^nx,1);
cols = zeros(nv*3^nx,1);
vals = zeros(nv*3^nx,1);
f = zeros(nv,ny);
kk = 1;

% Open waitbar
h = waitbar(0,'Computing matrices...');

% Loop on all alpha bases
for i = 1:nv
    
    alphai = alpha(:,i);
    
    if allzero(i)
        f(i,:) = 0;
    else
        alphaidx = alphai.*dx;
        alphaidx = alphaidx';
        
        f(i,:) = alphaidx*fval;
        
        % The matrix is symmetric, so compute only the upper half
        for j = i:nv
            
            alphaj = alpha(:,j);
            
            if allzero(j)
                
            else
                
                int = alphaidx*alphaj;
                
                % Fill the sparse matrix
                rows(kk) = i;
                cols(kk) = j;
                vals(kk) = int;
                kk = kk+1;
                if i ~= j
                    rows(kk) = j;
                    cols(kk) = i;
                    vals(kk) = int;
                    kk = kk+1;
                end
                
            end
            
        end
    end
    waitbar(i/nv);
end

% Transform matrix f in a unique column
f = f(:);

% Remove zero elements
idx = rows == 0;
rows(idx) = [];
cols(idx) = [];
vals(idx) = [];

% Create sparse matrix H
Hsingle = sparse(rows,cols,vals,nv,nv,numel(vals));

% Replicate the matrix for all the codomain dimensions
H = sparse([]);
for i = 1:ny
    H = blkdiag(H,Hsingle);
end

close(h);




