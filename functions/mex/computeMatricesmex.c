#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <math.h>

/* computeMatrices(P,np,nx,ny,nsamples
 *
 **/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    int i,j,k;
    int idx,idx2;
    int nsamplestot,msize,nv;
    int *rows, *cols, *part_i;
    double *vals, *f, dt, delta;

    /* Total number of samples */
    nsamplestot = pow(nsamples,nx);
    
    /* Number of vertices */
    nv = 1;
    for (i=0;i<nx,i++)
        nv*=(np[i]+1);
        
    /* Maximum number of non-zero elements in the matrix */
    msize = nv*pow(3,nx);
    
    /* Initialize arrays for the definition of the sparse matrix H*/
    rows = mxMalloc(msize*sizeof(int));
    cols = mxMalloc(msize*sizeof(int));
    vals = mxMalloc(msize*sizeof(double));
    
    /* Initialize matrix f */
    f = mxMalloc(nv*ny*sizeof(double));
    
    /* Initialize matrices */
    xmin_i = mxMalloc(nx*sizeof(double));
    xmax_i = mxMalloc(nx*sizeof(double));
    samples = mxMalloc(nx*nsamples*sizeof(double));
    part_i = mxMalloc(nx*sizeof(int));
    
    kk = 1;
    
    /* Loop on all alpha bases */
    for (i=0;i<nv;i++) {
        /* Indices of the subdivisions corresponding to the vertex of the
           i-th alpha basis */
        part_i = linear2matrix(i,np+1);
        
        /* Domain of the i-th alpha basis and samples on this domain */
        dt = 1;
        for (k=0;k<nx;k++) {
            
            if (part_i[k] == 1) {
                idx = getIndex(k,1,nx);
                xmin_i[k] = P[idx];
            }
            else {
                idx = getIndex(k,part_i[k]-1,nx);
                xmin_i[k] = P[idx];
            }
            
            if (part_i[k] == np[k]+1 {
                idx = getIndex(k,np[k]+1,nx);
                xmax_i[k] = P[idx];
            }
            else {
                idx = getIndex(k,part_i[k]+1,nx);
                xmax_i[k] = P[idx];
            }
            
            /* Create samples */
            delta = (xmax_i[k]-xmin_i[k])/(nsamples-1);
            idx = getIndex(k,0,nx);
            samples[idx] = xmin_i[k]; 
            for (j=1;j<nsamples;j++) {
                idx = getIndex(k,j,nx);
                idx2 = getIndex(k,j-1,nx);
                samples[idx] = samples[idx2]+delta;
            }
            idx = getIndex(k,1,nx);
            idx2 = getIndex(k,0,nx);
            dt *= (samples[idx]-samples[idx2]);
            
        }
        
        
    
        
    }
    
}


% Loop on all alpha bases
for i = 1:nv
    
    
   
   
    
    % Create samples in domain
    [x{1:nx}] = ndgrid(samples{:});
    X = zeros(nsamplestot,nx);
    for k = 1:nx
        X(:,k) = x{k}(:);
    end
    
    % Evaluate function to approximate
    if isa(fun,'MOBYFunction')
        fval = fun.eval(X);
    elseif isa(fun,'controller')
        NX = fun.getNumberOfStates();
        NP = fun.getNumberOfParameters();
        ND = fun.getNumberOfUnmeasurableInputs();
        
        XX = X(:,1:NX);
        PP = X(:,NX+1:NX+NP);
        DD = X(:,NX+NP+1:NX+NP+ND);
        XREF = X(:,NX+NP+ND+1:end);
        
        fval = fun.eval(XX,PP,DD,XREF);
    end
    
    % Index of not NaN values (i.e., points outside the domain)
    ninan = all(~isnan(fval),2);
    
    X = X(ninan,:);
    
    % Sample the i-th alpha bases in points X
    falpha_i = fpwas.alphaBasis(X,i);
    
    % Compute integral
    f(i,:) = dt*falpha_i'*fval(ninan,:);
    
    if any(ninan)
        
        % Loop on all alpha bases
        for j = 1:nv
            
            % The matrix is symmetric, so compute only the upper half
            if j >= i
                
                % Indices of the subdivisions corresponding to the vertex of the
                % j-th alpha basis
                if usemex
                    part_j = linear2matrixmex(j,np+1);
                else
                    part_j = linear2matrix(j,np+1);
                end
                % Distance between the subdivisions
                partdiff = part_i-part_j;
                
                % If the distance is greater than one, the domains of the i-th
                % and j-th alpha bases do not intersect, so H(i,j) = 0
                if any(abs(partdiff) > 1)
                    
                    % Do nothing...
                    
                else
                    
                    % Sample the j-th alpha basis in points X
                    falpha_j = fpwas.alphaBasis(X,j);
                    
                    % Compute the integral
                    int = dt*falpha_i'*falpha_j;
                    
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




