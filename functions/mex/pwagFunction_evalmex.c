/*
 * Acknowledgements
 *
 * Contributors:
 *
 * Alberto Oliveri (alberto.oliveri@unige.it)
 *
 * Copyright is with:
 *
 * Copyright (C) 2011 University of Genoa, Italy.
 *
 * -------------------------------------------------------------------------
 * Legal note:
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330,
 * Boston, MA  02111-1307  USA
 * -------------------------------------------------------------------------
 *
 */

#include "mex.h"
#include "matrix.h"
#include <stdio.h>

/* Given the indices (i,j) of the matrix row and column, and the number of
 * rows M of the matrix, this functions gives back the single index
 * idx = M*j+i */
int getIndex(i,j,M)
{
    int idx;
    idx = M*j+i;
    return idx;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    mwSize nedges, npts, ned, ndim;
    double *H, *K, *F, *G, *x, *u, *reg, *Hd, *Kd;
    double *indexes, *nindexes, *signs, *ifunctions;
    double  nr, nf;
    int ii, i, j, k, t, region, idx1, idx2;
    bool found, inside, outside;
    int ne, offset, nfuntot;
    double sum;
    double tol;
    
    
    /* Tolerance of 1e-5 */
    tol = 0.00001;
    
    /* Read inputs from mex interface */
    
    /* Number of functions */
    nf = *mxGetPr(prhs[0]);
    /* Number of regions */
    nr = *mxGetPr(prhs[1]);
    /* Matrix of coefficients H */
    H = mxGetPr(prhs[2]);
    /* Total number of edges */
    nedges = mxGetM(prhs[2]);
    /* Number of domain dimensions */
    ndim = mxGetN(prhs[2]);
    /* Number of coefficients K */
    K = mxGetPr(prhs[3]);
    /* Matrix of coefficients F */
    F = mxGetPr(prhs[4]);
    /* Total number of affine functions */
    nfuntot = mxGetM(prhs[4]);
    /* Matrix of coefficients G */
    G = mxGetPr(prhs[5]);
    /* Number of edges for each region */
    nindexes = mxGetPr(prhs[6]);
    /* Indices of the edges of each region */
    indexes = mxGetPr(prhs[7]);
    /* Indices of the affine functions for each region */
    ifunctions = mxGetPr(prhs[8]);
    /* Sign of the edges for each region */
    signs = mxGetPr(prhs[9]);
    /* H coefficients of the function domain */
    Hd = mxGetPr(prhs[10]);
    /* Number of edges of the domain polytope */
    ned = mxGetM(prhs[10]);   
    /* K coefficients of the function domain */
    Kd = mxGetPr(prhs[11]);
    /* Input points */
    x = mxGetPr(prhs[12]);
    /* Number of points */
    npts = mxGetM(prhs[12]);    
    
    /* Prepare outputs */
    plhs[0] = mxCreateDoubleMatrix(nf,npts,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,npts,mxREAL);
    
    u = mxGetPr(plhs[0]);
    reg = mxGetPr(plhs[1]);
       
    /* Loop on number of points */
    for(ii = 0; ii < npts; ii++) {
               
        /* Check if point is inside the domain */
        outside = 0;
        /* Loop on number of domain edges */
        for (i = 0; i < ned; i++) {
            sum = 0;
            /* Loop on number of dimensions */
            for (j = 0; j < ndim; j++) {
                idx1 = getIndex(i,j,ned);
                idx2 = getIndex(ii,j,npts);
                sum += Hd[idx1]*x[idx2];    /* --> Hd(i,j)*x(ii,j) */
            }
            /* Evaluate current edge */
            outside = sum > Kd[i]+tol;
            if (outside)
                break;
        }
        
        /* If the point is outside, set the value of the controller at 1e99,
         * and the region at 0*/
        if (outside) {
            for (k = 0; k < nf; k++) {
                idx1 = getIndex(k,ii,(int)nf);
                u[idx1] = 1.0e99;   /* --> u(k,ii) */
            }
            reg[ii] = 0;
        }
        
        /* Otherwise, look for the region */
        else {
            found = 0;
            i = 0;
            offset = 0;
            /* Loop on number of regions */
            while (!found && i<nr) {
                ne = nindexes[i];
                inside = 1;
                j = 0;
                /* Loop on number of edges of current region */
                while (inside && j<ne) {
                    sum = 0;
                    for (k = 0; k < ndim; k++) {
                        idx1 = getIndex(ii,k,npts);
                        idx2 = getIndex((int) indexes[offset+j]-1,k,nedges);
                        sum += signs[offset+j]*x[idx1]*H[idx2]; /* --> x(ii,k)*H(indexes(offset+j)-1,k) */
                    }
                    /* Evaluate inequality */
                    inside = sum <= signs[offset+j]*K[(int) indexes[offset+j]-1]+tol;
                    j++;
                }
                offset += ne;
                if (j == ne && inside)
                    found = 1;
                i++;
               
            }
            
            /* If the region has not been found, set u = 1e99 end reg = 0 */
            if (!found) {
                for (k = 0; k < nf; k++) {
                    idx1 = getIndex(k,ii,(int)nf);
                    u[idx1] = 1.0e99;   /* --> u(k,ii) */
                }
                reg[ii] = 0;
            }
            /* If the region has been found, compute controller */
            else {
                
                /* Index of region containing point x */
                region = i-1;
                
                offset = 0;
                t = 0;
                
                for (i = 0; i < nf; i++) {
                    sum = 0;
                    for (k = 0; k < ndim; k++) {
                        idx1 = getIndex((int) ifunctions[(int)nf*region+i]-1,k,nfuntot);
                        idx2 = getIndex(ii,k,npts);
                        sum += F[idx1]*x[idx2]; /* --> F(offset+ifunctions[((nu*nobj)*region+t)]-1,k)*x(ii,k) */
                    }
                    sum += G[(int) ifunctions[(int)nf*region+i]-1];
                    idx1 = getIndex(i,ii,(int)nf);
                    u[idx1] = sum;  /* --> u(i,ii) */
                }
                reg[ii] = region+1;
            }
        }
    }
}
