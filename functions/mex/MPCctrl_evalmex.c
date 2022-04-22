/* 
 Acknowledgements

 Contributors:
 
 * Alberto Oliveri (alberto.oliveri@unige.it)
 
 Copyright is with:
 
 * Copyright (C) 2011 University of Genoa, Italy.

-------------------------------------------------------------------------
 Legal note:
        This program is free software; you can redistribute it and/or
         modify it under the terms of the GNU General Public
         License as published by the Free Software Foundation; either
         version 2.1 of the License, or (at your option) any later version.

         This program is distributed in the hope that it will be useful,
         but WITHOUT ANY WARRANTY; without even the implied warranty of
         MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
         General Public License for more details.

         You should have received a copy of the GNU General Public
         License along with this library; if not, write to the 
         Free Software Foundation, Inc., 
         59 Temple Place, Suite 330, 
         Boston, MA  02111-1307  USA
-------------------------------------------------------------------------

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
    
    mwSize nedges, npts, ned, nif, ndim;
    double *H, *K, *F, *G, *x, *u, *reg, *Hd, *Kd;
    double *indexes, *nindexes, *nfunctions, *signs, *ifunctions;
    double  nobj, nr, nu;
    int ii, i, j, k, t, region, idx1, idx2;
    bool found, inside, outside;
    int ne, offset, nf, nnu, nfuntot;
    double sum;
    double tol;
    
    
    /* Tolerance of 1e-5 */
    tol = 0.00001;

    /* Read inputs from mex interface */
    nobj = *mxGetPr(prhs[0]);
    nu = *mxGetPr(prhs[1]);
    nr = *mxGetPr(prhs[2]);
    H = mxGetPr(prhs[3]);
    nedges = mxGetM(prhs[3]);
    ndim = mxGetN(prhs[3]);
    K = mxGetPr(prhs[4]);
    F = mxGetPr(prhs[5]);
    nfuntot = mxGetM(prhs[5]);
    G = mxGetPr(prhs[6]);
    nindexes = mxGetPr(prhs[7]);
    indexes = mxGetPr(prhs[8]);
    nfunctions = mxGetPr(prhs[9]);
    ifunctions = mxGetPr(prhs[10]);
    nif = mxGetM(prhs[10]);
    signs = mxGetPr(prhs[11]);
    Hd = mxGetPr(prhs[12]);
    ned = mxGetM(prhs[12]);   /* Number of edges of the domain polytope */
    Kd = mxGetPr(prhs[13]);
    x = mxGetPr(prhs[14]);
    npts = mxGetN(prhs[14]);   /* Number of points */
    
    nnu = nobj*nu;
    
    /* Prepare outputs */
    plhs[0] = mxCreateDoubleMatrix(nu*nobj,npts,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,npts,mxREAL);
    
    u = mxGetPr(plhs[0]);
    reg = mxGetPr(plhs[1]);
    
    /*printf("%f %f %f %f %f %f %f %f %f",indexes[0],indexes[1],indexes[2],indexes[3],indexes[4],indexes[5],indexes[6],indexes[7],indexes[8]);*/
    
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
                idx2 = getIndex(j,ii,ndim);    
                sum += Hd[idx1]*x[idx2];    /* --> Hd(i,j)*x(j,ii) */
            }
            /* Evaluate current edge */
            outside = sum > Kd[i];
            if (outside)
                break;
        }
        
        /* If the point is outside, set the value of the controller at 1e99,
         and the region at 0*/
        if (outside) {
            k = 0;
            for (i = 0; i < nobj; i++) {
                for (j = 0; j < nu; j++) {
                    idx1 = getIndex(k,ii,nnu);    
                    u[idx1] = 1.0e99;   /* --> u(k,ii) */
                    k++;
                }
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
                /*rintf("Region %d \n",i);*/
                ne = nindexes[i];
                inside = 1;
                j = 0;
                /* Loop on number of edges of current region */
                while (inside && j<ne) {
                    sum = 0;
                    for (k = 0; k < ndim; k++) {
                        idx1 = getIndex(k,ii,ndim);
                        idx2 = getIndex((int) indexes[offset+j]-1,k,nedges);
                        sum += signs[offset+j]*x[idx1]*H[idx2]; /* --> x(k,ii)*H(indexes(offset+j)-1,k) */
                        /*printf("%f x %f x %f +",signs[offset+j],x[idx1],H[idx2]);*/
                    }
                    /* Evaluate inequality */
                    inside = sum <= signs[offset+j]*K[(int) indexes[offset+j]-1]+tol;
                    /*printf("<= %f \n",signs[offset+j]*K[(int) indexes[offset+j]-1]+tol);*/
                    j++;
                }
                offset += ne;
                if (j == ne && inside)
                    found = 1;
                i++;
            }
            
            /* If the region has not been found, set u = 1e99 end reg = 0 */
            if (!found) {
                k = 0;
                for (i = 0; i < nobj; i++) {
                    for (j = 0; j < nu; j++) {
                        idx1 = getIndex(k,ii,nnu);
                        u[idx1] = 1.0e99;   /* --> u(k,ii) */
                        k++;
                    }
                }
                reg[ii] = 0;
            }
            /* If the region has been found, compute controller */
            else {
                
                /* Index of region containing point x */
                region = i-1;
                
                offset = 0;
                t = 0;
                /* Loop on number of controllers */
                for (i = 0; i < nobj; i++) {
                    /* Number of functions in current controller */
                    nf = nfunctions[i];
                    for (j = 0; j < nu; j++) {
                        sum = 0;
                        for (k = 0; k < ndim; k++) {
                            idx1 = getIndex(offset+(int) ifunctions[(int)(nu*nobj)*region+t]-1,k,nfuntot);
                            idx2 = getIndex(k,ii,ndim);
                            sum += F[idx1]*x[idx2]; /* --> F(offset+ifunctions[((nu*nobj)*region+t)]-1,k)*x(k,ii) */
                            /*printf("%f x %f +",F[idx1],x[idx2]);*/
                        }
                        sum += G[offset + (int) ifunctions[(int)(nu*nobj)*region+t]-1];
                        /*printf("%f \n ",G[offset + (int) ifunctions[(int) ((nu*nobj)*region+t)]-1]);*/
                        idx1 = getIndex(t,ii,nnu);
                        u[idx1] = sum;  /* --> u(t,ii) */
                        /*printf("%d -> %f \n ",idx1,u[idx1]);*/
                        t++;
                    }
                    offset += nf;
                }
                reg[ii] = region+1;
            }
        }
    }
}
