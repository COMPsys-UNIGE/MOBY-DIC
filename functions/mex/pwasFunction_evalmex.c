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

void bubblesort(double *x,double *xs,int n) {
    
    bool sorted, inside;
    int i;
    double tmp;
    
    /* Copy x into xs */
    for(i=0;i<n;i++)
        xs[i] = x[i];
    
    sorted = false;
    while(!sorted)
    {
        inside = false;
        for(i=0;i<n-1;i++)
        {
            if(xs[i+1] > xs[i])
            {
                tmp = xs[i];
                xs[i] = xs[i+1];
                xs[i+1] = tmp;
                inside = true;
            }
        }
        if(!inside)
            sorted = true;
    }
}

void findIndex(int *x, double *np, int ndim, int *addr) {
    
    int i,k,idx1;
    int coeff;
    
    for(i=0;i<ndim+1;i++) {
        coeff = 1;
        addr[i] = 0;
        for(k=0;k<ndim;k++) {
            idx1 = getIndex(k,i,ndim);
            addr[i] += coeff*x[idx1];
            coeff *= np[k]+1;
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    mwSize ndim,nf,nw,npts,npmax;
    double *P, *dP, *w, *np, *x, *xmin, *xmax, tol;
    double *xint, *xdec, *xdecs, *mu, *mu0, *u,sum;
    int *idx,*addr,*xs;
    bool outside;
    int idx1,idx2,idx3,ii,j,k;
    
    /* Tolerance of 1e-5 */
    tol = 0.00001;
    
    /* Read inputs from mex interface */
    
    /* Subdivisions per dimension */
    P = mxGetPr(prhs[0]);
    /* Maximum number of subdivisions per dimension */
    npmax = mxGetN(prhs[0]);
    /* Number of subdivisions per dimension */
    np = mxGetPr(prhs[1]);
    /* Number of domain dimensions */
    ndim = mxGetM(prhs[1]);
    /* Weights */
    w = mxGetPr(prhs[2]);
    /* Number of codomain dimensions */
    nf = mxGetN(prhs[2]);
    /* Number of weights */
    nw = mxGetM(prhs[2]);
    /* Minimum values of the function domain */
    xmin = mxGetPr(prhs[3]);
    /* Maximum values of the function domain */
    xmax = mxGetPr(prhs[4]);
    /* Input points */
    x = mxGetPr(prhs[5]);
    /* Number of points */
    npts = mxGetN(prhs[5]);
        
    /* Prepare outputs */
    plhs[0] = mxCreateDoubleMatrix(nf,npts,mxREAL);
    
    u = mxGetPr(plhs[0]);
    
    xint = mxMalloc(ndim*sizeof(double));
    xdec = mxMalloc(ndim*sizeof(double));
    xdecs = mxMalloc(ndim*sizeof(double));
    idx = mxMalloc(ndim*sizeof(int));
    addr = mxMalloc((ndim+1)*sizeof(int));
    xs = mxMalloc(ndim*(ndim+1)*sizeof(int));
    mu = mxMalloc((ndim+1)*sizeof(double));
    mu0 = mxMalloc((ndim+1)*sizeof(double));
    dP = mxMalloc(ndim*npmax*sizeof(double));
       
    /* Retrieve length of subdivisions */
    for(j=0; j<ndim; j++) {
        for(k=0; k<np[j]; k++) {
            idx1 = getIndex(j,k+1,ndim);
            idx2 = getIndex(j,k,ndim);
            dP[idx2] = P[idx1]-P[idx2];
        }
        idx2 = getIndex(j,(int)np[j],ndim);
        dP[idx2] = 1;
    }
    
    /* Loop on number of points */
    for(ii = 0; ii < npts; ii++) {
               
        /* Check if point is inside the domain */
        outside = 0;
        /* Loop on number of dimensions */
        for (j = 0; j < ndim; j++) {
            idx1 = getIndex(j,ii,ndim);
            /* If point iis outside bounds, set outside = 1 */
            if(x[idx1] < xmin[j] || x[idx1] > xmax[j])
            {
                outside = 1;
            }
        }
        
        
        /* If the point is outside, set the value of the controller at 1e99 */
        if (outside) {
            for (k = 0; k < nf; k++) {
                idx1 = getIndex(k,ii,(int)nf);
                u[idx1] = 1.0e99;   /* --> u(k,ii) */
            }
        }
        
        /* Otherwise, look for the region */
        else {
            
            /*Find bottom-left vertex */
            for (j = 0; j < ndim; j++) {
                idx1 = getIndex(j,ii,ndim);
                for (k = 0; k < np[j]; k++) {
                    idx2 = getIndex(j,k,ndim);
                    idx3 = getIndex(j,k+1,ndim);
                    /* x[j][ii] >= P[j][k] && x[j][ii] >= P[j][k+1] */
                    if(x[idx1] >= P[idx2] && x[idx1] < P[idx3]) {
                        idx[j] = k;
                        xint[j] = P[idx2];
                        xdec[j] = (x[idx1]-xint[j])/dP[idx2];
                    }
                }
                idx2 = getIndex(j,(int)np[j],ndim);
                
                if(x[idx1] == P[idx2]) {
                    idx[j] = np[j];
                    xint[j] = P[idx2];
                    xdec[j] = (x[idx1]-xint[j])/dP[idx2];
                }
            }
            
            /* Sort decimal parts */
            bubblesort(xdec,xdecs,ndim);
                        
            for(j=0;j<ndim;j++) {
                idx1 = getIndex(j,0,ndim);
                xs[idx1] = idx[j];
            }
            for(j=0;j<ndim;j++) {
                for(k=0;k<ndim;k++) {
                    idx1 = getIndex(k,j+1,ndim);
                    if(xdec[k] >= xdecs[j])
                        xs[idx1] = idx[k]+1;
                    else
                        xs[idx1] = idx[k];
                }
            }
            
            /* Find index of simplex vertices */
            findIndex(xs,np,ndim,addr);
            
            /* Compute coefficients mu */
            mu[0] = 1-xdecs[0];
            
            /* mu0[j] is true if mu[j] = 0 */
            mu0[0] = !mu[0];
            for(j=1;j<ndim;j++)
            {
                mu[j] = xdecs[j-1]-xdecs[j];
                mu0[j] = !mu[j];
            }
            mu[ndim] = xdecs[ndim-1];
            mu0[ndim] = !mu[ndim];
            
            /* Compute weighted sum */
            for(k=0;k<nf;k++) {
                sum = 0;
                for(j=0;j<ndim+1;j++) {
                    if(!mu0[j]) {
                        idx1 = getIndex(addr[j],k,nw);
                        sum += mu[j]*w[idx1];
                    }
                }
                idx2 = getIndex(k,ii,(int)nf);
                u[idx2] = sum;
                
            }
        }
        
        
    }
    
    mxFree(xint);
    mxFree(xdec);
    mxFree(xdecs);
    mxFree(idx);
    mxFree(addr);
    mxFree(xs);
    mxFree(mu);
    mxFree(mu0);
    mxFree(dP);
    
}

