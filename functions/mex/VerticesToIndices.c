/*
Acknowledgements
Contributors:

* Alberto Oliveri (alberto.oliveri@unige.it)

Copyright is with:

* Copyright (C) 2012 University of Genoa, Italy.


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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"


/****************************** PROTOTYPES **********************************/
int findIndex(double *v, double *P, int Ndim);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Given the indices (i,j) of the matrix row and column, and the number of
 * rows M of the matrix, this functions gives back the single index
 * idx = M*j+i */
int getIndex(i,j,M)
{
    int idx;
    idx = M*j+i;
    return idx;
}


int findIndex(double *v, double *P, int Ndim)
{
    int n,k, ind=0;
    int prod;
    
    ind=(int)v[0]+1;
    for(n=1;n<Ndim;n++)
    {
        prod=1;
        for(k=0;k<n;k++)
        {
            prod=prod*((int)P[k]+1);
        }
        ind=ind+((int)v[n])*prod;
    }
    
    return ind;
}

void findIndex(int *x, double *np, int ndim, int nV, int *addr) {
    
    int i,k,idx1;
    int coeff;
    
    for(i=0; i<nV; i++) {
        coeff = 1;
        addr[i] = 0;
        for(k = ndim-1; k>=0; k--) {
            idx1 = getIndex(k,i,ndim);
            addr[i] += coeff*x[idx1];
            coeff *= np[k]+1;
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, k;
	int Ndim, Nv, Nb;
	//double *P, *X, *out;
    double *P, *V;
    double *v;
    double *out;

    /* Vertices */
    V = mxGetPr(prhs[0]);
	/* Subdivisions per dimension */
    P = mxGetPr(prhs[1]);
    /* Number of subdivisions per dimension */
    np = mxGetPr(prhs[2]);
    /* Number of vertices */
    Nv = mxGetM(prhs[0]);
    /* Number of domain dimensions */
	Ndim = mxGetN(prhs[0]);
    /* Minimum values of the function domain */
    xmin = mxGetPr(prhs[3]);
    /* Maximum values of the function domain */
    xmax = mxGetPr(prhs[4]);
    
    /* Number of partition vertices */
    Nb=1;
	for (i=0; i<Ndim; i++)
	{
		Nb = Nb*(P[i]+ 1);
	}

    v = mxMalloc(ndim*sizeof(double));

    plhs[0] = mxCreateDoubleMatrix(Nv,1,mxREAL);
    addr = mxGetPr(plhs[0]);

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
            /* If the point lies exactly on the upper edge, then decrease
             * it slightly so that the value of the function can be computed */
            if(x[idx1] == xmax[j])
                    x[idx1] = x[idx1]-tol;
        }
        
        /* If the point is outside, set the index at 0 */
        if (outside) {
            addr[ii] = 0; 
        }
        
        /* Otherwise, look for the region */
        else {
            
            /*Find bottom-left vertex */
            for (j = 0; j < ndim; j++) {
                idx1 = getIndex(j,ii,ndim);
                for (k = 0; k < np[j]; k++) {
                    idx2 = getIndex(j,k,ndim);
                    idx3 = getIndex(j,k+1,ndim);
                    if(x[idx1] >= P[idx2] && x[idx1] < P[idx3]) {
                        idx[j] = k;
                    }
                }
            }
                        
            findIndex(idx,np,ndim,&addr[ii]);        
            
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    for(i=0; i<Nv; i++)
    {
        for(k=0; k<Ndim; k++)
        {
            idx = getIndex(i,k,Nv);
            v[k] = V[idx];  // V[i,k]
        }
        out[i] = (double)findIndex(v,P,Ndim);
        
        out[i] = (out[i]-1)*Nv+i+1;
    }
}

