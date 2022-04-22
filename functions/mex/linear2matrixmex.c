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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    mwSize n;
    double *siz,num,den,*matrixIdx,*linIdx;
    int i;
    
    /* Read inputs from mex interface */
    
    /* Linear index */
    linIdx = mxGetPr(prhs[0]);
    /* Matrix dimensions */
    siz = mxGetPr(prhs[1]);
    /* Number of dimensions*/
    n = mxGetM(prhs[1]);
    
    /* Prepare outputs */
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    
    matrixIdx = mxGetPr(plhs[0]);
    
    
    num = linIdx[0]-1;
    den = 1;
    for(i=0;i<n-1;i++)
        den*=siz[i];
    for(i=n-1;i>0;i--){
        matrixIdx[i] = (int)(num/den);
        num = num-matrixIdx[i]*den;
        den = den/siz[i-1];
    }
    matrixIdx[0] = (int)(num/den);
    for(i=0;i<n;i++)
        matrixIdx[i] = matrixIdx[i]+1;   
}

