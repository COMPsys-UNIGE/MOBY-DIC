/* Sort decimal parts */
            bubblesort(xdec,xdecs,ndim);
            
            printf("sorted\n");
                                   
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
            mu0[0] = !mu[0];
            for(j=1;j<ndim;j++)
            {
                mu[j] = xdecs[j-1]-xdecs[j];
                mu0[j] = !mu[j];
            }
            mu[ndim] = xdecs[ndim-1];
            mu0[ndim] = !mu[ndim];
            
            for(k=0;k<nf;k++) {
                sum = 0;
                for(j=0;j<ndim+1;j++) {
                    if(mu0[j]) {
                        idx1 = getIndex(addr[j],k,nw);
                        sum += mu[j]*w[idx1];
                    }
                }
                idx2 = getIndex(k,ii,(int)nf);
                u[idx2] = sum;

            }
        }