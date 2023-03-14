#ifndef _PWAS_FUN_H_
#define _PWAS_FUN_H_

#define nDim 4
#define nY 1
#define nWeight 2000

void getMu(float *z_sorted, float *mu);

void Sorting(float *number, int n);

void Transform(float *x, float *z);

void divideZ(float *Z, int *intZ, float *decZ);

void getSimplexVertices(int *intZ, float *decZ, float *sortedDecZ, int *addr);

void calculatePWAS(float *x, float *u);

float uFunction(int *vertexAddr, float *mu, int uIndex, float *x); 

#endif
