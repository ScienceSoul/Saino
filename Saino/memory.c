/*
 *  memory.c
 *  FlowErosModel alpha 3.
 *
 *  Created by Hakime SEDDIK on Thu Mar 11 2004.
 *  Copyright (c) 2004 Institute of Low Temperature Science. All rights reserved.
 *
 */

#include "memory.h"

#define FI_END 1
#define FREE_ARG char*

int *intvec(long nl, long nh)
{
	int *v;
	
	v = (int *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(int)));
	if (!v) errorfunct("intvec", "Allocation failure for the integer vector");
	return v-nl+FI_END;
	
}

unsigned long *ulongvec(long nl, long nh)
{
    unsigned long *v;
	
	v = (unsigned long *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(unsigned long)));
	if (!v) errorfunct("intvec", "Allocation failure for the long vector");
	return v-nl+FI_END;
    
}

double *doublevec(long nl, long nh)
{
	
	double *v;
	
	v = (double *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(double)));
	if (!v) errorfunct("doublevec", "Allocation failure for the double vector");
	
	return v-nl+FI_END;
	
}

float *floatvec(long nl, long nh)
{
	float *v;
	
	v = (float *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(float)));
	if (!v) errorfunct("floatvec", "Allocation failure for the float vector");
	
	return v-nl+FI_END;
	
}

double complex *cdoublevec(long nl, long nh)
{
	
	double complex *v;
	
	v = (double complex *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(double complex)));
	if (!v) errorfunct("cdoublevec", "Allocation failure for the complex double vector");
	
	return v-nl+FI_END;
	
}

double **doublematrix(long nrl, long nrh, long ncl, long nch)
{
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	double **m;
	
	m = (double **) malloc((size_t) ((nrow+FI_END)*sizeof(double*)));
	if (!m) errorfunct("doublematrix", "Allocation failure for the double matrix 1");
	m += FI_END;
	m -= nrl;
	
	m[nrl]=(double *) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(double)));
	if (!m[nrl]) errorfunct("doublematrix", "Allocation failure for the double matrix 2");
	m[nrl] += FI_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	return m;
	
}

float **floatmatrix(long nrl, long nrh, long ncl, long nch)
{
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	float **m;
	
	m = (float **) malloc((size_t) ((nrow+FI_END)*sizeof(float*)));
	if (!m) errorfunct("floatmatrix", "Allocation failure for the float matrix 1");
	m += FI_END;
	m -= nrl;
	
	m[nrl]=(float *) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(float)));
	if (!m[nrl]) errorfunct("floatmatrix", "Allocation failure for the float matrix 2");
	m[nrl] += FI_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	return m;
}

int **intmatrix(long nrl, long nrh, long ncl, long nch)
{
	
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	int **m;
	
	m = (int **) malloc((size_t) ((nrow+FI_END)*sizeof(int*)));
	if (!m) errorfunct("intmatrix", "Allocation failure for the integer matrix 1");
	m += FI_END;
	m -= nrl;
	
	m[nrl] = (int *) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(int)));
	if (!m[nrl]) errorfunct("intmatrix", "Allocation failure for the integer matrix 2");
	m[nrl] += FI_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	return m;
	
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh)
{
	long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ngh-ndl+1;
	double ***t;
	
	/*Allocate pointers to pointers to rows*/
	t = (double ***) malloc((size_t) ((nrow+FI_END)*sizeof(double**)));
	if (!t) errorfunct("d3tensor", "Allocation failure 1 in d3tensor()");
	t += FI_END;
	t -= nrl;
	
	/*Alocate pointers to rows and set pointers to them*/
	t[nrl] = (double **) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(double*)));
	if (!t[nrl]) errorfunct("d3tensor", "Allocation failure 2 in d3tensor()");
	t[nrl] += FI_END;
	t[nrl] -= ncl;
	
	/*Allocate rows and set pointers to them*/
	t[nrl][ncl] = (double *) malloc((size_t) ((nrow*ncol*ndep+FI_END)*sizeof(double)));
	if (!t[nrl][ncl]) errorfunct("d3tensor", "Allocation failure 3 in d3tensor()");
	t[nrl][ncl] += FI_END;
	t[nrl][ncl] -= ndl;
	
	for (j=ncl+1; j<=nch; j++) t[nrl][j] = t[nrl][j-1]+ndep;
	for (i=nrl+1; i<=nrh; i++) {
		t[i] = t[i-1]+ncol;
		t[i][ncl] = t[i-1][ncl]+ncol*ndep;
		for (j=ncl+1; j<=nch; j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*Return pointer to array of pointers to rows*/
	return t;	    
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh)
{
	long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ngh-ndl+1;
	float ***t;
	
	/*Allocate pointers to pointers to rows*/
	t = (float ***) malloc((size_t) ((nrow+FI_END)*sizeof(float**)));
	if (!t) errorfunct("d3tensor", "Allocation failure 1 in f3tensor()");
	t += FI_END;
	t -= nrl;
	
	/*Alocate pointers to rows and set pointers to them*/
	t[nrl] = (float **) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(float*)));
	if (!t[nrl]) errorfunct("d3tensor", "Allocation failure 2 in f3tensor()");
	t[nrl] += FI_END;
	t[nrl] -= ncl;
	
	/*Allocate rows and set pointers to them*/
	t[nrl][ncl] = (float *) malloc((size_t) ((nrow*ncol*ndep+FI_END)*sizeof(float)));
	if (!t[nrl][ncl]) errorfunct("d3tensor", "Allocation failure 3 in f3tensor()");
	t[nrl][ncl] += FI_END;
	t[nrl][ncl] -= ndl;
	
	for (j=ncl+1; j<=nch; j++) t[nrl][j] = t[nrl][j-1]+ndep;
	for (i=nrl+1; i<=nrh; i++) {
		t[i] = t[i-1]+ncol;
		t[i][ncl] = t[i-1][ncl]+ncol*ndep;
		for (j=ncl+1; j<=nch; j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*Return pointer to array of pointers to rows*/
	return t;	    
	
}

int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh)
{
	long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ngh-ndl+1;
	int ***t;
	
	/*Allocate pointers to pointers to rows*/
	t = (int ***) malloc((size_t) ((nrow+FI_END)*sizeof(int**)));
	if (!t) errorfunct("i3tensor", "Allocation failure 1 in i3tensor()");
	t += FI_END;
	t -= nrl;
	
	/*Alocate pointers to rows and set pointers to them*/
	t[nrl] = (int **) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(int*)));
	if (!t[nrl]) errorfunct("i3tensor", "Allocation failure 2 in i3tensor()");
	t[nrl] += FI_END;
	t[nrl] -= ncl;
	
	/*Allocate rows and set pointers to them*/
	t[nrl][ncl] = (int *) malloc((size_t) ((nrow*ncol*ndep+FI_END)*sizeof(int)));
	if (!t[nrl][ncl]) errorfunct("i3tensor", "Allocation failure 3 in i3tensor()");
	t[nrl][ncl] += FI_END;
	t[nrl][ncl] -= ndl;
	
	for (j=ncl+1; j<=nch; j++) t[nrl][j] = t[nrl][j-1]+ndep;
	for (i=nrl+1; i<=nrh; i++) {
		t[i] = t[i-1]+ncol;
		t[i][ncl] = t[i-1][ncl]+ncol*ndep;
		for (j=ncl+1; j<=nch; j++) t[i][j]=t[i][j-1]+ndep;
	}
	
	/*Return pointer to array of pointers to rows*/
	return t;	    
}

void free_dvector(double *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-FI_END));
}

void free_fvector(float *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-FI_END));
}

void free_ivector(int *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-FI_END));
}

void free_cdvector(double complex *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-FI_END));
}

void free_ulvector(unsigned long *v, long nl, long nh) 
{    
    free((FREE_ARG) (v+nl-FI_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{ 
	free((FREE_ARG) (m[nrl]+ncl-FI_END));
	free((FREE_ARG) (m+nrl-FI_END));
}

void free_fmatrix(float **m, long nrl, long nrh, long ncl, long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-FI_END));
	free((FREE_ARG) (m+nrl-FI_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-FI_END));
	free((FREE_ARG) (m+nrl-FI_END));
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
	/*Free a double d3tensor allocated by d3tensor*/
	free((FREE_ARG) (t[nrl][ncl]+ndl-FI_END));
	free((FREE_ARG) (t[nrl]+ncl-FI_END));
	free((FREE_ARG) (t+nrl-FI_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
	/*Free a double d3tensor allocated by d3tensor*/
	free((FREE_ARG) (t[nrl][ncl]+ndl-FI_END));
	free((FREE_ARG) (t[nrl]+ncl-FI_END));
	free((FREE_ARG) (t+nrl-FI_END));
}

void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
	/*Free a double d3tensor allocated by d3tensor*/
	free((FREE_ARG) (t[nrl][ncl]+ndl-FI_END));
	free((FREE_ARG) (t[nrl]+ncl-FI_END));
	free((FREE_ARG) (t+nrl-FI_END));
}

Nodes_t *NodeVec(long nl, long nh)
{
    Nodes_t *v;
    
    v = (Nodes_t *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(Nodes_t)));
    if (!v) errorfunct("NodeVec", "Allocation failure for the nodes vector");
    
    return v-nl+FI_END;
    
}

Element_t *ElementVec(long nl, long nh)
{
    Element_t *v;
    
    v = (Element_t *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(Element_t)));
    if (!v) errorfunct("ElementVec", "Allocation failure for the elements vector");
    
    return v-nl+FI_END;
}

BDElement_t *BDElementVec(long nl, long nh)
{
    BDElement_t *v;
    
    v = (BDElement_t *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(BDElement_t)));
    if (!v) errorfunct("BDElementVec", "Allocation failure for the boundary elements vector");
    
    return v-nl+FI_END;
}

void errorfunct(char head[], char message[], ...)
{   
    int np;
    va_list ap;
    
    np = 0;
    
    va_start(ap,message);
    np = va_arg(ap, int);
    va_end(ap);
    
    if (np > 0) {
        printf("%s: %s %d\n", head, message, np);
        exit(-1);
    } else {
        printf("%s: %s\n", head, message);
        exit(-1); 
    }

}

void warnfunct(char head[], char message[], ...) {
    
    int np;
    va_list ap;
    
    np = 0;
    
    va_start(ap,message);
    np = va_arg(ap, int);
    va_end(ap);

    if (np > 0) {
        printf("%s: %s %d\n", head, message, np);
    } else {
        printf("%s: %s\n", head, message);
    }
    
}



