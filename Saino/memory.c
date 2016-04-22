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

bool * __nonnull boolvec(long nl, long nh) {
    
    bool *v;
	
	v = (bool *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(bool)));
	if (!v) fatal("intvec", "Allocation failure for the bool vector.");
	return v-nl+FI_END;
}

int * __nonnull intvec(long nl, long nh)
{
	int *v;
	
	v = (int *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(int)));
	if (!v) fatal("intvec", "Allocation failure for the integer vector.");
	return v-nl+FI_END;
}

unsigned long * __nonnull ulongvec(long nl, long nh)
{
    unsigned long *v;
	
	v = (unsigned long *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(unsigned long)));
	if (!v) fatal("intvec", "Allocation failure for the long vector.");
	return v-nl+FI_END;
}

double * __nonnull doublevec(long nl, long nh)
{
	
	double *v;
	
	v = (double *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(double)));
	if (!v) fatal("doublevec", "Allocation failure for the double vector.");
	return v-nl+FI_END;
}

float * __nonnull floatvec(long nl, long nh)
{
	float *v;
	
	v = (float *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(float)));
	if (!v) fatal("floatvec", "Allocation failure for the float vector.");
	return v-nl+FI_END;
}

double complex * __nonnull cdoublevec(long nl, long nh)
{
	
	double complex *v;
	
	v = (double complex *)malloc((size_t) ((nh-nl+1+FI_END)*sizeof(double complex)));
	if (!v) fatal("cdoublevec", "Allocation failure for the complex double vector.");
	return v-nl+FI_END;
}

bool * __nonnull * __nonnull boolmatrix(long nrl, long nrh, long ncl, long nch) {
    
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	bool **m;
	
	m = (bool **) malloc((size_t) ((nrow+FI_END)*sizeof(bool*)));
	if (!m) fatal("intmatrix", "Allocation failure for the integer matrix 1.");
	m += FI_END;
	m -= nrl;
	
	m[nrl] = (bool *) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(bool)));
	if (!m[nrl]) fatal("intmatrix", "Allocation failure for the integer matrix 2.");
	m[nrl] += FI_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	return m;
}

int * __nonnull * __nonnull intmatrix(long nrl, long nrh, long ncl, long nch)
{
	
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	int **m;
	
	m = (int **) malloc((size_t) ((nrow+FI_END)*sizeof(int*)));
	if (!m) fatal("intmatrix", "Allocation failure for the integer matrix 1.");
	m += FI_END;
	m -= nrl;
	
	m[nrl] = (int *) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(int)));
	if (!m[nrl]) fatal("intmatrix", "Allocation failure for the integer matrix 2.");
	m[nrl] += FI_END;
	m[nrl] -= ncl;
	
	for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	return m;
}

double * __nonnull * __nonnull doublematrix(long nrl, long nrh, long ncl, long nch)
{
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	double **m;
	
	m = (double **) malloc((size_t) ((nrow+FI_END)*sizeof(double*)));
	if (!m) fatal("doublematrix", "Allocation failure for the double matrix 1.");
	m += FI_END;
	m -= nrl;
	
	m[nrl]=(double *) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(double)));
	if (!m[nrl]) fatal("doublematrix", "Allocation failure for the double matrix 2.");
	m[nrl] += FI_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	return m;
}

float * __nonnull * __nonnull floatmatrix(long nrl, long nrh, long ncl, long nch)
{
	long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	float **m;
	
	m = (float **) malloc((size_t) ((nrow+FI_END)*sizeof(float*)));
	if (!m) fatal("floatmatrix", "Allocation failure for the float matrix 1.");
	m += FI_END;
	m -= nrl;
	
	m[nrl]=(float *) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(float)));
	if (!m[nrl]) fatal("floatmatrix", "Allocation failure for the float matrix 2.");
	m[nrl] += FI_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	return m;
}

double complex * __nonnull * __nonnull cdoublematrix(long nrl, long nrh, long ncl, long nch) {
    
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	double complex **m;
	
	m = (double complex**) malloc((size_t) ((nrow+FI_END)*sizeof(double complex*)));
	if (!m) fatal("cdoublematrix", "Allocation failure for the double matrix 1.");
	m += FI_END;
	m -= nrl;
	
	m[nrl]=(double complex*) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(double complex)));
	if (!m[nrl]) fatal("cdoublematrix", "Allocation failure for the double matrix 2.");
	m[nrl] += FI_END;
	m[nrl] -= ncl;
	
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
	return m;
}

double * __nonnull * __nonnull * __nonnull d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh)
{
	long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ngh-ndl+1;
	double ***t;
	
	/*Allocate pointers to pointers to rows*/
	t = (double ***) malloc((size_t) ((nrow+FI_END)*sizeof(double**)));
	if (!t) fatal("d3tensor", "Allocation failure 1 in d3tensor().");
	t += FI_END;
	t -= nrl;
	
	/*Alocate pointers to rows and set pointers to them*/
	t[nrl] = (double **) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(double*)));
	if (!t[nrl]) fatal("d3tensor", "Allocation failure 2 in d3tensor().");
	t[nrl] += FI_END;
	t[nrl] -= ncl;
	
	/*Allocate rows and set pointers to them*/
	t[nrl][ncl] = (double *) malloc((size_t) ((nrow*ncol*ndep+FI_END)*sizeof(double)));
	if (!t[nrl][ncl]) fatal("d3tensor", "Allocation failure 3 in d3tensor().");
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

float * __nonnull * __nonnull * __nonnull f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh)
{
	long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ngh-ndl+1;
	float ***t;
	
	/*Allocate pointers to pointers to rows*/
	t = (float ***) malloc((size_t) ((nrow+FI_END)*sizeof(float**)));
	if (!t) fatal("d3tensor", "Allocation failure 1 in f3tensor().");
	t += FI_END;
	t -= nrl;
	
	/*Alocate pointers to rows and set pointers to them*/
	t[nrl] = (float **) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(float*)));
	if (!t[nrl]) fatal("d3tensor", "Allocation failure 2 in f3tensor().");
	t[nrl] += FI_END;
	t[nrl] -= ncl;
	
	/*Allocate rows and set pointers to them*/
	t[nrl][ncl] = (float *) malloc((size_t) ((nrow*ncol*ndep+FI_END)*sizeof(float)));
	if (!t[nrl][ncl]) fatal("d3tensor", "Allocation failure 3 in f3tensor().");
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

int * __nonnull * __nonnull * __nonnull i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh)
{
	long i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ngh-ndl+1;
	int ***t;
	
	/*Allocate pointers to pointers to rows*/
	t = (int ***) malloc((size_t) ((nrow+FI_END)*sizeof(int**)));
	if (!t) fatal("i3tensor", "Allocation failure 1 in i3tensor().");
	t += FI_END;
	t -= nrl;
	
	/*Alocate pointers to rows and set pointers to them*/
	t[nrl] = (int **) malloc((size_t) ((nrow*ncol+FI_END)*sizeof(int*)));
	if (!t[nrl]) fatal("i3tensor", "Allocation failure 2 in i3tensor().");
	t[nrl] += FI_END;
	t[nrl] -= ncl;
	
	/*Allocate rows and set pointers to them*/
	t[nrl][ncl] = (int *) malloc((size_t) ((nrow*ncol*ndep+FI_END)*sizeof(int)));
	if (!t[nrl][ncl]) fatal("i3tensor", "Allocation failure 3 in i3tensor().");
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

void free_bvector(bool * __nonnull v, long nl, long nh){
    
    free((FREE_ARG) (v+nl-FI_END));
}

void free_dvector(double * __nonnull v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-FI_END));
}

void free_fvector(float * __nonnull v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-FI_END));
}

void free_ivector(int * __nonnull v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-FI_END));
}

void free_cdvector(double complex * __nonnull v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-FI_END));
}

void free_ulvector(unsigned long * __nonnull v, long nl, long nh)
{    
    free((FREE_ARG) (v+nl-FI_END));
}

void free_bmatrix(bool * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch) {
    
    free((FREE_ARG) (m[nrl]+ncl-FI_END));
	free((FREE_ARG) (m+nrl-FI_END));
}

void free_dmatrix(double * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch)
{ 
	free((FREE_ARG) (m[nrl]+ncl-FI_END));
	free((FREE_ARG) (m+nrl-FI_END));
}

void free_fmatrix(float * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-FI_END));
	free((FREE_ARG) (m+nrl-FI_END));
}

void free_imatrix(int * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-FI_END));
	free((FREE_ARG) (m+nrl-FI_END));
}

void free_cdmatrix(double complex * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch) {
    
    free((FREE_ARG) (m[nrl]+ncl-FI_END));
	free((FREE_ARG) (m+nrl-FI_END));
}

void free_d3tensor(double * __nonnull * __nonnull * __nonnull t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
	/*Free a double d3tensor allocated by d3tensor*/
	free((FREE_ARG) (t[nrl][ncl]+ndl-FI_END));
	free((FREE_ARG) (t[nrl]+ncl-FI_END));
	free((FREE_ARG) (t+nrl-FI_END));
}

void free_f3tensor(float * __nonnull * __nonnull * __nonnull t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
	/*Free a double d3tensor allocated by d3tensor*/
	free((FREE_ARG) (t[nrl][ncl]+ndl-FI_END));
	free((FREE_ARG) (t[nrl]+ncl-FI_END));
	free((FREE_ARG) (t+nrl-FI_END));
}

void free_i3tensor(int * __nonnull * __nonnull * __nonnull t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
	/*Free a double d3tensor allocated by d3tensor*/
	free((FREE_ARG) (t[nrl][ncl]+ndl-FI_END));
	free((FREE_ARG) (t[nrl]+ncl-FI_END));
	free((FREE_ARG) (t+nrl-FI_END));
}

void __attribute__((overloadable)) fatal(char head[], char message[])
{
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "######################### SAINO Program ########################\n");
    fprintf(stderr, "##                    A FATAL ERROR occured                   ##\n");
    fprintf(stderr, "##        Please read the error message for diagnostic        ##\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s: %s\n", head, message);
    fprintf(stderr, "\n");
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "################################################################\n");
    exit(-1);
}

void __attribute__((overloadable)) fatal(char head[], char message[], int n)
{
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "######################### SAINO Program ########################\n");
    fprintf(stderr, "##                    A FATAL ERROR occured                   ##\n");
    fprintf(stderr, "##        Please read the error message for diagnostic        ##\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s: %s %d\n", head, message, n);
    fprintf(stderr, "\n");
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "################################################################\n");
    exit(-1);
}

void __attribute__((overloadable)) fatal(char head[], char message[], double n)
{
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "######################### SAINO Program ########################\n");
    fprintf(stderr, "##                    A FATAL ERROR occured                   ##\n");
    fprintf(stderr, "##        Please read the error message for diagnostic        ##\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s: %s %f\n", head, message, n);
    fprintf(stderr, "\n");
    fprintf(stderr, "################################################################\n");
    fprintf(stderr, "################################################################\n");
    exit(-1);
}

void __attribute__((overloadable)) warning(char head[], char message[])
{
    printf("%s: %s\n", head, message);
}

void __attribute__((overloadable)) warning(char head[], char message[], int n)
{
    printf("%s: %s %d\n", head, message, n);
}

void __attribute__((overloadable)) warning(char head[], char message[], double n)
{
    printf("%s: %s %f\n", head, message, n);
}
