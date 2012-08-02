/*
 *  memory.h
 *  FlowErosModel alpha 3.
 *
 *  Created by Hakime SEDDIK on Thu Mar 11 2004.
 *  Copyright (c) 2004 Institute of Low Temperature Science. All rights reserved.
 *
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>

#include "Constructors.h"

void errorfunct(char head[], char message[], ...);
void warnfunct(char head[], char message[], ...);
bool *boolvec(long nl, long nh);
int *intvec(long nl, long nh);
unsigned long *ulongvec(long nl, long nh);
double *doublevec(long nl, long nh);
float *floatvec(long nl, long nh);
double complex *cdoublevec(long nl, long nh);
bool **boolmatrix(long nrl, long nrh, long ncl, long nch);
double **doublematrix(long nrl, long nrh, long ncl, long nch);
float **floatmatrix(long nrl, long nrh, long ncl, long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh);
int **intmatrix(long nrl, long nrh, long ncl, long nch);
int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh);
void free_bvector(bool *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_fvector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cdvector(double complex *v, long nl, long nh);
void free_ulvector(unsigned long *v, long nl, long nh);
void free_bmatrix(bool **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_fmatrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);





