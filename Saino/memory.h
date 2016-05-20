/*
 *  memory.h
 *  Saino
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

void __attribute__((overloadable)) fatal(char head[]);
void __attribute__((overloadable)) fatal(char head[], char message[]);
void __attribute__((overloadable)) fatal(char head[], char message[], int n);
void __attribute__((overloadable)) fatal(char head[], char message[], double n);
void __attribute__((overloadable)) warning(char head[], char message[]);
void __attribute__((overloadable)) warning(char head[], char message[], int n);
void __attribute__((overloadable)) warning(char head[], char message[], double n);
bool * __nonnull boolvec(long nl, long nh);
int * __nonnull intvec(long nl, long nh);
unsigned long * __nonnull ulongvec(long nl, long nh);
double * __nonnull doublevec(long nl, long nh);
float * __nonnull floatvec(long nl, long nh);
double complex * __nonnull cdoublevec(long nl, long nh);
bool * __nonnull * __nonnull boolmatrix(long nrl, long nrh, long ncl, long nch);
int * __nonnull * __nonnull intmatrix(long nrl, long nrh, long ncl, long nch);
double * __nonnull * __nonnull doublematrix(long nrl, long nrh, long ncl, long nch);
float * __nonnull * __nonnull floatmatrix(long nrl, long nrh, long ncl, long nch);
double complex * __nonnull * __nonnull cdoublematrix(long nrl, long nrh, long ncl, long nch);
double * __nonnull * __nonnull * __nonnull d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh);
float * __nonnull * __nonnull * __nonnull f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh);
int * __nonnull * __nonnull * __nonnull i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh);
void free_bvector(bool * __nonnull v, long nl, long nh);
void free_dvector(double * __nonnull v, long nl, long nh);
void free_fvector(float * __nonnull v, long nl, long nh);
void free_ivector(int * __nonnull v, long nl, long nh);
void free_cdvector(double complex * __nonnull v, long nl, long nh);
void free_ulvector(unsigned long * __nonnull v, long nl, long nh);
void free_bmatrix(bool * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch);
void free_fmatrix(float * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch);
void free_cdmatrix(double complex * __nonnull * __nonnull m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(double * __nonnull * __nonnull * __nonnull t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_f3tensor(float * __nonnull * __nonnull * __nonnull t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_i3tensor(int * __nonnull * __nonnull * __nonnull t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);





