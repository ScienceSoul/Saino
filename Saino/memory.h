//===----------------------------------------------------------------------===//
//  memory.h
//  Saino
//
//  Created by Hakime SEDDIK on Thu Mar 11 2004.
//  Copyright (c) 2004 ScienceSoul. All rights reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include "Constructors.h"

void __attribute__((overloadable)) fatal(char head[_Nonnull]);
void __attribute__((overloadable)) fatal(char head[_Nonnull], char message[_Nonnull]);
void __attribute__((overloadable)) fatal(char head[_Nonnull], char message[_Nonnull], int n);
void __attribute__((overloadable)) fatal(char head[_Nonnull], char message[_Nonnull], double n);
void __attribute__((overloadable)) warning(char head[_Nonnull], char message[_Nonnull]);
void __attribute__((overloadable)) warning(char head[_Nonnull], char message[_Nonnull], int n);
void __attribute__((overloadable)) warning(char head[_Nonnull], char message[_Nonnull], double n);
bool * _Nonnull boolvec(long nl, long nh);
int * _Nonnull intvec(long nl, long nh);
unsigned long * _Nonnull ulongvec(long nl, long nh);
double * _Nonnull doublevec(long nl, long nh);
float * _Nonnull floatvec(long nl, long nh);
double complex * _Nonnull cdoublevec(long nl, long nh);
bool * _Nonnull * _Nonnull boolmatrix(long nrl, long nrh, long ncl, long nch);
int * _Nonnull * _Nonnull intmatrix(long nrl, long nrh, long ncl, long nch);
double * _Nonnull * _Nonnull doublematrix(long nrl, long nrh, long ncl, long nch);
float * _Nonnull * _Nonnull floatmatrix(long nrl, long nrh, long ncl, long nch);
double complex * _Nonnull * _Nonnull cdoublematrix(long nrl, long nrh, long ncl, long nch);
double * _Nonnull * _Nonnull * _Nonnull d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh);
float * _Nonnull * _Nonnull * _Nonnull f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh);
int * _Nonnull * _Nonnull * _Nonnull i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ngh);
void free_bvector(bool * _Nonnull v, long nl, long nh);
void free_dvector(double * _Nonnull v, long nl, long nh);
void free_fvector(float * _Nonnull v, long nl, long nh);
void free_ivector(int * _Nonnull v, long nl, long nh);
void free_cdvector(double complex * _Nonnull v, long nl, long nh);
void free_ulvector(unsigned long * _Nonnull v, long nl, long nh);
void free_bmatrix(bool * _Nonnull * _Nonnull m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double * _Nonnull * _Nonnull m, long nrl, long nrh, long ncl, long nch);
void free_fmatrix(float * _Nonnull * _Nonnull m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int * _Nonnull * _Nonnull m, long nrl, long nrh, long ncl, long nch);
void free_cdmatrix(double complex * _Nonnull * _Nonnull m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(double * _Nonnull * _Nonnull * _Nonnull t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_f3tensor(float * _Nonnull * _Nonnull * _Nonnull t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_i3tensor(int * _Nonnull * _Nonnull * _Nonnull t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);





