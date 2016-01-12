//
//  Utils.h
//  Saino
//
//  Created by Seddik hakime on 05/07/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>
#include <limits.h>
#include <time.h>

#include "memory.h"

int __attribute__((overloadable)) max(int x, int y);
int __attribute__((overloadable)) max(int x, int y, int z);
int __attribute__((overloadable)) max(int w, int x, int y, int z);
double __attribute__((overloadable)) max(double x, double y);
double __attribute__((overloadable)) max(double x, double y, double z);

int __attribute__((overloadable)) min(int x, int y);
int __attribute__((overloadable)) min(int x, int y, int z);
int __attribute__((overloadable)) min(int w, int x, int y, int z);
double __attribute__((overloadable)) min(double x, double y);
double __attribute__((overloadable)) min(double x, double y, double z);

int __attribute__((overloadable)) max_array(int * __nonnull a, int num_elements);
double __attribute__((overloadable)) max_array(double * __nonnull a, int num_elements);
int __attribute__((overloadable)) min_array(int * __nonnull a, int num_elements);
double __attribute__((overloadable)) min_array(double * __nonnull a, int num_elements);

bool __attribute__((overloadable)) all(int * __nonnull v, char mask, int val, int range);
bool __attribute__((overloadable)) all(float * __nonnull v, char mask, float val, int range);
bool __attribute__((overloadable)) all(double * __nonnull v, char mask, double val, int range);

bool __attribute__((overloadable)) all_in_range(int * __nonnull v, char mask, int val, int start, int range);

bool __attribute__((overloadable)) any(int * __nonnull v, char mask, int val, int range);
bool __attribute__((overloadable)) any(float * __nonnull v, char mask, float val, int range);
bool __attribute__((overloadable)) any(double * __nonnull v, char mask, double val, int range);

int __attribute__((overloadable)) count(int * __nonnull v, char mask, int val, int range);
int __attribute__((overloadable)) count(float * __nonnull v, char mask, float val, int range);
int __attribute__((overloadable)) count(double * __nonnull v, char mask, double val, int range);

void __attribute__((overloadable)) sort(unsigned long n, int * __nonnull arr);
void __attribute__((overloadable)) sort(unsigned long n, int * __nonnull arr, int * __nonnull brr);
void __attribute__((overloadable)) sort(unsigned long n, int * __nonnull arr, float * __nonnull brr);
void __attribute__((overloadable)) sort(unsigned long n, int * __nonnull arr, double * __nonnull brr);
void __attribute__((overloadable)) sort(unsigned long n, int * __nonnull arr, double complex * __nonnull brr); //TODO: Not implemented yet

void __attribute((overloadable)) reverse(int * __nonnull arr, size_t narr);
void __attribute((overloadable)) cshift(int * __nonnull arr, size_t narr, unsigned long shift);

char * __nullable dateAndTime(void);
void startAdvanceOutput(char * __nonnull solverName, char * __nonnull outputType);
void advanceOutput(int t, int n, double * __nullable dot_t, double * __nullable percent_t);