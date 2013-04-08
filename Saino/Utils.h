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

int __attribute__((overloadable)) min(int x, int y);
int __attribute__((overloadable)) min(int x, int y, int z);
int __attribute__((overloadable)) min(int w, int x, int y, int z);
double __attribute__((overloadable)) min(double x, double y);

int __attribute__((overloadable)) max_array(int *a, int num_elements);
double __attribute__((overloadable)) max_array(double *a, int num_elements);
int __attribute__((overloadable)) min_array(int *a, int num_elements);
double __attribute__((overloadable)) min_array(double *a, int num_elements);

void __attribute__((overloadable)) sort(unsigned long n, int *arr);
bool __attribute__((overloadable)) all(int *v, char mask, int val, int range);
bool __attribute__((overloadable)) all(float *v, char mask, float val, int range);
bool __attribute__((overloadable)) all(double *v, char mask, double val, int range);

bool __attribute__((overloadable)) all_in_range(int *v, char mask, int val, int start, int range);

bool __attribute__((overloadable)) any(int *v, char mask, int val, int range);
bool __attribute__((overloadable)) any(float *v, char mask, float val, int range);
bool __attribute__((overloadable)) any(double *v, char mask, double val, int range);

int __attribute__((overloadable)) count(int *v, char mask, int val, int range);
int __attribute__((overloadable)) count(float *v, char mask, float val, int range);
int __attribute__((overloadable)) count(double *v, char mask, double val, int range);

void __attribute__((overloadable)) sort(unsigned long n, int *arr, int *brr); 
void __attribute__((overloadable)) sort(unsigned long n, int *arr, float *brr);
void __attribute__((overloadable)) sort(unsigned long n, int *arr, double *brr);
void __attribute__((overloadable)) sort(unsigned long n, int *arr, double complex *brr);

void __attribute((overloadable)) reverse(int *arr, size_t narr);
void __attribute((overloadable)) cshift(int *arr, size_t narr, unsigned long shift);

char *dateAndTime(void);