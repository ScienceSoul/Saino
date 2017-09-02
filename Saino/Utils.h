//===----------------------------------------------------------------------===//
//  Utils.h
//  Saino
//
//  Created by Seddik hakime on 05/07/11.
//  Copyright (c) 2011 ScienceSoul. All rights reserved
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

int __attribute__((overloadable)) max_array(int * _Nonnull a, int num_elements);
double __attribute__((overloadable)) max_array(double * _Nonnull a, int num_elements);
int __attribute__((overloadable)) min_array(int * _Nonnull a, int num_elements);
double __attribute__((overloadable)) min_array(double * _Nonnull a, int num_elements);

bool __attribute__((overloadable)) all(int * _Nonnull v, char mask, int val, int range);
bool __attribute__((overloadable)) all(float * _Nonnull v, char mask, float val, int range);
bool __attribute__((overloadable)) all(double * _Nonnull v, char mask, double val, int range);

bool __attribute__((overloadable)) all_in_range(int * _Nonnull v, char mask, int val, int start, int range);

bool __attribute__((overloadable)) any(int * _Nonnull v, char mask, int val, int range);
bool __attribute__((overloadable)) any(float * _Nonnull v, char mask, float val, int range);
bool __attribute__((overloadable)) any(double * _Nonnull v, char mask, double val, int range);

int __attribute__((overloadable)) count(int * _Nonnull v, char mask, int val, int range);
int __attribute__((overloadable)) count(float * _Nonnull v, char mask, float val, int range);
int __attribute__((overloadable)) count(double * _Nonnull v, char mask, double val, int range);

void __attribute__((overloadable)) sort(unsigned long n, int * _Nonnull arr);
void __attribute__((overloadable)) sort(unsigned long n, int * _Nonnull arr, int * _Nonnull brr);
void __attribute__((overloadable)) sort(unsigned long n, int * _Nonnull arr, float * _Nonnull brr);
void __attribute__((overloadable)) sort(unsigned long n, int * _Nonnull arr, double * _Nonnull brr);
void __attribute__((overloadable)) sort(unsigned long n, int * _Nonnull arr, double complex * _Nonnull brr); //TODO: Not implemented yet

void __attribute((overloadable)) reverse(int * _Nonnull arr, size_t narr);
void __attribute((overloadable)) cshift(int * _Nonnull arr, size_t narr, unsigned long shift);

char * _Nullable dateAndTime(void);
void startAdvanceOutput(char * _Nonnull solverName, char * _Nonnull outputType);
void advanceOutput(int t, int n, double * _Nullable dot_t, double * _Nullable percent_t);

int increaseStackSize(void);
