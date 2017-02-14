//===----------------------------------------------------------------------===//
//  GPUUtils.h
//  Saino
//
//  Created by Hakime Seddik on 13/05/2016.
//  Copyright © 2016 ScienceSoul. All rights reserved.
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
#include <stdio.h>
#import <stdbool.h>
#include <stdarg.h>

// Defines what precision the GPU solver uses
// Value 0 is for double precsion (default mode)
// Value 1 is for single precision
static int precisionMode = 0;

#define createDeviceBuffer(t, ctx, options, size, hptr, error) ((t) == 1 ? \
                        clCreateBuffer(ctx, options, sizeof(cl_float)*size, hptr, &error) : \
                        clCreateBuffer(ctx, options, sizeof(cl_double)*size, hptr, &error) )

#define enqueueDeviceMapBuffer(t, queue, bf, block, map, offset, size, nevts, wl, evt, error) ((t) == 1 ? \
                clEnqueueMapBuffer(queue, bf, block, map, offset, sizeof(cl_float)*size, nevts, wl, evt, &error) : \
                clEnqueueMapBuffer(queue, bf, block, map, offset, sizeof(cl_double)*size, nevts, wl, evt, &error) )

#define setDeviceKernelArg(t, kern, idx, size, val) ((t) == 1 ? \
                clSetKernelArg(kern, idx, sizeof(cl_float)*size, &val) : \
                clSetKernelArg(kern, idx, sizeof(cl_double)*size, &val) )

#define mapmemcpy(t, target, source, size) ((t) == 1 ? \
                memcpy(target, source, size*sizeof(cl_float)) : \
                memcpy(target, source, size*sizeof(cl_double)) )

#define mapmemset(t, target, size) ((t) == 1 ? \
                memset(target, 0.0f, size*sizeof(cl_float)) : \
                memset(target, 0.0, size*sizeof(cl_double)) )

void setPrecision(bool single);
int precision(void);

typedef struct GlobalMemoryAllocationSize_t {
    
    int nb_char;
    int nb_uchar;
    int nb_short;
    int nb_ushort;
    int nb_int;
    int nb_uint;
    int nb_long;
    int nb_ulong;
    int nb_float;
    int nb_double;
} GlobalMemoryAllocationSize_t;

void initGlobalMemoryAllocation(GlobalMemoryAllocationSize_t * __nonnull allocationSize);
size_t computeGlobalMemoryAllocation(GlobalMemoryAllocationSize_t * __nonnull allocationSize);

void * __nonnull alloc_mem(int precision, size_t size);
int init_data(int precision, const char * __nonnull source_type, int count, ...);
int __attribute__((overloadable)) init_data_concat(const char * __nonnull source_type, float * __nonnull target, size_t size_target, bool reset, int count, ...);
int __attribute__((overloadable)) init_data_concat(const char * __nonnull source_type, double * __nonnull target, size_t size_target, bool reset, int count, ...);
