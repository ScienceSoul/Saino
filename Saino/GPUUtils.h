//
//  GPUUtils.h
//  Saino
//
//  Created by Hakime Seddik on 13/05/2016.
//  Copyright © 2016 Institute of Low Temperature Science. All rights reserved.
//


#include <stdio.h>
#import <stdbool.h>

// Defines what precision the GPU solver uses
// Value 0 is for double precsion (default mode)
// Value 1 is for single precision
static int precisionMode = 0;

#define createDeviceBuffer(t, ctx, options, size, hptr, error) ((t) == 1 ? \
                        clCreateBuffer(ctx, options, sizeof(cl_float)*size, (float *)hptr, &error) : \
                        clCreateBuffer(ctx, options, sizeof(cl_double)*size, hptr, &error) )

#define enqueueDeviceMapBuffer(t, queue, bf, block, map, offset, size, nevts, wl, evt, error) ((t) == 1 ? \
                clEnqueueMapBuffer(queue, bf, block, map, offset, sizeof(cl_float)*size, nevts, wl, evt, &error) : \
                clEnqueueMapBuffer(queue, bf, block, map, offset, sizeof(cl_double)*size, nevts, wl, evt, &error) )

#define setDeviceKernelArg(t, kern, idx, size, val) ((t) == 1 ? \
                clSetKernelArg(kern, idx, sizeof(cl_float)*size, (float *)&val) : \
                clSetKernelArg(kern, idx, sizeof(cl_double)*size, &val) )

void setPrecision(bool single);
int precision(void);