//===----------------------------------------------------------------------===//
//  GPUUtils.c
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

#include <OpenCL/OpenCL.h>
#include <malloc/malloc.h>
#include <string.h>
#include "GPUUtils.h"

static const char *type_int = "int";
static const char *type_float = "float";
static const char *type_double = "double";

static int indx = 0;
static size_t tot_size = 0;

void setPrecision(bool single) {
    
    if (single == true) {
        precisionMode = 1;
    }
}

int precision(void) {
    return precisionMode;
}

void initGlobalMemoryAllocation(GlobalMemoryAllocationSize_t * _Nonnull allocationSize) {

    *allocationSize = (GlobalMemoryAllocationSize_t){.nb_char=0, .nb_uchar=0, .nb_short=0, .nb_ushort=0,
                       .nb_int=0, .nb_uint=0, .nb_long=0, .nb_ulong=0, .nb_float=0, .nb_double=0};
}

size_t computeGlobalMemoryAllocation(GlobalMemoryAllocationSize_t * _Nonnull allocationSize) {
        
    return (sizeof(cl_char)*allocationSize->nb_char
            + sizeof(cl_uchar)*allocationSize->nb_uchar
            + sizeof(cl_short)*allocationSize->nb_short
            + sizeof(cl_ushort)*allocationSize->nb_ushort
            + sizeof(cl_int)*allocationSize->nb_int
            + sizeof(cl_uint)*allocationSize->nb_uint
            + sizeof(cl_long)*allocationSize->nb_long
            + sizeof(cl_ulong)*allocationSize->nb_ulong
            + sizeof(cl_float)*allocationSize->nb_float
            + sizeof(cl_double)*allocationSize->nb_double);
}

static bool check_type (const char *source_type) {
    
    bool correctType = false;
    
    if (strncmp(source_type, type_int, 6) == 0) {
        correctType = true;
    } else if (strncmp(source_type, type_float, 6) == 0) {
        correctType = true;
    } else if (strncmp(source_type, type_double, 6) == 0) {
        correctType = true;
    }
    
    return correctType;
}

static void check_size(va_list ap, int count, const char *source_type) {
    
    void *vpt;
    size_t size;
    
    for (int i=0; i<count; i+=2) {
        size = 0;
        if (strncmp(source_type, type_int, 6) == 0) {
            vpt = va_arg(ap, int *);
            size = va_arg(ap, size_t);
            tot_size = tot_size + size;
        } else if (strncmp(source_type, type_float, 6) == 0) {
            vpt = va_arg(ap, float *);
            size = va_arg(ap, size_t);
            tot_size = tot_size + size;
        }  else if (strncmp(source_type, type_double, 6) == 0) {
            vpt = va_arg(ap, double *);
            size = va_arg(ap, size_t);
            tot_size = tot_size + size;
        }
    }
}

void * _Nonnull alloc_mem(int precision, size_t size) {
    
    void *pt;
    
    if (precision == 1) {
        pt = malloc(sizeof(float) * size);
    } else {
        pt = malloc(sizeof(double) * size);
    }
    
    return pt;
}

/************************************************************************************************************
 
    Initialize a target with a source (both 1D arrays) using a precision mode (float or double) and a type
    describing the source (int, float, double). Several source-target pairs can be provided with a size,
    for example:
 
        init_data(.., int count, source1, target1, size1, source2, target2, size2, ...), 
 
    where count gives the total number of aguments following count.
 
 ***********************************************************************************************************/
int init_data(int precision, const char * _Nonnull source_type, int count, ...) {
    
    va_list ap;
    size_t size;
    int success = 0;
    
    if (check_type(source_type) == false) {
        fprintf(stdout, "init_data: incorrect type given to source array.\n");
        fprintf(stdout, "init_data: possible options are: int, float or double.\n");
        return -1;
    }
    
    va_start(ap, count);
    
    if (precision == 1) {
        float *target = NULL;
        for (int i=0; i<count; i+=3) {
            size = 0;
            if (strncmp(source_type, type_int, 6) == 0) {
                int *source = NULL;
                source = va_arg(ap, int *);
                target = va_arg(ap, float *);
                size = va_arg(ap, size_t);
                for (int j=0; j<size; j++) {
                    target[j] = (float)source[j];
                }
            } else if (strncmp(source_type, type_float, 6) == 0) {
                float *source = NULL;
                source = va_arg(ap, float *);
                target = va_arg(ap, float *);
                size = va_arg(ap, size_t);
                for (int j=0; j<size; j++) {
                    target[j] = source[j];
                }
            } else if (strncmp(source_type, type_double, 6) == 0) {
                double *source = NULL;
                source = va_arg(ap, double *);
                target = va_arg(ap, float *);
                size = va_arg(ap, size_t);
                for (int j=0; j<size; j++) {
                    target[j] = (float)source[j];
                }
            }
        }
    } else {
        double *target = NULL;
        for (int i=0; i<count; i+=3) {
            size = 0;
            if (strncmp(source_type, type_int, 6) == 0) {
                int *source = NULL;
                source = va_arg(ap, int *);
                target = va_arg(ap, double *);
                size = va_arg(ap, size_t);
                for (int j=0; j<size; j++) {
                    target[j] = (double)source[j];
                }
            } else if (strncmp(source_type, type_float, 6) == 0) {
                float *source = NULL;
                source = va_arg(ap, float *);
                target = va_arg(ap, double *);
                size = va_arg(ap, size_t);
                for (int j=0; j<size; j++) {
                    target[j] = (double)source[j];
                }
            } else if (strncmp(source_type, type_double, 6) == 0) {
                double *source = NULL;
                source = va_arg(ap, double *);
                target = va_arg(ap, double *);
                size = va_arg(ap, size_t);
                for (int j=0; j<size; j++) {
                    target[j] = source[j];
                }
            }
        }
    }
    
    va_end(ap);
    
    return success;
}

/************************************************************************************************************
 
    Initialize a target (float array) by concatenating several sources (all 1D arrays) using a type 
    describing the sources (int, float, double). Several sources can be provided and they will all be 
    concateneted to the same target with size size_target. The sources should also be provided with their
    respective size and the total size of all sources must not exceed the size of target, so for example:
 
        init_data_concat(.., size_target, reset, int count, source1, size1, source2, size2, ...),
 
    where size1+size2<=size_target and count gives the total number of aguments following count. 
    A reset option can be provided, if true, the indexing inside the target will be reset to 0.
 
 
***********************************************************************************************************/
int __attribute__((overloadable)) init_data_concat(const char * _Nonnull source_type, float * _Nonnull target, size_t size_target, bool reset, int count, ...) {
    
    va_list ap, ap2;
    size_t size;
    
    if (check_type(source_type) == false) {
        fprintf(stdout, "init_data_concat: incorrect type given to source array.\n");
        fprintf(stdout, "init_data_concat: possible options are: int, float or double.\n");
        return -1;
    }
    
    if (reset == true) {
        indx = 0;
        tot_size = 0;
    }
    
    va_start(ap, count);
    va_copy(ap2, ap);
    
    check_size(ap, count, source_type);
    if (tot_size > size_target) {
        fprintf(stdout, "init_data_concat: concatenated data bigger than the size of target.\n");
        return -1;
    }
    
    va_end(ap);
    
    for (int i=0; i<count; i+=2) {
        size = 0;
        if (strncmp(source_type, type_int, 6) == 0) {
            int *source;
            source = va_arg(ap2, int *);
            size = va_arg(ap2, size_t);
            for (int j=0; j<size; j++) {
                target[indx] = (float)source[j];
                indx++;
            }
        } else if (strncmp(source_type, type_float, 6) == 0) {
            float *source;
            source = va_arg(ap2, float *);
            size = va_arg(ap2, size_t);
            for (int j=0; j<size; j++) {
                target[indx] = source[j];
                indx++;
            }
        }  else if (strncmp(source_type, type_double, 6) == 0) {
            double *source;
            source = va_arg(ap2, double *);
            size = va_arg(ap2, size_t);
            for (int j=0; j<size; j++) {
                target[indx] = (float)source[j];
                indx++;
            }
        }
    }
    
    va_end(ap2);
    
    return 0;
}

/************************************************************************************************************
 
    Initialize a target (double array) by concatenating several sources (all 1D arrays) using a type
    describing the sources (int, float, double). Several sources can be provided and they will all be
    concateneted to the same target with size size_target. The sources should also be provided with their
    respective size and the total size of all sources must not exceed the size of target, so for example:
 
        init_data_concat(.., size_target, reset, int count, source1, size1, source2, size2, ...),
 
    where size1+size2<=size_target and count gives the total number of aguments following count.
    A reset option can be provided, if true, the indexing inside the target will be reset to 0.
 
***********************************************************************************************************/
int __attribute__((overloadable)) init_data_concat(const char * _Nonnull source_type, double * _Nonnull target, size_t size_target, bool reset, int count, ...) {
    
    va_list ap, ap2;
    size_t size;
    
    if (check_type(source_type) == false) {
        fprintf(stdout, "init_data_concat: incorrect type given to source array.\n");
        fprintf(stdout, "init_data_concat: possible options are: int, float or double.\n");
        return -1;
    }
    
    if (reset == true) {
        indx = 0;
        tot_size = 0;
    }
    
    va_start(ap, count);
    va_copy(ap2, ap);
    
    check_size(ap, count, source_type);
    if (tot_size > size_target) {
        fprintf(stdout, "init_data_concat: concatenated data bigger than the size of target.\n");
        return -1;
    }
    
    va_end(ap);
    
    for (int i=0; i<count; i+=2) {
        size = 0;
        if (strncmp(source_type, type_int, 6) == 0) {
            int *source;
            source = va_arg(ap2, int *);
            size = va_arg(ap2, size_t);
            for (int j=0; j<size; j++) {
                target[indx] = (double)source[j];
                indx++;
            }
        } else if (strncmp(source_type, type_float, 6) == 0) {
            float *source;
            source = va_arg(ap2, float *);
            size = va_arg(ap2, size_t);
            for (int j=0; j<size; j++) {
                target[indx] = (double)source[j];
                indx++;
            }
        }  else if (strncmp(source_type, type_double, 6) == 0) {
            double *source;
            source = va_arg(ap2, double *);
            size = va_arg(ap2, size_t);
            for (int j=0; j<size; j++) {
                target[indx] = source[j];
                indx++;
            }
        }
    }
    
    va_end(ap2);
    
    return 0;
}
