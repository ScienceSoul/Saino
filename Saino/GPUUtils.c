//
//  GPUUtils.c
//  Saino
//
//  Created by Hakime Seddik on 13/05/2016.
//  Copyright © 2016 Institute of Low Temperature Science. All rights reserved.
//

#include <OpenCL/OpenCL.h>
#include "GPUUtils.h"

void setPrecision(bool single) {
    
    if (single == true) {
        precisionMode = 1;
    }
}

int precision(void) {
    return precisionMode;
}

void initGlobalMemoryAllocation(GlobalMemoryAllocationSize_t * __nonnull allocationSize) {

    *allocationSize = (GlobalMemoryAllocationSize_t){.nb_char=0, .nb_uchar=0, .nb_short=0, .nb_ushort=0,
                       .nb_int=0, .nb_uint=0, .nb_long=0, .nb_ulong=0, .nb_float=0, .nb_double=0};
}

size_t computeGlobalMemoryAllocation(GlobalMemoryAllocationSize_t * __nonnull allocationSize) {
        
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
