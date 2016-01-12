//
//  OpenCLUtils.h
//  Saino
//
//  Created by Seddik hakime on 25/10/2013.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#include <OpenCL/OpenCL.h>

#ifndef Saino_OpenCLUtils_h
#define Saino_OpenCLUtils_h

int device_stats(cl_device_id __nonnull device_id);
int LoadFileIntoString(const char * __nonnull filename, char * __nonnull * __nullable text, size_t * __nonnull len);

#endif
