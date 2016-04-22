//
//  OpenCLUtils.h
//  Saino
//
//  Created by Seddik hakime on 25/10/2013.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#include <OpenCL/OpenCL.h>
#include <OpenGL/OpenGL.h>

cl_platform_id __nullable * __nullable find_platforms(cl_uint * __nullable numberofPlatforms, cl_uint * __nullable numberOfDevices);
cl_device_id __nullable find_single_device(void);
int device_stats(cl_device_id __nonnull device_id);
int LoadFileIntoString(const char * __nonnull filename, char * __nonnull * __nullable text, size_t * __nonnull len);

