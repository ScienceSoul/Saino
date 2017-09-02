//===----------------------------------------------------------------------===//
//  OpenCLUtils.h
//  Saino
//
//  Created by Seddik hakime on 25/10/2013.
//  Copyright (c) 2013 ScienceSoul. All rights reserved.
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
#include <OpenGL/OpenGL.h>

cl_platform_id _Nullable * _Nullable find_platforms(cl_uint * _Nullable numberofPlatforms, cl_uint * _Nullable numberOfDevices);
cl_device_id _Nullable find_single_device(void);
int device_info(cl_device_id _Nonnull device_id);
int device_stats(cl_device_id _Nonnull device_id);
int LoadFileIntoString(const char * _Nonnull filename, char * _Nonnull * _Nullable text, size_t * _Nonnull len);

