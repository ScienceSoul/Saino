//===----------------------------------------------------------------------===//
//  TimeProfile.c
//  Saino
//
//  Created by Hakime Seddik on 08/02/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
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

#include <stdio.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "TimeProfile.h"

static struct rusage usage;

static struct timeval tp;
static struct timezone tzp;

double cputime () {
    
    getrusage( RUSAGE_SELF, &usage );
    return (double) usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1.0e-6;
}

double realtime () {
    
    gettimeofday( &tp,&tzp );
    return (double) tp.tv_sec + tp.tv_usec*1.0e-6;
}

double cpumemory () {
    
    getrusage( RUSAGE_SELF, &usage );
    return (double) 1.0 * usage.ru_maxrss;
}

double machcore(uint64_t endTime, uint64_t startTime){
    
    uint64_t difference = endTime - startTime;
    static double conversion = 0.0;
    double value = 0.0;
    
    if( 0.0 == conversion )
    {
        mach_timebase_info_data_t info;
        kern_return_t err = mach_timebase_info( &info );
        
        if( 0 == err ){
            /* seconds */
            conversion = 1e-9 * (double) info.numer / (double) info.denom;
            /* nanoseconds */
            //conversion = (double) info.numer / (double) info.denom;
        }
    }
    
    value = conversion * (double) difference;
    
    return value;
}
