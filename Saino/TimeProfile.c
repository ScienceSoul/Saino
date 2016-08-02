//
//  TimeProfile.c
//  Saino
//
//  Created by Hakime Seddik on 08/02/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

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