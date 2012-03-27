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