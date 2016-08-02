//
//  TimeProfile.h
//  Saino
//
//  Created by Hakime Seddik on 08/02/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#include <time.h>
#include <mach/mach_time.h>

double cputime();
double realtime();
double cpumemory();
double machcore(uint64_t endTime, uint64_t startTime);
