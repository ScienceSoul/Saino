//
//  main.m
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import <Cocoa/Cocoa.h>

#import "FEMJob.h"
#import "Utils.h"
#import "TimeProfile.h"

int main(int argc, char *argv[])
{
    int initialize = 0;
    double cp, rt;
    
    FEMJob *job = [[FEMJob alloc] init];
    
    cp = cputime();
    rt = realtime();
    
    [job runWithInitialize:initialize];
    [job deallocation];
    
    // TODO: add support for parallel run
    NSLog(@"SAINO: JOB TOTAL TIME (CPU,REAL): %lf %lf\n", cputime()-cp, realtime()-rt);
    
    NSString *dateString = [NSString stringWithCString:dateAndTime() encoding:NSASCIIStringEncoding];
    NSLog(@"SAINO JOB FINISHED AT: %@\n", dateString);
    
#ifdef TEST // This is a TEMPORALLY hack so that an application test works.
    @autoreleasepool {
        
        return NSApplicationMain(argc, (const char **)argv);
    }
#endif
    
    return 0;
}
