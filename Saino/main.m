//
//  main.m
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import <Cocoa/Cocoa.h>

#import <SainoCore/SainoCore.h>

int main(int argc, char *argv[])
{
    
#ifdef TEST
    
    @autoreleasepool {
        
        return NSApplicationMain(argc, (const char **)argv);
    }
    
#else
    
    int initialize = 0;
    double cp, rt;
    
    FEMJob *job = [[FEMJob alloc] init];
    
    cp = cputime();
    rt = realtime();
    
    [job runWithInitialize:initialize];
    [job deallocation];
    
    // TODO: add support for parallel run
    fprintf(stdout, "SAINO: JOB TOTAL TIME (CPU,REAL): %lf %lf.\n", cputime()-cp, realtime()-rt);
    
    NSString *dateString = [NSString stringWithCString:dateAndTime() encoding:NSASCIIStringEncoding];
    fprintf(stdout, "SAINO JOB FINISHED AT: %s.\n", [dateString UTF8String]);
    
    return 0;
    
#endif
}
