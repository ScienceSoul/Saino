//
//  main.m
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import <Cocoa/Cocoa.h>

#import "FEMJob.h"

int main(int argc, char *argv[])
{
    
//    @autoreleasepool {
//        
//        return NSApplicationMain(argc, (const char **)argv);
//    
//    }
    
    int initialize = 0;
    
    NSLog(@"Test Saino starting...\n");
    FEMJob *job = [[FEMJob alloc] init];
    [job runWithInitialize:initialize];
    
    return 0;
}
