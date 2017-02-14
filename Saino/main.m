//===----------------------------------------------------------------------===//
//  main.m
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright (c) 2011 ScienceSoul. All rights reserved.
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
