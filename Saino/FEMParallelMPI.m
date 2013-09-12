//
//  FEMParallelMPI.m
//  Saino
//
//  Created by Hakime Seddik on 09/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMParallelMPI.h"

@implementation FEMParallelMPI

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
    }
    
    return self;
}

-(double)parallelReductionOfValue:(double)r operArg:(int *)oper_arg{
    
    double rsum;
    
    rsum = r;
    //TODO: add support for parallel run
    
    return rsum;
}

@end
