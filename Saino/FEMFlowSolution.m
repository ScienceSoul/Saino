//
//  FEMFlowSolution.m
//  Saino
//
//  Created by Seddik hakime on 15/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMFlowSolution.h"

@implementation FEMFlowSolution

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
    }
    
    return self;
}

-(void)deallocation:(FEMSolution *)solution {
    
}

-(void)solutionComputer:(FEMSolution *)solution model:(FEMModel *)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
}

@end
