//
//  FEMStressAnalysis.m
//  Saino
//
//  Created by Seddik hakime on 19/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMStressAnalysisSolution.h"

@implementation FEMStressAnalysisSolution

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
    }
    
    return self;
}

-(void)deallocation:(FEMSolution * __nonnull)solution {
    
}

-(void)solutionComputer:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
}

@end
