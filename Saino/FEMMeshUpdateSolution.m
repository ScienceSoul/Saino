//
//  FEMMeshUpdate.m
//  Saino
//
//  Created by Seddik hakime on 19/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMeshUpdateSolution.h"

@implementation FEMMeshUpdateSolution

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
    }
    
    return self;
}

-(void)fieldSolutionComputer:(FEMSolution *)solution model:(FEMModel *)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
}

@end
