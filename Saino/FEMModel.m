//
//  FEMModel.m
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMModel.h"

@implementation FEMModel

@synthesize dimension = _dimension;
@synthesize numberOfNodes = _numberOfNodes;
@synthesize numberOfBulkElements = _numberOfBulkElements;
@synthesize numberOfBoundaryElements = _numberOfBoundaryElements;
@synthesize numberOfBodyForces = _numberOfBodyForces;
@synthesize numberOfBoundaries = _numberOfBoundaries;
@synthesize variables = _variables;
@synthesize bodyForces = _bodyForces;
@synthesize boundaries = _boundaries;
@synthesize simulations = _simulations;

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
    }
    
    return self;
}

@end
