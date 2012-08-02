//
//  FEMVariable.m
//  Saino
//
//  Created by Seddik hakime on 27/07/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMVariable.h"

@implementation FEMVariable

@synthesize next = _next;
@synthesize  name = _name;
@synthesize nameLength = _nameLength;
@synthesize dofs = _dofs;
@synthesize nonLinConverged = _nonLinConverged;
@synthesize steadyConverged = _steadyConverged;
@synthesize nonLinIter = _nonLinIter;
@synthesize norm = _norm;
@synthesize prevNorm = _prevNorm;
@synthesize nonLinChange = _nonLinChange;
@synthesize steadyChange = _steadyChange;
@synthesize valid = _valid;
@synthesize output = _output;
@synthesize valuesChanged = _valuesChanged;
@synthesize secondary = _secondary;

- (id)init
{
    self = [super init];
    if (self) {
        // NOTE: The classe below needs to be assigned to allocated FEMVariable classes when used
        _next = nil;
        
        _norm = 0.0;
        _prevNorm = 0.0;
        _nonLinChange = 0.0;
        _steadyChange = 0.0;
        _nonLinConverged = -1;
        _steadyConverged = -1;
    }
    
    return self;
}

-(variableArraysContainer*)getContainers {
    
    return _containers;
}


@end
