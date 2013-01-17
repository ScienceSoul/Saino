//
//  FEMVariable.m
//  Saino
//
//  Created by Seddik hakime on 27/07/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMVariable.h"

@implementation FEMVariable

@synthesize  name = _name;
@synthesize primaryMesh = _primaryMesh;
@synthesize solution = _solution;
@synthesize nameLength = _nameLength;
@synthesize dofs = _dofs;
@synthesize nonLinConverged = _nonLinConverged;
@synthesize steadyConverged = _steadyConverged;
@synthesize nonLinIter = _nonLinIter;
@synthesize type = _type;
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
        
        _primaryMesh = nil;
        _solution = nil;
        _type = VARIABLE_ON_NODES;
        _norm = 0.0;
        _prevNorm = 0.0;
        _nonLinChange = 0.0;
        _steadyChange = 0.0;
        _nonLinConverged = -1;
        _steadyConverged = -1;
        
        _containers = (variableArraysContainer*)malloc(sizeof(variableArraysContainer));
        _containers->Perm = NULL;
        _containers->Values = NULL;
        _containers->PrevValues = NULL;
        _containers->PValues = NULL;
        _containers->NonLinValues = NULL;
        _containers->SteadyValues = NULL;
        _containers->EigenValues = NULL;
        _containers->EigenVectors = NULL;
        _containers->sizePerm = 0;
        _containers->sizeValues = 0;
        _containers->size1PrevValues = 0;
        _containers->size2PrevValues = 0;
        _containers->sizePValues = 0;
        _containers->sizeNonLinValues = 0;
        _containers->sizeSteadyValues = 0;
    }
    
    return self;
}

-(void)deallocation {
    free(_containers);
}

-(variableArraysContainer*)getContainers {
    
    return _containers;
}


@end
