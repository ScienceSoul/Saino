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
@synthesize numberOfBodies = _numberOfBodies;
@synthesize numberOfBodyForces = _numberOfBodyForces;
@synthesize numberOfBoundaries = _numberOfBoundaries;
@synthesize numberOfICs = _numberOfICs;
@synthesize numberOfSolutions = _numberOfSolutions;
@synthesize numberOfEquations = _numberOfEquations;
@synthesize mesh = _mesh;
@synthesize boundaryID = _boundaryID;
@synthesize solutions = _solutions;
@synthesize meshes = _meshes;
@synthesize bodyForces = _bodyForces;
@synthesize boundaries = _boundaries;
@synthesize equations = _equations;
@synthesize simulation = _simulation;
@synthesize variables = _variables;

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _numberOfNodes = 0;
        _numberOfBulkElements = 0;
        _numberOfBoundaryElements = 0;
        _numberOfBodies = 0;
        _numberOfBodyForces = 0;
        _numberOfBoundaries = 0;
        _numberOfICs = 0;
        _numberOfSolutions = 0;
        _mesh = nil;
        
        _containers = (modelArraysContainer*)malloc(sizeof(modelArraysContainer) * 1 );
        _containers->FreeSurfaceNodes = NULL;
        _containers->BoundaryCurvatures = NULL;
    }
    
    return self;
}

-(void)deallocation {
    free(_containers);
}

#pragma mark Elements getter

-(Element_t *)getCurrentElement {
    
    return _currentElement;
}

-(modelArraysContainer*)getContainers {
    
    return _containers;
}


@end
