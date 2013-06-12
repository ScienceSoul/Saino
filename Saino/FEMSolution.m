//
//  FEMSolution.m
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMSolution.h"
#import "FEMUtilities.h"
#import "Utils.h"


@implementation FEMSolution

#pragma mark Public methods...

@synthesize simulationID = _simulationID;
@synthesize coordinateSystemDimension  = _coordinateSystemDimension;
@synthesize timeOrder = _timeOrder;
@synthesize doneTime = _doneTime;
@synthesize order = _order;
@synthesize nOfEigenValues = _nOfEigenValues;
@synthesize solutionSolveWhen = _solutionSolveWhen;
@synthesize solutionMode = _solutionMode;
@synthesize numberOfActiveElements = _numberOfActiveElements;
@synthesize multigridLevel = _multigridLevel;
@synthesize multiGridTotal = _multiGridTotal;
@synthesize multiGridSweep = _multiGridSweep;
@synthesize alpha = _alpha;
@synthesize beta = _beta;
@synthesize dt = _dt;
@synthesize multigridSolution = _multigridSolution;
@synthesize multigridEqualPlit = _multigridEqualPlit;
@synthesize selector = _selector;
@synthesize plugInPrincipalClassInstance = _plugInPrincipalClassInstance;
@synthesize plugInName = _plugInName;
@synthesize matrix = _matrix;
@synthesize variable = _variable;
@synthesize mesh = _mesh;
@synthesize solutionInfo = _solutionInfo;
@synthesize normalTangentialName = _normalTangentialName;
@synthesize exportedVariables = _exportedVariables;
@synthesize valuesList = _valuesList;

#pragma mark Initializations

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _selector = NULL;
        _plugInPrincipalClassInstance = nil;
        
        _containers = (solutionArraysContainer*)malloc(sizeof(solutionArraysContainer));
        _containers->activeElements = NULL;
        _containers->defDofs = NULL;
        _containers->size1DefDofs = 0;
        _containers->size2DefDofs = 0;
    }
    
    return self;
}

-(void)deallocation {
    if (_containers->activeElements != NULL) {
        free_ivector(_containers->activeElements, 0, _containers->sizeActiveElements-1);
        _containers->activeElements = NULL;
    }
    if (_containers->defDofs != NULL) {
        free_imatrix(_containers->defDofs, 0, _containers->size1DefDofs-1, 0, _containers->size2DefDofs-1);
        _containers->defDofs = NULL;
    }
    free(_containers);
}

-(solutionArraysContainer*)getContainers {
    
    return _containers;
}

-(BOOL)instantiatePrincipalClassFromPlugIn:(NSBundle *)bundle {
    
    Class currPrincipalClass;
    FEMUtilities *utilities;
    
    utilities = [[FEMUtilities alloc] init];
    
    currPrincipalClass = [bundle principalClass];
    if (currPrincipalClass) {
        if ([utilities plugInClassIsValid:currPrincipalClass] == YES) {
            self.plugInPrincipalClassInstance = [[currPrincipalClass alloc] init];
            if (self.plugInPrincipalClassInstance) return YES;
        }
    }
    
    return NO;
}

@end
