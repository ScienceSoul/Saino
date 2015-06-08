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
@synthesize componentVariable = _componentVariable;
@synthesize secondary = _secondary;
@synthesize componentSecondaryVariable = _componentSecondaryVariable;

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
        _containers->ComponentValues = NULL;
        _containers->PrevValues = NULL;
        _containers->ComponentPrevValues = NULL;
        _containers->SecondaryToValues = NULL;
        _containers->ComponentSecondaryToValues = NULL;
        _containers->PValues = NULL;
        _containers->NonLinValues = NULL;
        _containers->SteadyValues = NULL;
        _containers->EigenValues = NULL;
        _containers->EigenVectors = NULL;
        _containers->ComponentEigenVectors = NULL;
        _containers->CValues = NULL;
        _containers->lowerLimitActive = NULL;
        _containers->upperLimitActive = NULL;
        _containers->sizePerm = 0;
        _containers->sizeValues = 0;
        _containers->sizeComponentValues = 0;
        _containers->size1PrevValues = 0;
        _containers->size2PrevValues = 0;
        _containers->size1ComponentPrevValues = 0;
        _containers->size2ComponentPrevValues = 0;
        _containers->sizeSecondaryToValues = 0;
        _containers->sizeComponentSecondaryToValues = 0;
        _containers->sizePValues = 0;
        _containers->sizeNonLinValues = 0;
        _containers->sizeSteadyValues = 0;
        _containers->sizeEigenValues = 0;
        _containers->size1EigenVectors = 0;
        _containers->size2EigenVectors = 0;
        _containers->size1ComponentEigenVectors = 0;
        _containers->size2ComponentEigenVectors = 0;
        _containers->sizeCValues = 0;
        _containers->sizeLowerLimitActive = 0;
        _containers->sizeUpperLimitActive = 0;
    }
    
    return self;
}

-(void)deallocation {
    int i;
    
    // When a variable is a component variable, a secondary or a component secondary variable, its containers dealocation
    // was not done during the deallocation of the mesh variables. Do it here.
    
    if (self.componentVariable == YES || self.secondary == YES || self.componentSecondaryVariable == YES) {
        
        _containers->Perm = NULL;
        
        if (_containers->ComponentValues != NULL) {
            free(_containers->ComponentValues);
            _containers->ComponentValues = NULL;
        }
        
        if (_containers->ComponentPrevValues != NULL) {
            for (i=0; _containers->size1ComponentPrevValues; i++) {
                free(_containers->ComponentPrevValues[i]);
            }
            free(_containers->ComponentPrevValues);
            _containers->ComponentPrevValues = NULL;
        }
        
        if (_containers->SecondaryToValues != NULL) {
            free(_containers->SecondaryToValues);
            _containers->SecondaryToValues = NULL;
        }
        
        if (_containers->ComponentSecondaryToValues != NULL) {
            free(_containers->ComponentSecondaryToValues);
            _containers->ComponentSecondaryToValues = NULL;
        }
        
        if (_containers->ComponentEigenVectors != NULL) {
            for (i=0; _containers->size1ComponentEigenVectors; i++) {
                free(_containers->ComponentEigenVectors[i]);
            }
            free(_containers->ComponentEigenVectors);
            _containers->ComponentEigenVectors = NULL;
        }
        free(_containers);
        _containers = NULL;
    }
}

/*
    Canonicalize the variable name.
    This means that it the variable name is of the form dummy[dummy1:n1 dummy2:n2] with n1 and n2 two integers,
    this method will return dummy (that is the standard naming form used by a solution variable).
    Otherwise it will return the name as it is.
*/
-(NSString *)canonicalizeName {
    
    NSRange ind = [self.name rangeOfString:@"["];
    if (ind.location != NSNotFound) {
        return [self.name substringToIndex:ind.location];
    } else {
        return self.name;
    }
}

-(variableArraysContainer*)getContainers {
    
    return _containers;
}


@end
