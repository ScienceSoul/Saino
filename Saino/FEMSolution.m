//===----------------------------------------------------------------------===//
//  FEMSolution.m
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

#import "FEMSolution.h"
#import "FEMUtilities.h"
#import "FEMFlowSolution.h"
#import "FEMMagneticInductionSolution.h"
#import "FEMStressAnalysisSolution.h"
#import "FEMMeshUpdateSolution.h"
#import "FEMHeatSolution.h"
#import "FEMHeatSolution_OpenCL.h"
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
@synthesize hasBuiltInSolution = _hasBuiltInSolution;
@synthesize builtInSolution = _builtInSolution;
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
        _hasBuiltInSolution = NO;
        _builtInSolution = nil;
        _plugInPrincipalClassInstance = nil;
        
        _plugInName = nil;
        _matrix = nil;
        _variable = nil;
        _mesh = nil;
        _normalTangentialName = nil;
        _exportedVariables = nil;
        
        _solutionInfo = [[NSMutableDictionary alloc] init];
        _valuesList = [[NSMutableArray alloc] init];
        
        _containers = (solutionArraysContainer*)malloc(sizeof(solutionArraysContainer));
        _containers->activeElements = NULL;
        _containers->defDofs = NULL;
        _containers->size1DefDofs = 0;
        _containers->size2DefDofs = 0;
    }
    
    return self;
}

-(void)deallocation {
    if (_builtInSolution != nil) {
        if ([(_solutionInfo)[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            FEMFlowSolution *solComputer = (FEMFlowSolution *)_builtInSolution;
            [solComputer deallocation:self];
        } else if ([(_solutionInfo)[@"equation"] isEqualToString:@"magnetic induction"] == YES) {
            FEMMagneticInductionSolution *solComputer = (FEMMagneticInductionSolution *)_builtInSolution;
            [solComputer deallocation:self];
        } else if ([(_solutionInfo)[@"equation"] isEqualToString:@"stress analysis"] == YES) {
            FEMStressAnalysisSolution *solComputer = (FEMStressAnalysisSolution *)_builtInSolution;
            [solComputer deallocation:self];
        } else if ([(_solutionInfo)[@"equation"] isEqualToString:@"mesh update"] == YES) {
            FEMMeshUpdateSolution *solComputer = (FEMMeshUpdateSolution *)_builtInSolution;
            [solComputer deallocation:self];
        } else if ([(_solutionInfo)[@"equation"] isEqualToString:@"heat equation"] == YES) {
            if ([(_solutionInfo)[@"parallel assembly"] boolValue] == YES) {
                FEMHeatSolution_OpenCL *solComputer = (FEMHeatSolution_OpenCL *)_builtInSolution;
                [solComputer deallocation:self];
            } else {
                FEMHeatSolution *solComputer = (FEMHeatSolution *)_builtInSolution;
                [solComputer deallocation:self];
            }
        }
    } else if (_plugInPrincipalClassInstance != nil) {
        id <SainoSolutionsComputer> instance = _plugInPrincipalClassInstance;
        [instance deallocation:self];
    }
    
    if (_containers->activeElements != NULL) {
        free_ivector(_containers->activeElements, 0, _containers->sizeActiveElements-1);
        _containers->activeElements = NULL;
    }
    if (_containers->defDofs != NULL) {
        free_imatrix(_containers->defDofs, 0, _containers->size1DefDofs-1, 0, _containers->size2DefDofs-1);
        _containers->defDofs = NULL;
    }
    free(_containers);
    _matrix = nil;
    _variable = nil;
    _mesh = nil;
    _normalTangentialName = nil;
    _exportedVariables = nil;
}

-(solutionArraysContainer * _Nonnull)getContainers {
    
    return _containers;
}

-(BOOL)instantiatePrincipalClassFromPlugIn:(NSBundle * _Nonnull)bundle {
    
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
