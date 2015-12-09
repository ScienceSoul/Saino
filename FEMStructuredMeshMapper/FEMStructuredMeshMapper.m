//
//  FEMStructuredMeshMapper.m
//  Saino
//
//  Created by Seddik hakime on 13/10/2015.
//  Copyright © 2015 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMStructuredMeshMapper.h"

@implementation FEMStructuredMeshMapper {
    
    int _nsize;
    int _pointerSize;
    int *_bottomPerm;
    int *_bottomPointer;
    int *_maskPerm;
    int *_topPerm;
    int *_topPointer;
    double *_bottomField;
    double *_coord;
    double *_field;
    double *_origCoord;
    double *_surface;
    double *_topField;
    BOOL _initialized;
    BOOL _maskExists;
    BOOL _visited;
}

- (id)init
{
    self = [super init];
    if (self) {
        _pointerSize = 0;
        _nsize = 0;
        _bottomPerm = NULL;
        _bottomPointer = NULL;
        _maskPerm = NULL;
        _topPerm = NULL;
        _topPointer = NULL;
        _bottomField = NULL;
        _coord = NULL;
        _field = NULL;
        _origCoord = NULL;
        _surface = NULL;
        _topField = NULL;
        _initialized = NO;
        _maskExists = NO;
        _visited = NO;
    }
    
    return self;
}

-(void)solutionComputer:(FEMSolution *)solution model:(FEMModel *)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
    int j, n, ibot, itop, bottomNode, topNode;
    double at0, at1, bottomVal, bottomValue0, heps, topVal, topValue0, xLoc, x0Bot, x0Loc, x0Top, wTop;
    NSString *varName;
    NSArray *bc;
    FEMVariable *var, *updateVar, *veloVar;
    Element_t *element;
    variableArraysContainer *meshUpdateContainers = NULL, *meshVeloContainers = NULL, *varContainers = NULL;
    listBuffer vector = { NULL, NULL, NULL, NULL, 0, 0, 0};
    BOOL displacementMode, found, gotUpdateVar, gotVeloVar;
    
    NSLog(@"FEMStructuredMeshMapper:solutionComputer: ----------------------------------------------\n");
    NSLog(@"FEMStructuredMeshMapper:solutionComputer: performing mapping on a structured mesh.\n");
    NSLog(@"FEMStructuredMeshMapper:solutionComputer: ----------------------------------------------\n");
    
    FEMCore *core = [FEMCore sharedCore];
    
    if (_initialized == NO) {
        FEMMeshUtils *meshUtilities = [[FEMMeshUtils alloc] init];
        _pointerSize = [meshUtilities detectExtrudedStructureMesh:solution.mesh solution:solution model:model externVariable:var needExternVariable:YES isTopActive:YES topNodePointer:_topPointer isBottomActive:YES bottomNodePointer:_bottomPointer isUpActive:NO upNodePointer:NULL isDownActive:NO downNodePointer:NULL numberOfLayers:NULL nodeLayer:NULL];
        varContainers = var.getContainers;
        _maskExists = (varContainers->Perm != NULL) ? YES : NO;
        if (_maskExists == YES) _maskPerm = varContainers->Perm;
        _coord = varContainers->Values;
        _nsize = varContainers->sizeValues;
        _initialized = YES;
        
        if ((solution.solutionInfo)[@"mesh update variable"] != nil) {
            _origCoord = doublevec(0, _nsize-1);
            memcpy(_origCoord, _coord, _nsize*sizeof(double));
        }
    }
    at0 = cputime();
    
    // End of initialization
    
    // Get either variable or constant values for top surface
    topNode = 0;
    topValue0 = 0.0;
    if ((solution.solutionInfo)[@"top surface level"] != nil) {
        topValue0 = [(solution.solutionInfo)[@"top surface level"] doubleValue];
        topNode = 1;
    } else {
        if ((solution.solutionInfo)[@"top surface variable name"] != nil) {
            varName = (solution.solutionInfo)[@"top surface variable name"];
            FEMUtilities *utilities = [[FEMUtilities alloc] init];
            varContainers = NULL;
            var = [utilities getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:nil info:&found];
            if (found == YES) {
                varContainers = var.getContainers;
                if (var.dofs != 1) {
                    errorfunct("FEMStructuredMeshMapper:solutionComputer", "Top surface should have only 1 dof.");
                } else {
                    if (varContainers != NULL) {
                        _topField = varContainers->Values;
                        _topPerm = varContainers->Perm;
                    }
                    topNode = 2;
                }
            } else {
                NSLog(@"FEMStructuredMeshMapper:solutionComputer: top surface variable is missing: %@.\n", varName);
                errorfunct("FEMStructuredMeshMapper:solutionComputer", "Program terminating now...");
            }
        }
    }
    
    if (topNode == 0) {
        FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
        if ([listUtilities listCheckPresentAnyBoundaryCondition:model name:@"top surface"] == YES) {
            topNode = 3;
            if (_field == NULL) {
                _field = doublevec(0, solution.mesh.maxElementNodes-1);
                _surface = doublevec(0, solution.mesh.maxElementNodes-1);
                memset(_field, 0.0, solution.mesh.maxElementNodes*sizeof(double) );
                memset(_surface, 0.0, solution.mesh.maxElementNodes*sizeof(double) );
            }
            for (int elem=0; elem<solution.mesh.numberOfBoundaryElements; elem++) {
                element = [core getBoundaryElement:solution atIndex:elem];
                bc = [core getBoundaryCondition:model forElement:element];
                if (bc != nil) {
                    n = element->Type.NumberOfNodes;
                    found = [listUtilities listGetReal:model inArray:bc forVariable:@"top surface" numberOfNodes:n indexes:element->NodeIndexes buffer:&vector minValue:NULL maxValue:NULL];
                    if (found == YES) memcpy(_surface, vector.vector, n*sizeof(double));
                    if (found == YES) {
                        for (int i=0; i<n; i++) {
                            _field[element->NodeIndexes[i]] = _surface[i];
                        }
                    }
                }
            }
        }
    }
    
    // Get either variable or constant values for bottom surface
    bottomNode = 0;
    bottomValue0 = 0.0;
    if ((solution.solutionInfo)[@"bottom surface level"] != nil) {
        bottomValue0 = [(solution.solutionInfo)[@"bottom surface level"] doubleValue];
        bottomNode = 1;
    } else {
        if ((solution.solutionInfo)[@"bottom surface variable name"] != nil) {
            varName = (solution.solutionInfo)[@"bottom surface variable name"];
            FEMUtilities *utilities = [[FEMUtilities alloc] init];
            varContainers = NULL;
            var = [utilities getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:nil info:&found];
            if (found == YES) {
                varContainers = var.getContainers;
                if (var.dofs != 1) {
                    errorfunct("FEMStructuredMeshMapper:solutionComputer", "Bottom surface should have only 1 dof.");
                } else {
                    if (varContainers != NULL) {
                        _bottomField = varContainers->Values;
                        _bottomPerm = varContainers->Perm;
                    }
                    bottomNode = 2;
                }
            } else {
                NSLog(@"FEMStructuredMeshMapper:solutionComputer: bottom surface variable is missing: %@.\n", varName);
                errorfunct("FEMStructuredMeshMapper:solutionComputer", "Program terminating now...");
            }
        }
    }
    
    if (bottomNode == 0) {
        FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
        if ([listUtilities listCheckPresentAnyBodyForce:model name:@"bottom surface"] == YES) {
            bottomNode = 3;
            if (_field == NULL) {
                _field = doublevec(0, solution.mesh.maxElementNodes-1);
                _surface = doublevec(0, solution.mesh.maxElementNodes-1);
                memset(_field, 0.0, solution.mesh.maxElementNodes*sizeof(double) );
                memset(_surface, 0.0, solution.mesh.maxElementNodes*sizeof(double) );
            }
            for (int elem=0; elem<solution.mesh.numberOfBoundaryElements; elem++) {
                element = [core getBoundaryElement:solution atIndex:elem];
                bc = [core getBoundaryCondition:model forElement:element];
                if (bc != nil) {
                    n = element->Type.NumberOfNodes;
                    found = [listUtilities listGetReal:model inArray:bc forVariable:@"bottom surface" numberOfNodes:n indexes:element->NodeIndexes buffer:&vector minValue:NULL maxValue:NULL];
                    if (found == YES) memcpy(_surface, vector.vector, n*sizeof(double));
                    if (found == YES) {
                        for (int i=0; i<n; i++) {
                            _field[element->NodeIndexes[i]] = _surface[i];
                        }
                    }
                }
            }
        }
    }
    
    // Get the velocity variable component
    gotVeloVar = NO;
    if ((solution.solutionInfo)[@"mesh velocity variable"] != nil) {
        varName = (solution.solutionInfo)[@"mesh velocity variable"];
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        veloVar = [utilities getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:nil info:&found];
        if (found == YES) {
            meshVeloContainers = veloVar.getContainers;
            if (veloVar.dofs == 1) {
                gotVeloVar = YES;
            } else {
                errorfunct("FEMStructuredMeshMapper:solutionComputer", "The size of mesh velocity must be one.");
            }
        } else {
            NSLog(@"FEMStructuredMeshMapper:solutionComputer: the variable does not exist: %@.\n", varName);
            errorfunct("FEMStructuredMeshMapper:solutionComputer", "Program terminating now...");
        }
    }
    
    // Get the mesh update variable component
    gotUpdateVar = NO;
    if ((solution.solutionInfo)[@"mesh update variable"] != nil) {
        varName = (solution.solutionInfo)[@"mesh update variable"];
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        updateVar = [utilities getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:nil info:&found];
        if (found == YES) {
            meshUpdateContainers = updateVar.getContainers;
            if (updateVar.dofs == 1) {
                gotUpdateVar = YES;
            } else {
                errorfunct("FEMStructuredMeshMapper:solutionComputer", "The size of mesh update must be one.");
            }
        } else {
            NSLog(@"FEMStructuredMeshMapper:solutionComputer: the variable does not exist: %@.\n", varName);
            errorfunct("FEMStructuredMeshMapper:solutionComputer", "Program terminating now...");
        }
    }
    
    if ((solution.solutionInfo)[@"displacement mode"] != nil) {
        displacementMode = [(solution.solutionInfo)[@"displacement mode"] boolValue];
    } displacementMode = NO;
    
    // Get the new mapping using linear interpolation from bottom to top
    if ((solution.solutionInfo)[@"minimum mesh height"] != nil) {
        heps = [(solution.solutionInfo)[@"minimum mesh height"] doubleValue];
    } heps = DBL_EPSILON;
    
    topVal = 0.0;
    bottomVal = 0.0;
    for (int i=0; i<_nsize; i++) {
        j = i;
        if (_maskExists) {
            j = _maskPerm[i];
            if (j < 0) continue;
        }
        itop = _topPointer[i];
        ibot = _bottomPointer[i];
        if (i == itop || i == ibot) continue;
        x0Top = _coord[itop];
        x0Bot = _coord[ibot];
        x0Loc = _coord[i];
        
        wTop = (x0Loc - x0Bot) / (x0Top - x0Bot);
        
        if (topNode == 1) {
            topVal = topValue0;
        } else if (topNode == 2) {
            topVal = _topField[_topPerm[itop]];
        } else if (topNode == 3) {
            topVal = _field [itop];
        } else {
            if (displacementMode == YES) {
                topVal = 0.0;
            } else {
                topVal = x0Top;
            }
        }
        
        if (bottomNode == 1) {
            bottomVal = bottomValue0;
        } else if (bottomNode == 2) {
            bottomVal = _bottomField[_bottomPerm[ibot]];
        } else if (bottomNode == 3) {
            bottomVal = _field[ibot];
        } else {
            if (displacementMode == YES) {
                bottomVal = 0.0;
            } else {
                bottomVal = x0Bot;
            }
        }
        
        if (topVal - bottomVal < heps) {
            Nodes_t *nodes = solution.mesh.getNodes;
            NSLog(@"FEMStructuredMeshMapper:solutionComputer: node %d, height %f, w %f.\n", i, _coord[i], wTop);
            NSLog(@"FEMStructuredMeshMapper:solutionComputer: position %f %f %f.\n", nodes->x[i], nodes->y[i], nodes->z[i]);
            NSLog(@"FEMStructuredMeshMapper:solutionComputer: topVal %f, botVal %f, dVal %f.\n", topVal, bottomVal, topVal-bottomVal);
            NSLog(@"FEMStructuredMeshMapper:solutionComputer: top and bottom get tangled: %f %f.\n", topVal, bottomVal);
            errorfunct("FEMStructuredMeshMapper:solutionComputer", "Program terminating now...");
        }
        
        xLoc = wTop * topVal + (1.0 - wTop) * bottomVal;
        
        if (displacementMode == YES) {
            if (gotVeloVar == YES) meshVeloContainers->Values[meshVeloContainers->Perm[i]] = xLoc / timeStep;
            _coord[i] = _coord[i] + xLoc;
        } else {
            if (gotVeloVar == YES) meshVeloContainers->Values[meshVeloContainers->Perm[i]] = (xLoc - _coord[i]) / timeStep;
            _coord[i] = xLoc;
        }
        if (gotUpdateVar == YES) meshUpdateContainers->Values[meshUpdateContainers->Perm[i]] = _coord[i] - _origCoord[i];
    }
    
    // Map to top and bottom themselves for very last
    for (int i=0; i<_nsize; i++) {
        if (_maskExists == YES) {
            if (_maskPerm[i] < 0) continue;
        }
        xLoc = 0.0;
        if (i == _topPointer[i]) {
            if (topNode == 1) {
                xLoc = topVal;
            } else if (topNode == 2) {
                xLoc = _topField[_topPerm[i]];
            } else if (topNode == 3) {
                xLoc = _field[i];
            }
        } else if (i == _bottomPointer[i]) {
            if (bottomNode == 1) {
                xLoc = bottomVal;
            } else if (bottomNode == 2) {
                xLoc = _bottomField[_bottomPerm[i]];
            } else if (bottomNode == 3) {
                xLoc = _field[i];
            }
        } else {
            continue;
        }
        
        if (displacementMode == YES) {
            if (gotVeloVar == YES) meshVeloContainers->Values[meshVeloContainers->Perm[i]] = xLoc / timeStep;
            _coord[i] = _coord[i] + xLoc;
        } else {
            if (gotVeloVar == YES) meshVeloContainers->Values[meshVeloContainers->Perm[i]] = (xLoc - _coord[i]) / timeStep;
            _coord[i] = xLoc;
        }
        if (gotUpdateVar == YES) meshUpdateContainers->Values[meshUpdateContainers->Perm[i]] = _coord[i] - _origCoord[i];
    }
    
    if (gotVeloVar == YES && _visited == NO) {
        if ((solution.solutionInfo)[@"mesh velocity first zero"] != nil) {
            if ([(solution.solutionInfo)[@"mesh velocity first zero"] boolValue] == YES)
                memset(meshVeloContainers->Values, 0.0, meshVeloContainers->sizeValues*sizeof(double) );
        }
    }
    at1 = cputime();
    
    NSLog(@"FEMStructuredMeshMapper:solutionComputer: active coordinate mapping time: %f.\n", at1-at0);
    
    _visited = YES;
    
    if (vector.vector != NULL) {
        free_dvector(vector.vector, 0, vector.m-1);
        vector.vector = NULL;
    }
}

-(void)deallocation:(FEMSolution *)solution {
    if (_origCoord != NULL) {
        free_dvector(_origCoord, 0, _nsize-1);
    }
    if (_field != NULL) {
        free_dvector(_field, 0, solution.mesh.maxElementNodes-1);
    }
    if (_surface != NULL) {
        free_dvector(_surface, 0, solution.mesh.maxElementNodes-1);
    }
    if (_bottomPointer != NULL) {
        free_ivector(_bottomPointer, 0, _pointerSize-1);
    }
    if (_topPointer != NULL) {
        free_ivector(_topPointer, 0, _pointerSize-1);
    }
}

@end
