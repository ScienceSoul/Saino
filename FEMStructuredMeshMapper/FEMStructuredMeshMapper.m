//===----------------------------------------------------------------------===//
//  FEMStructuredMeshMapper.m
//  Saino
//
//  Created by Seddik hakime on 13/10/2015.
//  Copyright © 2015 ScienceSoul. All rights reserved.
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

#import "FEMStructuredMeshMapper.h"

@implementation FEMStructuredMeshMapper {
    
    int _nsize;
    int * __nullable _bottomPerm;
    int * __nullable _bottomPointer;
    int * __nullable _downPointer;
    int * __nullable _maskPerm;
    int * __nullable _topPerm;
    int * __nullable _topPointer;
    int * __nullable _midPointer;
    int * __nullable _nodeLayer;
    int * __nullable _upPointer;
    double * __nullable _bottomField;
    double * __nullable _coord;
    double * __nullable _field;
    double * __nullable _origCoord;
    double * __nullable _topField;
    BOOL _initialized;
    BOOL _maskExists;
    BOOL _visited;
}

- (id)init
{
    self = [super init];
    if (self) {
        _nsize = 0;
        _bottomField = NULL;
        _bottomPerm = NULL;
        _bottomPointer = NULL;
        _coord = NULL;
        _downPointer = NULL;
        _field = NULL;
        _maskPerm = NULL;
        _midPointer = NULL;
        _nodeLayer = NULL;
        _origCoord = NULL;
        _topField = NULL;
        _topPerm = NULL;
        _topPointer = NULL;
        _upPointer = NULL;
        _initialized = NO;
        _maskExists = NO;
        _visited = NO;
    }
    
    return self;
}

/**************************************************************************
 
    This solution corresponds to Elmer from git on October 27 2015

**************************************************************************/

-(void)solutionComputer:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
    int j, n, ibot, imid, itop, bottomNode, topNode;
    int *tangleMaskPerm = NULL;
    double at0, at1, bottomVal, bottomValue0, heps, midVal, topVal, topValue0, xLoc, x0Bot, x0Loc, x0Mid, x0Top, wTop;
    double *tangleMask = NULL;
    NSString *varName;
    NSArray *bc;
    FEMVariable *var, *updateVar, *veloVar;
    Element_t *element;
    variableArraysContainer *maskVarContainers = NULL, *meshUpdateContainers = NULL, *meshVeloContainers = NULL, *varContainers = NULL;
    listBuffer vector = { NULL, NULL, NULL, NULL, 0, 0, 0};
    BOOL computeTangleMask = NO, displacementMode, found, gotUpdateVar, gotVeloVar, midLayerExists = NO;
    
    fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: ----------------------------------------------\n");
    fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: performing mapping on a structured mesh.\n");
    fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: ----------------------------------------------\n");
    
    FEMCore *core = [FEMCore sharedCore];
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    
    BOOL reinitialize;
    if (solution.solutionInfo[@"always detect structure"] != nil) {
        reinitialize = [solution.solutionInfo[@"always detect structure"] boolValue];
    } else reinitialize = NO;
    
    // Initialization
    if (_initialized == NO || reinitialize == YES) {
        if (_bottomPointer != NULL) free_ivector(_bottomPointer, 0, _nsize-1);
        if (_topPointer != NULL) free_ivector(_topPointer, 0, _nsize-1);
        FEMMeshUtils *meshUtilities = [[FEMMeshUtils alloc] init];
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        if (solution.solutionInfo[@"mapping mask variable"] != nil) {
            NSString *maskName = solution.solutionInfo[@"mapping mask variable"];
            FEMVariable *maskVar = [utilities getVariableFrom:solution.mesh.variables model:model name:maskName onlySearch:NULL maskName:nil info:&found];
            if (found == YES) {
                maskVarContainers = maskVar.getContainers;
                if (maskVarContainers != NULL) {
                    if (maskVarContainers->Perm != NULL) _maskExists = YES;
                }
            }
        }
        if (_maskExists == YES) {
            _nsize = maskVarContainers->sizePerm;
            fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: applying mask of size: %d.\n", _nsize);
        } else {
            _nsize = solution.mesh.numberOfNodes;
            fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: applying mask to the whole mesh.\n");
        }
        _topPointer = intvec(0, _nsize-1);
        _bottomPointer = intvec(0, _nsize-1);
        _upPointer = intvec(0, _nsize-1);
        _downPointer = intvec(0, _nsize-1);
        _midPointer = intvec(0, _nsize-1);
        var = [meshUtilities detectExtrudedStructureMesh:solution.mesh solution:solution model:model numberofNodes:_nsize ifMask:_maskExists mask:maskVarContainers isTopActive:YES topNodePointer:_topPointer isBottomActive:YES bottomNodePointer:_bottomPointer isUpActive:NO upNodePointer:_upPointer isDownActive:NO downNodePointer:_downPointer numberOfLayers:NULL isMidNode:YES midNodePointer:_midPointer midLayerExists:&midLayerExists isNodeLayer:NO nodeLayer:NULL];
        if (_maskExists == YES) _maskPerm = maskVarContainers->Perm;
        if (var == nil) {
            fprintf(stderr, "FEMStructuredMeshMapper:solutionComputer: problem in retrieving variable used for the mapping.\n");
            fatal("FEMStructuredMeshMapper:solutionComputer");
        } else varContainers = var.getContainers;
        _coord = varContainers->Values;
        
        if (_origCoord != NULL) free_dvector(_origCoord, 0, _nsize-1);
        _origCoord = doublevec(0, _nsize-1);
        
        _initialized = YES;
    }
    memcpy(_origCoord, _coord, _nsize*sizeof(double));
    at0 = cputime();
    
    // Detangling stuff
    double minHeight = 0.0;
    BOOL deTangle = NO;
    if (solution.solutionInfo[@"correct surface"] != nil) {
        deTangle = [solution.solutionInfo[@"correct surface"] boolValue];
    }
    if (deTangle == YES) {
        fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: correct surface in case of of intersecting upper and lower surfaces.\n");
        if (solution.solutionInfo[@"minimum height"] != nil) {
            minHeight = [solution.solutionInfo[@"minimum height"] doubleValue];
        } else if (solution.solutionInfo[@"minimum height"] == nil || minHeight <= 0.0) {
            fatal("EMStructuredMeshMapper:solutionComputer", "Minimum height either set to negative/zero or not found.");
        }
        NSString *tangleMaskVarName;
        if (solution.solutionInfo[@"correct surface mask"] != nil) {
            tangleMaskVarName = solution.solutionInfo[@"correct surface mask"];
            FEMUtilities *utilities = [[FEMUtilities alloc] init];
            FEMVariable *tangleMaskVar = [utilities getVariableFrom:solution.mesh.variables model:model name:tangleMaskVarName onlySearch:NULL maskName:nil info:&found];
            if (found == YES) {
                computeTangleMask = YES;
                if (tangleMaskVar.dofs != 1) fatal("EMStructuredMeshMapper:solutionComputer", "correct surface mask variable should have only 1 dof.");
                variableArraysContainer *tangleMaskVarContainers = tangleMaskVar.getContainers;
                tangleMask = tangleMaskVarContainers->Values;
                for (int i=0; i<tangleMaskVarContainers->sizeValues; i++) {
                    tangleMask[i] = 1.0;
                }
                tangleMaskPerm = tangleMaskVarContainers->Perm;
                fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: output of correct surface mask to: %s.\n", [tangleMaskVarName UTF8String]);
            } else {
                fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: ignoring the variable %s given as correct surface mask because it was not found.\n", [tangleMaskVarName UTF8String]);
            }
        }
    }
    
    // Get either variable or constant values for top surface
    topNode = 0;
    topValue0 = 0.0;
    if (solution.solutionInfo[@"top surface level"] != nil) {
        topValue0 = [solution.solutionInfo[@"top surface level"] doubleValue];
        topNode = 1;
    } else {
        if (solution.solutionInfo[@"top surface variable name"] != nil) {
            varName = solution.solutionInfo[@"top surface variable name"];
            FEMUtilities *utilities = [[FEMUtilities alloc] init];
            varContainers = NULL;
            var = [utilities getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:nil info:&found];
            if (found == YES) {
                if (var.dofs != 1) {
                    fatal("FEMStructuredMeshMapper:solutionComputer", "Top surface should have only 1 dof.");
                } else {
                    varContainers = var.getContainers;
                    if (varContainers != NULL) {
                        _topField = varContainers->Values;
                        _topPerm = varContainers->Perm;
                    }
                    topNode = 2;
                }
            } else {
                fprintf(stderr, "FEMStructuredMeshMapper:solutionComputer: top surface variable is missing: %s.\n", [varName UTF8String]);
                fatal("FEMStructuredMeshMapper:solutionComputer");
            }
        }
    }
    
    if (topNode == 0) {
        if ([listUtilities listCheckPresentAnyBoundaryCondition:model name:@"top surface"] == YES) {
            topNode = 3;
            if (reinitialize == YES) {
                if (_field != NULL) free_dvector(_field, 0, _nsize-1);
            }
            if (_field == NULL) {
                _field = doublevec(0, _nsize-1);
                memset(_field, 0.0, _nsize*sizeof(double) );
            }
            for (int elem=0; elem<solution.mesh.numberOfBoundaryElements; elem++) {
                element = [core getBoundaryElement:solution atIndex:elem];
                bc = [core getBoundaryCondition:model forElement:element];
                if (bc != nil) {
                    n = element->Type.NumberOfNodes;
                    found = [listUtilities listGetReal:model inArray:bc forVariable:@"top surface" numberOfNodes:n indexes:element->NodeIndexes buffer:&vector minValue:NULL maxValue:NULL];
                    if (found == YES) {
                        for (int i=0; i<n; i++) {
                            _field[element->NodeIndexes[i]] = vector.vector[i];
                        }
                    }
                }
            }
        }
    }
    
    // Get either variable or constant values for bottom surface
    bottomNode = 0;
    bottomValue0 = 0.0;
    if (solution.solutionInfo[@"bottom surface level"] != nil) {
        bottomValue0 = [solution.solutionInfo[@"bottom surface level"] doubleValue];
        bottomNode = 1;
    } else {
        if (solution.solutionInfo[@"bottom surface variable name"] != nil) {
            varName = solution.solutionInfo[@"bottom surface variable name"];
            FEMUtilities *utilities = [[FEMUtilities alloc] init];
            varContainers = NULL;
            var = [utilities getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:nil info:&found];
            if (found == YES) {
                if (var.dofs != 1) {
                    fatal("FEMStructuredMeshMapper:solutionComputer", "Bottom surface should have only 1 dof.");
                } else {
                    varContainers = var.getContainers;
                    if (varContainers != NULL) {
                        _bottomField = varContainers->Values;
                        _bottomPerm = varContainers->Perm;
                    }
                    bottomNode = 2;
                }
            } else {
                fprintf(stderr, "FEMStructuredMeshMapper:solutionComputer: bottom surface variable is missing: %s.\n", [varName UTF8String]);
                fatal("FEMStructuredMeshMapper:solutionComputer");
            }
        }
    }
    
    if (bottomNode == 0) {
        if ([listUtilities listCheckPresentAnyBoundaryCondition:model name:@"bottom surface"] == YES) {
            bottomNode = 3;
            if (_field == NULL) {
                _field = doublevec(0, _nsize-1);
                memset(_field, 0.0, _nsize*sizeof(double) );
            }
            for (int elem=0; elem<solution.mesh.numberOfBoundaryElements; elem++) {
                element = [core getBoundaryElement:solution atIndex:elem];
                bc = [core getBoundaryCondition:model forElement:element];
                if (bc != nil) {
                    n = element->Type.NumberOfNodes;
                    found = [listUtilities listGetReal:model inArray:bc forVariable:@"bottom surface" numberOfNodes:n indexes:element->NodeIndexes buffer:&vector minValue:NULL maxValue:NULL];
                    if (found == YES) {
                        for (int i=0; i<n; i++) {
                            _field[element->NodeIndexes[i]] = vector.vector[i];
                        }
                    }
                }
            }
        }
    }
    
    // Get either variable or constant values for mid surface
    if (midLayerExists == YES) {
        if (_field == NULL) {
            _field = doublevec(0, _nsize-1);
            memset(_field, 0.0, _nsize*sizeof(double) );
        }
        for (int elem=0; elem<solution.mesh.numberOfBoundaryElements; elem++) {
            element = [core getBoundaryElement:solution atIndex:elem];
            bc = [core getBoundaryCondition:model forElement:element];
            if (bc != nil) {
                n = element->Type.NumberOfNodes;
                found = [listUtilities listGetReal:model inArray:bc forVariable:@"mid surface" numberOfNodes:n indexes:element->NodeIndexes buffer:&vector minValue:NULL maxValue:NULL];
                if (found == YES) {
                    for (int i=0; i<n; i++) {
                        _field[element->NodeIndexes[i]] = vector.vector[i];
                    }
                }
            }
        }
    }
    
    // Get the velocity variable component
    gotVeloVar = NO;
    if (solution.solutionInfo[@"mesh velocity variable"] != nil) {
        varName = solution.solutionInfo[@"mesh velocity variable"];
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        veloVar = [utilities getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:nil info:&found];
        if (found == YES) {
            if (veloVar.dofs == 1) {
                meshVeloContainers = veloVar.getContainers;
                gotVeloVar = YES;
            } else {
                fatal("FEMStructuredMeshMapper:solutionComputer", "The size of mesh velocity must be one.");
            }
        } else {
            fprintf(stderr, "FEMStructuredMeshMapper:solutionComputer: variable does not exist: %s.\n", [varName UTF8String]);
            fatal("FEMStructuredMeshMapper:solutionComputer");
        }
    }
    
    // Get the mesh update variable component
    gotUpdateVar = NO;
    if (solution.solutionInfo[@"mesh update variable"] != nil) {
        varName = solution.solutionInfo[@"mesh update variable"];
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        updateVar = [utilities getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:nil info:&found];
        if (found == YES) {
            if (updateVar.dofs == 1) {
                meshUpdateContainers = updateVar.getContainers;
                gotUpdateVar = YES;
            } else {
                fatal("FEMStructuredMeshMapper:solutionComputer", "The size of mesh update must be one.");
            }
        } else {
            fprintf(stderr, "FEMStructuredMeshMapper:solutionComputer: variable does not exist: %s.\n", [varName UTF8String]);
            fatal("FEMStructuredMeshMapper:solutionComputer");
        }
    }
    
    if (solution.solutionInfo[@"displacement mode"] != nil) {
        displacementMode = [solution.solutionInfo[@"displacement mode"] boolValue];
    } displacementMode = NO;
    
    // Get the new mapping using linear interpolation from bottom to top
    if (solution.solutionInfo[@"minimum mesh height"] != nil) {
        heps = [solution.solutionInfo[@"minimum mesh height"] doubleValue];
    } heps = DBL_EPSILON;
    
    int tangleCount = 0;
    BOOL tangled = NO;
    topVal = 0.0;
    bottomVal = 0.0;
    imid = 0;
    for (int i=0; i<_nsize; i++) {
        j = i;
        if (_maskExists) {
            j = _maskPerm[i];
            if (j < 0) continue;
        }
        itop = _topPointer[i];
        ibot = _bottomPointer[i];
        if (midLayerExists == YES) imid = _midPointer[i];
        
        // Use the previous coordinates for determining the weights
        x0Top = _origCoord[itop];
        x0Bot = _origCoord[ibot];
        x0Loc = _origCoord[i];
        
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
        
        if (displacementMode == YES) {
            tangled = (topVal + x0Top < bottomVal + x0Bot + minHeight) ? YES : NO;
        } else {
            tangled = (topVal < bottomVal + minHeight) ? YES :  NO;
        }
        
        if (tangled == YES) {
            tangleCount++;
            if (deTangle == NO) {
                fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: mode: %d.\n", displacementMode);
                fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: top val: %f %f %f.\n", topVal, x0Top, topVal+x0Top);
                fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: bottom val: %f %f %f.\n", bottomVal, x0Bot, bottomVal+x0Bot);
                
                Nodes_t *nodes = solution.mesh.getNodes;
                //fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: node %d, height %f, w %f.\n", i, _coord[i], wTop);
                fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: position: %f %f %f.\n", nodes->x[i], nodes->y[i], nodes->z[i]);
                fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: topVal: %f, botVal: %f, dVal %f.\n", topVal, bottomVal, topVal-bottomVal);

            } else {
                if (displacementMode == YES) {
                    topVal = bottomVal + x0Bot + x0Top + minHeight;
                } else {
                    topVal = bottomVal + minHeight;
                }
                
                if (computeTangleMask == YES) tangleMask[tangleMaskPerm[i]] = -1.0;
                if (/* DISABLES CODE */ (NO)) {
                    fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: corrected negative height: %f = %f - %f. New upper value: %f.\n", topVal-bottomVal, topVal, bottomVal, _field[itop]);
                }
            }
        } else {
            if (computeTangleMask == YES) {
                if (tangleMask[tangleMaskPerm[itop]] == -1.0) {
                    tangleMask[tangleMaskPerm[i]] = -1.0;
                } else {
                    tangleMask[tangleMaskPerm[i]] = 1.0;
                }
            }
        }
        
        // New coordinates location
        if (midLayerExists == YES) {
            // With middle layer in two parts
            midVal = _field[imid];
            x0Mid = _origCoord[imid];
            if ((x0Top - x0Mid) * (x0Loc - x0Mid) > 0.0) {
                wTop = (x0Loc - x0Mid) / (x0Top - x0Mid);
                xLoc = wTop * topVal + (1.0 - wTop) *  midVal;
            } else {
                wTop = (x0Loc - x0Bot) / (x0Mid - x0Bot);
                xLoc = wTop * midVal + (1.0 - wTop) * bottomVal;
            }
        } else {
            // Otherwise in one part
            wTop = (x0Loc - x0Bot) / (x0Top - x0Bot);
            xLoc = wTop * topVal + (1.0 - wTop) * bottomVal;
        }
        
        if (displacementMode == YES) {
            if (gotVeloVar == YES) {
             if (meshVeloContainers->Perm[i] >= 0) meshVeloContainers->Values[meshVeloContainers->Perm[i]] = xLoc / timeStep;
            }
            _coord[i] = _coord[i] + xLoc;
        } else {
            if (gotVeloVar == YES) {
                if (meshVeloContainers->Perm[i] >= 0) meshVeloContainers->Values[meshVeloContainers->Perm[i]] = (xLoc - _origCoord[i]) / timeStep;
            }
            _coord[i] = xLoc;
        }
        if (gotUpdateVar == YES) meshUpdateContainers->Values[meshUpdateContainers->Perm[i]] = _coord[i] - _origCoord[i];
    }
    
    if (gotVeloVar == YES && _visited == NO) {
        if (solution.solutionInfo[@"mesh velocity first zero"] != nil) {
            if ([solution.solutionInfo[@"mesh velocity first zero"] boolValue] == YES) {
                if (gotVeloVar == YES) memset(meshVeloContainers->Values, 0.0, meshVeloContainers->sizeValues*sizeof(double) );
            }
        }
    }
    
    if (tangleCount > 0) {
        fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: there seems to be %d (out of %d) tangled nodes.\n", tangleCount, _nsize);
    }
    
    at1 = cputime();
    fprintf(stdout, "FEMStructuredMeshMapper:solutionComputer: active coordinate mapping time: %f.\n", at1-at0);
    
    _visited = YES;
    
    if (vector.vector != NULL) {
        free_dvector(vector.vector, 0, vector.m-1);
        vector.vector = NULL;
    }
}

-(void)deallocation:(FEMSolution * __nonnull)solution {
    if (_topPointer != NULL) {
        free_ivector(_topPointer, 0, _nsize-1);
    }
    if (_bottomPointer != NULL) {
        free_ivector(_bottomPointer, 0, _nsize-1);
    }
    if (_upPointer != NULL) {
        free_ivector(_upPointer, 0, _nsize-1);
    }
    if (_downPointer != NULL) {
        free_ivector(_downPointer, 0, _nsize-1);
    }
    if (_midPointer != NULL) {
        free_ivector(_midPointer, 0, _nsize-1);
    }
    if (_nodeLayer != NULL) {
        free_ivector(_nodeLayer, 0, _nsize-1);
    }
    if (_origCoord != NULL) {
        free_dvector(_origCoord, 0, _nsize-1);
    }
    if (_field != NULL) {
        free_dvector(_field, 0, _nsize-1);
    }
}

@end
