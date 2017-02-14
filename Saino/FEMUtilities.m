//===----------------------------------------------------------------------===//
//  FEMUtilities.m
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
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

#import "FEMUtilities.h"
#import "FEMCore.h"
#import "FEMPElementMaps.h"
#import "FEMProjector.h"
#import "FEMElementUtils.h"
#import "FEMMeshUtils.h"
#import "FEMMatrixCRS.h"
#import "FEMMatrixBand.h"
#import "FEMListUtilities.h"
#import "FEMInterpolation.h"
#import "FEMElementDescription.h"
#import "SainoSolutionsComputer.h"
#import "FEMInterpolateMeshToMesh.h"
#import "Utils.h"

@implementation FEMUtilities {
    char _backSlash;
}

#pragma mark Public methods

@synthesize ext = _ext;
@synthesize appSupportSubpath = _appSupportSubpath;

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _ext = @"bundle";
        _appSupportSubpath = @"Application Support/Saino/PlugIns";
        _backSlash = (char)92;
    }
    return self;
}

-(void)zeroTheNumberOfRows:(int)n inMatrix:(FEMMatrix * __nonnull)a {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    if (a.format == MATRIX_CRS) {
        crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix zeroRowInMatrix:a numberOfRows:n];
        
    } else if (a.format == MATRIX_LIST) {
        
        // TODO: implement the zeroRow method for list matrix.
        
    } else if (a.format == MATRIX_BAND || a.format == MATRIX_SBAND) {
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix zeroRowInMatrix:a numberOfRows:n];
    }
}

-(void)setMatrixElement:(FEMMatrix * __nonnull)a atIndex:(int)i andIndex:(int)j value:(double)value {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    if (a.format == MATRIX_CRS) {
        crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix setElementInMatrix:a row:i col:j value:value];
        
    } else if (a.format == MATRIX_LIST) {
        // TODO: implement the setMatrixElement method for list matrix.
        
    } else if (a.format == MATRIX_BAND || a.format == MATRIX_SBAND) {
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix setElementInMatrix:a row:i col:j value:value];
    }
}

-(int)initialPermutationInMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model solution:(FEMSolution *__nonnull)solution equation:(NSString * __nonnull)str permutation:(int * __nonnull)perm DGSolution:(BOOL * __nullable)dg globalBubbles:(BOOL * __nullable)gb {
    
    int i, j, t, n, e, edofs, fdofs, bdofs, ndofs;
    int indexes[128];
    int *defDofs, *edgeDofs = NULL, *faceDofs = NULL;
    BOOL foundDG, DG, GB, radiation, any;
    Element_t *element, *elements, *edges, *faces;
    FEMPElementMaps *elementMaps;
    FEMListUtilities *listUtilities;
    solutionArraysContainer *solContainers = NULL;
    
    int k = 0;
    edofs = mesh.maxEdgeDofs;
    fdofs = mesh.maxFaceDofs;
    bdofs = mesh.maxBdofs;
    
    elementMaps = [[FEMPElementMaps alloc] init];
    listUtilities = [FEMListUtilities sharedListUtilities];
    
    elements = mesh.getElements;
    edges = mesh.getEdges;
    faces = mesh.getFaces;
    
    GB = NO;
    if (gb != NULL) GB = *gb;
    
    DG = NO;
    if (dg != NULL) DG = *dg;
    
    foundDG = NO;
    if (DG == YES) {
        for (t=0; t<mesh.numberOfEdges; t++) {
            n = 0;
            element = edges[t].BoundaryInfo->Left;
            if (element != NULL) {
                if ([listUtilities checkElementEquation:model forElement:element andEquation:str] == YES) {
                    foundDG = (foundDG == YES || element->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<element->DGDOFs; j++) {
                        indexes[n] = element->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            element = edges[t].BoundaryInfo->Right;
            if (element != NULL) {
                if ([listUtilities checkElementEquation:model forElement:element andEquation:str] == YES) {
                    foundDG = (foundDG == YES || element->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<element->DGDOFs; j++) {
                        indexes[n] = element->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            for (i=0; i<n; i++) {
                j = indexes[i];
                if (perm[j] < 0 ) {
                    perm[j] = k;
                    k++;
                }
            }
        }
        
        for (t=0; t<mesh.numberOfFaces; t++) {
            n = 0;
            element = faces[t].BoundaryInfo->Left;
            if (element != NULL) {
                if ([listUtilities checkElementEquation:model forElement:element andEquation:str] == YES) {
                    foundDG = (foundDG == YES || element->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<element->DGDOFs; j++) {
                        indexes[n] = element->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            element = faces[t].BoundaryInfo->Right;
            if (element != NULL) {
                if ([listUtilities checkElementEquation:model forElement:element andEquation:str] == YES) {
                    foundDG = (foundDG == YES || element->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<element->DGDOFs; j++) {
                        indexes[n] = element->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            for (i=0; i<n; i++) {
                j = indexes[i];
                if (perm[j] < 0) {
                    perm[j] = k;
                    k++;
                }
            }
        }
        
        if (foundDG == YES) {
            return k; // Discontinuous galerkin!!!
        }
    }
    
    solContainers = solution.getContainers;
    defDofs = intvec(0, solContainers->size2DefDofs-1);
    
    any = NO;
    for (i=0; i<solContainers->size1DefDofs; i++) {
        if (solContainers->defDofs[i][5] >= 0) {
            any = YES;
            break;
        }
    }
    if (any == YES) {
        if (mesh.numberOfEdges > 0) {
            edgeDofs = intvec(0, mesh.numberOfEdges-1);
            memset( edgeDofs, 0, mesh.numberOfEdges*sizeof(int) );
        }
        
        if (mesh.numberOfFaces > 0) {
            faceDofs = intvec(0, mesh.numberOfFaces-1);
            memset( faceDofs, 0, mesh.numberOfFaces*sizeof(int) );
        }
        
        n = mesh.numberOfBulkElements + mesh.numberOfBoundaryElements;
        t = 0;
        while (t < n) {
            while (t < n) {
                if ([listUtilities checkElementEquation:model forElement:&elements[t] andEquation:str] == YES) break;
                t++;
            }
            if (t >= n) break;
            
            for (i=0; i<solContainers->size2DefDofs; i++) {
                defDofs[i] = solContainers->defDofs[elements[t].BodyID-1][i];
            }
            if (elements[t].EdgeIndexes != NULL) {
                for (i=0; i<elements[t].Type.NumberOfEdges; i++) {
                    j = elements[t].EdgeIndexes[i];
                    edgeDofs[j] = max(edgeDofs[j], [elementMaps getEdgeDofsForElement:&elements[t] polyDegree:defDofs[5]]);
                }
            }
            
            if (elements[t].FaceIndexes != NULL) {
                for (i=0; i<elements[t].Type.NumberOfFaces; i++) {
                    j = elements[t].FaceIndexes[i];
                    faceDofs[j] = max(faceDofs[j], [elementMaps getFaceDofsForElement:&elements[t] polyDegree:defDofs[5] faceNumber:i]);
                }
            }
            t++;
        }
    }
    
    n = mesh.numberOfBulkElements + mesh.numberOfBoundaryElements;
    t = 0;
    while (t < n) {
        while (t < n) {
            if ([listUtilities checkElementEquation:model forElement:&elements[t] andEquation:str] == YES) break;
            t++;
        }
        if (t >= n) break;
        
        for (i=0; i<solContainers->size2DefDofs; i++) {
            defDofs[i] = solContainers->defDofs[elements[t].BodyID-1][i];
        }
        ndofs = elements[t].NDOFs;
        if (defDofs[0] >= 0) ndofs = defDofs[0]*elements[t].Type.NumberOfNodes;
        for (i=0; i<ndofs; i++) {
            j = elements[t].NodeIndexes[i];
            if (perm[j] < 0) {
                perm[j] = k;
                k++;
            }
        }
        
        if (elements[t].EdgeIndexes != NULL) {
            for (i=0; i<elements[t].Type.NumberOfEdges; i++) {
                ndofs = 0;
                if (defDofs[1] >= 0) {
                    ndofs = defDofs[1];
                } else if (defDofs[5] >= 0) {
                    ndofs = edgeDofs[elements[t].EdgeIndexes[i]];
                    ndofs = max(edges[elements[t].EdgeIndexes[i]].BDOFs, ndofs);
                }
                
                for (e=0; e<ndofs; e++) {
                    j = mesh.numberOfNodes + edofs*elements[t].EdgeIndexes[i] + e;
                    if (perm[j] < 0) {
                        perm[j] = k;
                        k++;
                    }
                }
            }
        }
        
        if (elements[t].FaceIndexes != NULL) {
            for (i=0; i<elements[t].Type.NumberOfFaces; i++) {
                ndofs = 0;
                if (defDofs[2] >= 0) {
                    ndofs = defDofs[2];
                } else if (defDofs[5] >= 0) {
                    ndofs = faceDofs[elements[t].FaceIndexes[i]];
                    ndofs = max(faces[elements[t].FaceIndexes[i]].BDOFs, ndofs);
                }
                
                for (e=0; e<ndofs; e++) {
                    j = mesh.numberOfNodes + edofs*mesh.numberOfEdges + fdofs*elements[t].FaceIndexes[i] + e;
                    if (perm[j] < 0) {
                        perm[j] = k;
                        k++;
                    }
                }
            }
        }
        
        if (GB == YES || elements[t].BubbleIndexes != NULL) {
            ndofs = 0;
            if (defDofs[4] >= 0) {
                ndofs = defDofs[4];
            } else if (defDofs[5] >= 0) {
                ndofs = [elementMaps getBubbleDofsForElement:&elements[t] polyDegree:defDofs[5]];
                if (defDofs[5] == 0) ndofs = max(elements[t].BDOFs, ndofs);
            }
            
            for (i=0; i<ndofs; i++) {
                j = mesh.numberOfNodes + edofs*mesh.numberOfEdges + fdofs*mesh.numberOfFaces + elements[t].BubbleIndexes[i];
                if (perm[j] < 0) {
                    perm[j] = k;
                    k++;
                }
            }
        }
        t++;
    }
    
    radiation = [solution.solutionInfo[@"radiation solver"] boolValue];
    if (radiation == YES || [str isEqualToString:@"heat equation"]) {
        t = mesh.numberOfBulkElements;
        n = mesh.numberOfBulkElements + mesh.numberOfBoundaryElements;
        while (t < n) {
            if (elements[t].BoundaryInfo->GebhardtFactors != NULL) {
                for (i=0; i<elements[t].Type.NumberOfNodes; i++) {
                    j = elements[t].NodeIndexes[i];
                    if (perm[j] < 0) {
                        perm[j] = k;
                        k++;
                    }
                }
            }
            t++;
        }
    }
    
    t = mesh.numberOfBulkElements;
    n = mesh.numberOfBulkElements + mesh.numberOfBoundaryElements;
    while (t < n) {
        if (elements[t].Type.ElementCode == 102) {
            for (i=0; i<elements[t].Type.NumberOfNodes; i++) {
                j = elements[t].NodeIndexes[i];
                if (perm[j] < 0) {
                    perm[j] = k;
                    k++;
                }
            }
        }
        t++;
    }
    
    if (edgeDofs != NULL) {
        free_ivector(edgeDofs, 0, mesh.numberOfEdges-1);
    }
    if (faceDofs != NULL) {
        free_ivector(faceDofs, 0,  mesh.numberOfFaces-1);
    }
    free_ivector(defDofs, 0, solContainers->size2DefDofs-1);
    
    return k;
}

/********************************************************************************************
 
    Adds a new variable to the list of variables

*********************************************************************************************/
-(void)addVariableTo:(NSMutableArray * __nonnull)anArray mesh:(FEMMesh * __nullable)mesh solution:(FEMSolution * __nullable)solution name:(NSString * __nonnull)name dofs:(int )dofs container:(variableArraysContainer * __nonnull)aContainer component:(BOOL)component ifOutput:(BOOL * __nullable)output ifSecondary:(BOOL * __nullable)secondary type:(int * __nullable)aType {
    
    FEMVariable *newVariable;
    variableArraysContainer *varContainers = NULL;
    
    if (aContainer == NULL) {
        fprintf(stderr, "FEMUtilities:addVariableTo: container argument is not allocated.\n");
        fatal("FEMUtilities:addVariableTo");
    }
    
    for (FEMVariable *variable in anArray) {
        if ([variable.name isEqualToString:name] == YES) {
            return; // Variable alreadt exists, don't do anything;
        }
    }
    
    newVariable = [[FEMVariable alloc] init];
    varContainers = newVariable.getContainers;
    
    newVariable.name = name;
    
    newVariable.dofs = dofs;
    if (aContainer->Perm != NULL) {
        varContainers->Perm = aContainer->Perm;
        varContainers->sizePerm = aContainer->sizePerm;
    }
    
    newVariable.norm = 0.0;
    newVariable.prevNorm = 0.0;
    
    newVariable.componentVariable = NO;
    newVariable.secondary = NO;
    newVariable.componentSecondaryVariable = NO;
    
    if (component == YES) {
        // If component is true, then we are adding a variable which is itself
        // a component of another variable with dof > 1. In this case, we set
        // an array of pointers where each pointer points to the appropriate
        // location in the variable array to which the component variable belongs.
        
        if (aContainer->ComponentValues != NULL) {
            varContainers->ComponentValues = aContainer->ComponentValues;
            varContainers->sizeComponentValues = aContainer->sizeComponentValues;
            newVariable.componentVariable = YES;
        } else if (*secondary == YES && aContainer->ComponentSecondaryToValues != NULL) {
            varContainers->ComponentSecondaryToValues = aContainer->ComponentSecondaryToValues;
            varContainers->sizeComponentSecondaryToValues = aContainer->sizeComponentSecondaryToValues;
            newVariable.componentSecondaryVariable = YES;
        }
    } else {
        if (aContainer->Values != NULL) {
            varContainers->Values = aContainer->Values;
            varContainers->sizeValues = aContainer->sizeValues;
        } else if (secondary != NULL) {
            if (*secondary == YES && aContainer->SecondaryToValues != NULL) {
                varContainers->SecondaryToValues = aContainer->SecondaryToValues;
                varContainers->sizeSecondaryToValues = aContainer->sizeSecondaryToValues;
                newVariable.secondary = YES;
            }
        }
    }
    
    newVariable.nonLinChange = 0.0;
    newVariable.steadyChange = 0.0;
    newVariable.nonLinIter = 0;
    
    if (solution != nil) newVariable.solution = solution;
    if (mesh != nil) newVariable.primaryMesh = mesh;
    
    newVariable.valid = YES;
    newVariable.output = YES;
    newVariable.valuesChanged = YES;
    
    // Converged information undefined = -1 not 0, yes = 1
    newVariable.nonLinConverged = -1;
    newVariable.steadyConverged = -1;
    
    if (secondary != NULL) {
        fprintf(stdout, "FEMUtilities:addVariableTo: secondary: %s.", [name UTF8String]);
        newVariable.secondary = *secondary;
    }
    if (aType != NULL) {
        newVariable.type = *aType;
    }
    if (output != NULL) {
        newVariable.output = *output;
    }
    
    // Finally add the new variable to the array of variables
    [anArray addObject:newVariable];
}

-(void)addVectorTo:(NSMutableArray * __nonnull)anArray mesh:(FEMMesh * __nonnull)mesh solution:(FEMSolution * __nonnull)solution name:(NSString * __nonnull)name dofs:(int * __nullable)dofs container:(variableArraysContainer * __nonnull)aContainer ifOutput:(BOOL * __nullable)output ifSecondary:(BOOL * __nullable)secondary global:(BOOL * __nullable)global initValue:(double * __nullable)initValue {
    
    int i, j, k, ndofs, nsize;
    NSString *tmpName;
    variableArraysContainer *varContainers = NULL;
    
    if (dofs != NULL) {
        ndofs = *dofs;
    } else {
        ndofs = 1;
    }
    
    varContainers = allocateVariableContainer();
    
    if (aContainer->Values != NULL) {
        varContainers->Values = aContainer->Values;
        varContainers->sizeValues = aContainer->sizeValues;
    } else if (*secondary == YES && aContainer->SecondaryToValues != NULL) {
        // Some variable may have their values arrays pointing to a given column of
        // another variable PrevValues. This happens when the options "calculate velocity" or "calculate acceleration"
        // are given to a solution and when the solution time order >= 2. In this case, the solution variable creates
        // the variables:
        //   [solution.variable.name appendString:@"velocity"]
        //   [solution.variable.name appendString:@"acceleration"]
        // which values point to PrevValues columns.
        varContainers->SecondaryToValues = aContainer->SecondaryToValues;
        varContainers->sizeSecondaryToValues = aContainer->sizeSecondaryToValues;
    } else {
        if (aContainer->Perm != NULL) {
            nsize = max_array(aContainer->Perm, aContainer->sizePerm);
        } else if (global != NULL) {
            if (*global == YES) {
                nsize = 1;
            } else {
                nsize = mesh.numberOfNodes;
            }
        } else {
            nsize = mesh.numberOfNodes;
        }
        varContainers->Values = doublevec(0, (ndofs*nsize)-1);
        varContainers->sizeValues = (ndofs*nsize);
        memset( varContainers->Values, 0.0, (ndofs*nsize)*sizeof(double) );
    }
    
    if (initValue != NULL) {
        if (*secondary == YES && aContainer->SecondaryToValues != NULL) {
            for (i=0; i<varContainers->sizeSecondaryToValues; i++) {
                *(varContainers->SecondaryToValues[i]) = *initValue;
            }
        } else {
            for (i=0; i<varContainers->sizeValues; i++) {
                varContainers->Values[i] = *initValue;
            }            
        }
    }
    
    if (ndofs > 1) {
        for (i=1; i<=ndofs; i++) {
            tmpName = [self appendNameFromString:name component:&i];
            if (*secondary == YES && aContainer->SecondaryToValues != NULL) {
                varContainers->ComponentSecondaryToValues = malloc ( (aContainer->sizeSecondaryToValues/ndofs) * sizeof ( double * ));
                k = 0;
                 for (j=(i-1); j<aContainer->sizeSecondaryToValues; j+=ndofs) {
                     varContainers->ComponentSecondaryToValues[k] = aContainer->SecondaryToValues[j];
                     k++;
                 }
                [self addVariableTo:anArray mesh:mesh solution:solution name:tmpName dofs:1 container:varContainers component:YES ifOutput:output ifSecondary:secondary type:NULL];
            } else {
                varContainers->ComponentValues = malloc ( (varContainers->sizeValues/ndofs) * sizeof ( double * ));
                k = 0;
                for (j=(i-1); j<varContainers->sizeValues; j+=ndofs) {
                    varContainers->ComponentValues[k] = &varContainers->Values[j];
                    k++;
                }
                varContainers->sizeComponentValues = (varContainers->sizeValues/ndofs);
                [self addVariableTo:anArray mesh:mesh solution:solution name:tmpName dofs:1 container:varContainers component:YES ifOutput:output ifSecondary:secondary type:NULL];
            }
        }
    }
    
    [self addVariableTo:anArray mesh:mesh solution:solution name:name dofs:ndofs container:varContainers component:NO ifOutput:output ifSecondary:secondary type:NULL];
    free(varContainers);
}

-(FEMVariable * __nullable)getVariableFrom:(NSMutableArray * __nonnull)anArray model:(FEMModel * __nonnull)model name:(NSString * __nonnull)name onlySearch:(BOOL * __nullable)only maskName:(NSString * __nullable)maskName info:(BOOL * __nonnull)found {
    
    int i, j, k, l, n, dofs, aNumber;
    FEMVariable *var = nil, *pVar = nil, *tmp = nil;
    FEMListUtilities *listUtilities;
    FEMSolution *solution;
    FEMMesh *mesh, *currentMesh;
    FEMProjector *projector;
    variableArraysContainer *varContainers = NULL, *pvarContainers = NULL, *bufferContainers = NULL;
    NSString *searchedName, *tmpname;
    BOOL onlyThis, stat, globalBubbles, output, canonicalize=YES;
    
    NSRange ind = [name rangeOfString:@"["];
    if (ind.location != NSNotFound) canonicalize = NO;
    
    *found = NO;
    for (FEMVariable *variable in anArray) {
        if (canonicalize == YES) {
            searchedName = [variable canonicalizeName];
        } else {
            searchedName = variable.name;
        }
        if ([searchedName isEqualToString:name] == YES) {
            if (variable.valid == YES) {
                *found = YES;
                var = variable;
            }
        }
    }
    
    if (only != NULL) {
        if (*only == YES) {
            return var;
        }
    }
    
    for (FEMMesh *mesh in model.meshes) {
        if ([anArray isEqualToArray:mesh.variables] == NO) {
            onlyThis = YES;
            pVar = [self getVariableFrom:mesh.variables model:model name:name onlySearch:&onlyThis maskName:NULL info:&stat];
            if (pVar != nil) {
                if ([mesh isEqual:pVar.primaryMesh]) {
                    break;
                }
            }
        }
    }
    if (pVar == nil) return var;
    
    if (*found == NO ) {
        listUtilities = [FEMListUtilities sharedListUtilities];
        solution = (FEMSolution *)pVar.solution;
        if ( solution.solutionInfo[@"bubbles in global system"] != nil ) {
            globalBubbles = [solution.solutionInfo[@"bubbles in global system"] boolValue];
        } else {
            globalBubbles = YES;
        }
        
        mesh = (FEMMesh *)model.mesh;
        dofs = mesh.numberOfNodes * pVar.dofs;
        if (globalBubbles == YES) dofs = dofs + mesh.maxBdofs * mesh.numberOfBulkElements * pVar.dofs;
        
        var = [[FEMVariable alloc] init];
        varContainers = var.getContainers;
        varContainers->Values = doublevec(0, dofs-1);
        varContainers->sizeValues = dofs;
        memset( varContainers->Values, 0.0, dofs*sizeof(double) );
        varContainers->Perm = NULL;
        pvarContainers = pVar.getContainers;
        if (pvarContainers->Perm != NULL) {
            varContainers->Perm = intvec(0, (dofs/pVar.dofs)-1);
            varContainers->sizePerm = dofs/pVar.dofs;
            memset( varContainers->Perm, -1, (dofs/pVar.dofs)*sizeof(int) );
            
            n = [self initialPermutationInMesh:mesh model:model solution:solution equation:solution.solutionInfo[@"equation"] permutation:varContainers->Perm DGSolution:NULL globalBubbles:&globalBubbles];
            
            if (n == 0) n = mesh.numberOfNodes;
            
            if (n == mesh.numberOfNodes) {
                for (i=0; i<n; i++) {
                    varContainers->Perm[i] = i;
                }
            }
        }
        output = pVar.output;
        [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:pVar.name dofs:pVar.dofs container:varContainers component:NO ifOutput:&output ifSecondary:NULL type:NULL];
        
        varContainers->Perm = NULL;
        varContainers->Values = NULL;
        [var deallocation];
        
        onlyThis = YES;
        var = [self getVariableFrom:anArray model:model name:name onlySearch:&onlyThis maskName:NULL info:&stat];
        varContainers = var.getContainers;
        varContainers->PrevValues = NULL;
        if (pvarContainers->PrevValues != NULL) {
            varContainers->PrevValues = doublematrix(0, dofs-1, 0, pvarContainers->size2PrevValues-1);
            varContainers->size1PrevValues = dofs;
            varContainers->size2PrevValues = pvarContainers->size2PrevValues;
        }
        
        if ([pVar.name isEqualToString:@"flow solution"] == YES) {
            bufferContainers = allocateVariableContainer();
            bufferContainers->ComponentValues = malloc ( (dofs/pVar.dofs) * sizeof ( double * ));
            k = 0;
            for (i=0; i<varContainers->sizeValues; i+=pVar.dofs) {
                bufferContainers->ComponentValues[k] = &varContainers->Values[i];
                k++;
            }
            bufferContainers->sizeComponentValues = dofs/pVar.dofs;
            bufferContainers->Perm = varContainers->Perm;
            bufferContainers->sizePerm =varContainers->sizePerm;
            output = pVar.output;
            [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:@"velocity 1" dofs:1 container:bufferContainers component:YES ifOutput:&output ifSecondary:NULL type:NULL];
            
            free(bufferContainers);
            bufferContainers = NULL;
            onlyThis = YES;
            tmp = [self getVariableFrom:anArray model:model name:@"velocity 1" onlySearch:&onlyThis maskName:NULL info:&stat];
            bufferContainers = tmp.getContainers;
            if (varContainers->PrevValues != NULL) {
                // If the variable dof > 1, then we are dealing with a component variable. In this case,
                // we set a 2D array of pointers where each pointer points to the appropriate
                // location in the variable PrevValues array to which the component variable belongs.
                bufferContainers->ComponentPrevValues = malloc ( (dofs/pVar.dofs) * sizeof ( double ** ));
                for (i=0; i<(dofs/pVar.dofs); i++) {
                    bufferContainers->ComponentPrevValues[i] = malloc ( varContainers->size2PrevValues * sizeof ( double * ));
                }
                k = 0;
                for (i=0; i<varContainers->size1PrevValues; i+=pVar.dofs) {
                    for (j=0; j<varContainers->size2PrevValues; j++) {
                        bufferContainers->ComponentPrevValues[k][j] = &varContainers->PrevValues[i][j];
                    }
                    k++;
                }
                bufferContainers->size1ComponentPrevValues = dofs/pVar.dofs;
                bufferContainers->size2ComponentPrevValues = varContainers->size2PrevValues;
            }
            
            bufferContainers = NULL;
            bufferContainers = allocateVariableContainer();
            bufferContainers->ComponentValues = malloc ( (dofs/pVar.dofs) * sizeof ( double * ));
            k = 0;
            for (i=1; i<varContainers->sizeValues; i+=pVar.dofs) {
                bufferContainers->ComponentValues[k] = &varContainers->Values[i];
                k++;
            }
            bufferContainers->sizeComponentValues = dofs/pVar.dofs;
            bufferContainers->Perm = varContainers->Perm;
            bufferContainers->sizePerm = varContainers->sizePerm;
            [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:@"velocity 2" dofs:1 container:bufferContainers component:YES ifOutput:&output ifSecondary:NULL type:NULL];
            
            free(bufferContainers);
            bufferContainers = NULL;
            onlyThis = YES;
            tmp = [self getVariableFrom:anArray model:model name:@"velocity 2" onlySearch:&onlyThis maskName:NULL info:&stat];
            bufferContainers = tmp.getContainers;
            if (varContainers->PrevValues != NULL) {
                bufferContainers->ComponentPrevValues = malloc ( (dofs/pVar.dofs) * sizeof ( double ** ));
                for (i=0; i<(dofs/pVar.dofs); i++) {
                    bufferContainers->ComponentPrevValues[i] = malloc ( varContainers->size2PrevValues * sizeof ( double * ));
                }
                k = 0;
                for (i=1; i<varContainers->size1PrevValues; i+=pVar.dofs) {
                    for (j=0; j<varContainers->size2PrevValues; j++) {
                        bufferContainers->ComponentPrevValues[k][j] = &varContainers->PrevValues[i][j];
                    }
                    k++;
                }
                bufferContainers->size1ComponentPrevValues = dofs/pVar.dofs;
                bufferContainers->size2ComponentPrevValues = varContainers->size2PrevValues;
            }
            
            bufferContainers = NULL;
            bufferContainers = allocateVariableContainer();
            bufferContainers->ComponentValues = malloc ( (dofs/pVar.dofs) * sizeof ( double * ));
            k = 0;
            for (i=2; i<varContainers->sizeValues; i+=pVar.dofs) {
                bufferContainers->ComponentValues[k] = &varContainers->Values[i];
                k++;
            }
            bufferContainers->sizeComponentValues = dofs/pVar.dofs;
            bufferContainers->Perm = varContainers->Perm;
            bufferContainers->sizePerm = varContainers->sizePerm;
            if (pVar.dofs == 3) {
                [self addVariableTo:anArray mesh:pVar.primaryMesh solution:solution name:@"pressure" dofs:1 container:bufferContainers component:YES ifOutput:&output ifSecondary:NULL type:NULL];
            } else {
                [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:@"velocity 3" dofs:1 container:bufferContainers component:YES ifOutput:&output ifSecondary:NULL type:NULL];
                
                free(bufferContainers);
                bufferContainers = NULL;
                onlyThis = YES;
                tmp = [self getVariableFrom:anArray model:model name:@"velocity 3" onlySearch:&onlyThis maskName:NULL info:&stat];
                bufferContainers = tmp.getContainers;
                if (varContainers->PrevValues != NULL) {
                    bufferContainers->ComponentPrevValues = malloc ( (dofs/pVar.dofs) * sizeof ( double ** ));
                    for (i=0; i<(dofs/pVar.dofs); i++) {
                        bufferContainers->ComponentPrevValues[i] = malloc ( varContainers->size2PrevValues * sizeof ( double * ));
                    }
                    k = 0;
                    for (i=2; i<varContainers->size1PrevValues; i+=pVar.dofs) {
                        for (j=0; j<varContainers->size2PrevValues; j++) {
                            bufferContainers->ComponentPrevValues[k][j] = &varContainers->PrevValues[i][j];
                        }
                        k++;
                    }
                    bufferContainers->size1ComponentPrevValues = dofs/pVar.dofs;
                    bufferContainers->size2ComponentPrevValues = varContainers->size2PrevValues;
                }
                
                bufferContainers = NULL;
                bufferContainers = allocateVariableContainer();
                bufferContainers->ComponentValues = malloc ( (dofs/pVar.dofs) * sizeof ( double * ));
                k = 0;
                for (i=3; i<varContainers->sizeValues; i+=pVar.dofs) {
                    bufferContainers->ComponentValues[k] = &varContainers->Values[i];
                    k++;
                }
                bufferContainers->sizeComponentValues = dofs/pVar.dofs;
                bufferContainers->Perm = varContainers->Perm;
                bufferContainers->sizePerm = varContainers->sizePerm;
                [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:@"pressure" dofs:1 container:bufferContainers component:YES ifOutput:&output ifSecondary:NULL type:NULL];
            }
            
            free(bufferContainers);
            bufferContainers = NULL;
            onlyThis = YES;
            tmp = [self getVariableFrom:anArray model:model name:@"pressure" onlySearch:&onlyThis maskName:NULL info:&stat];
            bufferContainers = tmp.getContainers;
            if (varContainers->PrevValues != NULL) {
                bufferContainers->ComponentPrevValues = malloc ( (dofs/pVar.dofs) * sizeof ( double ** ));
                for (i=0; i<(dofs/pVar.dofs); i++) {
                    bufferContainers->ComponentPrevValues[i] = malloc ( varContainers->size2PrevValues * sizeof ( double * ));
                }
                k = 0;
                for (i=pVar.dofs-1; i<varContainers->size1PrevValues; i+=pVar.dofs) {
                    for (j=0; j<varContainers->size2PrevValues; j++) {
                        bufferContainers->ComponentPrevValues[k][j] = &varContainers->PrevValues[i][j];
                    }
                    k++;
                }
                bufferContainers->size1ComponentPrevValues = dofs/pVar.dofs;
                bufferContainers->size2ComponentPrevValues = varContainers->size2PrevValues;
            }
        } else {
            output = pVar.output;
            if (pVar.dofs > 1) {
                for (i=0; i<pVar.dofs; i++) {
                    bufferContainers = NULL;
                    bufferContainers = allocateVariableContainer();
                    bufferContainers->ComponentValues = malloc ( (dofs/pVar.dofs) * sizeof ( double * ));
                    k = 0;
                    for (j=i; j<varContainers->sizeValues; j+=pVar.dofs) {
                        bufferContainers->ComponentValues[k] = &varContainers->Values[j];
                        k++;
                    }
                    bufferContainers->sizeComponentValues = dofs/pVar.dofs;
                    bufferContainers->Perm = varContainers->Perm;
                    bufferContainers->sizePerm = varContainers->sizePerm;
                    aNumber = i+1;
                    tmpname = [self appendNameFromString:pVar.name component:&aNumber];
                    [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:tmpname dofs:1 container:bufferContainers component:YES ifOutput:&output ifSecondary:NULL type:NULL];
                    
                    free(bufferContainers);
                    bufferContainers = NULL;
                    onlyThis = YES;
                    tmp = [self getVariableFrom:anArray model:model name:tmpname onlySearch:&onlyThis maskName:NULL info:&stat];
                    bufferContainers = tmp.getContainers;
                    if (varContainers->PrevValues != NULL) {
                        bufferContainers->ComponentPrevValues = malloc ( (dofs/pVar.dofs) * sizeof ( double ** ));
                        for (j=0; j<(dofs/pVar.dofs); j++) {
                            bufferContainers->ComponentPrevValues[j] = malloc ( varContainers->size2PrevValues * sizeof ( double * ));
                        }
                        k = 0;
                        for (j=i; j<varContainers->size1PrevValues; j+=pVar.dofs) {
                            for (l=0; l<varContainers->size2PrevValues; l++) {
                                bufferContainers->ComponentPrevValues[k][l] = &varContainers->PrevValues[j][l];
                            }
                            k++;
                        }
                        bufferContainers->size1ComponentPrevValues = dofs/pVar.dofs;
                        bufferContainers->size2ComponentPrevValues = varContainers->size2PrevValues;
                    }
                }
            }
        }
        onlyThis = YES;
        var = [self getVariableFrom:anArray model:model name:name onlySearch:&onlyThis maskName:NULL info:&stat];
    }
    
    // Build a temporary variable array of variables to be interpolated
    NSMutableArray *tmpArryay;
    mesh = nil;
    mesh = (FEMMesh *)pVar.primaryMesh;
    if ([pVar.name isEqualToString:@"flow solution"]) {
        tmpArryay = [NSMutableArray arrayWithObjects:[self getVariableFrom:mesh.variables model:model name:@"velocity 1" onlySearch:NULL maskName:NULL info:&stat],
                     [self getVariableFrom:mesh.variables model:model name:@"velocity 2" onlySearch:NULL maskName:NULL info:&stat],
                     [self getVariableFrom:mesh.variables model:model name:@"velocity 3" onlySearch:NULL maskName:NULL info:&stat],
                     [self getVariableFrom:mesh.variables model:model name:@"pressure" onlySearch:NULL maskName:NULL info:&stat], nil];
    } else if (pVar.dofs > 1){
        for (i=0; i<pVar.dofs; i++) {
             aNumber = i+1;
            tmpname = [self appendNameFromString:pVar.name component:&aNumber];
            [tmpArryay addObject:[self getVariableFrom:mesh.variables model:model name:tmpname onlySearch:NULL maskName:NULL info:&stat]];
        }
    }
    
    //Interpolation
    currentMesh = (FEMMesh *)model.mesh;
    FEMInterpolateMeshToMesh *interpolateMesh = [[FEMInterpolateMeshToMesh alloc] init];
    if (maskName != NULL) {
        [interpolateMesh interpolateMesh:mesh toMesh:currentMesh oldVariables:tmpArryay newVariables:anArray model:model quadrantTree:NULL projector:nil mask:maskName unfoundNodes:NULL];
    } else {
        [interpolateMesh interpolateMesh:mesh toMesh:currentMesh oldVariables:tmpArryay newVariables:anArray model:model quadrantTree:NULL projector:projector mask:nil unfoundNodes:NULL];
    }
    
    onlyThis = YES;
    var = [self getVariableFrom:anArray model:model name:name onlySearch:&onlyThis maskName:NULL info:&stat];
    if ([var.name isEqualToString:@"flow solution"]) {
        tmp = [self getVariableFrom:anArray model:model name:@"velocity 1" onlySearch:&onlyThis maskName:NULL info:&stat];
        if (tmp != nil) {
            tmp.valid = YES;
            tmp.valuesChanged = YES;
        }
        
        tmp = [self getVariableFrom:anArray model:model name:@"velocity 2" onlySearch:&onlyThis maskName:NULL info:&stat];
        if (tmp != nil) {
            tmp.valid = YES;
            tmp.valuesChanged = YES;
        }
        
        if (var.dofs == 4) {
            tmp = nil;
            tmp = [self getVariableFrom:anArray model:model name:@"velocity 3" onlySearch:&onlyThis maskName:NULL info:&stat];
            if (tmp != nil) {
                tmp.valid = YES;
                tmp.valuesChanged = YES;
            }
        }
        
        tmp = [self getVariableFrom:anArray model:model name:@"pressure" onlySearch:&onlyThis maskName:NULL info:&stat];
        if (tmp != nil) {
            tmp.valid = YES;
            tmp.valuesChanged = YES;
        }
    } else if (var.dofs > 1){
        for (i=0; i<var.dofs; i++) {
            tmp = nil;
            aNumber = i+1;
            tmpname = [self appendNameFromString:pVar.name component:&aNumber];
            tmp = [self getVariableFrom:anArray model:model name:tmpname onlySearch:&onlyThis maskName:NULL info:&stat];
            if (tmp != nil) {
                tmp.valid = YES;
                tmp.valuesChanged = YES;
            }
        }
    }
    
    return var;
}

/*********************************************************************************************************
 Interpolate values in a curve given by linear table or splines.
*********************************************************************************************************/
-(double)interpolateCurveTvalues:(double * __nonnull)tValues fValues:(double * __nonnull)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double * __nullable)cubicCoeff {
    
    int i;
    double f;
    BOOL cubic;
    double tval[2], fval[2], coeffval[2];
    
    // This is a misuse of the interpolation in case of standard
    // dependency of type y = a*x
    if (n == 1) {
        return fValues[0] * t;
    }
    
    for (i=0; i<n; i++) {
        if (tValues[i] >= t) break;
    }
    if (i > n-1) i = n-1;
    if (i < 1) i = 1;
    
    cubic = (cubicCoeff != NULL) ? YES : NO;
    cubic = (cubic == YES && t >= tValues[0] && t <=tValues[n-1]) ? YES : NO;
    
    if (cubic == YES) {
        tval[0] = tValues[i-1];
        tval[1] = tValues[i];
        fval[0] = fValues[i-1];
        fval[1] = fValues[i];
        coeffval[0] = cubicCoeff[i-1];
        coeffval[1] = cubicCoeff[i];
        f = [self cublicSplineX:tval Y:fval R:coeffval T:t];
    } else {
        f = (t - tValues[i-1]) / (tValues[i] - tValues[i-1]);
        f = (1.0 - f)*fValues[i-1] + f*fValues[i];
    }
    return f;
}

/*********************************************************************************************************
    Derivative a curve given by linear table or splines.
*********************************************************************************************************/
-(double)derivateCurveTvalues:(double * __nonnull)tValues fValues:(double * __nonnull)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double * __nullable)cubicCoeff {
    
    int i;
    double f;
    BOOL cubic;
    double tval[2], fval[2], coeffval[2];
    
    for (i=0; i<n; i++) {
        if (tValues[i] >= t) break;
    }
    
    if (i > n-1) i = n-1;
    if (i < 1) i = 1;
    
    cubic = (cubicCoeff != NULL) ? YES : NO;
    cubic = (cubic == YES && t >= tValues[0] && t <= tValues[n-1]) ? YES : NO;
    
    if (cubic == YES) {
        tval[0] = tValues[i-1];
        tval[1] = tValues[i];
        fval[0] = fValues[i-1];
        fval[1] = fValues[i];
        coeffval[0] = cubicCoeff[i-1];
        coeffval[1] = cubicCoeff[i];
        f = [self cublicSplineX:tval Y:fval R:coeffval T:t];
    } else {
        f = (fValues[i] - fValues[i-1]) / (tValues[i] - tValues[i-1]);
    }
    return f;
}

-(void)solveLinearSystem2x2:(double[][2])a afterSolve:(double * __nonnull)x rightHandSide:(double * __nonnull)b {
    
    double detA;
    
    detA = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    
    if (detA == 0.0) {
        fatal("FEMUtilities:solveLinearSystem2x2", "Singular matrix, bad...");
        return;
    }
    
    detA = 1.0 / detA;
    x[0] = detA * ( a[1][1] * b[0] - a[0][1] * b[1] );
    x[1] = detA * ( a[0][0] * b[1] - a[1][0] * b[0] );
}

-(void)solveLinearSystem3x3:(double[][3])a afterSolve:(double * __nonnull)x rightHandSide:(double * __nonnull)b {
    
    double c[2][2], y[2], g[2], s, t, q;
    
    if ( (fabs(a[0][0]) > fabs(a[0][1])) && (fabs(a[0][0]) > fabs(a[0][2])) ) {
        q = 1.0 / a[0][0];
        s = q * a[1][0];
        t = q * a[2][0];
        c[0][0] = a[1][1] - s * a[0][1];
        c[0][1] = a[1][2] - s * a[0][2];
        c[1][0] = a[2][1] - t * a[0][1];
        c[1][1] = a[2][2] - t * a[0][2];
        
        g[0] = b[1] - s * b[0];
        g[1] = b[2] - t * b[0];
        [self solveLinearSystem2x2:c afterSolve:y rightHandSide:g];
        
        x[1] = y[0];
        x[2] = y[1];
        x[0] = q * ( b[0] - a[0][1] * x[1] - a[0][2] * x[2] );
    } else if (fabs(a[0][1]) > fabs(a[0][2])) {
        q = 1.0 / a[0][1];
        s = q * a[1][1];
        t = q * a[2][1];
        c[0][0] = a[1][0] - s * a[0][0];
        c[0][1] = a[1][2] - s * a[0][2];
        c[1][0] = a[2][0] - t * a[0][0];
        c[1][1] = a[2][2] - t * a[0][2];
        
        g[0] = b[1] - s * b[0];
        g[1] = b[2] - t * b[0];
        [self solveLinearSystem2x2:c afterSolve:y rightHandSide:g];
        
        x[0] = y[0];
        x[2] = y[1];
        x[1] = q * ( b[0] - a[0][0] * x[0] - a[0][2] * x[2] );
    } else {
        q = 1.0 / a[0][2];
        s = q * a[1][2];
        t = q * a[2][2];
        c[0][0] = a[1][0] - s * a[0][0];
        c[0][1] = a[1][1] - s * a[0][1];
        c[1][0] = a[2][0] - t * a[0][0];
        c[1][1] = a[2][1] - t * a[0][1];
        
        g[0] = b[1] - s * b[0];
        g[1] = b[2] - t * b[0];
        [self solveLinearSystem2x2:c afterSolve:y rightHandSide:g];
        
        x[0] = y[0];
        x[1] = y[1];
        x[2] = q * ( b[0] - a[0][0] * x[0] - a[0][1] * x[1] );
    }
}

-(FEMMatrix * __nonnull)meshProjectorMesh1:(FEMMesh * __nonnull)mesh1 mesh2:(FEMMesh * __nonnull)mesh2 model:(FEMModel * __nonnull)model useQuadrantTree:(BOOL * __nullable)quadrantTree transpose:(BOOL * __nullable)trans {
    
    FEMProjector *projector;
    FEMMatrix *projectorMatrix;
    
    FEMInterpolateMeshToMesh *interpolateMesh = [[FEMInterpolateMeshToMesh alloc] init];
    
    projector = [[FEMProjector alloc] init];
    if (quadrantTree != NULL) {
        [interpolateMesh interpolateQMesh:mesh1 toMesh:mesh2 oldVariables:nil newVariables:nil model:model quadrantTree:quadrantTree projector:projector mask:nil nodesPresent:NULL newMaskPerm:NULL];
    } else {
        [interpolateMesh interpolateQMesh:mesh1 toMesh:mesh2 oldVariables:nil newVariables:nil model:model quadrantTree:NULL projector:projector mask:nil nodesPresent:NULL newMaskPerm:NULL];
    }
    
    projectorMatrix = projector.matrix.copy;
    if (trans != NULL) {
        if (*trans == YES) {
            projectorMatrix = projector.tMatrix.copy;
        }
    }
    
    return projectorMatrix;
}

-(double)cublicSplineX:(double * __nonnull)x Y:(double * __nonnull)y R:(double * __nonnull)r T:(double)t {
    
    double s, a, b, c, h, lt;
    
    h = x[1] - x[0];
    a = -2.0 * (y[1] - y[0]) + (r[0] + r[1]) * h;
    b = 3.0 * (y[1] - y[0]) - (2.0*r[0] + r[1]) * h;
    c = r[0] * h;
    
    lt = (t - x[0]) / h;
    s = ( (a*lt + b) * lt + c ) / h;
    
    return s;
}

-(NSString * __nullable)appendNameFromString:(NSString * __nonnull)string component:(int * __nullable)component {
    
    int i, j, dofs, comp, dofsTot;
    NSString *str;
    NSMutableString *str1;
    NSRange ind, ind1;
    
    if (string == nil) return nil;
    
    ind = [string rangeOfString:@"["];
    
    comp = 0;
    if (component != NULL) comp = *component;
    
    if (ind.location == NSNotFound) {
        str1 = [NSMutableString stringWithString:string];
        if (comp > 0) {
            [str1 appendString:@" "];
            [str1 appendString:[NSString stringWithFormat:@"%d",comp]];
            str = [NSString stringWithString:str1];
        }
    } else {
        j = 0;
        dofsTot = 0;
        while (1) {
            str = [string substringFromIndex:j];
            if (j == 0) {
                ind = [string rangeOfString:@"["];
            } else {
                ind.location = -1;
            }
            ind1 = [str rangeOfString:@":"];
            if (ind1.location == NSNotFound) fatal("FEMUtilities:appendNameFromString", "Missing separator ':' in variable definition using '[ ]' syntax.");
            if (j == 0) {
                if (ind1.location < ind.location) fatal("FEMUtilities:appendNameFromString", "Syntax error in variable definition.");
            } else {
                if (ind1.location == 0) fatal("FEMUtilities:appendNameFromString", "Syntax error in variable definition.");
            }
            i = (int)ind1.location + 1;
            while (1) {
                if ([str characterAtIndex:i] == ' ' || [str characterAtIndex:i] == ']') break;
                i++;
            }
            dofs = [[[str substringFromIndex:ind1.location+1] substringToIndex:i-(ind1.location+1)] intValue];
            dofsTot = dofsTot + dofs;
            if (dofsTot >= *component) break;
            while ([str characterAtIndex:i] == ' ') {
                i++;
            }
            j = i;
        }
        str = [str substringFromIndex:ind.location+1];
        ind1 =  [str rangeOfString:@":"];
        str = [str substringToIndex:ind1.location];
        if (dofs > 1) {
            dofs = *component - dofsTot + dofs;
            str1 = [NSMutableString stringWithString:str];
            [str1 appendString:@" "];
            [str1 appendString:[NSString stringWithFormat:@"%d",dofs]];
            str = [NSString stringWithString:str1];
        }
    }
    
    return str;
}

-(NSString * __nonnull)appendNameFromVariable:(FEMVariable * __nonnull)variable component:(int * __nullable)component {
    
    NSString *str;
    NSMutableString *str1;
    
    if ([variable.name isEqualToString:@"flow solution"]) {
        str = @"flow solution";
        if (component == NULL) return str;
        if (*component == variable.dofs) {
            str = @"pressure";
            return str;
        } else {
            str1 = [NSMutableString stringWithString:@"velocity"];
            [str1 appendString:@" "];
            [str1 appendString:[NSString stringWithFormat:@"%d",*component]];
            str = [NSString stringWithString:str1];
        }
    } else {
        str = [self appendNameFromString:variable.name component:component];
    }
    
    return str;
}

/******************************************************************************************
    Given the body of a keyword find the first free keyword in a dictionary.
    This is intended use of this is in the solution construction to declare exported 
    variables without the risk of running over some existing ones.
******************************************************************************************/
-(NSString * __nonnull)nextFreeKeyword:(NSString * __nonnull)keyword0 dictionary:(NSMutableDictionary * __nonnull)dictionary {
    
    int i;
    NSMutableString *str;
    
    for (i=1; i<=9999; i++) {
        @autoreleasepool {
            str = [NSMutableString stringWithString:keyword0];
            [str appendString:@" "];
            [str appendString:[NSString stringWithFormat:@"%d",i]];
            if (dictionary[str] == nil) break;
        }
    }
    
    return str;
}

/******************************************************
    Check the feasibility of solution options
******************************************************/
-(void)checkOptionsInSolution:(FEMSolution * __nonnull)solution {
    
    NSString *str;
    
    if (solution.solutionInfo[@"linear system solver"] != nil) {
        str = solution.solutionInfo[@"linear system solver"];
        
        if ([str isEqualToString:@"direct"] == YES) {
            if (solution.solutionInfo[@"linear system direct method"] != nil) {
                str = solution.solutionInfo[@"linear system direct method"];
                
                //TODO: app support for parallel run where in case of a parallel run
                // the direct solver must be MUMPS
                if ([str isEqualToString:@"mumps"] == YES) fatal("FEMUtilities:checkOptionsInSolution", "Currenly no serial version of the mumps solver implemented.");
                
                if ([str isEqualToString:@"banded"] == YES) {
                    
                } else if ([str isEqualToString:@"umfpack"] == YES || [str isEqualToString:@"big umfpack"] == YES) {
                    
                } else if ([str isEqualToString:@"mumps"] == YES) {
                    
                } else if ([str isEqualToString:@"superlu"] == YES) {
                    fatal("FEMUtilities:checkOptionsInSolution", "SuperLU solver is not available.");
                } else if ([str isEqualToString:@"pardiso"] == YES) {
                    fatal("FEMUtilities:checkOptionsInSolution", "Pardiso solver is not available.");
                } else if ([str isEqualToString:@"cholmod"] == YES || [str isEqualToString:@"spqr"] == YES) {
                    fatal("FEMUtilities:checkOptionsInSolution", "Cholmod solver is not available.");
                } else {
                    fatal("FEMUtilities:checkOptionsInSolution", "Unknown direct solver method.");
                }
            } else {
                //TODO: add support for parallel run since then it should be mumps by default.
                fprintf(stdout, "FEMUtilities:checkOptionsInSolution: setting the linear system direct method to %s.\n", [str UTF8String]);
                [solution.solutionInfo setObject:@"umfpack" forKey:@"linear system direct method"];
            }
        }
    } // If "linear system solver" is not given, it will be set by default to iterative later in the processing
}


/********************************************************************************************************
    Add the generic components to each solution
    A few solutions are for historical reasons given a special treatment
********************************************************************************************************/
-(void)addEquationBasicsToSolution:(FEMSolution * __nonnull)solution name:(NSString * __nonnull)name model:(FEMModel * __nonnull)model transient:(BOOL)transient {
    
    int i, j, k, l, maxDim, minVal, maxVal, dofs, ndeg, maxNDOFs, maxDGDOFs, maxEDOFs, maxFDOFs, maxBDOFs,
        matrixFormat, nrows, nsize, type;
    int *perm;
    double initValue;
    NSString *eq, *str, *varName, *tmpName;
    NSMutableString *string;
    BOOL isAssemblySolution, isCoupledSolution, isBlockSolution, variableGlobal, variableOutput, found, stat,
         globalBubbles, bandwidthOptimize, discontinuousGalerkin;
    NSRange range;
    Element_t *elements = NULL, *edges = NULL, *faces = NULL;
    FEMVariable *variable = nil, *newVariable = nil;
    FEMListUtilities *listUtilities;
    FEMElementUtils *elementUtils;
    variableArraysContainer *bufferContainers = NULL, *varContainers = NULL;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    
    // If there is a matrix level "Flux Corrected Transport" then it's required
    // to use global matrices for time integration
    if ([solution.solutionInfo[@"linear system fct"] boolValue] == YES) {
        [listUtilities addLogicalInClassList:solution theVariable:@"use global mass matrix" withValue:YES];
    }
    
    if (solution.solutionInfo[@"equation"] != nil) {
        eq = solution.solutionInfo[@"equation"];
        elements = solution.mesh.getElements;
        maxDim = 0;
        for (i=0; i<solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements; i++) {
            if ([listUtilities checkElementEquation:model forElement:&elements[i] andEquation:eq]) {
                maxDim = max(elements[i].Type.dimension, maxDim);
            }
        }
        [solution.solutionInfo setObject:@(maxDim) forKey:@"active mesh dimension"];
    }
    
    // TODO: Here Elmer calls a procedure post-fixed with "_Init" for a given solver.
    // Do we need to do that?
    
    solution.solutionMode = SOLUTION_MODE_DEFAULT;
    if ([solution.solutionInfo[@"auxiliary solution"] boolValue] == YES) solution.solutionMode = SOLUTION_MODE_AUXILIARY;
    if ([solution.solutionInfo[@"coupled solution"] boolValue] == YES) solution.solutionMode = SOLUTION_MODE_COUPLED;
    if ([solution.solutionInfo[@"block solution"] boolValue] == YES) solution.solutionMode = SOLUTION_MODE_BLOCK;
    if ([solution.solutionInfo[@"assembly solution"] boolValue] == YES) solution.solutionMode = SOLUTION_MODE_ASSEMBLY;
    
    if (solution.solutionMode == SOLUTION_MODE_DEFAULT) {
        if (solution.solutionInfo[@"equation"] != nil) {
            eq = solution.solutionInfo[@"equation"];
        } else {
            solution.solutionMode = SOLUTION_MODE_AUXILIARY;
        }
    }
    
    isCoupledSolution = (solution.solutionMode == SOLUTION_MODE_COUPLED) ? YES : NO;
    isBlockSolution = (solution.solutionMode == SOLUTION_MODE_BLOCK) ? YES : NO;
    isAssemblySolution = (solution.solutionMode == SOLUTION_MODE_ASSEMBLY) ? YES : NO;
    isAssemblySolution = (isAssemblySolution == YES || (isCoupledSolution == YES && (solution.plugInPrincipalClassInstance == nil && solution.hasBuiltInSolution == NO)) || (isBlockSolution == YES && (solution.plugInPrincipalClassInstance == nil && solution.hasBuiltInSolution == NO))) ? YES : NO;
    
    // Default order of equation
    solution.order = 1;
    solution.timeOrder = 1;
    
    found = NO;
    if (transient == YES) {
        if (solution.solutionInfo[@"time stepping method"] != nil) {
            str = solution.solutionInfo[@"time stepping method"];
            found = YES;
        } else {
            str = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"time stepping method" info:&found];
            if (found == YES) {
                [solution.solutionInfo setObject:str forKey:@"time stepping method"];
            }
        }
        
        if (found == YES) {
            if ([str isEqualToString:@"bdf"] == YES) {
                if (solution.solutionInfo[@"bdf order"] != nil) {
                    solution.order = [solution.solutionInfo[@"bdf order"] intValue];
                    if (solution.order < 1) solution.order = 1;
                    if (solution.order > 5) solution.order = 5;
                } else {
                    minVal = 1;
                    maxVal = 5;
                    solution.order = [listUtilities listGetInteger:model inArray:model.simulation.valuesList forVariable:@"bdf order" info:&found minValue:&minVal maxValue:&maxVal];
                }
                if (found == NO) {
                    solution.order = 2;
                    fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: bdf order set by default to 2.\n");
                }
            } else if ([str isEqualToString:@"runge-kutta"]) {
                minVal = 2;
                maxVal = 4;
                solution.order = [listUtilities listGetInteger:model inArray:model.simulation.valuesList forVariable:@"runge-kutta order" info:&found minValue:&minVal maxValue:&maxVal];
                if (found == NO) solution.order = 2;
            }
        } else {
            fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: time stepping method set by default to implicit Euler.\n");
            [solution.solutionInfo setObject:@"implicit euler" forKey:@"time stepping method"];
        }
    }
    
    dofs = model.dimension;
    initValue = 0.0;
    
    // These are historical solutions that may be built-in on some MDFs
    // Special strategies are used for them.
    // TODO: we may migrate them to plug-ins?
    
    if ([name isEqualToString:@"navier-stokes"] == YES) {
        if (solution.solutionInfo[@"variable"] == nil) {
            dofs = model.dimension;
            if (model.coordinates == cylindric_symmetric) dofs++;
            if (dofs == 3) {
                [solution.solutionInfo setObject:@"flow solution[velocity:3 pressure:1]" forKey:@"variable"];
            } else {
                [solution.solutionInfo setObject:@"flow solution[velocity:2 pressure:1]" forKey:@"variable"];
            }
        }
        initValue = 1.0e-6;
    } else if ([name isEqualToString:@"magnetic induction"] == YES) {
        if (solution.solutionInfo[@"variable"] == nil) {
            [solution.solutionInfo setObject:@"-dofs 3 magnetic field" forKey:@"variable"];
            [solution.solutionInfo setObject:@"electric current[electric current:3]"
                                      forKey:[self nextFreeKeyword:@"exported variable" dictionary:solution.solutionInfo]];
        }
    } else if ([name isEqualToString:@"stress analysis"] == YES) {
        if (solution.solutionInfo[@"variable"] == nil) {
            if (dofs == 2) {
                [solution.solutionInfo setObject:@"-dofs 2 displacement" forKey:@"variable"];
            } else {
                [solution.solutionInfo setObject:@"-dofs 3 displacement" forKey:@"variable"];
            }
        }
    } else if ([name isEqualToString:@"mesh update"] == YES) {
        if (solution.solutionInfo[@"variable"] == nil) {
            if (dofs == 2) {
                [solution.solutionInfo setObject:@"-dofs 2 mesh update" forKey:@"variable"];
            } else {
                [solution.solutionInfo setObject:@"-dofs 3 mesh update" forKey:@"variable"];
            }
        }
        
        if (transient == YES) {
            if (dofs == 2) {
                [solution.solutionInfo setObject:@"-dofs 2 mesh velocity"
                                          forKey:[self nextFreeKeyword:@"exported variable" dictionary:solution.solutionInfo]];
            } else {
                [solution.solutionInfo setObject:@"-dofs 3 mesh velocity"
                                          forKey:[self nextFreeKeyword:@"exported variable" dictionary:solution.solutionInfo]];
            }
        }
    } else if ([name isEqualToString:@"heat equation"] == YES) {
        if (solution.solutionInfo[@"variable"] == nil) {
            [solution.solutionInfo setObject:@"temperature" forKey:@"variable"];
        }
        
        if (solution.solutionInfo[@"radiation solver"] == nil) {
            [solution.solutionInfo setObject:@YES forKey:@"radiation solver"];
        }
    }
    
    // We allocate memory for the string to hold the name of the plug-in only if we find that we are
    // using a plug-in when we parse the MDF. So if we get nil, we are not using a plug-in.
    if (solution.plugInName != nil) { // We are woking with a solution computer provided by a plug-in
        // Load the plug-in bundle
        BOOL useAppSupportPath = NO;
        NSBundle *solutionBundle = [self loadBundle:solution.plugInName useApplicationSupportPath:&useAppSupportPath];
        if (solutionBundle == nil) {
            fprintf(stderr, "FEMUtilities:addEquationBasicsToSolution: error loading plug-in bundle.\n");
            fatal("FEMUtilities:addEquationBasicsToSolution");
        }
        // Instantiate the plug-in principal class
        if ([solution instantiatePrincipalClassFromPlugIn:solutionBundle] == NO) {
            fprintf(stderr, "FEMUtilities:addEquationBasicsToSolution: error instanciating plug-in principal class.\n");
            fatal("FEMUtilities:addEquationBasicsToSolution");
        }
    } else { // We are woking with a built-in solution computer
        solution.hasBuiltInSolution = YES;
    }
    
    //Initialize and get the variable
    solution.timeOrder = 0;
    solution.matrix = nil;
    
    if (solution.solutionInfo[@"variable"] == nil) {
        //Variable does not exist
        variable = [[FEMVariable alloc] init];
        solution.variable = variable;
    } else if (isCoupledSolution == YES && (solution.plugInPrincipalClassInstance == nil && solution.hasBuiltInSolution == NO)) {
        // Coupled solver may inherit the matrix only if procedure is given
        
    } else if (isBlockSolution == YES && (solution.plugInPrincipalClassInstance == nil && solution.hasBuiltInSolution == NO)) {
        // Block solver may inherit the matrix only if procedure is given

    } else {
        varName = solution.solutionInfo[@"variable"];
        
        // It may be a normal field variable or a global (0D) variable
        variableGlobal = [solution.solutionInfo[@"variable global"] boolValue];
        
        if (solution.solutionInfo[@"variable output"] != nil) {
            variableOutput = [solution.solutionInfo[@"variable output"] boolValue];
        } else variableOutput = YES;
        
        if (solution.solutionInfo[@"variable dofs"] != nil) {
            dofs = [solution.solutionInfo[@"variable dofs"] intValue];
            if (dofs < 1) dofs = 1;
        } else {
            j = 0;
            dofs = 0;
            str = [NSString stringWithString:varName];
            while (1) {
                str = [str substringFromIndex:j];
                range = [str rangeOfString:@":"];
                if (range.location == NSNotFound) break;
                i = (int)range.location + 1;
                while (1) {
                    if ([str characterAtIndex:i] == ' ' || [str characterAtIndex:i] == ']') break;
                    i++;
                }
                k = [[[str substringFromIndex:range.location+1] substringToIndex:i-(range.location+1)] intValue];
                // If we get 0 then we are dealing with an invalid dof number
                if (k == 0) {
                    fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: invalid dof in variable definiton.\n");
                    fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: the incorrect value was: %s.\n",
                          [[[str substringFromIndex:range.location+1] substringToIndex:i-(range.location+1)] UTF8String]);
                }
                dofs = dofs + k;
                j = i+1;
            }
        }
        
        range = [varName rangeOfString:@"-"];
        if (range.location != NSNotFound) {
            if ([[[varName substringFromIndex:range.location] substringToIndex:(range.location+9)-range.location] isEqualToString:@"-nooutput"] == YES) {
                variableOutput = NO;
                i = (int)(range.location+9);
                while ([varName characterAtIndex:i] == ' ') {
                    i++;
                }
                varName = [varName substringFromIndex:i];
            } else if ([[[varName substringFromIndex:range.location] substringToIndex:(range.location+7)-range.location] isEqualToString:@"-global"] == YES) {
                variableGlobal = YES;
                i = (int)(range.location+7);
                while ([varName characterAtIndex:i] == ' ') {
                    i++;
                }
                varName = [varName substringFromIndex:i];
            } else if ([[[varName substringFromIndex:range.location] substringToIndex:(range.location+5)-range.location] isEqualToString:@"-dofs"] == YES) {
                i = (int)(range.location+5);
                while ([varName characterAtIndex:i] == ' ') {
                    i++;
                }
                j = i;
                while (1) {
                    if ([varName characterAtIndex:j] == ' ') break;
                    j++;
                }
                dofs = [[[varName substringFromIndex:i] substringToIndex:j-i] intValue];
                // If we get 0 then we are dealing with an invalid dof number
                if (dofs == 0) {
                    fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: invalid dof in variable definiton.\n");
                    fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: the incorrect value was: %s.\n", [[[varName substringFromIndex:i] substringToIndex:j-i] UTF8String]);
                }
                while ([varName characterAtIndex:j] == ' ') {
                    j++;
                }
                varName = [varName substringFromIndex:j];
            }
        }
        if (dofs == 0) dofs = 1;
        
        // If the variable is global then it has nothing to do with the mesh
        // and it may be simply allocated
        if (variableGlobal == YES) {
            solution.solutionMode = SOLUTION_MODE_GLOBAL;
            bufferContainers = allocateVariableContainer();
            bufferContainers->Values = doublevec(0, dofs-1);
            memset( bufferContainers->Values, 0.0, dofs*sizeof(double) );
            bufferContainers->sizeValues = dofs;
            bufferContainers->Perm = NULL;
            [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:varName dofs:dofs container:bufferContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
            solution.variable = [self getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:NULL info:&found];
            if (dofs > 1) {
                for (i=1; i<=dofs; i++) {
                    tmpName = [self appendNameFromString:varName component:&i];
                    bufferContainers->ComponentValues = malloc ( 1 * sizeof ( double * ));
                    bufferContainers->ComponentValues[0] = &bufferContainers->Values[i-1];
                    bufferContainers->sizeComponentValues = 1;
                    [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:tmpName dofs:1 container:bufferContainers component:YES ifOutput:NULL ifSecondary:NULL type:NULL];
                }
            }
            free(bufferContainers);
        } else {
            // If variable is a field variable create a permutation and matrix related to it
            if (solution.solutionInfo[@"equation"] != nil) {
                eq = solution.solutionInfo[@"equation"];
            } else {
                fatal("FEMUtilities:addEquationBasicsToSolution", "Variable exists but the equation is not defined.");
            }
            found = NO;
            for (FEMEquation *equation in model.equations) {
                if ([listUtilities listGetLogical:model inArray:equation.valuesList forVariable:eq info:&stat] == YES) {
                    found = YES;
                    break;
                }
            }
            if (found == NO) {
                fprintf(stderr, "FEMUtilities:addEquationBasicsToSolution: variable %s exists but it's not associated to any equation.\n", [varName UTF8String]);
                fatal("FEMUtilities:addEquationBasicsToSolution");
            }
            
            // Compute the size of the permutation vector
            ndeg = 0;
            if (YES) {
                eq = solution.solutionInfo[@"equation"];
                maxNDOFs = 0;
                maxDGDOFs = 0;
                elements = solution.mesh.getElements;
                faces = solution.mesh.getFaces;
                edges = solution.mesh.getEdges;
                for (i=0; i<solution.mesh.numberOfBulkElements; i++) {
                    maxDGDOFs = max(maxDGDOFs, elements[i].DGDOFs);
                    maxNDOFs = max(maxNDOFs, elements[i].NDOFs);
                }
                
                maxEDOFs = 0;
                for (i=0; i<solution.mesh.numberOfEdges; i++) {
                    maxEDOFs = max(maxEDOFs, edges[i].BDOFs);
                }
                
                maxFDOFs = 0;
                for (i=0; i<solution.mesh.numberOfFaces; i++) {
                    maxFDOFs = max(maxFDOFs, faces[i].BDOFs);
                }
                
                maxBDOFs = 0;
                for (i=0; i<solution.mesh.numberOfBulkElements; i++) {
                    maxBDOFs = max(maxBDOFs, elements[i].BDOFs);
                }
                
                if (solution.solutionInfo[@"bubbles in global system"] != nil) {
                    globalBubbles = [solution.solutionInfo[@"bubbles in global system"] boolValue];
                } else globalBubbles = YES;
                
                ndeg = ndeg + solution.mesh.numberOfNodes;
                if (maxEDOFs > 0) ndeg = ndeg + maxEDOFs * solution.mesh.numberOfEdges;
                if (maxFDOFs > 0) ndeg = ndeg + maxFDOFs * solution.mesh.numberOfFaces;
                if (globalBubbles == YES) ndeg = ndeg + maxBDOFs * solution.mesh.numberOfBulkElements;
                if (solution.solutionInfo[@"discontinuous galerkin"] != nil) {
                    if ([solution.solutionInfo[@"discontinuous galerkin"] boolValue] == YES) {
                        ndeg = max( ndeg, maxDGDOFs * (solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements) );
                    }
                }
            }
            
            if (solution.solutionInfo[@"radiation solution"] != nil) {
                if ([solution.solutionInfo[@"radiation solution"] boolValue] == YES) {
                    //TODO: Need to implement this
                }
            }
            
            if (solution.solutionInfo[@"optimize bandwidth"] != nil) {
                bandwidthOptimize = [solution.solutionInfo[@"optimize bandwidth"] boolValue];
            } else bandwidthOptimize = YES;
            [self checkOptionsInSolution:solution];
            
            perm = intvec(0, ndeg-1);
            memset( perm, -1, ndeg*sizeof(int) );
            
            elementUtils = [[FEMElementUtils alloc] init];
            matrixFormat = MATRIX_CRS;
            if (solution.solutionInfo[@"discontinuous galerkin"] != nil) {
                discontinuousGalerkin = [solution.solutionInfo[@"discontinuous galerkin"] boolValue];
            } else discontinuousGalerkin = NO;
            solution.matrix = [elementUtils createMatrixInModel:model forSolution:solution mesh:solution.mesh dofs:dofs permutation:perm sizeOfPermutation:ndeg matrixFormat:matrixFormat optimizeBandwidth:bandwidthOptimize equationName:eq discontinuousGalerkinSolution:&discontinuousGalerkin globalBubbles:&globalBubbles nodalDofsOnly:NULL projectorDofs:NULL];
            nrows = dofs * ndeg;
            if (solution.matrix != nil) nrows = solution.matrix.numberOfRows;
            
            // Basically the solution could be matrix free but still the matrix is used
            // here temporarily since it's needed when making the permutation vector
            if (solution.solutionInfo[@"no matrix"] != nil) {
                if ([solution.solutionInfo[@"no matrix"] boolValue] == YES) {
                    solution.solutionMode = SOLUTION_MODE_MATRIXFREE;
                    [solution.matrix deallocation];
                    solution.matrix = nil;
                }
            }
            
            if (nrows > 0) {
                bufferContainers = allocateVariableContainer();
                bufferContainers->Values = doublevec(0, nrows-1);
                for (i=0; i<nrows; i++) {
                    bufferContainers->Values[i] = initValue;
                }
                bufferContainers->sizeValues = nrows;
                bufferContainers->Perm = perm;
                bufferContainers->sizePerm = ndeg;
                [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:varName dofs:dofs container:bufferContainers component:NO ifOutput:&variableOutput ifSecondary:NULL type:NULL];
                solution.variable = [self getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:NULL info:&found];
                                
                if (dofs > 1) {
                    for (i=1; i<=dofs; i++) {
                        tmpName = [self appendNameFromString:varName component:&i];
                        bufferContainers->ComponentValues = malloc ( (nrows/dofs) * sizeof ( double * ));
                        k = 0;
                        for (j=(i-1); j<=nrows-dofs+(i-1); j+=dofs) {
                            bufferContainers->ComponentValues[k] = &bufferContainers->Values[j];
                            k++;
                        }
                        bufferContainers->sizeComponentValues = nrows/dofs;
                        [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:tmpName dofs:1 container:bufferContainers component:YES ifOutput:NULL ifSecondary:NULL type:NULL];
                    }                    
                }
                free(bufferContainers);
            }
            
            //TODO: add support for parallel run here
            if (solution.solutionInfo[@"discontinuous galerkin"] != nil) {
                if ([solution.solutionInfo[@"discontinuous galerkin"] boolValue] == YES) solution.variable.type = VARIABLE_ON_NODES_ON_ELEMENTS;
            }
        }
    }
    
    // Add the exported variables which are typically auxiliary variables derived from
    // the solution without their own matrix equation 
    l = 0;
    while (1) {
        l++;
        str = [self appendNameFromString:@"exported variable" component:&l];
        if (solution.solutionInfo[str] == nil) break;
        varName = solution.solutionInfo[str];
        
        string = (NSMutableString *)[self appendNameFromString:@"exported variable" component:&l];
        [string appendString:@" output"];
        if (solution.solutionInfo[string] != nil) {
            variableOutput = [solution.solutionInfo[string] boolValue];
        } else variableOutput = YES;
        
        string = (NSMutableString *)[self appendNameFromString:@"exported variable" component:&l];
        [string appendString:@" dofs"];
        if (solution.solutionInfo[string] != nil) {
            dofs = [solution.solutionInfo[string] intValue];
        } else {
            j = 0;
            dofs = 0;
            str = [NSString stringWithString:varName];
            while (1) {
                str = [str substringFromIndex:j];
                range = [str rangeOfString:@":"];
                if (range.location == NSNotFound) break;
                i = (int)range.location + 1;
                while (1) {
                    if ([str characterAtIndex:i] == ' ' || [str characterAtIndex:i] == ']') break;
                    i++;
                }
                k = [[[str substringFromIndex:range.location+1] substringToIndex:i-(range.location+1)] intValue];
                // If we get 0 then we are dealing with an invalid dof number
                if (k == 0) {
                    fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: invalid dof in variable definiton.\n");
                    fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: the incorrect value was: %s.\n",
                          [[[str substringFromIndex:range.location+1] substringToIndex:i-(range.location+1)] UTF8String]);
                }
                dofs = dofs + k;
                j = i+1;
            }
        }
        
        variableOutput = YES;
        variableGlobal = NO;
        
        range = [varName rangeOfString:@"-"];
        if (range.location != NSNotFound) {
            if ([[[varName substringFromIndex:range.location] substringToIndex:(range.location+9)-range.location] isEqualToString:@"-nooutput"] == YES) {
                variableOutput = NO;
                i = (int)(range.location+9);
                while ([varName characterAtIndex:i] == ' ') {
                    i++;
                }
                varName = [varName substringFromIndex:i];
            } else if ([[[varName substringFromIndex:range.location] substringToIndex:(range.location+7)-range.location] isEqualToString:@"-global"] == YES) {
                variableGlobal = YES;
                i = (int)(range.location+7);
                while ([varName characterAtIndex:i] == ' ') {
                    i++;
                }
                varName = [varName substringFromIndex:i];
            } else if ([[[varName substringFromIndex:range.location] substringToIndex:(range.location+5)-range.location] isEqualToString:@"-dofs"] == YES) {
                i = (int)(range.location+5);
                while ([varName characterAtIndex:i] == ' ') {
                    i++;
                }
                j = i;
                while (1) {
                    if ([varName characterAtIndex:j] == ' ') break;
                    j++;
                }
                dofs = [[[varName substringFromIndex:i] substringToIndex:j-i] intValue];
                // If we get 0 then we are dealing with an invalid dof number
                if (dofs == 0) {
                    fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: invalid dof in variable definiton.\n");
                    fprintf(stdout, "FEMUtilities:addEquationBasicsToSolution: the incorrect value was: %s.\n", [[[varName substringFromIndex:i] substringToIndex:j-i] UTF8String]);
                }
                while ([varName characterAtIndex:j] == ' ') {
                    j++;
                }
                varName = [varName substringFromIndex:j];
            }
        }
        if (dofs == 0) dofs = 1;

        newVariable = [self getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:NULL info:&found];
        
        if (newVariable == nil) {
            if (variableGlobal == YES) {
                nsize = dofs;
                perm = NULL;
            } else {
                varContainers = solution.variable.getContainers;
                nsize = dofs * varContainers->sizeValues / solution.variable.dofs;
                perm = varContainers->Perm;
            }
            
            if (perm != NULL) {
                bufferContainers = allocateVariableContainer();
                bufferContainers->Values = doublevec(0, nsize-1);
                memset( bufferContainers->Values, 0.0, nsize*sizeof(double) );
                bufferContainers->sizeValues = nsize;
                bufferContainers->Perm = perm;
                bufferContainers->sizePerm = varContainers->sizePerm;
                type = solution.variable.type;
                [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:varName dofs:dofs container:bufferContainers component:NO ifOutput:&variableOutput ifSecondary:NULL type:&type];
            } else {
                bufferContainers = allocateVariableContainer();
                bufferContainers->Values = doublevec(0, nsize-1);
                memset( bufferContainers->Values, 0.0, nsize*sizeof(double) );
                bufferContainers->sizeValues = nsize;
                type = solution.variable.type;
                [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:varName dofs:dofs container:bufferContainers component:NO ifOutput:&variableOutput ifSecondary:NULL type:&type];
            }
            
            if (dofs > 1 && variableOutput == NO) {
                for (i=1; i<=dofs; i++) {
                    tmpName = [self appendNameFromString:varName component:&i];
                    bufferContainers->ComponentValues = malloc ( (nsize/dofs) * sizeof ( double * ));
                    k = 0;
                    for (j=(i-1); j<=nsize-dofs+(i-1); j+=dofs) {
                        bufferContainers->ComponentValues[k] = & bufferContainers->Values[j];
                        k++;
                    }
                    bufferContainers->sizeComponentValues = nsize/dofs;
                    type = solution.variable.type;
                    if (perm != NULL) {
                        [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:tmpName dofs:1 container:bufferContainers component:YES ifOutput:&variableOutput ifSecondary:NULL type:&type];
                    } else {
                        [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:tmpName dofs:1 container:bufferContainers component:YES ifOutput:&variableOutput ifSecondary:NULL type:&type];
                    }
                }
            }
            free(bufferContainers);
        }
    }
    
    // Check for special solutions to be obtained only at
    // a certain instances during simulation
    solution.solutionSolveWhen = SOLUTION_SOLVE_ALWAYS;
    
    if (solution.solutionInfo[@"invoke solution computer"] != nil) {
        if ([solution.solutionInfo[@"invoke solution computer"] isEqualToString:@"never"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_NEVER;
        } else if ([solution.solutionInfo[@"invoke solution computer"] isEqualToString:@"always"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_ALWAYS;
        } else if ([solution.solutionInfo[@"invoke solution computer"] isEqualToString:@"after all"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_ALL;
        } else if ([solution.solutionInfo[@"invoke solution computer"] isEqualToString:@"before all"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_ALL;
        } else if ([solution.solutionInfo[@"invoke solution computer"] isEqualToString:@"before time step"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_TIME;
        } else if ([solution.solutionInfo[@"invoke solution computer"] isEqualToString:@"after time step"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_TIME;
        } else if ([solution.solutionInfo[@"invoke solution computer"] isEqualToString:@"before saving"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_SAVE;
        } else if ([solution.solutionInfo[@"invoke solution computer"] isEqualToString:@"after saving"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_SAVE;
        } else {
            solution.solutionSolveWhen = SOLUTION_SOLVE_ALWAYS;
        }
    }
    
    if (solution.solutionInfo[@"before all"] != nil) {
        if ([solution.solutionInfo[@"before all"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_ALL;
        }
    } else if (solution.solutionInfo[@"before simulation"] != nil) {
        if ([solution.solutionInfo[@"before simulation"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_ALL;
        }
    } else if (solution.solutionInfo[@"after all"] != nil) {
        if ([solution.solutionInfo[@"after all"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_ALL;
        }
    } else if (solution.solutionInfo[@"after simulation"] != nil) {
        if ([solution.solutionInfo[@"after simulation"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_ALL;
        }
    } else if (solution.solutionInfo[@"before time step"] != nil) {
        if ([solution.solutionInfo[@"before time step"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_TIME;
        }
    } else if (solution.solutionInfo[@"after time step"] != nil) {
        if ([solution.solutionInfo[@"after time step"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_TIME;
        }
    } else if (solution.solutionInfo[@"before saving"] != nil) {
        if ([solution.solutionInfo[@"before saving"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_SAVE;
        }
    } else if (solution.solutionInfo[@"after saving"] != nil) {
        if ([solution.solutionInfo[@"after saving"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_SAVE;
        }
    }
    
    // TODO: Do we need to add support for LinSolve?
}

/********************************************************************************************************
    Add information that is typically only needed if there's a matrix equation to work with.
    This should be called only after both the solution vector and matrix have been created.
********************************************************************************************************/
-(void)addEquationToSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model transient:(BOOL)transient {
    
    int i, j, k, l, dofs, nrows, n, mgLevels;
    double *sol;
    NSString *tmpName, *methodName;
    NSMutableString *varName, *str;
    BOOL found, onlySearch, secondary, variableOutput, harmonicAnal, eigenAnal, complexFlag, multigridActive, mgAlgebraic;
    FEMCore *core;
    FEMVariable *var;
    FEMMesh *newMesh, *oldMesh;
    FEMMatrix *oldMatrix, *newMatrix, *saveMatrix;
    FEMListUtilities *listUtilities;
    FEMMeshUtils *meshUtilities;
    variableArraysContainer *bufferContainers = NULL, *variableContainers = NULL, *varContainers = NULL;
    matrixArraysContainer *matContainers = NULL;
    listBuffer freqv = { NULL, NULL, NULL, NULL, 0, 0, 0 };
    
    solution.doneTime = 0;
    if (solution.variable == nil) return;
    variableContainers = solution.variable.getContainers;
    if (variableContainers->Values == NULL) return;
    
    // If soft limiters are applied then also loads must be calculated
    if ([solution.solutionInfo[@"apply limiter"] boolValue] == YES) {
        [solution.solutionInfo setObject:@YES forKey:@"calculate loads"];
    }
    
    // Create the variable needed for this computation of nodal loads: r = b-Ax
    if ([solution.solutionInfo[@"calculate loads"] boolValue] == YES) {
        varName = [NSMutableString stringWithString:[solution.variable canonicalizeName]];
        [varName appendString:@" loads"];
        var = [self getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:NULL info:&found];
        if (var == nil) {
            sol = doublevec(0, variableContainers->sizeValues-1);
            dofs = solution.variable.dofs;
            memset( sol, 0.0, variableContainers->sizeValues*sizeof(double) );
            nrows = variableContainers->sizeValues;
            variableOutput = solution.variable.output;
            
            bufferContainers = allocateVariableContainer();
            bufferContainers->Values = sol;
            bufferContainers->sizeValues = variableContainers->sizeValues;
            bufferContainers->Perm = variableContainers->Perm;
            bufferContainers->sizePerm = variableContainers->sizePerm;
            [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:varName dofs:dofs container:bufferContainers component:NO ifOutput:&variableOutput ifSecondary:NULL type:NULL];
            
            if (dofs > 1) {
                for (i=1; i<=dofs; i++) {
                    tmpName = [self appendNameFromString:varName component:&i];
                    bufferContainers->ComponentValues = malloc ( (variableContainers->sizeValues/dofs) * sizeof ( double * ));
                    k = 0;
                    for (j=(i-1); j<=nrows-dofs+(i-1); j+=dofs) {
                        bufferContainers->ComponentValues[k] = &sol[j];
                        k++;
                    }
                    bufferContainers->sizeComponentValues = variableContainers->sizeValues/dofs;
                     [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:tmpName dofs:1 container:bufferContainers component:YES ifOutput:NULL ifSecondary:NULL type:NULL];
                }
            }
            sol = NULL;
            free(bufferContainers);
        }
    }
    
    solution.nOfEigenValues = 0;
    solution.multigridLevel = 1;
    solution.multiGridTotal = 0;
    solution.multigridSolution = NO;
    solution.multigridEqualPlit = NO;
    
    harmonicAnal = [solution.solutionInfo[@"harmonic analysis"] boolValue];
    
    if (solution.matrix != nil) {
        matContainers = solution.matrix.getContainers;
        if (matContainers->RHS == NULL) {
            matContainers->RHS = doublevec(0, solution.matrix.numberOfRows-1);
            matContainers->sizeRHS = solution.matrix.numberOfRows;
            memset( matContainers->RHS, 0.0, solution.matrix.numberOfRows*sizeof(double) );
            
            if (harmonicAnal == YES) {
                matContainers->RHS_im = doublevec(0, solution.matrix.numberOfRows-1);
                matContainers->sizeRHS_im = solution.matrix.numberOfRows;
                memset( matContainers->RHS_im, 0.0, solution.matrix.numberOfRows*sizeof(double) );
            }
        }
    }
    
    eigenAnal = [solution.solutionInfo[@"eigen analysis"] boolValue];
    
    if (transient == YES && eigenAnal == NO && harmonicAnal == NO) {
        solution.timeOrder = 1;
        if (solution.solutionInfo[@"time derivative order"] != nil) {
            k = [solution.solutionInfo[@"time derivative order"] intValue];
            if (k < 0) k = 0;
            if (k > 2) k = 2;
            solution.timeOrder = min(max(1, k), 2);
        }
        
        if (solution.matrix != nil) {
            matContainers->Force = doublematrix(0, solution.matrix.numberOfRows-1, 0, (solution.timeOrder+1)-1);
            matContainers->size1force = solution.matrix.numberOfRows;
            matContainers->size2Force = solution.timeOrder + 1;
            memset( *matContainers->Force, 0.0, (solution.matrix.numberOfRows*(solution.timeOrder+1))*sizeof(double) );
        }
        
        if (variableContainers->PrevValues == NULL) {
            if (solution.timeOrder == 2) {
                variableContainers->PrevValues = doublematrix(0, variableContainers->sizeValues-1, 0, 4);
                variableContainers->size1PrevValues = variableContainers->sizeValues;
                variableContainers->size2PrevValues = 5;
            } else if (solution.order > solution.timeOrder) {
                variableContainers->PrevValues = doublematrix(0, variableContainers->sizeValues-1, 0, solution.order-1);
                variableContainers->size1PrevValues = variableContainers->sizeValues;
                variableContainers->size2PrevValues = solution.order;
            } else {
                variableContainers->PrevValues = doublematrix(0, variableContainers->sizeValues-1, 0, solution.timeOrder-1);
                variableContainers->size1PrevValues = variableContainers->sizeValues;
                variableContainers->size2PrevValues = solution.timeOrder;

            }
            memset( *variableContainers->PrevValues, 0.0, (variableContainers->size1PrevValues*variableContainers->size2PrevValues)*sizeof(double) );
            
            if (solution.variable.dofs > 1) {
                if ([[solution.variable canonicalizeName] isEqualToString:@"flow solution"]) {
                    onlySearch = YES;
                    for (k=1; k<=solution.variable.dofs-1; k++) {
                        str = [NSMutableString stringWithString:@"velocity "];
                        [str appendString:[NSString stringWithFormat:@"%d",k]];
                        var = [self getVariableFrom:solution.mesh.variables model:model name:str onlySearch:&onlySearch maskName:NULL info:&found];
                        varContainers = var.getContainers;
                        varContainers->ComponentPrevValues = malloc ( (variableContainers->size1PrevValues/solution.variable.dofs) * sizeof ( double ** ));
                        for (i=0; i<(variableContainers->size1PrevValues/solution.variable.dofs); i++) {
                            varContainers->ComponentPrevValues[i] = malloc ( variableContainers->size2PrevValues * sizeof ( double * ));
                        }
                        l = 0;
                        for (i=k-1; i<variableContainers->size1PrevValues; i+=solution.variable.dofs) {
                            for (j=0; j<variableContainers->size2PrevValues; j++) {
                                varContainers->ComponentPrevValues[l][j] = &variableContainers->PrevValues[i][j];
                            }
                            l++;
                        }
                        varContainers->size1ComponentPrevValues = (variableContainers->size1PrevValues/solution.variable.dofs);
                        varContainers->size2ComponentPrevValues = variableContainers->size2PrevValues;
                    }
                    varContainers = NULL;
                    var = [self getVariableFrom:solution.mesh.variables model:model name:@"pressure" onlySearch:&onlySearch maskName:NULL info:&found];
                    varContainers = var.getContainers;
                    varContainers->ComponentPrevValues = malloc ( (variableContainers->size1PrevValues/solution.variable.dofs) * sizeof ( double ** ));
                    for (i=0; i<(variableContainers->size1PrevValues/solution.variable.dofs); i++) {
                        varContainers->ComponentPrevValues[i] = malloc ( variableContainers->size2PrevValues * sizeof ( double * ));
                    }
                    l = 0;
                    for (i=3; i<variableContainers->size1PrevValues; i+=solution.variable.dofs) {
                        for (j=0; j<variableContainers->size2PrevValues; j++) {
                            varContainers->ComponentPrevValues[l][j] = &variableContainers->PrevValues[i][j];
                        }
                        l++;
                    }
                    varContainers->size1ComponentPrevValues = (variableContainers->size1PrevValues/solution.variable.dofs);
                    varContainers->size2ComponentPrevValues = variableContainers->size2PrevValues;
                } else {
                    onlySearch = YES;
                    for (k=1; k<=solution.variable.dofs; k++) {
                        str = (NSMutableString *)[self appendNameFromString:solution.variable.name component:&k];
                        var = [self getVariableFrom:solution.mesh.variables model:model name:str onlySearch:&onlySearch maskName:NULL info:&found];
                        varContainers = var.getContainers;
                        varContainers->ComponentPrevValues = malloc ( (variableContainers->size1PrevValues/solution.variable.dofs) * sizeof ( double ** ));
                        for (i=0; i<(variableContainers->size1PrevValues/solution.variable.dofs); i++) {
                            varContainers->ComponentPrevValues[i] = malloc ( variableContainers->size2PrevValues * sizeof ( double * ));
                        }
                        l = 0;
                        for (i=k-1; i<variableContainers->size1PrevValues; i+=solution.variable.dofs) {
                            for (j=0; j<variableContainers->size2PrevValues; j++) {
                                varContainers->ComponentPrevValues[l][j] = &variableContainers->PrevValues[i][j];
                            }
                            l++;
                        }
                        varContainers->size1ComponentPrevValues = (variableContainers->size1PrevValues/solution.variable.dofs);
                        varContainers->size2ComponentPrevValues = variableContainers->size2PrevValues;
                    }
                }
            }
        }
        
        if ([solution.solutionInfo[@"calculate velocity"] boolValue] == YES || [solution.solutionInfo[@"nonlinear calculate velocity"] boolValue] == YES) {
            if (solution.timeOrder < 1) {
                fprintf(stdout, "FEMUtilities:addEquationToSolution: velocity computation implemented only for time-dependent equations.\n");
            } else if (solution.timeOrder == 1) {
                str = [NSMutableString stringWithString:solution.variable.name];
                [str appendString:@" velocity"];
                bufferContainers = allocateVariableContainer();
                bufferContainers->Perm = variableContainers->Perm;
                bufferContainers->sizePerm = variableContainers->sizePerm;
                dofs = solution.variable.dofs;
                [self addVectorTo:solution.mesh.variables mesh:solution.mesh solution:solution name:str dofs:&dofs container:bufferContainers ifOutput:NULL ifSecondary:NULL global:NULL initValue:NULL];
                free(bufferContainers);
            } else if (solution.timeOrder >= 2) {
                str = [NSMutableString stringWithString:solution.variable.name];
                [str appendString:@" velocity"];
                bufferContainers = allocateVariableContainer();
                bufferContainers->Perm = variableContainers->Perm;
                bufferContainers->sizePerm = variableContainers->sizePerm;
                bufferContainers->SecondaryToValues = malloc ( variableContainers->size1PrevValues * sizeof ( double * ));
                bufferContainers->sizeSecondaryToValues = variableContainers->size1PrevValues;
                for (i=0; i<variableContainers->size1PrevValues; i++) {
                    bufferContainers->SecondaryToValues[i] = &variableContainers->PrevValues[i][0];
                }
                dofs = solution.variable.dofs;
                secondary = YES;
                [self addVectorTo:solution.mesh.variables mesh:solution.mesh solution:solution name:str dofs:&dofs container:bufferContainers ifOutput:NULL ifSecondary:&secondary global:NULL initValue:NULL];
                free(bufferContainers);
            }
        }
        
        if ([solution.solutionInfo[@"calculate acceleration"] boolValue] == YES) {
            if (solution.timeOrder == 1) {
                fprintf(stdout, "FEMUtilities:addEquationToSolution: acceleration computation implemented only for 2nd order time equations\n");
            } else if (solution.timeOrder >= 2) {
                str = [NSMutableString stringWithString:solution.variable.name];
                [str appendString:@" acceleration"];
                bufferContainers = allocateVariableContainer();
                bufferContainers->Perm = variableContainers->Perm;
                bufferContainers->sizePerm = variableContainers->sizePerm;
                bufferContainers->SecondaryToValues = malloc ( variableContainers->size1PrevValues * sizeof ( double * ));
                bufferContainers->sizeSecondaryToValues = variableContainers->size1PrevValues;
                for (i=0; i<variableContainers->size1PrevValues; i++) {
                    bufferContainers->SecondaryToValues[i] = &variableContainers->PrevValues[i][1];
                }
                dofs = solution.variable.dofs;
                secondary = YES;
                [self addVectorTo:solution.mesh.variables mesh:solution.mesh solution:solution name:str dofs:&dofs container:bufferContainers ifOutput:NULL ifSecondary:&secondary global:NULL initValue:NULL];
                free(bufferContainers);
            }
        }
    } else {
        solution.timeOrder = 0;
        
        if ([solution.solutionInfo[@"calculate derivative"] boolValue] == YES) {
            str = [NSMutableString stringWithString:solution.variable.name];
            [str appendString:@" derivative"];
            bufferContainers = allocateVariableContainer();
            bufferContainers->Perm = variableContainers->Perm;
            bufferContainers->sizePerm = variableContainers->sizePerm;
            dofs = solution.variable.dofs;
            [self addVectorTo:solution.mesh.variables mesh:solution.mesh solution:solution name:str dofs:&dofs container:bufferContainers ifOutput:NULL ifSecondary:NULL global:NULL initValue:NULL];
            free(bufferContainers);
        }
        
        if (eigenAnal == YES) {
            if (solution.solutionInfo[@"eigen system complex"] != nil) {
                complexFlag = [solution.solutionInfo[@"eigen system complex"] boolValue];
            } else complexFlag = NO;
            
            if (solution.solutionInfo[@"eigen system values"] != nil) {
                n = [solution.solutionInfo[@"eigen system values"] intValue];
                if (n > 0) {
                    solution.nOfEigenValues = n;
                    if (variableContainers->EigenValues == NULL) {
                        variableContainers->EigenValues = cdoublevec(0, n-1);
                        variableContainers->sizeEigenValues = n;
                        variableContainers->EigenVectors = cdoublematrix(0, n-1, 0, variableContainers->sizeValues-1);
                        variableContainers->size1EigenVectors = n;
                        variableContainers->size2EigenVectors = variableContainers->sizeValues;
                        
                        memset( variableContainers->EigenValues, 0.0, n*sizeof(double complex) );
                        memset( *variableContainers->EigenVectors, 0.0, (n*variableContainers->sizeValues)*sizeof(double complex) );
                        
                        for (k=1; k<=solution.variable.dofs; k++) {
                            str = (NSMutableString *)[self appendNameFromString:solution.variable.name component:&k];
                            onlySearch = YES;
                            var = [self getVariableFrom:solution.mesh.variables model:model name:str onlySearch:&onlySearch maskName:NULL info:&found];
                            if (var != nil) {
                                varContainers = var.getContainers;
                                varContainers->EigenValues = variableContainers->EigenValues;
                                if (complexFlag == YES) {
                                    varContainers->EigenVectors = variableContainers->EigenVectors;
                                } else {
                                    varContainers->ComponentEigenVectors = malloc ( variableContainers->size1EigenVectors * sizeof ( double complex ** ));
                                    for (i=0; i<variableContainers->size1EigenVectors; i++) {
                                        varContainers->ComponentEigenVectors[i] = malloc ( (variableContainers->size2EigenVectors/solution.variable.dofs) * sizeof ( double complex * ));
                                    }
                                    l = 0;
                                    for (i=0; i<variableContainers->size1EigenVectors; i++) {
                                        for (j=k-1; j<variableContainers->size2EigenVectors; j+=solution.variable.dofs) {
                                            varContainers->ComponentEigenVectors[i][l] = &variableContainers->EigenVectors[i][j];
                                            l++;
                                        }
                                    }
                                    varContainers->size1ComponentEigenVectors = variableContainers->size1EigenVectors;
                                    varContainers->size2ComponentEigenVectors = (variableContainers->size2EigenVectors/solution.variable.dofs);
                                }
                            }
                        }
                    }
                    matContainers = solution.matrix.getContainers;
                    matContainers->MassValues = doublevec(0, matContainers->sizeValues-1);
                    matContainers->sizeMassValues = matContainers->sizeValues;
                    memset( matContainers->MassValues, 0.0, matContainers->sizeMassValues*sizeof(double) );
                }
            }
        } else if (harmonicAnal == YES) {
            n = 0;
            if (solution.solutionInfo[@"harmonic system values"] != nil) {
                n = [solution.solutionInfo[@"harmonic system values"] intValue];
            }
            if (n > 1) {
                listUtilities = [FEMListUtilities sharedListUtilities];
                found = [listUtilities listGetConstRealArray:model inArray:solution.valuesList forVariable:@"frequency" buffer:&freqv];
                if (found == YES) {
                    if (freqv.m < n) {
                        fatal("FEMUtilities:addEquationToSolution", "The solution option 'frequency' must be at least the same size as the Harmonic system values.");
                    }
                    if (freqv.matrix != NULL) {
                        free_dmatrix(freqv.matrix, 0, freqv.m-1, 0, freqv.n-1);
                        freqv.matrix = NULL;
                    }
                } else {
                    fatal("FEMUtilities:addEquationToSolution", "The solution option 'fequency' must be given for harmonic analysis.");
                }
            } else {
                n = 1;
            }
            
            solution.nOfEigenValues = n;
            if (variableContainers->EigenValues == NULL) {
                variableContainers->EigenValues = cdoublevec(0, n-1);
                variableContainers->sizeEigenValues = n;
                variableContainers->EigenVectors = cdoublematrix(0, n-1, 0, variableContainers->sizeValues-1);
                variableContainers->size1EigenVectors = n;
                variableContainers->size2EigenVectors = variableContainers->sizeValues;
                
                memset( variableContainers->EigenValues, 0.0, n*sizeof(double complex) );
                memset( *variableContainers->EigenVectors, 0.0, (n*variableContainers->sizeValues)*sizeof(double complex) );
                
                for (k=1; k<=solution.variable.dofs; k++) {
                    str = (NSMutableString *)[self appendNameFromString:solution.variable.name component:&k];
                    onlySearch = YES;
                    var = [self getVariableFrom:solution.mesh.variables model:model name:str onlySearch:&onlySearch maskName:NULL info:&found];
                    if (var != nil) {
                        varContainers = var.getContainers;
                        varContainers->EigenValues = variableContainers->EigenValues;
                        varContainers->ComponentEigenVectors = malloc ( variableContainers->size1EigenVectors * sizeof ( double complex ** ));
                        for (i=0; i<variableContainers->size1EigenVectors; i++) {
                            varContainers->ComponentEigenVectors[i] = malloc ( (variableContainers->size2EigenVectors/solution.variable.dofs) * sizeof ( double complex * ));
                        }
                        l = 0;
                        for (i=0; i<variableContainers->size1EigenVectors; i++) {
                            for (j=k-1; j<variableContainers->size2EigenVectors; j+=solution.variable.dofs) {
                                varContainers->ComponentEigenVectors[i][l] = &variableContainers->EigenVectors[i][j];
                                l++;
                            }
                        }
                        varContainers->size1ComponentEigenVectors = variableContainers->size1EigenVectors;
                        varContainers->size2ComponentEigenVectors = (variableContainers->size2EigenVectors/solution.variable.dofs);
                    }
                }
            }
            matContainers = solution.matrix.getContainers;
            matContainers->MassValues = doublevec(0, matContainers->sizeValues-1);
            matContainers->sizeMassValues = matContainers->sizeValues;
            memset( matContainers->MassValues, 0.0, matContainers->sizeMassValues*sizeof(double) );
        }
    }
    
    if (solution.matrix != nil) {
        if (solution.solutionInfo[@"linear system symmetric"] != nil) {
            solution.matrix.symmetric = [solution.solutionInfo[@"linear system symmetric"] boolValue];
        }
        
        if (solution.solutionInfo[@"lumped mass matrix"] != nil) {
            solution.matrix.lumped = [solution.solutionInfo[@"lumped mass matrix"] boolValue];
        }
        
        multigridActive = ([solution.solutionInfo[@"linear system solver"] isEqualToString:@"multigrid"] == YES
                           || [solution.solutionInfo[@"linear system preconditioning"] isEqualToString:@"multigrid"] == YES) ? YES : NO;
        
        // Check for multigrid solver
        if (multigridActive == YES) {
            meshUtilities = [[FEMMeshUtils alloc] init];
            
            // Multigrid may be either solver or preconditioner, is it solver?
            if (solution.solutionInfo[@"linear system solver"] != nil) {
                solution.multigridSolution = ([solution.solutionInfo[@"linear system solver"] isEqualToString:@"multigrid"] == YES) ? YES : NO;
            }
            
            // There are four different methods: algebraic, cluster, p and geometric
            if (solution.solutionInfo[@"mg method"] != nil) {
                methodName = solution.solutionInfo[@"mg method"];
                mgAlgebraic = ([methodName isEqualToString:@"algebraic"] == YES || [methodName isEqualToString:@"cluster"] == YES
                               || [methodName isEqualToString:@"p"] == YES) ? YES : NO;
            } else {
                mgAlgebraic = ([solution.solutionInfo[@"mg algebraic"] boolValue] == YES || [solution.solutionInfo[@"mg cluster"] boolValue] == YES
                               || [solution.solutionInfo[@"mg pelement"] boolValue] == YES) ? YES : NO;
            }
            
            if (solution.solutionInfo[@"mg levels"] != nil) {
                mgLevels = [solution.solutionInfo[@"mg levels"] intValue];
                if (mgLevels < 1) mgLevels = 1;
            } else {
                if (solution.solutionInfo[@"multigrid levels"] != nil) {
                    mgLevels = [solution.solutionInfo[@"multigrid levels"] intValue];
                    if (mgLevels < 1) mgLevels = 1;
                } else {
                    if (mgAlgebraic == YES) {
                        mgLevels = 10;
                    } else {
                        fatal("FEMUtilities:addEquationToSolution", "'mg levels' must be defined for geometric multigrid.");
                    }
                }
            }
            solution.multiGridTotal = mgLevels;
            
            // In case of geometric multigrid, make the hierarchical meshes
            if (mgAlgebraic == NO) {
                // Check if h/2 splitting of mesh requested
                if (solution.solutionInfo[@"mg equal split"] != nil) {
                    solution.multigridEqualPlit = [solution.solutionInfo[@"mg equal split"] boolValue];
                }
                
                if (solution.multigridEqualPlit == YES) {
                    // TODO: add support for parallel run, here initialization of a parallel matrix
                    solution.multigridLevel = 1;
                    
                    while (solution.multigridLevel < solution.multiGridTotal) {
                        if (solution.mesh.child != nil) {
                            newMesh = solution.mesh.child;
                            
                            oldMesh = solution.mesh;
                            oldMatrix = solution.matrix;
                            [meshUtilities updateMesh:newMesh inSolution:solution model:model];
                            solution.mesh.changed = NO;
                        } else {
                            newMesh = [meshUtilities splitMeshEqual:solution.mesh model:model nodal:NULL sizeNodal:NULL];
                            newMesh.next = model.meshes;
                            NSMutableArray *meshArray = [[NSMutableArray alloc] init];
                            [meshArray addObject:newMesh];
                            model.meshes = meshArray;
                            
                            oldMesh = solution.mesh;
                            oldMatrix = solution.matrix;
                            [meshUtilities updateMesh:newMesh inSolution:solution model:model];
                            solution.mesh.changed = NO;
                            
                            newMesh.parent = oldMesh;
                            oldMesh.child = newMesh;
                            [newMesh.name setString:oldMesh.name];
                        }
                        newMesh.outputActive = YES;
                        oldMesh.outputActive = NO;
                        
                        newMatrix = solution.matrix;
                        newMatrix.parent = oldMatrix;
                        oldMatrix.child = newMatrix;
                        // TODO: add support for parallel run, here initialization of a parallel matrix
                        solution.multigridLevel++;
                    }
                } else {
                    // TODO: add support for parallel run, here initialization of a parallel matrix
                    oldMesh = solution.mesh;
                    var = solution.variable;
                    saveMatrix = solution.matrix;
                    
                    solution.multigridLevel = 1;
                    while (solution.multigridLevel < solution.multiGridTotal) {
                        if (solution.mesh.parent != nil) {
                            newMesh = solution.mesh.parent;
                            oldMatrix = solution.matrix;
                            [meshUtilities updateMesh:newMesh inSolution:solution model:model];
                            solution.mesh.changed = NO;
                            newMatrix = solution.matrix;
                            newMatrix.child = oldMatrix;
                            oldMatrix.parent = newMatrix;
                            // TODO: add support for parallel run, here initialization of a parallel matrix
                        }
                        solution.multigridLevel++;
                    }
                    solution.mesh = oldMesh;
                    solution.variable = var;
                    solution.matrix = saveMatrix;
                    [meshUtilities setCurrentMesh:solution.mesh inModel:model];
                }
            }
            [meshUtilities setStabilizationParametersInMesh:solution.mesh model:model];
        }
        
        // Set the default verbosity of the iterative solvers accordingly with the global verbosity
        if (solution.solutionInfo[@"linear system residual output"] == nil) {
            core = [FEMCore sharedCore];
            k = 1;
            if ([core.outputLevelMask[4] boolValue] == NO) {
                k = 0;
            } else if ([core.outputLevelMask[5] boolValue] == NO) {
                k = 10;
            }
            if (k != 1) {
                [solution.solutionInfo setObject:@(k) forKey:@"linear system residual output"];
            }
        }
    }
}

-(BOOL)isFileNameQualified:(NSString * __nonnull)file {
    
    NSRange range;
    
    if (file == nil) return NO;
    
    range = [file rangeOfString:@":"];
    return (range.location != NSNotFound || [file characterAtIndex:0] == '/' || [file characterAtIndex:0] == _backSlash) ? YES : NO;
}

/****************************************************************
    Given the filename0 (and suffix0) find the 1st free filename
    that does not exist in the current working directory
****************************************************************/
-(NSMutableString * __nullable)nextFreeFileName:(NSString * __nonnull)fileName0 suffix:(NSString * __nullable)suffix0 lastExisting:(BOOL * __nullable)lastExisting {
    
    int no;
    NSRange range;
    NSString *prefix, *suffix;
    NSMutableString *str, *fileName, *prevFileName;
    NSFileManager *fileManager;
    
    if (fileName0 == nil) return nil;
    
    range = [fileName0 rangeOfString:@"." options:NSBackwardsSearch];
    if (range.location != NSNotFound) {
        prefix = [fileName0 substringToIndex:range.location];
        suffix = [fileName0 substringFromIndex:range.location];
    } else {
        prefix = [NSString stringWithString:prefix];
        if (suffix0 != nil) {
            str = [NSMutableString stringWithString:@"."];
            [str appendString:suffix0];
            suffix = [NSString stringWithString:str];
        } else {
            suffix = @".dat";
        }
    }
    
    prevFileName = [NSMutableString stringWithString:@""];
    fileName = [NSMutableString stringWithString:@""];
    fileManager = [NSFileManager defaultManager];
    
    for (no=1; no<=9999; no++) {
        if (no > 0) [prevFileName setString:fileName];
        [fileName setString:prefix];
        [fileName appendString:[NSString stringWithFormat:@"%d",no]];
        [fileName appendString:suffix];
        if ([fileManager fileExistsAtPath:fileName] == NO) break;
    }
    
    if (lastExisting != NULL) {
        if (*lastExisting == YES) [fileName setString:prevFileName];
    }
    
    return fileName;
}

-(NSBundle * __nullable)loadBundle:(NSString * __nonnull)bundleName useApplicationSupportPath:(BOOL * __nullable)useApplicationSupportPath {
    
    NSBundle *currentBundle;
    NSArray *librarySearchPaths;
    NSMutableArray *bundleSearchPaths = [NSMutableArray array];
    BOOL found = NO, searchApplicationSupportPath;
    
    if (useApplicationSupportPath != NULL) {
        searchApplicationSupportPath = *useApplicationSupportPath;
    } else searchApplicationSupportPath = NO;
    
    if (searchApplicationSupportPath == YES) {
        // NSSearchPathForDirectoriesInDomains will return:
        //      /Users/seddikhakime/Library
        //      /Library
        //      /Network/Library
        librarySearchPaths = NSSearchPathForDirectoriesInDomains(NSLibraryDirectory, NSAllDomainsMask - NSSystemDomainMask, YES);
        for (NSString *currPath in librarySearchPaths) {
            [bundleSearchPaths addObject:[currPath stringByAppendingPathComponent:self.appSupportSubpath]];
        }
        
        NSString *builtInPlugInsPath = [[NSBundle mainBundle] builtInPlugInsPath];
        if (builtInPlugInsPath != nil) [bundleSearchPaths addObject:builtInPlugInsPath];
        
        NSString *searchedSolutionBundle = [bundleName stringByAppendingPathExtension:self.ext];
        
        for (NSString *currPath in bundleSearchPaths) {
            NSDirectoryEnumerator *bundleEnum;
            bundleEnum = [[NSFileManager defaultManager] enumeratorAtPath:currPath];
            if (bundleEnum) {
                for (NSString *currBundlePath in bundleEnum) {
                    if ([currBundlePath isEqualToString:searchedSolutionBundle] == YES) {
                        currentBundle = [NSBundle bundleWithPath:[currPath stringByAppendingPathComponent:currBundlePath]];
                        found = YES;
                        break;
                    }
                }
            }
            if (found == YES) break;
        }
    } else {
        currentBundle = [NSBundle bundleWithPath:[bundleName stringByAppendingPathExtension:self.ext]];
    }
    
    return currentBundle;
}

-(BOOL)plugInClassIsValid:(Class __nonnull)plugInClass {
    
    if ([plugInClass conformsToProtocol:@protocol(SainoSolutionsComputer)] == YES) {
        if ([plugInClass instancesRespondToSelector:@selector(solutionComputer:model:timeStep:transientSimulation:)] == YES &&
                [plugInClass instancesRespondToSelector:@selector(deallocation:)] == YES) {
            return YES;
        }
    }
    return NO;
}

-(void)generateColor:(RGBColors * __nonnull)color {
    
    int red = arc4random() % 256;
    int green = arc4random() % 256;
    int blue = arc4random() % 256;
    
    color->red = (float)red / 255.0;
    color->green = (float)green / 255.0;
    color->blue = (float)blue / 255.0;
}

@end
