//
//  FEMUtilities.m
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMUtilities.h"
#import "FEMKernel.h"
#import "FEMEquation.h"
#import "FEMPElementMaps.h"
#import "FEMProjector.h"
#import "FEMElementUtils.h"
#import "FEMMeshUtils.h"
#import "FEMMatrixCRS.h"
#import "FEMMatrixBand.h"
#import "FEMListUtilities.h"
#import "FEMInterpolation.h"
#import "FEMElementDescription.h"
#import "Utils.h"

@interface FEMUtilities ()
-(void)FEMUtils_applyProjector:(NSMutableArray *)variables model:(FEMModel *)aModel fromMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh projector:(FEMProjector *)projector;
-(void)FEMUtils_interpolateQMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh oldVariables:(NSMutableArray *)oldVar newVariables:(NSMutableArray *)newVar model:(FEMModel *)aModel quadrantTree:(BOOL *)useQuandrant withProjector:(BOOL)withProj projector:(FEMProjector *)projector mask:(NSString *)maskName nodesPresent:(BOOL *)foundNodes;
-(void)FEMUtils_interpolateMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh oldVariables:(NSMutableArray *)oldVar newVariables:(NSMutableArray *)newVar model:(FEMModel *)aModel quadrantTree:(BOOL *)useQuandrant withProjector:(BOOL)withProj projector:(FEMProjector *)projector mask:(NSString *)maskName;
@end

@implementation FEMUtilities {
    char _backSlash;
}

#pragma Private methods

-(void)FEMUtils_applyProjector:(NSMutableArray *)variables model:(FEMModel *)aModel fromMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh projector:(FEMProjector *)projector {
    
    int i, j;
    variableArraysContainer *varContainers = NULL, *oldContainers = NULL, *newContainers = NULL;
    FEMVariable *oldSol, *newSol;
    FEMMatrixCRS *crsMatrix;
    double *bf1, *bf2;
    BOOL only, stat;
    
    crsMatrix = [[FEMMatrixCRS alloc] init];
    
    for (FEMVariable *variable in variables) {
        varContainers = variable.getContainers;
        
        if (varContainers->sizeValues == variable.dofs) continue;
        if (variable.secondary == YES) continue;
        
        if (variable.dofs == 1 && [variable.name isEqualToString:@"coordinate"] == YES) {
            only = YES;
            oldSol = [self getVariableFrom:oldMesh.variables model:aModel name:variable.name onlySearch:&only maskName:NULL info:&stat];
            newSol = [self getVariableFrom:newMesh.variables model:aModel name:variable.name onlySearch:&only maskName:NULL info:&stat];
            
            if (newSol == nil) continue;
            
            oldContainers = oldSol.getContainers;
            newContainers = newSol.getContainers;
            
            [crsMatrix applyProjector:projector.matrix values:oldContainers->Values permutation:oldContainers->Perm values:newContainers->Values permutation:newContainers->Perm transpose:NULL];
            
            
            if (oldContainers->PrevValues != NULL) {
                bf1 = doublevec(0, oldContainers->size1PrevValues-1);
                bf2 = doublevec(0, newContainers->size1PrevValues-1);
                for (i=0; i<oldContainers->size2PrevValues; i++) {
                    for (j=0; j<oldContainers->size1PrevValues; j++) {
                        bf1[j] = oldContainers->PrevValues[j][i];
                    }
                    for (j=0; j<newContainers->size1PrevValues; j++) {
                        bf2[j] = newContainers->PrevValues[j][i];
                    }
                    [crsMatrix applyProjector:projector.matrix values:bf1 permutation:oldContainers->Perm values:bf2 permutation:newContainers->Perm transpose:NULL];
                    for (j=0; j<oldContainers->size1PrevValues; j++) {
                        newContainers->PrevValues[j][i] = bf2[j];
                    }
                }
                free_dvector(bf1, 0, oldContainers->size1PrevValues-1);
                free_dvector(bf2, 0, newContainers->size1PrevValues-1);
            }
        }
    }
}

/**************************************************************************************************
 
 Interpolates velues of all variables from a mesh associated with the old model to the mesh
 of the new model.
 
 (FEMMesh *)oldMesh         ->  Old mesh class
 (FEMMesh *)newMesh         ->  New mesh class
 (FEMVariable *)oldVar      ->  Old model variable class
 (FEMVariable *)newVar      ->  New model variable class
 (BOOL *)useQuandrant       ->  If true use the RootQuadrant of the old mesh in interpolation
 (FEMProjector *)projector  ->  Use projector between meshes for interpolation if available
 (NSString *)maskName       ->  Mask the old variable set by the given mask name when trying
 to define the interpolation
 (BOOL *)foundNodes         ->  List of nodes where the interpolation was a success
 
 ***************************************************************************************************/
-(void)FEMUtils_interpolateQMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh oldVariables:(NSMutableArray *)oldVar newVariables:(NSMutableArray *)newVar model:(FEMModel *)aModel quadrantTree:(BOOL *)useQuandrant withProjector:(BOOL)withProj projector:(FEMProjector *)projector mask:(NSString *)maskName nodesPresent:(BOOL *)foundNodes {
    
    int dim, epsTries;
    int i, j, k, l, n, bfId, qTreeFails, totFails;
    Nodes_t *elementNodes;
    double point[3], localCoordinates[3];
    FEMVariable *oldSol, *newSol;
    double *newValues, *elementValues;
    Element_t *element, *oldElements;
    Nodes_t *oldNodes, *newNodes;
    Quadrant_t *leafQuadrant, *rootQuadrant, *oldMeshQuadrant;
    double *basis;
    double boundingBox[6], *valuesPtr, u, v, w;
    BOOL useQTree, tryQtree, found, searchOnly;
    int *rows, *cols;
    int *rInd;
    BOOL maskExists, all;
    double eps1 = 0.1, eps2, epsGlobal, epsLocal, max, numericEps = 1.0e-12;
    double *values, *localU, *localV, *localW;
    FEMListUtilities *listUtilities;
    FEMInterpolation *interpolation;
    FEMElementDescription *elementDescription;
    variableArraysContainer *varContainers = NULL, *newVarContainers = NULL, *oldVarContainers = NULL;
    matrixArraysContainer *matContainers = NULL, *tmatContainers = NULL;
    
    typedef struct Epntr_t {
        Element_t *element;
    } Epntr_t;
    Epntr_t *elementPointers;
    
    listUtilities = [[FEMListUtilities alloc] init];
    interpolation = [[FEMInterpolation alloc] init];
    elementDescription = [[FEMElementDescription alloc] init];
    
    // If projector argument given, search for existing projector matrix or
    // generate new projector if not already there
    if (withProj == YES) {
        for (FEMProjector *projector in newMesh.projectors) {
            if ([projector.mesh isEqualTo:oldMesh]) {
                if (oldVar != nil) [self FEMUtils_applyProjector:oldVar model:aModel fromMesh:oldMesh toMesh:newMesh projector:projector];
                break;
            }
        }
        
        n = newMesh.numberOfNodes;
        localU = doublevec(0, n-1);
        localV = doublevec(0, n-1);
        localW = doublevec(0, n-1);
        elementPointers = (Epntr_t*) malloc( sizeof(Epntr_t) * n );
        for (i=0; i<n; i++) {
            elementPointers[i].element = NULL;
        }
    }
    
    // Check if using spatial division hierarchy for the search
    rootQuadrant = oldMesh.getQuadrant;
    dim = aModel.dimension;
    
    if (useQuandrant == NULL) {
        useQTree = YES;
    } else {
        useQTree = *useQuandrant;
    }
    
    if (useQTree == YES) {
        if (rootQuadrant == NULL) {
            oldNodes = oldMesh.getNodes;
            
            boundingBox[0] = min_array(oldNodes->x, oldMesh.numberOfNodes);
            boundingBox[1] = min_array(oldNodes->y, oldMesh.numberOfNodes);
            boundingBox[2] = min_array(oldNodes->z, oldMesh.numberOfNodes);
            boundingBox[3] = max_array(oldNodes->x, oldMesh.numberOfNodes);
            boundingBox[4] = max_array(oldNodes->y, oldMesh.numberOfNodes);
            boundingBox[5] = max_array(oldNodes->z, oldMesh.numberOfNodes);
            
            max = -HUGE_VAL;
            k = 3;
            valuesPtr = &boundingBox[0];
            for (i=0; i<3; i++) {
                if ( (*(valuesPtr + k) - *(valuesPtr + i)) > max ) {
                    max = *(valuesPtr + k) - *(valuesPtr + i);
                }
                k++;
            }
            eps2 = 0.1 * max;
            for (i=0; i<3; i++) {
                *(valuesPtr + i) = *(valuesPtr + i) - eps2;
            }
            for (i=3; i<6; i++) {
               *(valuesPtr + i) = *(valuesPtr + i) + eps2;
            }
            
            oldMeshQuadrant = oldMesh.getQuadrant;
            [interpolation buildQuadrantTreeForMesh:oldMesh model:aModel boundingBox:boundingBox rootQuadrant:oldMeshQuadrant];
            rootQuadrant = oldMeshQuadrant;
        }
    }
    
    // If we use mask
    maskExists = (maskName != nil) ? YES : NO;
    
    n = oldMesh.maxElementNodes;
    elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(elementNodes);
    elementNodes->x = doublevec(0, n-1);
    elementNodes->y = doublevec(0, n-1);
    elementNodes->z = doublevec(0, n-1);
    
    elementValues = doublevec(0, n-1);
    newValues = doublevec(0, newMesh.numberOfNodes-1);
    
    epsGlobal = [listUtilities listGetConstReal:aModel inArray:aModel.simulation.valuesList forVariable:@"interpolation global epsilon" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) epsGlobal = 2.0e-10;
    epsLocal = [listUtilities listGetConstReal:aModel inArray:aModel.simulation.valuesList forVariable:@"interpolation local epsilon" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) epsLocal = 1.0e-10;
    epsTries = [listUtilities listGetConstReal:aModel inArray:aModel.simulation.valuesList forVariable:@"interpolation maximum iterations" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) epsTries = 12;
    
    qTreeFails = 0;
    totFails = 0;
    
    newNodes = newMesh.getNodes;
    oldElements = oldMesh.getElements;
    
    searchOnly = YES;
    
    // Loop over all nodes in the new mesh
    for (i=0; i<newMesh.numberOfNodes; i++) {
        
        element = NULL;
        point[0] = newNodes->x[i];
        point[1] = newNodes->y[i];
        point[2] = newNodes->z[i];
        
        // Find in which old mesh bulk element the point belongs to
        found = NO;
        tryQtree = (rootQuadrant != NULL && useQTree == YES) ? YES : NO;
        
        if (tryQtree) {
            // Find the last existing quadrant that the point belongs to
            leafQuadrant = NULL;
            [interpolation findLeafElementsForPoint:point dimension:dim rootQuadrant:rootQuadrant leafQuadrant:leafQuadrant];
            
            if (leafQuadrant != NULL) {
                // Go through the bulk elements in the last child quadrant only.
                // Try to find matching element with progressively sloppier tests.
                // All at most 100% of slack
                eps1 = epsGlobal;
                eps2 = epsLocal;
                
                for (j=0; j<epsTries; j++) {
                    for (k=0; k<leafQuadrant->nElementsInQuadrant; k++) {
                        
                        element = &oldElements[leafQuadrant->elements[k]];
                        
                        if (maskExists == YES) {
                            if ((aModel.bodies)[oldElements[leafQuadrant->elements[k]].BodyID-1][@"body force"] == nil) continue;
                            bfId = [(aModel.bodies)[oldElements[leafQuadrant->elements[k]].BodyID][@"body force"] intValue];
                            if ([listUtilities listCheckPresentVariable:maskName inArray:(aModel.bodyForces)[bfId-1]] == NO) continue;
                        }
                        
                        n = oldElements[leafQuadrant->elements[k]].Type.NumberOfNodes;
                        for (l=0; l<n; l++) {
                            elementNodes->x[l] = oldNodes->x[oldElements[leafQuadrant->elements[k]].NodeIndexes[l]];
                            elementNodes->y[l] = oldNodes->y[oldElements[leafQuadrant->elements[k]].NodeIndexes[l]];
                            elementNodes->z[l] = oldNodes->z[oldElements[leafQuadrant->elements[k]].NodeIndexes[l]];
                        }
                        found = [interpolation isPointInElement:&oldElements[leafQuadrant->elements[k]] elementNodes:elementNodes point:point localCoordinates:localCoordinates globalEpsilon:&eps1 localEpsilon:&eps2 numericEpsilon:&numericEps globalDistance:NULL localDistance:NULL model:aModel];
                        if (found == YES) break;
                    }
                    if (found == YES) break;
                    
                    eps1 = 10.0 * eps1;
                    eps2 = 10.0 * eps2;
                    if (eps1 > 1.0) break;
                }
            }
        }
        
        if (tryQtree == NO || (found == NO && foundNodes == NULL)) {
            
            // Go through all old mesh bulk elements
            for (k=0; k<oldMesh.numberOfBulkElements; k++) {
                element = &oldElements[k];
                n = oldElements[k].Type.NumberOfNodes;
                
                for (l=0; l<n; l++) {
                    elementNodes->x[l] = oldNodes->x[oldElements[k].NodeIndexes[l]];
                    elementNodes->y[l] = oldNodes->y[oldElements[k].NodeIndexes[l]];
                    elementNodes->z[l] = oldNodes->z[oldElements[k].NodeIndexes[l]];
                }
                
                found = [interpolation isPointInElement:&oldElements[k] elementNodes:elementNodes point:point localCoordinates:localCoordinates globalEpsilon:NULL localEpsilon:NULL numericEpsilon:NULL globalDistance:NULL localDistance:NULL model:aModel];
                if (found == YES) {
                    if (tryQtree == YES) qTreeFails++;
                    break;
                }
            }
        }
        
        if (found == NO) {
            element = NULL;
            if (foundNodes == NULL) {
                NSLog(@"FEMUtils_interpolateQMesh: point %d was not found in any of the elements!\n", i);
                totFails++;
            }
            continue;
        }
        if (foundNodes != NULL) foundNodes[i] = YES;
        
        // Found element in oldModel
        if (withProj == YES) {
            elementPointers[i].element = element;
            localU[i] = localCoordinates[0];
            localV[i] = localCoordinates[1];
            localW[i] = localCoordinates[2];
        }
        
        if (oldVar == nil || withProj == YES) continue;
        
        // Go through all variables to be interpolated
        for (FEMVariable *variable in oldVar) {
            
            varContainers = variable.getContainers;
            if (varContainers->sizeValues == variable.dofs) continue;
            
            if (variable.secondary == YES) continue;
            
            newSol = nil;
            oldSol = nil;
            if (variable.dofs == 1 && [variable.name isEqualToString:@"coordinate"] == NO) {
                
                // Interpolate variable at point in element
                newSol = [self getVariableFrom:newMesh.variables model:aModel name:variable.name onlySearch:&searchOnly maskName:NULL info:&found];
                if (newSol == nil) continue;
                
                oldSol = [self getVariableFrom:oldMesh.variables model:aModel name:variable.name onlySearch:&searchOnly maskName:NULL info:&found];
                
                newVarContainers = newSol.getContainers;
                oldVarContainers = oldSol.getContainers;
                
                // Check that the node was found in the old mesh
                if (element != NULL) {
                    //Check for rounding errors
                    all = YES;
                    for (k=0; k<element->Type.NumberOfNodes; k++) {
                        if (oldVarContainers->Perm[element->NodeIndexes[k]] > 0) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) {
                        if (newVarContainers->Perm[i] != 0.0) {
                            for (j=0; j<n; j++) {
                                elementValues[j] = oldVarContainers->Values[oldVarContainers->Perm[element->NodeIndexes[j]]];
                            }
                            newVarContainers->Values[newVarContainers->Perm[i]] = [elementDescription interpolateInElement:element nodalValues:elementValues evaluatedAt:localCoordinates[0] andAt:localCoordinates[1] andAt:localCoordinates[2] withBasis:NULL];
                            
                            if (oldVarContainers->PrevValues != NULL) {
                                for (j=0; j<oldVarContainers->size2PrevValues; j++) {
                                    for (k=0; k<n; k++) {
                                        elementValues[k] = oldVarContainers->PrevValues[oldVarContainers->Perm[element->NodeIndexes[k]]][j];
                                    }
                                    newVarContainers->PrevValues[newVarContainers->Perm[i]][j] = [elementDescription interpolateInElement:element nodalValues:elementValues evaluatedAt:localCoordinates[0] andAt:localCoordinates[1] andAt:localCoordinates[2] withBasis:NULL];
                                }
                            }
                        }
                    }
                } else {
                    // TODO: Not sure what this is fore!
                    if (newVarContainers->Perm[i] != 0) newValues[newVarContainers->Perm[i]] = 0.0;
                }
            }
        }
    }
    
    if (foundNodes == NULL) {
        if (qTreeFails > 0) {
            NSLog(@"FEMUtils_interpolateQMesh: number of points not found in quadtree: %d", qTreeFails);
            if (totFails == 0) {
                NSLog(@"FEMUtils_interpolateQMesh: all nodes still found by N^2 dummy search\n");
            }
        }
        if (totFails > 0) {
            NSLog(@"FEMUtils_interpolateQMesh: number of points not found: %d\n", totFails);
        }
    }
    
    // Construct mesh projector if required. Next time around will use
    // the existing projector to interpolate values
    if (withProj == YES) {
        n = newMesh.numberOfNodes;
        
        rows = intvec(0, (n+1)-1);
        rows[0] = 0;
        for (i=1; i<n+1; i++) {
            if (elementPointers[i-1].element != NULL) {
                rows[i] = rows[i-1] + elementPointers[i-1].element->Type.NumberOfNodes;
            } else {
                rows[i] = rows[i-1]+1;
            }
        }
        
        cols = intvec(0, rows[n]-1);
        values = doublevec(0, rows[n]-1);
        memset( cols, 0, (rows[n]*sizeof(cols)) );
        memset( values, 0.0, (rows[n]*sizeof(values)) );

        projector = [[FEMProjector alloc] init];
        projector.matrix.numberOfRows = n;
        matContainers = projector.matrix.getContainers;
        matContainers->Rows = rows;
        matContainers->Cols = cols;
        matContainers->Values = values;
        
        projector.mesh = oldMesh;
        [newMesh.projectors addObject:projector];
        
        rInd = intvec(0, oldMesh.numberOfNodes-1);
        basis = doublevec(0, oldMesh.maxElementNodes-1);
        memset( rInd, 0, (oldMesh.numberOfNodes*sizeof(rInd)) );
        
        found = YES;
        for (i=0; i<n; i++) {
            element = elementPointers[i].element;
            if (element == NULL) {
                found = NO;
                cols[rows[i]] = 0;
                values[rows[i]] = 0.0;
                continue;
            }
            
            k = element->Type.NumberOfNodes;
            for (j=0; j<k; j++) {
                rInd[element->NodeIndexes[k]] = rInd[element->NodeIndexes[k]] + 1;
            }
            
            u = localU[i];
            v = localV[i];
            w = localW[i];
            
            for (j=0; j<k; j++) {
                basis[j] = 0.0;
            }
            for (j=0; j<k; j++) {
                l = rows[i] + j;
                cols[l] = element->NodeIndexes[j];
                basis[j] = 1.0;
                values[l] = [elementDescription interpolateInElement:element nodalValues:basis evaluatedAt:u andAt:v andAt:w withBasis:NULL];
                basis[j] = 0.0;
            }
        }
        free_dvector(basis, 0, oldMesh.maxElementNodes);
        free(elementPointers);
        free_dvector(localU, 0, newMesh.numberOfNodes-1);
        free_dvector(localV, 0, newMesh.numberOfNodes-1);
        free_dvector(localW, 0, newMesh.numberOfNodes-1);
        
        rows = NULL;
        cols = NULL;
        values = NULL;
        // Store also the transpose of the projector
        if (found == YES) {
            n = oldMesh.numberOfNodes;
            rows = intvec(0, (n+1)-1);
            rows[0] = 0;
             for (i=1; i<n+1; i++) {
                 rows[i] = rows[i-1] + rInd[i-1];
             }
            cols = intvec(0, rows[n]-1);
            values = doublevec(0, rows[n]-1);
            projector.tMatrix.numberOfRows = n;
            tmatContainers = projector.tMatrix.getContainers;
            tmatContainers->Rows = rows;
            tmatContainers->Cols = cols;
            tmatContainers->Values = values;
            
            memset( rInd, 0, (oldMesh.numberOfNodes*sizeof(rInd)) );
            for (i=0; i<projector.matrix.numberOfRows; i++) {
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    k = matContainers->Cols[j];
                    l = rows[l] + rInd[k];
                    rInd[k] = rInd[k] + 1;
                    cols[l] = i;
                    values[l] = matContainers->Values[j];
                }
            }
        }
        
        free_ivector(rInd, 0,  oldMesh.numberOfNodes-1);
        if (oldVar != nil) [self FEMUtils_applyProjector:oldVar model:aModel fromMesh:oldMesh toMesh:newMesh projector:projector];
    }
    
    [elementDescription deallocation];
    
    free_dvector(elementNodes->x, 0, oldMesh.maxElementNodes-1);
    free_dvector(elementNodes->y, 0, oldMesh.maxElementNodes-1);
    free_dvector(elementNodes->z, 0, oldMesh.maxElementNodes-1);
    free(elementNodes);
    
    free_dvector(elementValues, 0, oldMesh.maxElementNodes-1);
    free_dvector(newValues, 0,  newMesh.numberOfNodes-1);
}

-(void)FEMUtils_interpolateMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh oldVariables:(NSMutableArray *)oldVar newVariables:(NSMutableArray *)newVar model:(FEMModel *)aModel quadrantTree:(BOOL *)useQuandrant withProjector:(BOOL)withProj projector:(FEMProjector *)projector mask:(NSString *)maskName {
    
    [self FEMUtils_interpolateQMesh:oldMesh toMesh:newMesh oldVariables:oldVar newVariables:newVar model:aModel quadrantTree:useQuandrant withProjector:withProj projector:projector mask:maskName nodesPresent:NULL];
    
    // TODO: The rest of this method deals with a parallel mesh
}

#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _backSlash = (char)92;
    }
    return self;
}

-(FEMMatrix *)allocateMatrix {
    
    FEMMatrix *matrix;
    
    matrix = [[FEMMatrix alloc] init];
    
    return matrix;
}

-(void)zeroTheNumberOfRows:(int)n inMatrix:(FEMMatrix *)a {
    
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

-(void)setMatrixElement:(FEMMatrix *)a atIndex:(int)i andIndex:(int)j value:(double)value {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    if (a.format == MATRIX_CRS) {
        crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix setMatrixElementInMatrix:a atIndex:i andIndex:j value:value];
        
    } else if (a.format == MATRIX_LIST) {
        // TODO: implement the setMatrixElement method for list matrix.
        
    } else if (a.format == MATRIX_BAND || a.format == MATRIX_SBAND) {
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix setMatrixElementInMatrix:a atIndex:i andIndex:j value:value];
    }
}

/******************************************************************************************************
 
    Check if given element belongs to a body for which given equation should be solved

******************************************************************************************************/
-(BOOL)checkEquationForElement:(Element_t *)element model:(FEMModel *)aModel equation:(NSString *)str {
    
    int k, bodyID;
    BOOL flag, found;
    FEMListUtilities *listUtilities;
    FEMEquation *equation;
    
    listUtilities = [[FEMListUtilities alloc] init];
    
    flag = NO;
    bodyID = element->BodyID;
    if (bodyID > 0 && bodyID <= aModel.numberOfBodies) {
        k = [(aModel.bodies)[bodyID][@"equation"] intValue];
        if (k < 1 || k > aModel.numberOfEquations) errorfunct("checkEquationForElement", "Equation number not in range");
        equation = (aModel.equations)[k];
        if (k > 0) flag = [listUtilities listGetLogical:aModel inArray:equation.valuesList forVariable:str info:&found];
    }
    return flag;
}

-(int)initialPermutationInMesh:(FEMMesh *)aMesh model:(FEMModel *)aModel solution:(FEMSolution *)aSolution equation:(NSString *)str permutation:(int *)perm DGSolution:(BOOL *)dg globalBubbles:(BOOL *)gb {
    
    int i, j, t, n, e, edofs, fdofs, bdofs, ndofs;
    int indexes[128];
    int *defDofs, *edgeDofs, *faceDofs;
    int k;
    BOOL foundDG, DG, GB, radiation, any;
    Element_t *element, *elements, *edges, *faces;
    FEMPElementMaps *elementMaps;
    FEMListUtilities *listUtilities;
    solutionArraysContainer *solContainers = NULL;
    
    k = 0;
    edofs = aMesh.maxEdgeDofs;
    fdofs = aMesh.maxFaceDofs;
    bdofs = aMesh.maxBdofs;
    
    elementMaps = [[FEMPElementMaps alloc] init];
    listUtilities = [[FEMListUtilities alloc] init];
    
    elements = aMesh.getElements;
    edges = aMesh.getEdges;
    faces = aMesh.getFaces;
    
    GB = NO;
    if (gb != NULL) GB = *gb;
    
    DG = NO;
    if (dg != NULL) DG = *dg;
    
    foundDG = NO;
    if (DG == YES) {
        for (t=0; t<aMesh.numberOfEdges; t++) {
            n = 0;
            element = edges[t].BoundaryInfo->Left;
            if (element != NULL) {
                if ([self checkEquationForElement:element model:aModel equation:str] == YES) {
                    foundDG = (foundDG == YES || element->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<element->DGDOFs; j++) {
                        indexes[n] = element->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            element = edges[t].BoundaryInfo->Right;
            if (element != NULL) {
                if ([self checkEquationForElement:element model:aModel equation:str] == YES) {
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
        
        for (t=0; t<aMesh.numberOfFaces; t++) {
            n = 0;
            element = faces[t].BoundaryInfo->Left;
            if (element != NULL) {
                if ([self checkEquationForElement:element model:aModel equation:str] == YES) {
                    foundDG = (foundDG == YES || element->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<element->DGDOFs; j++) {
                        indexes[n] = element->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            element = faces[t].BoundaryInfo->Right;
            if (element != NULL) {
                if ([self checkEquationForElement:element model:aModel equation:str] == YES) {
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
    
    solContainers = aSolution.getContainers;
    defDofs = intvec(0, solContainers->size2DefDofs-1);
    
    any = NO;
    for (i=0; i<solContainers->size1DefDofs; i++) {
        if (solContainers->defDofs[i][5] >= 0) {
            any = YES;
            break;
        }
    }
    if (any == YES) {
        if (aMesh.numberOfEdges > 0) {
            edgeDofs = intvec(0, aMesh.numberOfEdges-1);
            memset( edgeDofs, 0, (aMesh.numberOfEdges*sizeof(edgeDofs)) );
        }
        
        if (aMesh.numberOfFaces > 0) {
            faceDofs = intvec(0, aMesh.numberOfFaces-1);
            memset( faceDofs, 0, (aMesh.numberOfFaces*sizeof(faceDofs)) );
        }
        
        n = aMesh.numberOfBulkElements + aMesh.numberOfBoundaryElements;
        t = 0;
        while (t < n) {
            while (t < n) {
                if ([self checkEquationForElement:&elements[t] model:aModel equation:str] == YES) break;
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
    
    n = aMesh.numberOfBulkElements + aMesh.numberOfBoundaryElements;
    t = 0;
    while (t < n) {
        while (t < n) {
            if ([self checkEquationForElement:&elements[t] model:aModel equation:str] == YES) break;
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
                    j = aMesh.numberOfNodes + edofs*elements[t].EdgeIndexes[i] + e;
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
                    j = aMesh.numberOfNodes + edofs*aMesh.numberOfEdges + fdofs*elements[t].FaceIndexes[i] + e;
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
                j = aMesh.numberOfNodes + edofs*aMesh.numberOfEdges + fdofs*aMesh.numberOfFaces + elements[t].BubbleIndexes[i];
                if (perm[j] < 0) {
                    perm[j] = k;
                    k++;
                }
            }
        }
        t++;
    }
    
    radiation = [(aSolution.solutionInfo)[@"radiation solver"] boolValue];
    if (radiation == YES || [str isEqualToString:@"heat equation"]) {
        t = aMesh.numberOfBulkElements;
        n = aMesh.numberOfBulkElements + aMesh.numberOfBoundaryElements;
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
    
    t = aMesh.numberOfBulkElements;
    n = aMesh.numberOfBulkElements + aMesh.numberOfBoundaryElements;
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
        free_ivector(edgeDofs, 0, aMesh.numberOfEdges-1);
        edgeDofs = NULL;
    }
    if (faceDofs != NULL) {
        free_ivector(faceDofs, 0,  aMesh.numberOfFaces-1);
        faceDofs = NULL;
    }
    free_ivector(defDofs, 0, solContainers->size2DefDofs-1);
    
    return k;
}

/********************************************************************************************
 
    Adds a new variable to the list of variables

*********************************************************************************************/
-(void)addVariableTo:(NSMutableArray *)anArray mesh:(FEMMesh *)aMesh solution:(FEMSolution *)aSolution name:(NSString *)name dofs:(int )dofs container:(variableArraysContainer *)aContainer component:(BOOL)component ifOutput:(BOOL *)output ifSecondary:(BOOL *)secondary type:(int *)aType {
    
    FEMVariable *newVariable;
    variableArraysContainer *varContainers = NULL;
    
    for (FEMVariable *variable in anArray) {
        if ([variable.name isEqualToString:name] == YES) {
            return; // Variable alreadt exists, don't do anything;
        }
    }
    
    newVariable = [[FEMVariable alloc] init];
    varContainers = newVariable.getContainers;
    
    newVariable.dofs = dofs;
    if (aContainer->Perm != NULL) {
        varContainers->Perm = aContainer->Perm;
        varContainers->sizePerm = aContainer->sizePerm;
    }
    
    newVariable.norm = 0.0;
    newVariable.prevNorm = 0.0;
    
    if (component == YES) {
        // If component is true, then we are adding a variable which is itself
        // a component of another variable with dof > 1. In this case, we set
        // an array of pointers where each pointer points to the appropriate
        // location in the variable array to which the component variable belongs.
        
        if (varContainers->ComponentValues != NULL) {
            varContainers->ComponentValues = aContainer->ComponentValues;
            varContainers->sizeComponentValues = aContainer->sizeComponentValues;
        } else if (*secondary == YES && varContainers->ComponentSecondaryToValues != NULL) {
            varContainers->ComponentSecondaryToValues = aContainer->ComponentSecondaryToValues;
            varContainers->sizeComponentSecondaryToValues = aContainer->sizeComponentSecondaryToValues;
        }
    } else {
        if (varContainers->Values != NULL) {
            varContainers->Values = aContainer->Values;
            varContainers->sizeValues = aContainer->sizeValues;
        } else if (*secondary == YES && varContainers->SecondaryToValues != NULL) {
            varContainers->SecondaryToValues = aContainer->SecondaryToValues;
            varContainers->sizeSecondaryToValues = aContainer->sizeSecondaryToValues;
        }
    }
    
    newVariable.nonLinChange = 0.0;
    newVariable.steadyChange = 0.0;
    newVariable.nonLinIter = 0;
    
    newVariable.solution = aSolution;
    newVariable.primaryMesh = aMesh;
    
    newVariable.valid = YES;
    newVariable.output = YES;
    newVariable.secondary = NO;
    newVariable.valuesChanged = YES;
    
    // Converged information undefined = -1 not 0, yes = 1
    newVariable.nonLinConverged = -1;
    newVariable.steadyConverged = -1;
    
    if (secondary != NULL) {
        NSLog(@"Secondary: %@", name);
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

-(void)addVectorTo:(NSMutableArray *)anArray mesh:(FEMMesh *)aMesh solution:(FEMSolution *)aSolution name:(NSString *)name dofs:(int *)dofs container:(variableArraysContainer *)aContainer ifOutput:(BOOL *)output ifSecondary:(BOOL *)secondary global:(BOOL *)global initValue:(double *)initValue {
    
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
        //   "solution variable" velocity
        //   "solution variable" acceleration
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
                nsize = aMesh.numberOfNodes;
            }
        } else {
            nsize = aMesh.numberOfNodes;
        }
        varContainers->Values = doublevec(0, (ndofs*nsize)-1);
        varContainers->sizeValues = (ndofs*nsize);
        memset( varContainers->Values, 0.0, ((ndofs*nsize)*sizeof(varContainers->Values)) );
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
                 for (j=(i-1); i<aContainer->sizeSecondaryToValues; j+=ndofs) {
                     varContainers->ComponentSecondaryToValues[k] = aContainer->SecondaryToValues[j];
                     k++;
                 }
                [self addVariableTo:anArray mesh:aMesh solution:aSolution name:tmpName dofs:1 container:varContainers component:YES ifOutput:output ifSecondary:secondary type:NULL];
            } else {
                varContainers->ComponentValues = malloc ( (varContainers->sizeValues/ndofs) * sizeof ( double * ));
                k = 0;
                for (j=(i-1); i<varContainers->sizeValues; j+=ndofs) {
                    varContainers->ComponentValues[k] = &varContainers->Values[j];
                    k++;
                }
                varContainers->sizeComponentValues = (varContainers->sizeValues/ndofs);
                [self addVariableTo:anArray mesh:aMesh solution:aSolution name:tmpName dofs:1 container:varContainers component:YES ifOutput:output ifSecondary:secondary type:NULL];
            }
        }
    }
    
    [self addVariableTo:anArray mesh:aMesh solution:aSolution name:name dofs:ndofs container:varContainers component:NO ifOutput:output ifSecondary:secondary type:NULL];
    free(varContainers);
}

-(FEMVariable *)getVariableFrom:(NSMutableArray *)anArray model:(FEMModel *)aModel name:(NSString *)name onlySearch:(BOOL *)only maskName:(NSString *)maskName info:(BOOL *)found {
    
    int i, j, k, l, n, dofs, aNumber;
    FEMVariable *var = nil, *pVar = nil, *tmp = nil;
    FEMListUtilities *listUtilities;
    FEMSolution *solution;
    FEMMesh *mesh, *currentMesh;
    FEMProjector *projector;
    variableArraysContainer *varContainers = NULL, *pvarContainers = NULL, *bufferContainers = NULL;
    NSString *tmpname;
    BOOL onlyThis, stat, globalBubbles, output;
    
    *found = NO;
    for (FEMVariable *variable in anArray) {
     
        if ([variable.name isEqualToString:name] == YES) {
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
    
    for (FEMMesh *mesh in aModel.meshes) {
        if ([anArray isEqualToArray:mesh.variables] == NO) {
            onlyThis = YES;
            pVar = [self getVariableFrom:mesh.variables model:aModel name:name onlySearch:&onlyThis maskName:NULL info:&stat];
            if (pVar != nil) {
                if ([mesh isEqual:pVar.primaryMesh]) {
                    break;
                }
            }
        }
       
    }
    if (pVar == nil) return var;
    
    if (*found == NO ) {
        listUtilities = [[FEMListUtilities alloc] init];
        solution = (FEMSolution *)pVar.solution;
        if ( (solution.solutionInfo)[@"bubbles in global system"] != nil ) {
            globalBubbles = [(solution.solutionInfo)[@"bubbles in global system"] boolValue];
        } else {
            globalBubbles = YES;
        }
        
        mesh = (FEMMesh *)aModel.mesh;
        dofs = mesh.numberOfNodes * pVar.dofs;
        if (globalBubbles == YES) dofs = dofs + mesh.maxBdofs * mesh.numberOfBulkElements * pVar.dofs;
        
        var = [[FEMVariable alloc] init];
        varContainers = var.getContainers;
        varContainers->Values = doublevec(0, dofs-1);
        varContainers->sizeValues = dofs;
        memset( varContainers->Values, 0.0, (dofs*sizeof(varContainers->Values)) );
        varContainers->Perm = NULL;
        pvarContainers = pVar.getContainers;
        if (pvarContainers->Perm != NULL) {
            varContainers->Perm = intvec(0, (dofs/pVar.dofs)-1);
            varContainers->sizePerm = dofs/pVar.dofs;
            memset( varContainers->Perm, -1, ((dofs/pVar.dofs)*sizeof(varContainers->Perm)) );
            
            n = [self initialPermutationInMesh:mesh model:aModel solution:solution equation:(solution.solutionInfo)[@"equation"] permutation:varContainers->Perm DGSolution:NULL globalBubbles:&globalBubbles];
            
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
        var = [self getVariableFrom:anArray model:aModel name:name onlySearch:&onlyThis maskName:NULL info:&stat];
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
            tmp = [self getVariableFrom:anArray model:aModel name:@"velocity 1" onlySearch:&onlyThis maskName:NULL info:&stat];
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
            tmp = [self getVariableFrom:anArray model:aModel name:@"velocity 2" onlySearch:&onlyThis maskName:NULL info:&stat];
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
                tmp = [self getVariableFrom:anArray model:aModel name:@"velocity 3" onlySearch:&onlyThis maskName:NULL info:&stat];
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
            tmp = [self getVariableFrom:anArray model:aModel name:@"pressure" onlySearch:&onlyThis maskName:NULL info:&stat];
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
                    tmp = [self getVariableFrom:anArray model:aModel name:tmpname onlySearch:&onlyThis maskName:NULL info:&stat];
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
        var = [self getVariableFrom:anArray model:aModel name:name onlySearch:&onlyThis maskName:NULL info:&stat];
    }
    
    // Build a temporary variable array of variables to be interpolated
    NSMutableArray *tmpArryay;
    mesh = nil;
    mesh = (FEMMesh *)pVar.primaryMesh;
    if ([pVar.name isEqualToString:@"flow solution"]) {
        tmpArryay = [NSMutableArray arrayWithObjects:[self getVariableFrom:mesh.variables model:aModel name:@"velocity 1" onlySearch:NULL maskName:NULL info:&stat],
                     [self getVariableFrom:mesh.variables model:aModel name:@"velocity 2" onlySearch:NULL maskName:NULL info:&stat],
                     [self getVariableFrom:mesh.variables model:aModel name:@"velocity 3" onlySearch:NULL maskName:NULL info:&stat],
                     [self getVariableFrom:mesh.variables model:aModel name:@"pressure" onlySearch:NULL maskName:NULL info:&stat], nil];
    } else if (pVar.dofs > 1){
        for (i=0; i<pVar.dofs; i++) {
             aNumber = i+1;
            tmpname = [self appendNameFromString:pVar.name component:&aNumber];
            [tmpArryay addObject:[self getVariableFrom:mesh.variables model:aModel name:tmpname onlySearch:NULL maskName:NULL info:&stat]];
        }
    }
    
    //Interpolation
    currentMesh = (FEMMesh *)aModel.mesh;
    if (maskName != NULL) {
        [self FEMUtils_interpolateMesh:mesh toMesh:currentMesh oldVariables:tmpArryay newVariables:anArray model:aModel quadrantTree:NULL withProjector:NO projector:nil mask:maskName];
    } else {
        [self FEMUtils_interpolateMesh:mesh toMesh:currentMesh oldVariables:tmpArryay newVariables:anArray model:aModel quadrantTree:NULL withProjector:YES projector:projector mask:nil];
    }
    
    onlyThis = YES;
    var = [self getVariableFrom:anArray model:aModel name:name onlySearch:&onlyThis maskName:NULL info:&stat];
    if ([var.name isEqualToString:@"flow solution"]) {
        tmp = [self getVariableFrom:anArray model:aModel name:@"velocity 1" onlySearch:&onlyThis maskName:NULL info:&stat];
        if (tmp != nil) {
            tmp.valid = YES;
            tmp.valuesChanged = YES;
        }
        
        tmp = [self getVariableFrom:anArray model:aModel name:@"velocity 2" onlySearch:&onlyThis maskName:NULL info:&stat];
        if (tmp != nil) {
            tmp.valid = YES;
            tmp.valuesChanged = YES;
        }
        
        if (var.dofs == 4) {
            tmp = nil;
            tmp = [self getVariableFrom:anArray model:aModel name:@"velocity 3" onlySearch:&onlyThis maskName:NULL info:&stat];
            if (tmp != nil) {
                tmp.valid = YES;
                tmp.valuesChanged = YES;
            }
        }
        
        tmp = [self getVariableFrom:anArray model:aModel name:@"pressure" onlySearch:&onlyThis maskName:NULL info:&stat];
        if (tmp != nil) {
            tmp.valid = YES;
            tmp.valuesChanged = YES;
        }
    } else if (var.dofs > 1){
        for (i=0; i<var.dofs; i++) {
            tmp = nil;
            aNumber = i+1;
            tmpname = [self appendNameFromString:pVar.name component:&aNumber];
            tmp = [self getVariableFrom:anArray model:aModel name:tmpname onlySearch:&onlyThis maskName:NULL info:&stat];
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
-(double)interpolateCurveTvalues:(double *)tValues fValues:(double *)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double *)cubicCoeff {
    
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
    
    cubic = NO;
    if (cubicCoeff != NULL) cubic = YES;
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

-(void)solveLinearSystem2x2:(double **)a afterSolve:(double *)x rightHandSide:(double *)b {
    
    double detA;
    
    detA = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    
    if (detA == 0.0) {
        errorfunct("solveLinearSystem2x2", "Singular matrix, bad!!!");
        return;
    }
    
    detA = 1.0 / detA;
    x[0] = detA * ( a[1][1] * b[0] - a[0][1] * b[1] );
    x[1] = detA * ( a[0][0] * b[1] - a[1][0] * b[0] );
    
}

-(void)solveLinearSystem3x3:(double **)a afterSolve:(double *)x rightHandSide:(double *)b {
    
    double **c, *y, *g, s, t, q;
    
    c = doublematrix(0, 1, 0, 1);
    y = doublevec(0, 1);
    g = doublevec(0, 1);
    
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
    
    free_dmatrix(c, 0, 1, 0, 1);
    free_dvector(y, 0, 1);
    free_dvector(g, 0, 1);
    
}

-(FEMMatrix *)meshProjector:(FEMMesh *)mesh1 secondmesh:(FEMMesh *)mesh2 model:(FEMModel *)aModel useQuadrantTree:(BOOL *)quadrantTree transpose:(BOOL *)trans {
    
    FEMProjector *projector;
    FEMMatrix *projectorMatrix;
    
    if (quadrantTree != NULL) {
        [self FEMUtils_interpolateQMesh:mesh1 toMesh:mesh2 oldVariables:nil newVariables:nil model:aModel quadrantTree:quadrantTree withProjector:YES projector:projector mask:nil nodesPresent:NULL];
    } else {
        [self FEMUtils_interpolateQMesh:mesh1 toMesh:mesh2 oldVariables:nil newVariables:nil model:aModel quadrantTree:NULL withProjector:YES projector:projector mask:nil nodesPresent:NULL];
    }
    
    projectorMatrix = projector.matrix;
    if (trans != NULL) {
        if (*trans == YES) {
            projectorMatrix = projector.tMatrix;
        }
    }
    
    return projectorMatrix;
}

-(double)cublicSplineX:(double *)x Y:(double *)y R:(double *)r T:(double)t {
    
    double s, a, b, c, h, lt;
    
    h = x[1] - x[0];
    a = -2.0 * (y[1] - y[0]) + (r[0] + r[1]) * h;
    b = 3.0 * (y[1] - y[0]) - (2.0*r[0] + r[1]) * h;
    c = r[0] * h;
    
    lt = (t - x[0]) / h;
    s = ( (a*lt + b) * lt + c ) / h;
    
    return s;
}

-(NSString *)appendNameFromString:(NSString *)string component:(int *)component {
    
    int i, j, dofs, comp, dofsTot;
    NSString *str;
    NSMutableString *str1;
    NSRange ind, ind1;
    
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
            if (ind1.location == NSNotFound) errorfunct("appendNameFromString", "Missing separator ':' in variable definition using '[ ]' syntax.");
            if (j == 0) {
                if (ind1.location < ind.location) errorfunct("appendNameFromString", "Syntax error in variable definition.");
            } else {
                if (ind1.location == 0) errorfunct("appendNameFromString", "Syntax error in variable definition.");
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
        if (dofs > 0) {
            dofs = *component - dofsTot + dofs;
            str1 = [NSMutableString stringWithString:str];
            [str1 appendString:@" "];
            [str1 appendString:[NSString stringWithFormat:@"%d",dofs]];
            str = [NSString stringWithString:str1];
        }
    }
    
    return str;
}

-(NSString *)appendNameFromVariable:(FEMVariable *)variable component:(int *)component {
    
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
-(NSString *)nextFreeKeyword:(NSString *)keyword0 dictionary:(NSMutableDictionary *)dictionary {
    
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
 Check the feasability of solution options
******************************************************/
-(void)checkOptionsInSolution:(FEMSolution *)solution {
    
    NSString *str;
    
    if ((solution.solutionInfo)[@"linear system solver"] != nil) {
        str = (solution.solutionInfo)[@"linear system solver"];
        
        if ([str isEqualToString:@"direct"] == YES) {
            if ((solution.solutionInfo)[@"linear system direct method"] != nil) {
                str = (solution.solutionInfo)[@"linear system direct method"];
                
                //TODO: app support for parallel run where in case of a parallel run
                // the direct must be MUMPS
                if ([str isEqualToString:@"mumps"] == YES) errorfunct("checkOptionsInSolution", "Currenly no serial version of the mumps solver implemented,");
                
                if ([str isEqualToString:@"banded"] == YES) {
                    
                } else if ([str isEqualToString:@"umfpack"] == YES || [str isEqualToString:@"big umfpack"] == YES) {
                    
                } else if ([str isEqualToString:@"mumps"] == YES) {
                    
                } else if ([str isEqualToString:@"superlu"] == YES) {
                    errorfunct("checkOptionsInSolution", "SuperLU solver is not available.");
                } else if ([str isEqualToString:@"pardiso"] == YES) {
                    errorfunct("checkOptionsInSolution", "Pardiso solver is not available.");
                } else if ([str isEqualToString:@"cholmod"] == YES || [str isEqualToString:@"spqr"] == YES) {
                    errorfunct("checkOptionsInSolution", "Cholmod solver is not available.");
                } else {
                    errorfunct("checkOptionsInSolution", "Unknown direct solver method.");
                }
            } else {
                //TODO: add support for parallel run since then it should be mumps by default.
                str = @"umfpack";
                NSLog(@"checkOptionsInSolution: setting the linear system direct method to %@\n", str);
                [solution.solutionInfo setObject:str forKey:@"linear system direct method"];
            }
        }
    } // If "linear system solver" is not given, it will be set by default to iterative later in the processing
}


/********************************************************************************************************
    Add the generic components to each solution
    A few solutions are for historical reasons given a special treatment
********************************************************************************************************/
-(void)addEquationBasicsToSolution:(FEMSolution *)solution name:(NSString *)name model:(FEMModel *)model transient:(BOOL)transient {
    
    int i, j, k, l, maxDim, minVal, maxVal, dofs, ndeg, maxNDOFs, maxDGDOFs, maxEDOFs, maxFDOFs, maxBDOFs,
        matrixFormat, nrows, nsize, type;
    int *perm;
    double initValue, *sol;
    NSString *eq, *str, *varName, *tmpName;
    NSMutableString *string;
    NSNumber *yesNumber;
    BOOL isAssemblySolution, isCoupledSolution, isBlockSolution, variableGlobal, variableOutput, found, stat,
         globalBubbles, bandwidthOptimize, discontinuousGalerkin;
    NSRange range;
    Element_t *elements = NULL, *edges = NULL, *faces = NULL;
    FEMVariable *variable = nil, *newVariable = nil;
    FEMListUtilities *listUtilities;
    FEMElementUtils *elementUtils;
    variableArraysContainer *bufferContainers = NULL, *varContainers = NULL;
    
    listUtilities = [[FEMListUtilities alloc] init];
    
    // If there is a matrix level "Flux Corrected Transport" then it's required
    // to use global matrices for time integration
    if ([(solution.solutionInfo)[@"linear system fct"] boolValue] == YES) {
        [listUtilities addLogicalInClassList:solution theVariable:@"use global mass matrix" withValue:YES];
    }
    
    if ((solution.solutionInfo)[@"equation"] != nil) {
        eq = (solution.solutionInfo)[@"equation"];
        elements = solution.mesh.getElements;
        maxDim = 0;
        for (i=0; i<solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements; i++) {
            if ([listUtilities checkElementEquation:model forElement:&elements[i] andEquation:eq]) {
                maxDim = max(elements[i].Type.dimension, maxDim);
            }
        }
        [listUtilities addIntegerInClassList:solution theVariable:@"active mesh dimension" withValue:maxDim];
    }
    
    // TODO: Here Elmer calls a procedure post-fixed with "_Init" for a given solver.
    // Do we need to do that?
    
    solution.solutionMode = SOLUTION_MODE_DEFAULT;
    if ([(solution.solutionInfo)[@"auxiliary solution"] boolValue] == YES) solution.solutionMode = SOLUTION_MODE_AUXILIARY;
    if ([(solution.solutionInfo)[@"coupled solution"] boolValue] == YES) solution.solutionMode = SOLUTION_MODE_COUPLED;
    if ([(solution.solutionInfo)[@"block solution"] boolValue] == YES) solution.solutionMode = SOLUTION_MODE_BLOCK;
    if ([(solution.solutionInfo)[@"assembly solution"] boolValue] == YES) solution.solutionMode = SOLUTION_MODE_ASSEMBLY;
    
    if (solution.solutionMode == SOLUTION_MODE_DEFAULT) {
        if ((solution.solutionInfo)[@"equation"] != nil) {
            eq = (solution.solutionInfo)[@"equation"];
        } else {
            solution.solutionMode = SOLUTION_MODE_AUXILIARY;
        }
    }
    
    isCoupledSolution = (solution.solutionMode == SOLUTION_MODE_COUPLED) ? YES : NO;
    isBlockSolution = (solution.solutionMode == SOLUTION_MODE_BLOCK) ? YES : NO;
    isAssemblySolution = (solution.solutionMode == SOLUTION_MODE_ASSEMBLY) ? YES : NO;
    isAssemblySolution = (isAssemblySolution == YES || (isCoupledSolution == YES && (solution.plugInPrincipalClassInstance == nil || solution.selector == NULL)) || (isBlockSolution == YES && (solution.plugInPrincipalClassInstance == nil || solution.selector == NULL))) ? YES : NO;
    
    // Default order of equation
    solution.order = 1;
    solution.timeOrder = 1;
    
    found = NO;
    if (transient == YES) {
        if ((solution.solutionInfo)[@"time stepping method"] != nil) {
            str = (solution.solutionInfo)[@"time stepping method"];
            found = YES;
        } else {
            str = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"time stepping method" info:&found];
            if (found == YES) {
                [solution.solutionInfo setObject:str forKey:@"time stepping method"];
            }
        }
        
        if (found == YES) {
            if ([str isEqualToString:@"bdf"] == YES) {
                if ((solution.solutionInfo)[@"bdf order"] != nil) {
                    solution.order = [(solution.solutionInfo)[@"bdf order"] intValue];
                    if (solution.order < 1) solution.order = 1;
                    if (solution.order > 5) solution.order = 5;
                } else {
                    minVal = 1;
                    maxVal = 5;
                    solution.order = [listUtilities listGetInteger:model inArray:model.simulation.valuesList forVariable:@"bdf order" info:&found minValue:&minVal maxValue:&maxVal];
                }
                if (found == NO) {
                    solution.order = 2;
                    NSLog(@"addEquationBasicsToSolution: BDF order set by default to 2.\n");
                }
            } else if ([str isEqualToString:@"runge-kutta"]) {
                minVal = 2;
                maxVal = 4;
                solution.order = [listUtilities listGetInteger:model inArray:model.simulation.valuesList forVariable:@"runge-kutta order" info:&found minValue:&minVal maxValue:&maxVal];
                if (found == NO) solution.order = 2;
            }
        } else {
            NSLog(@"addEquationBasicsToSolution: time stepping method set by default to implicit Euler\n");
            [solution.solutionInfo setObject:@"implicit euler" forKey:@"time stepping method"];
        }
    }
    
    dofs = model.dimension;
    initValue = 0.0;
    
    // These are historical solutions that may be built-in on some MDFs
    // Special strategies are used for them.
    // TODO: we may migrate them to plug-ins?
    
    if ([name isEqualToString:@"navier-stokes"] == YES) {
        if ((solution.solutionInfo)[@"variable"] == nil) {
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
        if ((solution.solutionInfo)[@"variable"] == nil) {
            [solution.solutionInfo setObject:@"-dofs 3 magnetic field" forKey:@"variable"];
            [solution.solutionInfo setObject:@"electric current[electric current:3]"
                                      forKey:[self nextFreeKeyword:@"exported variable" dictionary:solution.solutionInfo]];
        }
    } else if ([name isEqualToString:@"stress analysis"] == YES) {
        if ((solution.solutionInfo)[@"variable"] == nil) {
            if (dofs == 2) {
                [solution.solutionInfo setObject:@"-dofs 2 displacement" forKey:@"variable"];
            } else {
                [solution.solutionInfo setObject:@"-dofs 3 displacement" forKey:@"variable"];
            }
        }
    } else if ([name isEqualToString:@"mesh update"] == YES) {
        if ((solution.solutionInfo)[@"variable"] == nil) {
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
        if ((solution.solutionInfo)[@"variable"] == nil) {
            [solution.solutionInfo setObject:@"temperature" forKey:@"variable"];
        }
        
        if ((solution.solutionInfo)[@"radiation solver"] == nil) {
            yesNumber = @YES;
            [solution.solutionInfo setObject:yesNumber forKey:@"radiation solver"];
        }
    }
    
    // We allocate memory for the string to hold the name of the plug-in only if we find that we are
    // using a plug-in when we parse the MDF. So if we get nil, we are not using a plug-in.
    if (solution.plugInName != nil) {
        // TODO: This is here where we should load the solution plug-in by its name (if present) and
        // create an instance of the plug-in principal class. We should also check here if
        // the plug-in class conforms to the protocol for solving equations.

    } else { // We are woking with a built-in solution computer
        solution.selector = @selector(fieldSolutionComputer::::);
    }
    
    //Initialize and get the variable
    solution.timeOrder = 0;
    solution.matrix = nil;
    
    if ((solution.solutionInfo)[@"variable"] == nil) {
        //Variable does not exist
        variable = [[FEMVariable alloc] init];
        solution.variable = variable;
    } else if (isCoupledSolution == YES && (solution.plugInPrincipalClassInstance == nil || solution.selector == NULL)) {
        
    } else if (isBlockSolution == YES && (solution.plugInPrincipalClassInstance == nil || solution.selector == NULL)) {
        
    } else {
        varName = (solution.solutionInfo)[@"variable"];
        
        // It may be a normal field variable or a global (0D) variable
        variableGlobal = [(solution.solutionInfo)[@"variable global"] boolValue];
        
        if ((solution.solutionInfo)[@"variable output"] != nil) {
            variableOutput = [(solution.solutionInfo)[@"variable output"] boolValue];
        } else variableOutput = YES;
        
        if ((solution.solutionInfo)[@"variable dofs"] != nil) {
            dofs = [(solution.solutionInfo)[@"variable dofs"] intValue];
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
                    NSLog(@"addEquationBasicsToSolution: invalid dof in variable definiton.\n");
                    NSLog(@"addEquationBasicsToSolution: the incorrect value was: %@\n",
                          [[str substringFromIndex:range.location+1] substringToIndex:i-(range.location+1)]);
                }
                dofs = dofs + k;
                j = i+1;
            }
        }
        
        if ([varName characterAtIndex:0] == '-') {
            if ([[[varName substringFromIndex:0] substringToIndex:10] isEqualToString:@"-nooutput "] == YES) {
                variableOutput = NO;
                varName = [varName substringFromIndex:10];
            }
            if ([[[varName substringFromIndex:0] substringToIndex:8] isEqualToString:@"-global "] == YES) {
                variableGlobal = YES;
                varName = [varName substringFromIndex:8];
            }
            if ([[[varName substringFromIndex:0] substringToIndex:6] isEqualToString:@"-dofs "] == YES) {
                i = 5;
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
                    NSLog(@"addEquationBasicsToSolution: invalid dof in variable definiton.\n");
                    NSLog(@"addEquationBasicsToSolution: the incorrect value was: %@\n", [[varName substringFromIndex:i] substringToIndex:j-i]);
                }
                varName = [varName substringFromIndex:j+1];
            }
        }
        if (dofs == 0) dofs = 1;
        
        // If the variable is global then it has nothing to do with the mesh
        // and it may be simply allocated
        if (variableGlobal == YES) {
            solution.solutionMode = SOLUTION_MODE_GLOBAL;
            sol = doublevec(0, dofs-1);
            memset( sol, 0.0, (dofs*sizeof(sol)) );
            
            bufferContainers = allocateVariableContainer();
            bufferContainers->Values = sol;
            bufferContainers->sizeValues = dofs;
            bufferContainers->Perm = NULL;
            [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:varName dofs:dofs container:bufferContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
            solution.variable = [self getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:NULL info:&found];
            if (dofs > 1) {
                for (i=1; i<=dofs; i++) {
                    tmpName = [self appendNameFromString:varName component:&i];
                    bufferContainers->ComponentValues = malloc ( 1 * sizeof ( double * ));
                    bufferContainers->ComponentValues[0] = &sol[i-1];
                    bufferContainers->sizeComponentValues = 1;
                    [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:tmpName dofs:1 container:bufferContainers component:YES ifOutput:NULL ifSecondary:NULL type:NULL];
                }
            }
            free(bufferContainers);
        } else {
            // If variable is a field variable create a permutation and matrix related to it
            if ((solution.solutionInfo)[@"equation"] != nil) {
                eq = (solution.solutionInfo)[@"equation"];
            } else {
                errorfunct("addEquationBasicsToSolution", "Variable exists but the equation is not defined.");
            }
            found = NO;
            for (FEMEquation *equation in model.equations) {
                if ([listUtilities listGetLogical:model inArray:equation.valuesList forVariable:eq info:&stat] == YES) {
                    found = YES;
                    break;
                }
            }
            if (found == NO) {
                NSLog(@"addEquationBasicsToSolution: variable %@ exists but it's not associated to any equation\n", varName);
                errorfunct("addEquationBasicsToSolution", "Program terminating now...");
            }
            
            // Compute the size of the permutation vector
            ndeg = 0;
            if (YES) {
                eq = (solution.solutionInfo)[@"equation"];
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
                
                if ((solution.solutionInfo)[@"bubbles in global system"] != nil) {
                    globalBubbles = [(solution.solutionInfo)[@"bubbles in global system"] boolValue];
                } else globalBubbles = YES;
                
                ndeg = ndeg + solution.mesh.numberOfNodes;
                if (maxEDOFs > 0) ndeg = ndeg + maxEDOFs * solution.mesh.numberOfEdges;
                if (maxFDOFs > 0) ndeg = ndeg + maxFDOFs * solution.mesh.numberOfFaces;
                if (globalBubbles == YES) ndeg = ndeg + maxBDOFs * solution.mesh.numberOfBulkElements;
                if ((solution.solutionInfo)[@"discontinuous galerkin"] != nil) {
                    if ([(solution.solutionInfo)[@"discontinuous galerkin"] boolValue] == YES) {
                        ndeg = max( ndeg, maxDGDOFs * (solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements) );
                    }
                }
            }
            
            if ((solution.solutionInfo)[@"radiation solution"] != nil) {
                if ([(solution.solutionInfo)[@"radiation solution"] boolValue] == YES) {
                    //TODO: Need to implement this
                }
            }
            
            if ((solution.solutionInfo)[@"optimize bandwidth"] != nil) {
                bandwidthOptimize = [(solution.solutionInfo)[@"optimize bandwidth"] boolValue];
            } else bandwidthOptimize = YES;
            [self checkOptionsInSolution:solution];
            
            perm = intvec(0, ndeg-1);
            memset( perm, 0, (ndeg*sizeof(perm)) );
            
            elementUtils = [[FEMElementUtils alloc] init];
            matrixFormat = MATRIX_CRS;
            if ((solution.solutionInfo)[@"discontinuous galerkin"] != nil) {
                discontinuousGalerkin = [(solution.solutionInfo)[@"discontinuous galerkin"] boolValue];
            } else discontinuousGalerkin = NO;
            solution.matrix = [elementUtils createMatrixInModel:model forSolution:solution mesh:solution.mesh dofs:dofs permutation:perm sizeOfPermutation:ndeg matrixFormat:matrixFormat optimizeBandwidth:bandwidthOptimize equationName:eq discontinuousGalerkinSolution:&discontinuousGalerkin globalBubbles:&globalBubbles];
            nrows = dofs * ndeg;
            if (solution.matrix != nil) nrows = solution.matrix.numberOfRows;
            
            // Basically the solution could be matrix free but still the matrix is used
            // here temporarily since it's needed when making the permutation vector
            if ((solution.solutionInfo)[@"no matrix"] != nil) {
                if ([(solution.solutionInfo)[@"no matrix"] boolValue] == YES) {
                    solution.solutionMode = SOLUTION_MODE_MATRIXFREE;
                    [solution.matrix deallocation];
                    solution.matrix = nil;
                }
            }
            
            if (nrows > 0) {
                sol = doublevec(0, nrows-1);
                for (i=0; i<nrows; i++) {
                    sol[i] = initValue;
                }
                
                bufferContainers = allocateVariableContainer();
                bufferContainers->Values = sol;
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
                        for (j=(i-1); i<=nrows-dofs+(i-1); j+=dofs) {
                            bufferContainers->ComponentValues[k] = &sol[j];
                            k++;
                        }
                        bufferContainers->sizeComponentValues = nrows/dofs;
                        [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:tmpName dofs:1 container:bufferContainers component:YES ifOutput:NULL ifSecondary:NULL type:NULL];
                    }                    
                }
                free(bufferContainers);
            }
            
            //TODO: add support for parallel run here
            if ((solution.solutionInfo)[@"discontinuous galerkin"] != nil) {
                if ([(solution.solutionInfo)[@"discontinuous galerkin"] boolValue] == YES) solution.variable.type = VARIABLE_ON_NODES_ON_ELEMENTS;
            }
        }
    }
    
    // Add the exported variables which are typically auxiliary variables derived from
    // the solution without their own matrix equation 
    l = 0;
    while (1) {
        l++;
        str = [self appendNameFromString:@"exported variable" component:&l];
        if ((solution.solutionInfo)[str] == nil) break;
        varName = (solution.solutionInfo)[str];
        
        string = (NSMutableString *)[self appendNameFromString:@"exported variable" component:&l];
        [string appendString:@" output"];
        if ((solution.solutionInfo)[string] != nil) {
            variableOutput = [(solution.solutionInfo)[string] boolValue];
        } else variableOutput = YES;
        
        string = (NSMutableString *)[self appendNameFromString:@"exported variable" component:&l];
        [string appendString:@" dofs"];
        if ((solution.solutionInfo)[string] != nil) {
            dofs = [(solution.solutionInfo)[string] intValue];
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
                    NSLog(@"addEquationBasicsToSolution: invalid dof in variable definiton.\n");
                    NSLog(@"addEquationBasicsToSolution: the incorrect value was: %@\n",
                          [[str substringFromIndex:range.location+1] substringToIndex:i-(range.location+1)]);
                }
                dofs = dofs + k;
                j = i+1;
            }
        }
        
        variableOutput = YES;
        variableGlobal = NO;
        
        if ([varName characterAtIndex:0] == '-') {
            if ([[[varName substringFromIndex:0] substringToIndex:10] isEqualToString:@"-nooutput "] == YES) {
                variableOutput = NO;
                varName = [varName substringFromIndex:10];
            }
            if ([[[varName substringFromIndex:0] substringToIndex:8] isEqualToString:@"-global "] == YES) {
                variableGlobal = YES;
                varName = [varName substringFromIndex:8];
            }
            if ([[[varName substringFromIndex:0] substringToIndex:6] isEqualToString:@"-dofs "] == YES) {
                i = 5;
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
                    NSLog(@"addEquationBasicsToSolution: invalid dof in variable definiton.\n");
                    NSLog(@"addEquationBasicsToSolution: the incorrect value was: %@\n", [[varName substringFromIndex:i] substringToIndex:j-i]);
                }
                varName = [varName substringFromIndex:j+1];
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
            
            sol = doublevec(0, nsize-1);
            memset( sol, 0.0, (nsize*sizeof(sol)) );
            if (perm != NULL) {
                bufferContainers = allocateVariableContainer();
                bufferContainers->Values = sol;
                bufferContainers->sizeValues = nsize;
                bufferContainers->Perm = perm;
                bufferContainers->sizePerm = varContainers->sizePerm;
                type = solution.variable.type;
                [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:varName dofs:dofs container:bufferContainers component:NO ifOutput:&variableOutput ifSecondary:NULL type:&type];
            } else {
                bufferContainers = allocateVariableContainer();
                bufferContainers->Values = sol;
                bufferContainers->sizeValues = nsize;
                type = solution.variable.type;
                [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:varName dofs:dofs container:bufferContainers component:NO ifOutput:&variableOutput ifSecondary:NULL type:&type];
            }
            
            if (dofs > 1 && variableOutput == NO) {
                for (i=1; i<=dofs; i++) {
                    tmpName = [self appendNameFromString:varName component:&i];
                    bufferContainers->ComponentValues = malloc ( (nsize/dofs) * sizeof ( double * ));
                    k = 0;
                    for (j=(i-1); i<=nsize-dofs+(i-1); j+=dofs) {
                        bufferContainers->ComponentValues[k] = &sol[j];
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
    
    if ((solution.solutionInfo)[@"solve equation"] != nil) {
        if ([(solution.solutionInfo)[@"solve equation"] isEqualToString:@"never"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_NEVER;
        } else if ([(solution.solutionInfo)[@"solve equation"] isEqualToString:@"always"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_ALWAYS;
        } else if ([(solution.solutionInfo)[@"solve equation"] isEqualToString:@"after all"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_ALL;
        } else if ([(solution.solutionInfo)[@"solve equation"] isEqualToString:@"before all"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_ALL;
        } else if ([(solution.solutionInfo)[@"solve equation"] isEqualToString:@"before time step"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_TIME;
        } else if ([(solution.solutionInfo)[@"solve equation"] isEqualToString:@"after time step"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_TIME;
        } else if ([(solution.solutionInfo)[@"solve equation"] isEqualToString:@"before saving"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_SAVE;
        } else if ([(solution.solutionInfo)[@"solve equation"] isEqualToString:@"after saving"] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_SAVE;
        } else {
            solution.solutionSolveWhen = SOLUTION_SOLVE_ALWAYS;
        }
    }
    
    if ((solution.solutionInfo)[@"before all"] != nil) {
        if ([(solution.solutionInfo)[@"before all"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_ALL;
        }
    } else if ((solution.solutionInfo)[@"before simulation"] != nil) {
        if ([(solution.solutionInfo)[@"before simulation"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_ALL;
        }
    } else if ((solution.solutionInfo)[@"after all"] != nil) {
        if ([(solution.solutionInfo)[@"after all"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_ALL;
        }
    } else if ((solution.solutionInfo)[@"after simulation"] != nil) {
        if ([(solution.solutionInfo)[@"after simulation"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_ALL;
        }
    } else if ((solution.solutionInfo)[@"before time step"] != nil) {
        if ([(solution.solutionInfo)[@"before time step"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_TIME;
        }
    } else if ((solution.solutionInfo)[@"after time step"] != nil) {
        if ([(solution.solutionInfo)[@"after time step"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_TIME;
        }
    } else if ((solution.solutionInfo)[@"before saving"] != nil) {
        if ([(solution.solutionInfo)[@"before saving"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AHEAD_SAVE;
        }
    } else if ((solution.solutionInfo)[@"after saving"] != nil) {
        if ([(solution.solutionInfo)[@"after saving"] boolValue] == YES) {
            solution.solutionSolveWhen = SOLUTION_SOLVE_AFTER_SAVE;
        }
    }
    
    // TODO: Do we need to add support for LinSolve?
}

/********************************************************************************************************
    Add information that is typically only needed if there's a matrix equation to work with.
    This should be called only after both the solution vector and matrix have been created.
********************************************************************************************************/
-(void)addEquationToSolution:(FEMSolution *)solution model:(FEMModel *)model transient:(BOOL)transient {
    
    int i, j, k, l, dofs, nrows, n, mgLevels;
    double *sol;
    NSString *tmpName, *methodName;
    NSMutableString *varName, *str;
    NSNumber *aNumber;
    BOOL found, onlySearch, secondary, variableOutput, harmonicAnal, eigenAnal, complexFlag, multigridActive, mgAlgebraic;
    FEMKernel *kernel;
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
    if ([(solution.solutionInfo)[@"apply limiter"] boolValue] == YES) {
        aNumber = @YES;
        [solution.solutionInfo setObject:aNumber forKey:@"calculate loads"];
    }
    
    // Create the variable needed for this computation of nodal loads: r = b-Ax
    if ([(solution.solutionInfo)[@"calculate loads"] boolValue] == YES) {
        varName = [NSMutableString stringWithString:solution.variable.name];
        [varName appendString:@" loads"];
        var = [self getVariableFrom:solution.mesh.variables model:model name:varName onlySearch:NULL maskName:NULL info:&found];
        if (var == nil) {
            sol = doublevec(0, variableContainers->sizeValues-1);
            dofs = solution.variable.dofs;
            memset( sol, 0.0, (variableContainers->sizeValues*sizeof(sol)) );
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
                    for (j=(i-1); i<=nrows-dofs+(i-1); j+=dofs) {
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
    
    harmonicAnal = [(solution.solutionInfo)[@"harmonic analysis"] boolValue];
    
    if (solution.matrix != nil) {
        matContainers = solution.matrix.getContainers;
        if (matContainers->RHS == NULL) {
            matContainers->RHS = doublevec(0, solution.matrix.numberOfRows-1);
            matContainers->sizeRHS = solution.matrix.numberOfRows;
            memset( matContainers->RHS, 0.0, (solution.matrix.numberOfRows*sizeof(matContainers->RHS)) );
            
            if (harmonicAnal == YES) {
                matContainers->RHS_im = doublevec(0, solution.matrix.numberOfRows-1);
                matContainers->sizeRHS_im = solution.matrix.numberOfRows;
                memset( matContainers->RHS_im, 0.0, (solution.matrix.numberOfRows*sizeof(matContainers->RHS_im)) );
            }
        }
    }
    
    eigenAnal = [(solution.solutionInfo)[@"eigen analysis"] boolValue];
    
    if (transient == YES && eigenAnal == NO && harmonicAnal == NO) {
        solution.timeOrder = 1;
        if ((solution.solutionInfo)[@"time derivative order"] != nil) {
            k = [(solution.solutionInfo)[@"time derivative order"] intValue];
            if (k < 0) k = 0;
            if (k > 2) k = 2;
            solution.timeOrder = min(max(1, k), 2);
        }
        
        if (solution.matrix != nil) {
            matContainers->Force = doublematrix(0, solution.matrix.numberOfRows-1, 0, (solution.timeOrder+1)-1);
            matContainers->size1force = solution.matrix.numberOfRows;
            matContainers->size2Force = solution.timeOrder + 1;
            for (i=0; i<solution.matrix.numberOfRows; i++) {
                for (j=0; j<(solution.timeOrder+1); j++) {
                    matContainers->Force[i][j] = 0.0;
                }
            }
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
            for (i=0; i<variableContainers->size1PrevValues; i++) {
                for (j=0; j<variableContainers->size2PrevValues; j++) {
                    variableContainers->PrevValues[i][j] = 0.0;
                }
            }
            
            if (solution.variable.dofs > 1) {
                if ([solution.variable.name isEqualToString:@"flow solution"]) {
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
        
        if ([(solution.solutionInfo)[@"calculate velocity"] boolValue] == YES || [(solution.solutionInfo)[@"nonlinear calculate velocity"] boolValue] == YES) {
            if (solution.timeOrder < 1) {
                NSLog(@"addEquationToSolution: velocity computation implemented only for time-dependent equations\n");
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
        
        if ([(solution.solutionInfo)[@"calculate acceleration"] boolValue] == YES) {
            if (solution.timeOrder == 1) {
                NSLog(@"addEquationToSolution: acceleration computation implemented only for 2nd order time equations\n");
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
                secondary = YES;
                [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:str dofs:solution.variable.dofs container:bufferContainers component:NO ifOutput:NULL ifSecondary:&secondary type:NULL];
                
                for (i=1; i<=solution.variable.dofs; i++) {
                    str = [NSMutableString stringWithString:solution.variable.name];
                    [str appendString:@" acceleration"];
                    [str appendString:@" "];
                    [str appendString:[NSString stringWithFormat:@"%d",i]];
                    varContainers->ComponentSecondaryToValues = malloc ( (variableContainers->size1PrevValues/solution.variable.dofs) * sizeof ( double * ));
                    k = 0;
                    for (j=(i-1); i<variableContainers->size1PrevValues; j+=solution.variable.dofs) {
                         bufferContainers->ComponentSecondaryToValues[k] =  &variableContainers->PrevValues[j][1];
                        k++;
                    }
                    [self addVariableTo:solution.mesh.variables mesh:solution.mesh solution:solution name:str dofs:1 container:bufferContainers component:YES ifOutput:NULL ifSecondary:&secondary type:NULL];
                    
                }
                free(bufferContainers);
            }
        }
    } else {
        solution.timeOrder = 0;
        
        if ([(solution.solutionInfo)[@"calculate derivative"] boolValue] == YES) {
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
            if ((solution.solutionInfo)[@"eigen system complex"] != nil) {
                complexFlag = [(solution.solutionInfo)[@"eigen system complex"] boolValue];
            } else complexFlag = NO;
            
            if ((solution.solutionInfo)[@"eigen system values"] != nil) {
                n = [(solution.solutionInfo)[@"eigen system values"] intValue];
                if (n > 0) {
                    solution.nOfEigenValues = n;
                    if (variableContainers->EigenValues == NULL) {
                        variableContainers->EigenValues = cdoublevec(0, n-1);
                        variableContainers->sizeEigenValues = n;
                        variableContainers->EigenVectors = cdoublematrix(0, n-1, 0, variableContainers->sizeValues-1);
                        variableContainers->size1EigenVectors = n;
                        variableContainers->size2EigenVectors = variableContainers->sizeValues;
                        
                        for (i=0; i<n; i++) {
                            variableContainers->EigenValues[i] = 0.0;
                            for (j=0; j<variableContainers->sizeValues; j++) {
                                variableContainers->EigenVectors[i][j] = 0.0;
                            }
                        }
                        
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
                    memset( matContainers->MassValues, 0.0, (matContainers->sizeMassValues*sizeof(matContainers->MassValues)) );
                }
            }
        } else if (harmonicAnal == YES) {
            n = 0;
            if ((solution.solutionInfo)[@"harmonic system values"] != nil) {
                n = [(solution.solutionInfo)[@"harmonic system values"] intValue];
            }
            if (n > 1) {
                listUtilities = [[FEMListUtilities alloc] init];
                found = [listUtilities listGetConstRealArray:model inArray:solution.valuesList forVariable:@"frequency" buffer:&freqv];
                if (found == YES) {
                    if (freqv.m < n) {
                        errorfunct("addEquationToSolution", "The solution option 'frequency' must be at least the same size as the Harmonic system values.");
                    }
                } else {
                    errorfunct("addEquationToSolution", "The solution option 'fequency' must be given for harmonic analysis.");
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
                
                for (i=0; i<n; i++) {
                    variableContainers->EigenValues[i] = 0.0;
                    for (j=0; j<variableContainers->sizeValues; j++) {
                        variableContainers->EigenVectors[i][j] = 0.0;
                    }
                }
                
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
            memset( matContainers->MassValues, 0.0, (matContainers->sizeMassValues*sizeof(matContainers->MassValues)) );
        }
    }
    
    if (solution.matrix != nil) {
        if ((solution.solutionInfo)[@"linear system symmetric"] != nil) {
            solution.matrix.symmetric = [(solution.solutionInfo)[@"linear system symmetric"] boolValue];
        }
        
        if ((solution.solutionInfo)[@"lumped mass matrix"] != nil) {
            solution.matrix.lumped = [(solution.solutionInfo)[@"lumped mass matrix"] boolValue];
        }
        
        multigridActive = ([(solution.solutionInfo)[@"linear system solver"] isEqualToString:@"multigrid"] == YES
                           || [(solution.solutionInfo)[@"linear system preconditioning"] isEqualToString:@"multigrid"] == YES) ? YES : NO;
        
        // Check for multigrid solver
        if (multigridActive == YES) {
            meshUtilities = [[FEMMeshUtils alloc] init];
            
            // Multigrid may be either solver or preconditioner, is it solver?
            if ((solution.solutionInfo)[@"linear system solver"] != nil) {
                solution.multigridSolution = ([(solution.solutionInfo)[@"linear system solver"] isEqualToString:@"multigrid"] == YES) ? YES : NO;
            }
            
            // There are four different methods: algebraic, cluster, p and geometric
            if ((solution.solutionInfo)[@"mg method"] != nil) {
                methodName = (solution.solutionInfo)[@"mg method"];
                mgAlgebraic = ([methodName isEqualToString:@"algebraic"] == YES || [methodName isEqualToString:@"cluster"] == YES
                               || [methodName isEqualToString:@"p"] == YES) ? YES : NO;
            } else {
                mgAlgebraic = ([(solution.solutionInfo)[@"mg algebraic"] boolValue] == YES || [(solution.solutionInfo)[@"mg cluster"] boolValue] == YES
                               || [(solution.solutionInfo)[@"mg pelement"] boolValue] == YES) ? YES : NO;
            }
            
            if ((solution.solutionInfo)[@"mg levels"] != nil) {
                mgLevels = [(solution.solutionInfo)[@"mg levels"] intValue];
                if (mgLevels < 1) mgLevels = 1;
            } else {
                if ((solution.solutionInfo)[@"multigrid levels"] != nil) {
                    mgLevels = [(solution.solutionInfo)[@"multigrid levels"] intValue];
                    if (mgLevels < 1) mgLevels = 1;
                } else {
                    if (mgAlgebraic == YES) {
                        mgLevels = 10;
                    } else {
                        errorfunct("addEquationToSolution", "'mg levels' must be defined for geometric multigrid.");
                    }
                }
            }
            solution.multiGridTotal = mgLevels;
            
            // In case of geometric multigrid, make the hierarchical meshes
            if (mgAlgebraic == NO) {
                // Check if h/2 splitting of mesh requested
                if ((solution.solutionInfo)[@"mg equal split"] != nil) {
                    solution.multigridEqualPlit = [(solution.solutionInfo)[@"mg equal split"] boolValue];
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
            [meshUtilities SetStabilizationParametersInMesh:solution.mesh model:model];
        }
        
        // Set the default verbosity of the iterative solvers accordingly with the global verbosity
        if ((solution.solutionInfo)[@"linear system residual output"] == nil) {
            kernel = [FEMKernel sharedKernel];
            k = 1;
            if ([kernel.outputLevelMask[4] boolValue] == NO) {
                k = 0;
            } else if ([kernel.outputLevelMask[5] boolValue] == NO) {
                k = 10;
            }
            if (k != 1) {
                aNumber = @(k);
                [solution.solutionInfo setObject:aNumber forKey:@"linear system residual output"];
            }
        }
    }
}

-(BOOL)isFileNameQualified:(NSString *)file {
    
    NSRange range;
    
    range = [file rangeOfString:@":"];
    return (range.location != NSNotFound || [file characterAtIndex:0] == '/' || [file characterAtIndex:0] == _backSlash);
}

-(NSMutableString *)nextFreeFileName:(NSString *)fileName0 suffix:(NSString *)suffix0 lastExisting:(BOOL *)lastExisting {
    
    int no;
    NSRange range;
    NSString *prefix, *suffix;
    NSMutableString *str, *fileName, *prevFileName;
    NSFileManager *fileManager;
    
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

@end
