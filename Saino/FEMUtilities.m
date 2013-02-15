//
//  FEMUtilities.m
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMUtilities.h"
#import "FEMEquation.h"
#import "FEMPElementMaps.h"
#import "FEMProjector.h"

@interface FEMUtilities ()
-(void)FEMUtils_applyProjector:(NSMutableArray *)variables model:(FEMModel *)aModel fromMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh projector:(FEMProjector *)projector;
-(void)FEMUtils_interpolateQMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh oldVariables:(NSMutableArray *)oldVar newVariables:(NSMutableArray *)newVar model:(FEMModel *)aModel quadrantTree:(BOOL *)useQuandrant withProjector:(BOOL)withProj projector:(FEMProjector *)projector mask:(NSString *)maskName nodesPresent:(BOOL *)foundNodes;
-(void)FEMUtils_interpolateMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh oldVariables:(NSMutableArray *)oldVar newVariables:(NSMutableArray *)newVar model:(FEMModel *)aModel quadrantTree:(BOOL *)useQuandrant withProjector:(BOOL)withProj projector:(FEMProjector *)projector mask:(NSString *)maskName;
@end

@implementation FEMUtilities

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
-(void)addVariableTo:(NSMutableArray *)anArray mesh:(FEMMesh *)aMesh solution:(FEMSolution *)aSolution name:(NSString *)name dofs:(int )dofs container:(variableArraysContainer *)aContainer ifOutput:(BOOL *)output ifSecondary:(BOOL *)secondary type:(int *)aType {
    
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
    
    varContainers->Values = aContainer->Values;
    varContainers->sizeValues = aContainer->sizeValues;
    
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

-(FEMVariable *)getVariableFrom:(NSMutableArray *)anArray model:(FEMModel *)aModel name:(NSString *)name onlySearch:(BOOL *)only maskName:(NSString *)maskName info:(BOOL *)found {
    
    int i, j, k, l, n, dofs;
    FEMVariable *var = nil, *pVar = nil, *tmp = nil;
    FEMListUtilities *listUtilities;
    FEMSolution *solution;
    FEMMesh *mesh, *currentMesh;
    FEMProjector *projector;
    variableArraysContainer *varContainers = NULL, *pvarContainers = NULL, *bufferContainers = NULL;
    NSMutableString *tmpname;
    NSNumber *aNumber;
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
        [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:pVar.name dofs:pVar.dofs container:varContainers ifOutput:&output ifSecondary:NULL type:NULL];
        
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
            bufferContainers = (variableArraysContainer*)malloc(sizeof(variableArraysContainer) * 1 );
            bufferContainers->Values = doublevec(0, (dofs/pVar.dofs)-1);
            k = 0;
            for (i=0; i<varContainers->sizeValues; i+=pVar.dofs) {
                bufferContainers->Values[k] = varContainers->Values[i];
                k++;
            }
            bufferContainers->sizeValues = dofs/pVar.dofs;
            bufferContainers->Perm = varContainers->Perm;
            bufferContainers->sizePerm = dofs/pVar.dofs;
            output = pVar.output;
            [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:@"velocity 1" dofs:1 container:bufferContainers ifOutput:&output ifSecondary:NULL type:NULL];
            
            free(bufferContainers);
            bufferContainers = NULL;
            onlyThis = YES;
            tmp = [self getVariableFrom:anArray model:aModel name:@"velocity 1" onlySearch:&onlyThis maskName:NULL info:&stat];
            bufferContainers = tmp.getContainers;
            if (varContainers->PrevValues != NULL) {
                bufferContainers->PrevValues = doublematrix(0,  (dofs/pVar.dofs)-1, 0, varContainers->size2PrevValues-1);
                k = 0;
                for (i=0; i<varContainers->sizeValues; i+=pVar.dofs) {
                    for (j=0; j<varContainers->size2PrevValues; j++) {
                        bufferContainers->PrevValues[k][j] = varContainers->PrevValues[i][j];
                    }
                    k++;
                }
                bufferContainers->size1PrevValues = dofs/pVar.dofs;
                bufferContainers->size2PrevValues = varContainers->size2PrevValues;
            }
            
            bufferContainers = NULL;
            bufferContainers = (variableArraysContainer*)malloc(sizeof(variableArraysContainer) * 1 );
            bufferContainers->Values = doublevec(0, (dofs/pVar.dofs)-1);
            k = 0;
            for (i=1; i<varContainers->sizeValues; i+=pVar.dofs) {
                bufferContainers->Values[k] = varContainers->Values[i];
                k++;
            }
            bufferContainers->sizeValues = dofs/pVar.dofs;
            bufferContainers->Perm = varContainers->Perm;
            bufferContainers->sizePerm = dofs/pVar.dofs;
            [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:@"velocity 2" dofs:1 container:bufferContainers ifOutput:&output ifSecondary:NULL type:NULL];
            
            free(bufferContainers);
            bufferContainers = NULL;
            onlyThis = YES;
            tmp = [self getVariableFrom:anArray model:aModel name:@"velocity 2" onlySearch:&onlyThis maskName:NULL info:&stat];
            bufferContainers = tmp.getContainers;
            if (varContainers->PrevValues != NULL) {
                bufferContainers->PrevValues = doublematrix(0,  (dofs/pVar.dofs)-1, 0, varContainers->size2PrevValues-1);
                k = 0;
                for (i=1; i<varContainers->sizeValues; i+=pVar.dofs) {
                    for (j=0; j<varContainers->size2PrevValues; j++) {
                        bufferContainers->PrevValues[k][j] = varContainers->PrevValues[i][j];
                    }
                    k++;
                }
                bufferContainers->size1PrevValues = dofs/pVar.dofs;
                bufferContainers->size2PrevValues = varContainers->size2PrevValues;
            }
            
            bufferContainers = NULL;
            bufferContainers = (variableArraysContainer*)malloc(sizeof(variableArraysContainer) * 1 );
            bufferContainers->Values = doublevec(0, (dofs/pVar.dofs)-1);
            k = 0;
            for (i=2; i<varContainers->sizeValues; i+=pVar.dofs) {
                bufferContainers->Values[k] = varContainers->Values[i];
                k++;
            }
            bufferContainers->sizeValues = dofs/pVar.dofs;
            bufferContainers->Perm = varContainers->Perm;
            bufferContainers->sizePerm = dofs/pVar.dofs;
            if (pVar.dofs == 3) {
                [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:@"pressure" dofs:1 container:bufferContainers ifOutput:&output ifSecondary:NULL type:NULL];
            } else {
                [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:@"velocity 3" dofs:1 container:bufferContainers ifOutput:&output ifSecondary:NULL type:NULL];
                
                free(bufferContainers);
                bufferContainers = NULL;
                onlyThis = YES;
                tmp = [self getVariableFrom:anArray model:aModel name:@"velocity 3" onlySearch:&onlyThis maskName:NULL info:&stat];
                bufferContainers = tmp.getContainers;
                if (varContainers->PrevValues != NULL) {
                    bufferContainers->PrevValues = doublematrix(0,  (dofs/pVar.dofs)-1, 0, varContainers->size2PrevValues-1);
                    k = 0;
                    for (i=2; i<varContainers->sizeValues; i+=pVar.dofs) {
                        for (j=0; j<varContainers->size2PrevValues; j++) {
                            bufferContainers->PrevValues[k][j] = varContainers->PrevValues[i][j];
                        }
                        k++;
                    }
                    bufferContainers->size1PrevValues = dofs/pVar.dofs;
                    bufferContainers->size2PrevValues = varContainers->size2PrevValues;
                }
                
                bufferContainers = NULL;
                bufferContainers = (variableArraysContainer*)malloc(sizeof(variableArraysContainer) * 1 );
                bufferContainers->Values = doublevec(0, (dofs/pVar.dofs)-1);
                k = 0;
                for (i=3; i<varContainers->sizeValues; i+=pVar.dofs) {
                    bufferContainers->Values[k] = varContainers->Values[i];
                    k++;
                }
                bufferContainers->sizeValues = dofs/pVar.dofs;
                bufferContainers->Perm = varContainers->Perm;
                bufferContainers->sizePerm = dofs/pVar.dofs;
                [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:@"pressure" dofs:1 container:bufferContainers ifOutput:&output ifSecondary:NULL type:NULL];
            }
            
            free(bufferContainers);
            bufferContainers = NULL;
            onlyThis = YES;
            tmp = [self getVariableFrom:anArray model:aModel name:@"pressure" onlySearch:&onlyThis maskName:NULL info:&stat];
            bufferContainers = tmp.getContainers;
            if (varContainers->PrevValues != NULL) {
                bufferContainers->PrevValues = doublematrix(0,  (dofs/pVar.dofs)-1, 0, varContainers->size2PrevValues-1);
                k = 0;
                for (i=pVar.dofs-1; i<varContainers->sizeValues; i+=pVar.dofs) {
                    for (j=0; j<varContainers->size2PrevValues; j++) {
                        bufferContainers->PrevValues[k][j] = varContainers->PrevValues[i][j];
                    }
                    k++;
                }
                bufferContainers->size1PrevValues = dofs/pVar.dofs;
                bufferContainers->size2PrevValues = varContainers->size2PrevValues;
            }
        } else {
            output = pVar.output;
            if (pVar.dofs > 1) {
                for (i=0; i<pVar.dofs; i++) {
                    bufferContainers = NULL;
                    bufferContainers = (variableArraysContainer*)malloc(sizeof(variableArraysContainer) * 1 );
                    bufferContainers->Values = doublevec(0, (dofs/pVar.dofs)-1);
                    k = 0;
                    for (j=i; j<varContainers->sizeValues; j+=pVar.dofs) {
                        bufferContainers->Values[k] = varContainers->Values[j];
                        k++;
                    }
                    bufferContainers->sizeValues = dofs/pVar.dofs;
                    bufferContainers->Perm = varContainers->Perm;
                    bufferContainers->sizePerm = dofs/pVar.dofs;
                    tmpname = [NSMutableString stringWithString:pVar.name];
                    aNumber = [NSNumber numberWithInt:i+1];
                    [tmpname appendString:@" "];
                    [tmpname appendString:[aNumber stringValue]];
                    [self addVariableTo:anArray mesh:pVar.primaryMesh solution:pVar.solution name:tmpname dofs:1 container:bufferContainers ifOutput:&output ifSecondary:NULL type:NULL];
                    
                    free(bufferContainers);
                    bufferContainers = NULL;
                    onlyThis = YES;
                    tmp = [self getVariableFrom:anArray model:aModel name:tmpname onlySearch:&onlyThis maskName:NULL info:&stat];
                    bufferContainers = tmp.getContainers;
                    if (varContainers->PrevValues != NULL) {
                        bufferContainers->PrevValues = doublematrix(0,  (dofs/pVar.dofs)-1, 0, varContainers->size2PrevValues-1);
                        k = 0;
                        for (j=i; j<varContainers->sizeValues; j+=pVar.dofs) {
                            for (l=0; l<varContainers->size2PrevValues; l++) {
                                bufferContainers->PrevValues[k][l] = varContainers->PrevValues[j][l];
                            }
                            k++;
                        }
                        bufferContainers->size1PrevValues = dofs/pVar.dofs;
                        bufferContainers->size2PrevValues = varContainers->size2PrevValues;
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
            tmpname = [NSMutableString stringWithString:pVar.name];
            aNumber = [NSNumber numberWithInt:i+1];
            [tmpname appendString:@" "];
            [tmpname appendString:[aNumber stringValue]];
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
            tmpname = [NSMutableString stringWithString:pVar.name];
            aNumber = [NSNumber numberWithInt:i+1];
            [tmpname appendString:@" "];
            [tmpname appendString:[aNumber stringValue]];
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
    
    int comp;
    NSString *str;
    NSMutableString *str1;
    NSRange location;
    
    location = [string rangeOfString:@"["];
    
    comp = 0;
    if (component != NULL) comp = *component;
    
    if (location.location == NSNotFound) {
        str1 = [NSMutableString stringWithString:string];
        if (comp > 0) {
            [str1 appendString:@" "];
            [str1 appendString:[NSString stringWithFormat:@"%d",comp]];
            str = [NSString stringWithString:str1];
        }
    } else {
        //TODO: Implement this particular part
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

@end
