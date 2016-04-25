//
//  FEMInterpolateMeshToMesh.m
//  Saino
//
//  Created by Seddik hakime on 12/11/2015.
//  Copyright © 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMInterpolateMeshToMesh.h"
#import "FEMCore.h"
#import "FEMSolution.h"
#import "FEMNumericIntegration.h"
#import "FEMMeshUtils.h"
#import "FEMUtilities.h"
#import "FEMInterpolation.h"
#import "FEMElementDescription.h"
#import "FEMListMatrix.h"

@interface FEMInterpolateMeshToMesh ()
-(void)FEMInterpolateMeshToMesh_applyProjector:(NSMutableArray * __nonnull)variables model:(FEMModel * __nonnull)model fromMesh:(FEMMesh * __nonnull)oldMesh toMesh:(FEMMesh * __nonnull)newMesh projector:(FEMProjector * __nonnull)projector utilities:(FEMUtilities * __nonnull)utilities;
@end

@implementation FEMInterpolateMeshToMesh

/************************************************************
 
    Method corresponds to Elmer from git on October 27 2015
 
************************************************************/
-(void)FEMInterpolateMeshToMesh_applyProjector:(NSMutableArray * __nonnull)variables model:(FEMModel * __nonnull)model fromMesh:(FEMMesh * __nonnull)oldMesh toMesh:(FEMMesh * __nonnull)newMesh projector:(FEMProjector * __nonnull)projector utilities:(FEMUtilities * __nonnull)utilities {
    
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
        if (variable.isSecondary == YES) continue;
        
        NSString *string = [variable.name substringToIndex:10];
        if (variable.dofs == 1 && [string isEqualToString:@"coordinate"] == NO) {
            only = YES;
            oldSol = [utilities getVariableFrom:oldMesh.variables model:model name:variable.name onlySearch:&only maskName:NULL info:&stat];
            newSol = [utilities getVariableFrom:newMesh.variables model:model name:variable.name onlySearch:&only maskName:NULL info:&stat];
            
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
                    for (j=0; j<newContainers->size1PrevValues; j++) {
                        newContainers->PrevValues[j][i] = bf2[j];
                    }
                }
                free_dvector(bf1, 0, oldContainers->size1PrevValues-1);
                free_dvector(bf2, 0, newContainers->size1PrevValues-1);
            }
        }
    }
}

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
    }
    
    return self;
}

/********************************************************************************************************
 
    Interpolates values of all variables from a mesh associated with the old model to the mesh
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
    (int *)newMaskPerm         ->  Mask the new variable set by the given maskName when trying to define
                                   the interpolation
 
 
    Method corresponds to Elmer from git on October 27 2015

********************************************************************************************************/
-(void)interpolateQMesh:(FEMMesh * __nonnull)oldMesh toMesh:(FEMMesh * __nonnull)newMesh oldVariables:(NSMutableArray * __nullable)oldVar newVariables:(NSMutableArray * __nullable)newVar model:(FEMModel * __nonnull)model quadrantTree:(BOOL * __nullable)useQuandrant projector:(FEMProjector * __nullable)projector mask:(NSString * __nullable)maskName nodesPresent:(BOOL * __nullable)nodesPresent newMaskPerm:(int * __nullable)newMaskPerm {
    
    int i, j, k, l=-1, n, bfId, qTreeFails, np, nrow, totFails;
    int dim, epsTries, passiveCoordinate;
    Nodes_t *elementNodes = NULL;
    double point[3], localCoordinates[3];
    FEMVariable *oldSol, *newSol;
    double *elementValues, *newValues, **rotBasis = NULL, **wBasis = NULL;
    Element_t *element = NULL, *oldElements = NULL;
    Nodes_t *oldNodes= NULL, *newNodes = NULL, *nodes = NULL;
    Quadrant_t *leafQuadrant = NULL, *rootQuadrant = NULL;
    double boundingBox[6], *valuesPtr = NULL, u, v, w;
    BOOL edgeBasis, found, piolaT, searchOnly, stat, tryLinear, tryQtree, useQTree;
    int *rows = NULL, *cols = NULL;
    int *rInd = NULL;
    BOOL maskExists, all;
    double eps1 = 0.1, eps2, epsGlobal, epsLocal, maxVal, numericEps;
    double *values = NULL, *localU = NULL, *localV = NULL, *localW = NULL;
    variableArraysContainer *varContainers = NULL, *newVarContainers = NULL, *oldVarContainers = NULL;
    matrixArraysContainer *matContainers = NULL, *tmatContainers = NULL;
    FEMNumericIntegration *integration;
    
    typedef struct Epntr_t {
        Element_t *element;
    } Epntr_t;
    Epntr_t *elementPointers = NULL;
    
    FEMCore *core = [FEMCore sharedCore];
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    FEMInterpolation *interpolation = [[FEMInterpolation alloc] init];
    FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    FEMPElementMaps *elementMaps = [[FEMPElementMaps alloc] init];
    
    //TODO: add support for parallel run
    BOOL parallel = NO;
    
    // If projector argument given, search for existing projector matrix or
    // generate new projector if not already there
    if (projector != nil) {
        for (FEMProjector *projector in newMesh.projectors) {
            if ([projector.mesh isEqualTo:oldMesh]) {
                if (oldVar != nil) [self FEMInterpolateMeshToMesh_applyProjector:oldVar model:model fromMesh:oldMesh toMesh:newMesh projector:projector utilities:utilities];
                return;
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
    dim = model.dimension;
    
    if (useQuandrant == NULL) {
        useQTree = YES;
    } else {
        useQTree = *useQuandrant;
    }
    
    if (useQTree == YES) {
        if (rootQuadrant == NULL) {
            oldNodes = oldMesh.getNodes;
            
            vDSP_minvD(oldNodes->x, 1, &boundingBox[0], oldMesh.numberOfNodes);
            vDSP_minvD(oldNodes->y, 1, &boundingBox[1], oldMesh.numberOfNodes);
            vDSP_minvD(oldNodes->z, 1, &boundingBox[2], oldMesh.numberOfNodes);
            vDSP_maxvD(oldNodes->x, 1, &boundingBox[3], oldMesh.numberOfNodes);
            vDSP_maxvD(oldNodes->y, 1, &boundingBox[4], oldMesh.numberOfNodes);
            vDSP_maxvD(oldNodes->z, 1, &boundingBox[5], oldMesh.numberOfNodes);
            
            maxVal = -HUGE_VAL;
            k = 3;
            valuesPtr = &boundingBox[0];
            for (i=0; i<3; i++) {
                if ( (*(valuesPtr + k) - *(valuesPtr + i)) > maxVal ) {
                    maxVal = *(valuesPtr + k) - *(valuesPtr + i);
                }
                k++;
            }
            eps2 = 0.1 * maxVal;
            for (i=0; i<3; i++) {
                *(valuesPtr + i) = *(valuesPtr + i) - eps2;
            }
            for (i=3; i<6; i++) {
                *(valuesPtr + i) = *(valuesPtr + i) + eps2;
            }
            
             // Create mother of all quadrants
            Quadrant_t *quadrant = (Quadrant_t*) malloc(sizeof(Quadrant_t));
            [interpolation buildQuadrantTreeForMesh:oldMesh model:model boundingBox:boundingBox rootQuadrant:quadrant];
            [oldMesh assignQuadrant:quadrant];
            rootQuadrant = oldMesh.getQuadrant;
        }
    }
    
    // If we use a mask
    maskExists = (maskName != nil) ? YES : NO;
    
    n = oldMesh.maxElementNodes;
    elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(elementNodes);
    elementNodes->x = doublevec(0, n-1);
    elementNodes->y = doublevec(0, n-1);
    elementNodes->z = doublevec(0, n-1);
    
    elementValues = doublevec(0, n-1);
    newValues = doublevec(0, newMesh.numberOfNodes-1);
    
    epsGlobal = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"interpolation global epsilon" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) epsGlobal = 2.0e-10;
    epsLocal = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"interpolation local epsilon" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) epsLocal = 1.0e-10;
    epsTries = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"interpolation maximum iterations" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) epsTries = 12;
    
    numericEps = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"interpolation numeric epsilon" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) numericEps = 1.0e-10;
    
    FEMSolution *solution = (FEMSolution *)model.solution;
    passiveCoordinate = 0;
    if (solution != nil) {
        if ((solution.solutionInfo)[@"interpolation passive coordinate"] != nil) {
            passiveCoordinate = [(solution.solutionInfo)[@"interpolation passive coordinate"] intValue];
        }
    }
    
    qTreeFails = 0;
    totFails = 0;
    
    edgeBasis = NO;
    if (solution != nil) {
        if ((solution.solutionInfo)[@"edge basis"] != nil) {
            edgeBasis = [(solution.solutionInfo)[@"edge basis"] boolValue];
        }
    }
    
    piolaT = NO;
    if (edgeBasis == YES) {
        if ((solution.solutionInfo)[@"use piola transform"] != nil) {
            piolaT = [(solution.solutionInfo)[@"use piola transform"] boolValue];
        }
    }
    
    tryLinear = [listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"try linear search if qtree fails" info:&found];
    if (found == NO) tryLinear = YES;
    
    newNodes = newMesh.getNodes;
    oldElements = oldMesh.getElements;
    
    searchOnly = YES;
    int foundCnt = 0;
    
    // Loop over all nodes in the new mesh
    for (i=0; i<newMesh.numberOfNodes; i++) {
        
        element = NULL;
        // Only get the variable for the requested nodes
        if (newMaskPerm != NULL) {
            if (newMaskPerm[i] < 0) continue;
        }
        
        point[0] = newNodes->x[i];
        point[1] = newNodes->y[i];
        point[2] = newNodes->z[i];
        
        if (passiveCoordinate != 0) {
            point[passiveCoordinate-1] = 0.0;
        }
        
        // Find in which old mesh bulk element the point belongs to
        found = NO;
        tryQtree = (rootQuadrant != NULL && useQTree == YES) ? YES : NO;
        
        if (tryQtree) {
            // Find the last existing quadrant that the point belongs to
            leafQuadrant = NULL;
            leafQuadrant = [interpolation findLeafElementsForPoint:point dimension:dim rootQuadrant:rootQuadrant];
            
            if (leafQuadrant != NULL) {
                FEMBodyForce *bodyForceAtID;
                // Go through the bulk elements in the last child quadrant only.
                // Try to find matching element with progressively sloppier tests.
                // Allow at most 100% of slack
                eps1 = epsGlobal;
                eps2 = epsLocal;
                
                for (j=0; j<epsTries; j++) {
                    for (k=0; k<leafQuadrant->nElementsInQuadrant; k++) {
                        
                        element = &oldElements[leafQuadrant->elements[k]];
                        
                        if (maskExists == YES) {
                            if ((model.bodies)[oldElements[leafQuadrant->elements[k]].BodyID-1][@"body force"] == nil) continue;
                            bfId = [(model.bodies)[oldElements[leafQuadrant->elements[k]].BodyID][@"body force"] intValue];
                            bodyForceAtID = (model.bodyForces)[bfId-1];
                            if ([listUtilities listCheckPresentVariable:maskName inArray:bodyForceAtID.valuesList] == NO) continue;
                        }
                        
                        n = oldElements[leafQuadrant->elements[k]].Type.NumberOfNodes;
                        for (l=0; l<n; l++) {
                            elementNodes->x[l] = oldNodes->x[oldElements[leafQuadrant->elements[k]].NodeIndexes[l]];
                            elementNodes->y[l] = oldNodes->y[oldElements[leafQuadrant->elements[k]].NodeIndexes[l]];
                            elementNodes->z[l] = oldNodes->z[oldElements[leafQuadrant->elements[k]].NodeIndexes[l]];
                        }
                        found = [interpolation isPointInElement:&oldElements[leafQuadrant->elements[k]] elementNodes:elementNodes point:point localCoordinates:localCoordinates globalEpsilon:&eps1 localEpsilon:&eps2 numericEpsilon:&numericEps globalDistance:NULL localDistance:NULL model:model elementDescription:elementDescription elementMaps:elementMaps edgeBasis:&piolaT];
                        if (found == YES) break;
                    }
                    if (found == YES) break;
                    
                    eps1 = 10.0 * eps1;
                    eps2 = 10.0 * eps2;
                    if (eps1 > 1.0) break;
                }
            }
        }
        
        if (tryQtree == NO || (found == NO && parallel == NO && tryLinear == YES)) {
            
            // Go through all old mesh bulk elements
            for (k=0; k<oldMesh.numberOfBulkElements; k++) {
                element = &oldElements[k];
                n = oldElements[k].Type.NumberOfNodes;
                
                for (l=0; l<n; l++) {
                    elementNodes->x[l] = oldNodes->x[oldElements[k].NodeIndexes[l]];
                    elementNodes->y[l] = oldNodes->y[oldElements[k].NodeIndexes[l]];
                    elementNodes->z[l] = oldNodes->z[oldElements[k].NodeIndexes[l]];
                }
                
                found = [interpolation isPointInElement:&oldElements[k] elementNodes:elementNodes point:point localCoordinates:localCoordinates globalEpsilon:NULL localEpsilon:NULL numericEpsilon:NULL globalDistance:NULL localDistance:NULL model:model elementDescription:elementDescription elementMaps:elementMaps edgeBasis:NULL];
                if (found == YES) {
                    if (tryQtree == YES) qTreeFails++;
                    break;
                }
            }
        }
        
        if (found == NO) {
            element = NULL;
            if (parallel == NO) {
                fprintf(stdout, "FEMInterpolateMeshToMesh:interpolateQMesh: point %d was not found in any of the elements.\n", i);
                totFails++;
            }
            continue;
        }
        if (nodesPresent != NULL) nodesPresent[i] = YES;
        
        // Found element in oldModel
        if (projector != nil) {
            foundCnt++;
            elementPointers[i].element = element;
            localU[i] = localCoordinates[0];
            localV[i] = localCoordinates[1];
            localW[i] = localCoordinates[2];
        }
        
        if (oldVar == nil || projector != nil) continue;
        
        // Go through all variables to be interpolated
        for (FEMVariable *variable in oldVar) {
            
            varContainers = variable.getContainers;
            if (varContainers->sizeValues == variable.dofs) continue;
            
            if (variable.isSecondary == YES) continue;
            
            newSol = nil;
            oldSol = nil;
            NSString *string = [variable.name substringToIndex:10];
            if (variable.dofs == 1 && [string isEqualToString:@"coordinate"] == NO) {
                
                // Interpolate variable at point in element
                newSol = [utilities getVariableFrom:newMesh.variables model:model name:variable.name onlySearch:&searchOnly maskName:NULL info:&found];
                if (newSol == nil) continue;
                
                oldSol = [utilities getVariableFrom:oldMesh.variables model:model name:variable.name onlySearch:&searchOnly maskName:NULL info:&found];
                
                newVarContainers = newSol.getContainers;
                oldVarContainers = oldSol.getContainers;
                
                // Check that the node was found in the old mesh
                if (element != NULL) {
                    //Check for rounding errors
                    all = YES;
                    for (k=0; k<element->Type.NumberOfNodes; k++) {
                        if (oldVarContainers->Perm[element->NodeIndexes[k]] >= 0) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) {
                        if (newVarContainers->Perm[i] >= 0.0) {
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
                    if (newVarContainers->Perm[i] >= 0) newValues[newVarContainers->Perm[i]] = 0.0;
                }
            }
        }
    }
    
    if (parallel == NO) {
        if (qTreeFails > 0) {
            fprintf(stdout, "FEMInterpolateMeshToMesh:interpolateQMesh: number of points not found in quadtree: %d.", qTreeFails);
            if (totFails == 0) {
                fprintf(stdout, "FEMInterpolateMeshToMesh:interpolateQMesh: all nodes still found by N^2 dummy search.\n");
            }
        }
        if (totFails == 0) {
            fprintf(stdout, "FEMInterpolateMeshToMesh:interpolateQMesh: found all nodes in the target nesh.\n");
        } else {
            fprintf(stdout, "FEMInterpolateMeshToMesh:interpolateQMesh: points not found: %d (found %d).\n", totFails, newMesh.numberOfNodes-totFails);
        }
    }
    
    // Construct mesh projector if required. Next time around will use
    // the existing projector to interpolate values
    if (projector != nil) {
        n = newMesh.numberOfNodes;
        int *indexes = intvec(0, 99);
        
        // The critical value of basis function that is accepted to the
        // projector. Note that the sum of weights is one, so we
        // know the scale for this one
        double eps_basis = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"interpolation basis epsilon" info:&found minValue:NULL maxValue:NULL];
        if (found == NO) eps_basis = 0.0e-12;
        
        rows = intvec(0, (n+1)-1);
        rows[0] = 0;
        BOOL projectorAllocated = NO;
        
    label:
        nrow = 0;
        for (i=0; i<n; i++) {
            element = NULL;
            if (edgeBasis == YES && oldMesh.parent != nil) {
                Element_t *faces = oldMesh.parent.getFaces;
                element = &faces[elementPointers[i].element->ElementIndex-1];
                if (faces[elementPointers[i].element->ElementIndex-1].BoundaryInfo != NULL) {
                    Element_t *parent = faces[elementPointers[i].element->ElementIndex-1].BoundaryInfo->Left;
                    if (parent != NULL) {
                        k = faces[elementPointers[i].element->ElementIndex-1].Type.NumberOfNodes;
                        np = parent->Type.NumberOfNodes;
                    }
                }
            } else {
                element = elementPointers[i].element;
            }
            found = (element != NULL) ? YES : NO;
            
            if (found == NO) {
                // It seems unnevessary to make a matrix entry in case no target element is found
                if (/* DISABLES CODE */ (NO)) {
                    if (projectorAllocated == YES) {
                        cols[nrow] = 0;
                        values[nrow] = 0.0;
                    }
                    nrow++;;
                }
            } else {
                u = localU[i];
                v = localV[i];
                w = localW[i];
                if (nodes == NULL) {
                    nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
                    initNodes(nodes);
                }
                if (edgeBasis == YES) {
                    [core getNodes:solution model:model inElement:element resultNodes:nodes numberOfNodes:NULL mesh:nil];
                } else {
                    [core getNodes:solution model:model inElement:element resultNodes:nodes numberOfNodes:NULL mesh:oldMesh];
                }
                
                BOOL notDG = (solution != nil) ? YES : NO;
                memset( indexes, -1, 100*sizeof(int) );
                k = [core getElementDofsSolution:solution model:model forElement:element atIndexes:indexes disableDiscontinuousGalerkin:&notDG];
                
                np = [core getNumberOfNodesForElement:element];
                BOOL any = NO;
                for (j=0; j<np; j++) {
                    if (indexes[j] > element->NodeIndexes[j]) {
                        any = YES;
                        break;
                    }
                }
                if (any == YES) np = 0;
                
                if (integration == nil) {
                    integration = [[FEMNumericIntegration alloc] init];
                    if ([integration allocation:oldMesh] == NO) fatal("FEMInterpolateMeshToMesh:interpolateQMesh", "Allocation error in FEMNumericIntegration.");
                }
                
                if (edgeBasis == YES) {
                    wBasis = doublematrix(0, k-1, 0, 2);
                    rotBasis = doublematrix(0, k-1, 0, 2);
                    if (piolaT == YES) {
                        //TODO: Needs to be implemented because this requires changes in FEMNumericIntegration
                        // In Elmer: ElementInfo(Element,Nodes,u,v,w,detJ,Vals,EdgeBasis=WBasis )
                    } else {
                        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:oldMesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                        stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:oldMesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:u withBubbles:NO basisDegree:NULL];
                        [elementDescription getEdgeBasisElement:element wBasis:wBasis rotWBasis:rotBasis basis:integration.basis dBasisdx:integration.basisFirstDerivative];
                        
                    }
                } else {
                    stat = [integration setBasisForElement:element elementNodes:nodes inMesh:oldMesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                }
                
                double rowSum = 0.0;
                for (j=0; j<k; j++) {
                    if (j < np) rowSum = rowSum + integration.basis[j];
                    if (projectorAllocated == NO) {
                        if (edgeBasis == NO || (edgeBasis == YES && j < np)) {
                            nrow++;
                        } else {
                            nrow = nrow + 3;
                        }
                    }
                }
                
                if (projectorAllocated == YES) {
                    for (j=0; j<k; j++) {
                        if (edgeBasis == NO) rInd[indexes[j]] = rInd[indexes[j]] + 1;
                        
                        // Always normalize the weights to one!
                        if (edgeBasis == NO || (edgeBasis == YES && j < np)) {
                            cols[nrow] = indexes[j];
                            values[nrow] = integration.basis[j] / rowSum;
                            nrow++;
                        } else {
                            cols[nrow] = -indexes[j];
                            values[nrow] = wBasis[j-np][0];
                            nrow++;
                            cols[nrow] = -indexes[j];
                            values[nrow] = wBasis[j-np][1];
                            nrow++;
                            cols[nrow] = -indexes[j];
                            values[nrow] = wBasis[j-np][2];
                            nrow++;
                        }
                    }
                }
                if (edgeBasis == YES) {
                    free_dmatrix(wBasis, 0, k-1, 0, 2);
                    free_dmatrix(rotBasis, 0, k-1, 0, 2);
                }
            }
            rows[i+1] = nrow;
        }
        
        free_dvector(nodes->x, 0, nodes->numberOfNodes-1);
        free_dvector(nodes->y, 0, nodes->numberOfNodes-1);
        free_dvector(nodes->z, 0, nodes->numberOfNodes-1);
        free(nodes);
        nodes = NULL;
        [integration deallocation:oldMesh];
        integration = nil;
        
        if (projectorAllocated == NO) {
            cols = intvec(0, rows[n]-1);
            values = doublevec(0, rows[n]-1);
            memset( cols, 0, rows[n]*sizeof(int) );
            memset( values, 0.0, rows[n]*sizeof(double) );

            projector.matrix.numberOfRows = n;
            matContainers = projector.matrix.getContainers;
            matContainers->Rows = rows;
            matContainers->sizeRows = n+1;
            matContainers->Cols = cols;
            matContainers->sizeCols = rows[n];
            matContainers->Values = values;
            matContainers->sizeValues = rows[n];
            
            projector.mesh = oldMesh;
            [newMesh.projectors addObject:projector];
            
            if (edgeBasis == NO) {
                rInd = intvec(0, oldMesh.numberOfNodes-1);
                memset( rInd, 0, oldMesh.numberOfNodes*sizeof(int) );
            }
            
            projectorAllocated = YES;
            goto label;
        }
        
        free_ivector(indexes, 0, 99);
        free(elementPointers);
        free_dvector(localU, 0, newMesh.numberOfNodes-1);
        free_dvector(localV, 0, newMesh.numberOfNodes-1);
        free_dvector(localW, 0, newMesh.numberOfNodes-1);
        
        // Store also the transpose of the projector
        if (edgeBasis == NO) {
            if (found == YES) {
                n = oldMesh.numberOfNodes;
                // Needed for some matrices
                matContainers = projector.matrix.getContainers;
                n = max(n, max_array(matContainers->Cols, matContainers->sizeCols));
                
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
                tmatContainers->sizeRows = (n+1);
                tmatContainers->Cols = cols;
                tmatContainers->sizeCols = rows[n];
                tmatContainers->Values = values;
                tmatContainers->sizeValues = rows[n];
                
                memset( rInd, 0, oldMesh.numberOfNodes*sizeof(int) );
                for (i=0; i<projector.matrix.numberOfRows; i++) {
                    for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                        k = matContainers->Cols[j];
                        l = rows[k] + rInd[k];
                        rInd[k] = rInd[k] + 1;
                        cols[l] = i;
                        values[l] = matContainers->Values[j];
                    }
                }
            }
            free_ivector(rInd, 0, oldMesh.numberOfNodes-1);
        }
        
        if (oldVar != nil) [self FEMInterpolateMeshToMesh_applyProjector:oldVar model:model fromMesh:oldMesh toMesh:newMesh projector:projector utilities:utilities];
    }
    
    free_dvector(elementNodes->x, 0, oldMesh.maxElementNodes-1);
    free_dvector(elementNodes->y, 0, oldMesh.maxElementNodes-1);
    free_dvector(elementNodes->z, 0, oldMesh.maxElementNodes-1);
    free(elementNodes);
    
    free_dvector(elementValues, 0, oldMesh.maxElementNodes-1);
    free_dvector(newValues, 0,  newMesh.numberOfNodes-1);
    [elementMaps deallocation];
}

/*******************************************************************************
 
    Map results from mesh to mesh. The from-Mesh is stored in an octree from
    which it is relatively fast to find the to-nodes. When the node is found
    interpolation is performed. Optionally there may be an existing projector
    that speeds up the interpolation.
 
    Method corresponds to Elmer from git on October 27 2015

*******************************************************************************/
-(void)interpolateMesh:(FEMMesh * __nonnull)oldMesh toMesh:(FEMMesh * __nonnull)newMesh oldVariables:(NSMutableArray * __nullable)oldVar newVariables:(NSMutableArray * __nullable)newVar model:(FEMModel * __nonnull)model quadrantTree:(BOOL * __nullable)useQuandrant projector:(FEMProjector * __nullable)projector mask:(NSString * __nullable)maskName unfoundNodes:(BOOL * __nullable)unfoundNodes {
    
    BOOL *foundNodes = (BOOL*)malloc(sizeof(BOOL) * newMesh.numberOfNodes );
    memset( foundNodes, 0, newMesh.numberOfNodes*sizeof(BOOL) );
    
    [self interpolateQMesh:oldMesh toMesh:newMesh oldVariables:oldVar newVariables:newVar model:model quadrantTree:useQuandrant projector:projector mask:maskName nodesPresent:foundNodes newMaskPerm:NULL];
    
    // unFoundNodes should allocated by the caller
    if (unfoundNodes != NULL) {
        for (int i=0; i<newMesh.numberOfNodes; i++) {
            unfoundNodes[i] = !foundNodes[i];
        }
    }
    free(foundNodes);
    
    // TODO: Add support for parallel run
    // The rest of this method deals with a parallel mesh
}

/***********************************************************************************
 
    Create a projector for mapping between interfaces using the Galerkin method
    A temporal mesh structure with a node for each Gaussian integration point is
    created. The this projector matrix is transferred to a projector on the nodal
    coordinates.
 
    Method corresponds to Elmer from git on October 27 2015
 
***********************************************************************************/
-(FEMMatrix * __nonnull)weightedProjectorMesh2:(FEMMesh * __nonnull)bMesh2 mesh1:(FEMMesh * __nonnull)bMesh1 inversePermutation2:(int * __nonnull)invPerm2 sizeInversePermutation2:(int)sizeInversePermutation2 inversePermutation1:(int * __nonnull)invPerm1 sizeInversePermutation1:(int)sizeInversePermutation1 useQuadrantTree:(BOOL)useQuadrantTree repeating:(BOOL)repeating antiRepeating:(BOOL)antiRepeating periodicScale:(double)periodicScale nodalJump:(BOOL)nodalJump model:(FEMModel * __nonnull)model {
    
    int n, np, eqindsize, *indexes, pp, qq;
    double **rotWBasis = NULL, u, v, w, **wBasis = NULL, vq[3], wq[3];
    FEMMatrix *gaussProjector, *projector;
    Element_t *bMesh1Elements = NULL, *bMesh2Elements = NULL, *element = NULL;
    Nodes_t *nodes = NULL;
    GaussIntegrationPoints *IP = NULL;
    BOOL axisSym, stat;
    
    fprintf(stdout, "FEMInterpolateMeshToMesh:weightedProjectorMesh2: creating Galerkin projector between two boundaries.\n");
    
    FEMCore *core = [FEMCore sharedCore];
    FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    
    indexes = intvec(0, 99);
    
    Nodes_t *gaussNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    Nodes_t *elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    Nodes_t *realNodes = bMesh1.getNodes;
    int nbNodes = bMesh1.numberOfNodes;
    
    BOOL edgeBasis = NO;
    FEMSolution *solution = (FEMSolution *)model.solution;
    if (solution != nil) {
        if ((solution.solutionInfo)[@"edge basis"] != nil) {
            edgeBasis = [(solution.solutionInfo)[@"edge basis"] boolValue];
        }
    }
    
    BOOL piolaT = NO;
    if (edgeBasis == YES) {
        if ((solution.solutionInfo)[@"use piola transform"] != nil) {
            piolaT = [(solution.solutionInfo)[@"use piola transform"] boolValue];
        }
        fprintf(stdout, "FEMInterpolateMeshToMesh:weightedProjectorMesh2: accounting for edge elements in projector.\n");
    }
    
    int relOrder = 0;
    if (solution != nil) {
        if((solution.solutionInfo)[@"projector relative integration order"] != nil) {
            relOrder = [(solution.solutionInfo)[@"projector relative integration order"] intValue];
        }
    }
    if (relOrder > 1) relOrder = 1;
    if (relOrder < -1) relOrder = -1;
    
    // Calculate the total number of Gaussian integration points
    // and allocate space for the node structures
    int nbGaussPoints = 0;
    bMesh1Elements = bMesh1.getElements;
    bMesh2Elements = bMesh2.getElements;
    FEMNumericIntegration *integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:bMesh1] == NO) fatal("FEMInterpolateMeshToMesh:weightedProjectorMesh2", "Allocation error in FEMNumericIntegration.");
    for (int i=0; i<bMesh1.numberOfBulkElements; i++) {
        //TODO: the following call is not complete, so it won't work.
        // We need to update the interface of GaussQuadrature so that it supports edge basis
        // In Elmer: IntegStuff = GaussPoints( Element, RelOrder=RelOrder, EdgeBasis=PiolaT )
        IP = GaussQuadrature(&bMesh1Elements[i], NULL, &relOrder);
        nbGaussPoints = nbGaussPoints + IP->n;
    }
    fprintf(stdout, "FEMInterpolateMeshToMesh:weightedProjectorMesh2: number of nodes and Gauss points: %d %d.\n", nbNodes, nbGaussPoints);
    
    gaussNodes->x = doublevec(0, nbGaussPoints-1);
    gaussNodes->y = doublevec(0, nbGaussPoints-1);
    gaussNodes->z = doublevec(0, nbGaussPoints-1);
    int sizeGaussNodes = nbGaussPoints;
    
    // Change the local coordinates of the bMesh2 to match to corresponding faces
    if (edgeBasis == YES) {
        int *xPerm = intvec(0, bMesh2.parent.numberOfNodes-1);
        memset(xPerm, 0, bMesh2.parent.numberOfNodes*sizeof(int) );
        for (int i=0; i<sizeInversePermutation2; i++) {
            xPerm[invPerm2[i]] = i;
        }
        
        Element_t *faces = bMesh2.parent.getFaces;
        for (int i=0; i<bMesh2.numberOfBulkElements; i++) {
            for (int j=0; j<bMesh2Elements[i].Type.NumberOfNodes; j++) {
                bMesh2Elements[i].NodeIndexes[j] = xPerm[faces[bMesh2Elements[i].ElementIndex-1].NodeIndexes[j]];
            }
        }
    }
    
    axisSym = NO;
    if (model.coordinates == axis_symmetric || model.coordinates == cylindric_symmetric) {
        if (nodalJump == YES) {
            axisSym = YES;
        } else if (solution != nil) {
            if ((solution.solutionInfo)[@"projector metrics"] != nil) {
                axisSym = [(solution.solutionInfo)[@"projector metrics"] boolValue];
            }
        }
        if (axisSym == YES) fprintf(stdout, "FEMInterpolateMeshToMesh:weightedProjectorMesh2: projector will be weighted for axi-symmetry.\n");
    }
    
    int totSize = bMesh1.parent.numberOfNodes + bMesh1.parent.numberOfEdges;
    if (solution != nil) {
        totSize = solution.matrix.numberOfRows / solution.variable.dofs;
    }
    
    if (edgeBasis == YES) {
        elementNodes->x = doublevec(0, bMesh1.parent.maxElementDofs-1);
        elementNodes->y = doublevec(0, bMesh1.parent.maxElementDofs-1);
        elementNodes->z = doublevec(0, bMesh1.parent.maxElementDofs-1);
        wBasis = doublematrix(0, bMesh1.parent.maxElementDofs-1, 0, 2);
        rotWBasis = doublematrix(0, bMesh1.parent.maxElementDofs-1, 0, 2);
        eqindsize = totSize;
    } else {
        elementNodes->x = doublevec(0, bMesh1.maxElementDofs-1);
        elementNodes->y = doublevec(0, bMesh1.maxElementDofs-1);
        elementNodes->z = doublevec(0, bMesh1.maxElementDofs-1);
        eqindsize = bMesh1.numberOfNodes;
    }
    
    eqindsize = 0;
    for (int i=0; i<bMesh1.numberOfBulkElements; i++) {
        if (edgeBasis == YES) {
            Element_t *faces = bMesh1.parent.getFaces;
            n = [core getElementDofsSolution:solution model:model forElement:&faces[bMesh1Elements[i].ElementIndex-1] atIndexes:indexes disableDiscontinuousGalerkin:NULL];
            np = [core getNumberOfNodesForElement:&faces[bMesh1Elements[i].ElementIndex-1]];
        } else {
            n = bMesh1Elements[i].Type.NumberOfNodes;
            np = n;
            memcpy(indexes, bMesh1Elements[i].NodeIndexes, n*sizeof(int));
        }
        eqindsize = max(eqindsize, max_array(indexes, n));
    }
    
    // Create the nodal coordinates for all Gaussian integration points
    nbGaussPoints = 0;
    for (int i=0; i<bMesh1.numberOfBulkElements; i++) {
        n = bMesh1Elements[i].Type.NumberOfNodes;
        for (int j=0; j<n; j++) {
            elementNodes->x[j] = realNodes->x[bMesh1Elements[i].NodeIndexes[j]];
            elementNodes->y[j] = realNodes->y[bMesh1Elements[i].NodeIndexes[j]];
            elementNodes->z[j] = realNodes->z[bMesh1Elements[i].NodeIndexes[j]];
        }
        
        //TODO: the following call is not complete, so it won't work.
        // We need to update the interface of GaussQuadrature so that it supports edge basis
        // In Elmer: IntegStuff = GaussPoints( Element, RelOrder=RelOrder, EdgeBasis=PiolaT )
        IP = GaussQuadrature(&bMesh1Elements[i], NULL, &relOrder);
        for (int j=0; j<IP->n; j++) {
            u = IP->u[j];
            v = IP->v[j];
            w = IP->w[j];
            if (piolaT == YES) {
                //TODO: the following call is not complete, so it won't work.
                // We need to update the FEMNumericalIntegration class methods so that they support edge basis
                // In Elmer: ElementInfo( Element, ElementNodes, u,v,w, detJ, Basis, EdgeBasis=WBasis )
            } else {
                stat = [integration setBasisForElement:&bMesh1Elements[i] elementNodes:elementNodes inMesh:bMesh1 firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
            }
            gaussNodes->x[nbGaussPoints] = cblas_ddot(n, integration.basis, 1, elementNodes->x, 1);
            gaussNodes->y[nbGaussPoints] = cblas_ddot(n, integration.basis, 1, elementNodes->y, 1);
            gaussNodes->z[nbGaussPoints] = cblas_ddot(n, integration.basis, 1, elementNodes->z, 1);
            nbGaussPoints++;
        }
    }
    
    Nodes_t *bMesh1Nodes = bMesh1.getNodes;
    bMesh1Nodes = gaussNodes;
    bMesh1.numberOfNodes = nbGaussPoints;
    
    // Create the mirror node flag and map the nodes of Mesh1 to be
    // in the interval of Mesh2
    BOOL *mirrorNode = NULL;
    if (repeating == YES) {
        if (antiRepeating == YES) {
            mirrorNode = (BOOL*)malloc(sizeof(BOOL) * bMesh1.numberOfNodes );
            memset(mirrorNode, 0, bMesh1.numberOfNodes*sizeof(BOOL) );
        }
        FEMMeshUtils *meshUtils = [[FEMMeshUtils alloc] init];
        int size = bMesh1.numberOfNodes;
        [meshUtils preRotationalProjectorMesh1:bMesh1 mesh2:bMesh2 mirrorNode:mirrorNode sizeMirrorNode:&size];
    }
    
    // Create the projector for Gaussian integration points
    gaussProjector = [utilities meshProjectorMesh1:bMesh2 mesh2:bMesh1 model:model useQuadrantTree:&useQuadrantTree transpose:NULL];
    matrixArraysContainer *gaussProjectorContainers = gaussProjector.getContainers;
    int *rows = gaussProjectorContainers->Rows;
    int *cols = gaussProjectorContainers->Cols;
    double *values = gaussProjectorContainers->Values;
    
    // If there are mirror nodes, change the sign
    if (antiRepeating == YES) {
        FEMMeshUtils *meshUtils = [[FEMMeshUtils alloc] init];
        int size = bMesh1.numberOfNodes;
        [meshUtils postRotationalProjector:gaussProjector mirrorNode:mirrorNode sizeMirrorNode:&size];
        if (mirrorNode != NULL) free(mirrorNode);
    }
    
    // Transfer the projector on the Gaussian points to that on
    // nodal points utilizing the flexibility of the list matrix
    projector = [[FEMMatrix alloc] init];
    projector.format = MATRIX_LIST;
    projector.projectorType = PROJECTOR_TYPE_GALERKIN;
    
    int *equind = intvec(0, eqindsize-1);
    memset(equind, 0, eqindsize*sizeof(int) );
    
    matrixArraysContainer *matrixContainers = projector.getContainers;
    matrixContainers->InvPerm = intvec(0, eqindsize-1);
    matrixContainers->sizeInvPerm = eqindsize;
    memset(matrixContainers->InvPerm, 0, eqindsize*sizeof(int) );
    
    FEMListMatrix *listMatrix = [[FEMListMatrix alloc] init];
    
    int ind = 0;
    int noGaussPoints = 0;
    double val, x, weight;
    BOOL any;
    for (int i=0; i<bMesh1.numberOfBulkElements; i++) {
        if (edgeBasis == YES) {
            Element_t *faces = bMesh1.parent.getFaces;
            element = &faces[bMesh1Elements[i].ElementIndex-1];
            n = [core getElementDofsSolution:solution model:model forElement:&faces[bMesh1Elements[i].ElementIndex-1] atIndexes:indexes disableDiscontinuousGalerkin:NULL];
            np = [core getNumberOfNodesForElement:&faces[bMesh1Elements[i].ElementIndex-1]];
            any = NO;
            for (int j=0; j<np; j++) {
                if (indexes[j] > faces[bMesh1Elements[i].ElementIndex-1].NodeIndexes[j]) {
                    any = YES;
                    break;
                }
            }
            if (any == YES) np = 0;
            if (nodes == NULL) {
                nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
                initNodes(nodes);
            }
            [core getNodes:solution model:model inElement:&faces[bMesh1Elements[i].ElementIndex-1] resultNodes:nodes numberOfNodes:NULL mesh:nil];
        } else {
            element = &bMesh1Elements[i];
            n = bMesh1Elements[i].Type.NumberOfNodes;
            np = n;
            memcpy(indexes, bMesh1Elements[i].NodeIndexes, n*sizeof(int));
            for (int j=0; j<n; j++) {
                elementNodes->x[j] = realNodes->x[indexes[j]];
                elementNodes->y[j] = realNodes->y[indexes[j]];
                elementNodes->z[j] = realNodes->z[indexes[j]];
            }
        }
        
        //TODO: the following call is not complete, so it won't work.
        // We need to update the interface of GaussQuadrature so that it supports edge basis
        // In Elmer: IntegStuff = GaussPoints( Element, RelOrder=RelOrder, EdgeBasis=PiolaT )
        IP = GaussQuadrature(&bMesh1Elements[i], NULL, &relOrder);
        for (int j=0; j<IP->n; j++) {
            u = IP->u[j];
            v = IP->v[j];
            w = IP->w[j];
            
            if (edgeBasis == YES) {
                if (piolaT == YES) {
                    // We need to update the FEMNumericalIntegration class methods so that they support edge basis
                    // In Elmer: ElementInfo( Element, Nodes, u,v,w, detJ, Basis, EdgeBasis=WBasis )
                } else {
                    stat = [integration setBasisForElement:element elementNodes:nodes inMesh:bMesh1 firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                    stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:bMesh1 firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                    stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:bMesh1 firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
                    [elementDescription getEdgeBasisElement:element wBasis:wBasis rotWBasis:rotWBasis basis:integration.basis dBasisdx:integration.basisFirstDerivative];
                }
            } else {
                stat = [integration setBasisForElement:element elementNodes:elementNodes inMesh:bMesh1 firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                stat = [integration setMetricDeterminantForElement:element elementNodes:elementNodes inMesh:bMesh1 firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
            }
            
            // Modify weight so that the projector is consistent with the coordinate system
            weight = integration.metricDeterminant * IP->s[j];
            if (axisSym == YES) {
                if (edgeBasis == YES) {
                    x = cblas_ddot(np, integration.basis, 1, nodes->x, 1);
                } else {
                    x = cblas_ddot(np, integration.basis, 1, elementNodes->x, 1);
                }
                weight = weight * x;
            }
            
            // Do the numbering of new dofs
            // This needs to be done here because the nodal jump
            // needs the index related to (p, q) pair
            for (int p=0; p<np; p++) {
                if (equind[indexes[p]] < 0) {
                    equind[indexes[p]] = ind;
                    if (edgeBasis == YES) {
                        matrixContainers->InvPerm[ind] = indexes[p];
                    } else {
                        matrixContainers->InvPerm[ind] = invPerm1[indexes[p]];
                    }
                    ind++;
                }
            }
            
            for (int p=0; p<np; p++) {
                val = weight * integration.basis[p];
                for (int q=0; q<np; q++) {
                    qq = indexes[q];
                    if (edgeBasis == NO) qq = invPerm1[qq];
                    [listMatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:equind[indexes[p]] andIndex:qq value:integration.basis[q]*val setValue:NULL];
                }
                for (int q=rows[noGaussPoints]; q<=rows[noGaussPoints+1]-1; q++) {
                    qq = cols[q];
                    if (qq < 0) break;
                    if (edgeBasis == NO) qq = invPerm2[qq];
                    [listMatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:equind[indexes[p]] andIndex:qq value:-periodicScale*values[q]*val setValue:NULL];
                }
            }
            
            if (edgeBasis == YES) {
                for (int p=np; p<n; p++) {
                    pp = p - np;
                    for (int k=0; k<3; k++) {
                        wq[k] = wBasis[pp][k];
                    }
                    
                    if (equind[indexes[p]] < 0) {
                        equind[indexes[p]] = ind;
                        matrixContainers->InvPerm[ind] = indexes[p];
                        ind++;
                    }
                    for (int q=np; q<n; q++) {
                        qq = q -  np;
                        for (int k=0; k<3; k++) {
                            vq[k] = wBasis[qq][k];
                        }
                        [listMatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:equind[indexes[p]] andIndex:indexes[q] value:weight*cblas_ddot(3, vq, 1, wq, 1) setValue:NULL];
                    }
                    for (int q=rows[noGaussPoints]+np; q<=rows[noGaussPoints+1]-1; q+=3) {
                        if (cols[q] > 0) fatal("FEMInterpolateMeshToMesh:weightedProjectorMesh2", "Error in this method.");
                        vq[0] = values[q];
                        vq[1] = values[q+1];
                        vq[2] = values[q+2];
                        [listMatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:equind[indexes[p]] andIndex:-cols[q] value:-periodicScale*weight*cblas_ddot(3, vq, 1, wq, 1) setValue:NULL];
                    }
                }
            }
            noGaussPoints++;
        }
    }
    
    [listMatrix convertToCRSMatrix:projector];
    
    bMesh1Nodes = realNodes;
    bMesh1.numberOfNodes = nbNodes;
    
    [integration deallocation:bMesh1];
    free_ivector(indexes, 0, 99);
    
    if (edgeBasis == YES) {
        free_dvector(elementNodes->x, 0, bMesh1.parent.maxElementDofs-1);
        free_dvector(elementNodes->y, 0, bMesh1.parent.maxElementDofs-1);
        free_dvector(elementNodes->z, 0, bMesh1.parent.maxElementDofs-1);
        free(elementNodes);
        free_dmatrix(wBasis, 0, bMesh1.parent.maxElementDofs-1, 0, 2);
        free_dmatrix(rotWBasis, 0, bMesh1.parent.maxElementDofs-1, 0, 2);
    } else {
        free_dvector(elementNodes->x, 0, bMesh1.maxElementDofs-1);
        free_dvector(elementNodes->y, 0, bMesh1.maxElementDofs-1);
        free_dvector(elementNodes->z, 0, bMesh1.maxElementDofs-1);
        free(elementNodes);
    }
    
    free_dvector(gaussNodes->x, 0, sizeGaussNodes-1);
    free_dvector(gaussNodes->y, 0, sizeGaussNodes-1);
    free_dvector(gaussNodes->z, 0, sizeGaussNodes-1);
    
    if (nodes != NULL) {
        free_dvector(nodes->x, 0, nodes->numberOfNodes-1);
        free_dvector(nodes->y, 0, nodes->numberOfNodes-1);
        free_dvector(nodes->z, 0, nodes->numberOfNodes-1);
        free(nodes);
    }
    
    return projector;
}

@end