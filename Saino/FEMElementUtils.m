//===----------------------------------------------------------------------===//
//  FEMElementUtils.m
//  Saino
//
//  Created by Seddik hakime on 28/12/12.
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

#import "FEMElementUtils.h"
#import "FEMListUtilities.h"
#import "FEMUtilities.h"
#import "FEMListMatrix.h"
#import "FEMBandwidthOptimize.h"
#import "FEMMatrixCRS.h"
#import "FEMMatrixBand.h"
#import "FEMNumericIntegration.h"
#import "FEMCoordinateSystems.h"
#import "GaussIntegration.h"
#import "Utils.h"

@interface FEMElementUtils ()
-(ListMatrix_t * _Nonnull)FEMElementUtils_makeListMatrixInModel:(FEMModel * _Nonnull)model solution:(FEMSolution * _Nonnull)solution mesh:(FEMMesh * _Nonnull)mesh reorder:(int * _Nonnull)reorder sizeOfReorder:(int)sizeOfReorder localNodes:(int)localNodes equation:(NSString * _Nullable)equation dgSolver:(BOOL * _Nonnull)dgSolver globalBubbles:(BOOL * _Nonnull)globalBubbles nodalDofsOnly:(BOOL * _Nullable)nodalDofsOnly projectorDofs:(BOOL * _Nullable)projectorDofs;
-(void)FEMElementUtils_initializeMatrix:(FEMMatrix * _Nonnull)matrix size:(int)n list:(ListMatrix_t * _Nonnull)list reorder:(int * _Nonnull)reorder invInitialReorder:(int * _Nonnull)invInitialReorder dofs:(int)dofs;
@end

@implementation FEMElementUtils

#pragma mark Private methods

/************************************************************************************************************************
 
    Create a list matrix given the mesh, the active domains and the element type related to the solver. The list
    matrix is flexible since it can account for any entries. Also constraints and periodic BCs may give rise to 
    entries in the list matrix topology.
 
    Method corresponds to Elmer from git on October 27 2015
 
************************************************************************************************************************/
-(ListMatrix_t * _Nonnull)FEMElementUtils_makeListMatrixInModel:(FEMModel * _Nonnull)model solution:(FEMSolution * _Nonnull)solution mesh:(FEMMesh * _Nonnull)mesh reorder:(int * _Nonnull)reorder sizeOfReorder:(int)sizeOfReorder localNodes:(int)localNodes equation:(NSString * _Nullable)equation dgSolver:(BOOL * _Nonnull)dgSolver globalBubbles:(BOOL * _Nonnull)globalBubbles nodalDofsOnly:(BOOL * _Nullable)nodalDofsOnly projectorDofs:(BOOL * _Nullable)projectorDofs {
    
    int t, i, j, k, l, m, k1, k2, n, edofs, fdofs, bdofs;
    int indexSize, numberOfFactors;
    int invPerm[localNodes], *indexes = NULL;
    BOOL doProjectors, foundDG, found, radiation;
    FEMMatrix *projector;
    FEMListUtilities *listUtilities;
    FEMListMatrix *listMatrix;
    ListMatrix_t *list = NULL;
    ListMatrixEntry_t *cList, *lptr;
    Element_t *elements, *element, *edges, *faces;
    matrixArraysContainer *matContainers = NULL;
    modelArraysContainer *modelContainers = NULL;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    
    listMatrix = [[FEMListMatrix alloc] init];
    list = [listMatrix allocateMatrix:localNodes];
    
    bdofs = mesh.maxBdofs;
    edofs = mesh.maxEdgeDofs;
    fdofs = mesh.maxFaceDofs;
    
    if (nodalDofsOnly != NULL) {
        if (*nodalDofsOnly == YES) {
            edofs = 0;
            fdofs = 0;
        }
    }
    
    if (projectorDofs != NULL) {
        doProjectors = *projectorDofs;
    } else doProjectors = YES;
    
    indexSize = 128;
    indexes = intvec(0, indexSize-1);
    
    elements = mesh.getElements;
    faces = mesh.getFaces;
    edges = mesh.getEdges;
    
    if (edofs > 0 && edges == NULL) {
        fprintf(stdout, "FEMElementUtils:FEMElementUtils_makeListMatrixInModel: edge dofs requested but no edges exist in mesh.\n");
        edofs = 0;
    }
    if (fdofs > 0 && faces == NULL) {
        fprintf(stdout, "FEMElementUtils:FEMElementUtils_makeListMatrixInModel: face dofs requested but no faces exist in mesh.\n");
        fdofs = 0;
    }
    
    // Create the permutation for the Discontinuous Galerkin solver
    foundDG = NO;
    if (*dgSolver == YES && equation != nil) {
        for (t=0; t<mesh.numberOfEdges; t++) {
            n = 0;
            if (edges[t].BoundaryInfo->Left != NULL) {
                if ([listUtilities checkElementEquation:model forElement:edges[t].BoundaryInfo->Left andEquation:equation] == YES) {
                    foundDG = (foundDG == YES || edges[t].BoundaryInfo->Left->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<edges[t].BoundaryInfo->Left->DGDOFs; j++) {
                        indexes[n] = edges[t].BoundaryInfo->Left->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            if (edges[t].BoundaryInfo->Right != NULL) {
                if ([listUtilities checkElementEquation:model forElement:edges[t].BoundaryInfo->Right andEquation:equation] == YES) {
                    foundDG = (foundDG == YES || edges[t].BoundaryInfo->Right->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<edges[t].BoundaryInfo->Right->DGDOFs; j++) {
                        indexes[n] = edges[t].BoundaryInfo->Right->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            for (i=0; i<n; i++) {
                k1 = reorder[indexes[i]];
                if (k1 < 0) continue;
                for (j=0; j<n; j++) {
                    k2 = reorder[indexes[j]];
                    if (k2 < 0) continue;
                    lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:k1 andIndex:k2];
                }
            }
        }
        
        for (t=0; t<mesh.numberOfFaces; t++) {
            n = 0;
            if (faces[t].BoundaryInfo->Left != NULL) {
                if ([listUtilities checkElementEquation:model forElement:faces[t].BoundaryInfo->Left andEquation:equation] == YES) {
                    foundDG = (foundDG == YES || faces[t].BoundaryInfo->Left->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<faces[t].BoundaryInfo->Left->DGDOFs; j++) {
                        indexes[n] = faces[t].BoundaryInfo->Left->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            if (faces[t].BoundaryInfo->Right != NULL) {
                if ([listUtilities checkElementEquation:model forElement:faces[t].BoundaryInfo->Right andEquation:equation] == NO) {
                    foundDG = (foundDG == YES || faces[t].BoundaryInfo->Right->DGDOFs > 0) ? YES : NO;
                    for (j=0; j<faces[t].BoundaryInfo->Right->DGDOFs; j++) {
                        indexes[n] = faces[t].BoundaryInfo->Right->DGIndexes[j];
                        n++;
                    }
                }
            }
            
            for (i=0; i<n; i++) {
                k1 = reorder[indexes[i]];
                if (k1 < 0) continue;
                for (j=0; j<n; j++) {
                    k2 = reorder[indexes[j]];
                    if (k2 < 0) continue;
                    lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:k1 andIndex:k2];
                }
            }
        }
    }
    
    if (foundDG == NO) {
        t = 0;
        while (t < mesh.numberOfBulkElements + mesh.numberOfBoundaryElements) {
            element = &elements[t];
            
            if (equation != nil) {
                while (t < mesh.numberOfBulkElements + mesh.numberOfBoundaryElements) {
                    element = &elements[t];
                    if ([listUtilities checkElementEquation:model forElement:element andEquation:equation] == YES) break;
                    t++;
                }
                if (t >= mesh.numberOfBulkElements + mesh.numberOfBoundaryElements) break;
            }
            
            n = element->NDOFs + element->Type.NumberOfEdges * edofs + element->Type.NumberOfFaces * fdofs;
            
            if (*globalBubbles == YES) n = n + element->BDOFs;
            
            if (n > indexSize) {
                if (indexes != NULL) free_ivector(indexes, 0, indexSize-1);
                indexSize = n;
                indexes = intvec(0, n-1);
            }
            
            n = 0;
            for (i=0; i<element->NDOFs; i++) {
                indexes[n] = element->NodeIndexes[i];
                n++;
            }
            
            if (edges != NULL) {
                if (element->EdgeIndexes != NULL) {
                    for (j=0; j<element->Type.NumberOfEdges; j++) {
                        for (i=0; i<edges[element->EdgeIndexes[j]].BDOFs; i++) {
                            indexes[n] = edofs * element->EdgeIndexes[j] + i + mesh.numberOfNodes;
                            n++;
                        }
                    }
                }
            }
            
            if (faces != NULL) {
                if (element->FaceIndexes != NULL) {
                    for (j=0; j<element->Type.NumberOfFaces; j++) {
                        for (i=0; i<faces[element->FaceIndexes[j]].BDOFs; i++) {
                            indexes[n] = fdofs * element->FaceIndexes[j] + i + mesh.numberOfNodes + edofs * mesh.numberOfEdges;
                            n++;
                        }
                    }
                }
            }
            
            if (*globalBubbles == YES && element->BubbleIndexes != NULL) {
                for (i=0; i<element->BDOFs; i++) {
                    indexes[n] = fdofs * mesh.numberOfFaces + mesh.numberOfNodes + edofs * mesh.numberOfEdges + element->BubbleIndexes[i];
                    n++;
                }
            }
            for (i=0; i<n; i++) {
                k1 = reorder[indexes[i]];
                if (k1 < 0) continue;
                for (j=0; j<n; j++) {
                    k2 = reorder[indexes[j]];
                    if (k2 < 0) continue;
                    lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:k1 andIndex:k2];
                }
            }
            t++;
        }
        
        if (indexes != NULL) free_ivector(indexes, 0, indexSize-1);
        
        // Diffuse gray radiation condition
        radiation = [listUtilities listGetLogical:model inArray:solution.valuesList forVariable:@"radiation solver" info:&found];
        if (found == NO && equation != nil) radiation = (radiation == YES || [equation isEqualToString:@"heat equation"] == YES) ? YES : NO;
        
        if (radiation == YES) {
            for (i=mesh.numberOfBulkElements; i<mesh.numberOfBulkElements + mesh.numberOfBoundaryElements; i++) {
                
                if (elements[i].BoundaryInfo->GebhardtFactors != NULL) {
                    for (j=0; j<elements[i].Type.NumberOfNodes; j++) {
                        k1 = reorder[elements[i].NodeIndexes[j]];
                        numberOfFactors = elements[i].BoundaryInfo->GebhardtFactors->NumberOfImplicitFactors;
                        
                        for (n=0; n<numberOfFactors; n++) {
                            for (k=0; k<elements[elements[i].BoundaryInfo->GebhardtFactors->Elements[n]].Type.NumberOfNodes; k++) {
                                k2 = reorder[elements[elements[i].BoundaryInfo->GebhardtFactors->Elements[n]].NodeIndexes[k]];
                                lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:k1 andIndex:k2];
                            }
                        }
                    }
                }
            }
        }
        
        for (i=mesh.numberOfBulkElements; i<mesh.numberOfBulkElements + mesh.numberOfBoundaryElements; i++) {
            if (elements[i].Type.ElementCode < 102 || elements[i].Type.ElementCode >= 200) continue;
            
            k1 = reorder[elements[i].NodeIndexes[0]];
            if (k1 >= 0) {
                for (k=0; k<elements[i].Type.NumberOfNodes; k++) {
                    k2 = reorder[elements[i].NodeIndexes[k]];
                    if (k2 >= 0) {
                        lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:k1 andIndex:k2];
                        lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:k2 andIndex:k1];
                    }
                }
            }
            
            // This is a connection element, make a matrix connection for that
            if (elements[i].Type.ElementCode == 102) {
                k2 = reorder[elements[i].NodeIndexes[1]];
                if (k2 >= 0) lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:k2 andIndex:k2];
            }
        }
        
        
        // Add connection from projectors. These are only needed if the projector is treated
        // implicitly. For explicit projectors or when using Lagrange coefficients, the
        // connections are not needed.
        if (doProjectors == YES) {
            int bd = 0;
            for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                bd++;
                projector = boundaryCondition.pMatrix;
                if (projector == nil) continue;
                
                if ([listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"periodic bc explicit" info:&found] == YES) continue;
                if ([listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"periodic bc use lagrange coefficient" info:&found] == YES) continue;
                
                fprintf(stdout, "FEMElementUtils:FEMElementUtils_makeListMatrixInModel: adding matrix topology for BC: %d.", bd);
                
                matContainers = projector.getContainers;
                for (i=0; i<projector.numberOfRows; i++) {
                    k = reorder[matContainers->InvPerm[i]];
                    if (k >= 0) {
                        for (l=matContainers->Rows[i]; l<=matContainers->Rows[i+1]-1; l++) {
                            if (matContainers->Cols[l] < 0) continue;
                            m = reorder[matContainers->Cols[l]];
                            if (m >= 0) {
                                lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:k andIndex:m];
                                lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:m andIndex:k]; // Keep structure symmetric
                                cList = list[k].Head;
                                while (cList != NULL) {
                                    lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:m andIndex:cList->Index];
                                    lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:cList->Index andIndex:m]; // Keep structure symmetric
                                    cList = cList->Next;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    k = 0;
    for (i=0; i<sizeOfReorder; i++) {
        if (reorder[i] >= 0) {
            invPerm[reorder[i]] = k;
            k++;
        }
    }
    
    model.totalMatrixElements = 0;
    modelContainers = model.getContainers;
    for (i=0; i<localNodes; i++) {
        modelContainers->rowNonZeros[invPerm[i]] = list[i].Degree;
        model.totalMatrixElements = model.totalMatrixElements + list[i].Degree;
    }
    return list;
}

/***************************************************************************************************************
    Initialize a CRS format matrix to the effect that it will be ready to accept values when the method
    glueLocalMatrix (defined in the FEMMatrixCRS class) is called (build up the index tables of 
    a CRS format matrix).
****************************************************************************************************************/
-(void)FEMElementUtils_initializeMatrix:(FEMMatrix * _Nonnull)matrix size:(int)n list:(ListMatrix_t * _Nonnull)list reorder:(int * _Nonnull)reorder invInitialReorder:(int * _Nonnull)invInitialReorder dofs:(int)dofs {
    
    int i, j, k, l, m, k1, k2;
    ListMatrixEntry_t *cList;
    FEMMatrixCRS *crsMatrix;
    matrixArraysContainer *matrixContainers = NULL;

    crsMatrix = [[FEMMatrixCRS alloc] init];
    matrixContainers = matrix.getContainers;
    
    for (i=0; i<n; i++) {
        for (l=0; l<dofs; l++) {
            cList = list[i].Head;
            j = reorder[invInitialReorder[i]];
            k1 = dofs * j + l;
            k2 = matrixContainers->Rows[k1] - 1;
            while (cList != NULL) {
                k = reorder[invInitialReorder[cList->Index]];
                k = dofs * k;
                for (m=k; m<k+dofs; m++) {
                    k2 = k2 + 1;
                    matrixContainers->Cols[k2] = m;
                }
                cList = cList->Next;
            }
        }
    }
    
    if (matrix.format == MATRIX_CRS) [crsMatrix sortMatrix:matrix alsoValues:NULL];
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

/********************************************************************************************
 
    Creates and return a matrix
 
    Method corresponds to Elmer from git on October 27 2015
 
********************************************************************************************/
-(FEMMatrix * _Nonnull)createMatrixInModel:(FEMModel * _Nonnull)model forSolution:(FEMSolution * _Nonnull)solution mesh:(FEMMesh * _Nonnull)mesh dofs:(int)dofs permutation:(int * _Nonnull)perm sizeOfPermutation:(int)permSize matrixFormat:(int)matrixFormat optimizeBandwidth:(BOOL)optimizeBandwidth equationName:(NSString * _Nullable)equation discontinuousGalerkinSolution:(BOOL * _Nullable)dgSolution globalBubbles:(BOOL * _Nullable)gbBubbles nodalDofsOnly:(BOOL * _Nullable)nodalDofsOnly projectorDofs:(BOOL * _Nullable)projectorDofs {
    
    int i, j, k, l, m, n, p, k1, edofs, fdofs, bdofs, cols;
    int *invInitialReorder;
    BOOL dg, gb, useOptimized, found, all;
    NSMutableString *str1, *str2;
    Element_t *elements, *edges, *faces;
    modelArraysContainer *modelContainers = NULL;
    ListMatrix_t *listMatrix;
    FEMMatrix *matrix, *constraint;
    FEMBandwidthOptimize *optimizeBW;
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    matrixArraysContainer *matContainers;
    listBuffer ivals = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    FEMListUtilities *listUtils = [FEMListUtilities sharedListUtilities];
    FEMUtilities *utils = [[FEMUtilities alloc] init];
    
    modelContainers = model.getContainers;
    
    dg = NO;
    if (dgSolution != NULL) dg = *dgSolution;
    
    gb = NO;
    if (gbBubbles != NULL) gb = *gbBubbles;
    
    elements = mesh.getElements;
    edges = mesh.getEdges;
    faces = mesh.getFaces;
    
    edofs = 0;
    for (i=0; i<mesh.numberOfEdges; i++) {
        edofs = max(edofs, edges[i].BDOFs);
    }
    mesh.maxEdgeDofs = edofs;
    
    fdofs = 0;
    for (i=0; i<mesh.numberOfFaces; i++) {
        fdofs = max(fdofs, faces[i].BDOFs);
    }
    mesh.maxFaceDofs = fdofs;
    
    bdofs = 0;
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        bdofs = max(bdofs, elements[i].BDOFs);
    }
    mesh.maxBdofs = bdofs;
    
    memset( perm, -1, permSize*sizeof(int) );
    
    if (equation != nil) {
        k = [utils initialPermutationInMesh:mesh model:model solution:solution equation:equation permutation:perm DGSolution:&dg globalBubbles:&gb];
        if (k <= 0) {
            return nil;
        }
    } else {
        k = permSize;
    }
    
    if (k == permSize) {
        if (nodalDofsOnly != NULL) {
            if (*nodalDofsOnly == YES) k = mesh.numberOfNodes;
        }
        for (i=0; i<k; i++) {
            perm[i] = i;
        }
    }
    
    invInitialReorder = intvec(0, k-1);
    memset( invInitialReorder, -1, k*sizeof(int) );
    for (i=0; i<permSize; i++) {
        if (perm[i] >= 0 ) invInitialReorder[perm[i]] = i;
    }
    
    useOptimized = [listUtils listGetLogical:model inArray:solution.valuesList forVariable:@"optimize bandwidth use always" info:&found];
    
    // Check if matrix class really need to be created
    if ([listUtils listGetLogical:model inArray:solution.valuesList forVariable:@"no matrix" info:&found] == YES) return nil;
    
    //Compute matrix structure and do bandwidth optimization if requested
    modelContainers->rowNonZeros = intvec(0, k-1);
    modelContainers->sizeRowNonZeros = k;
    memset(modelContainers->rowNonZeros , 0, k*sizeof(int) );
    listMatrix = NULL;
    
    optimizeBW = [[FEMBandwidthOptimize alloc] init];
    
    if (equation != nil) {
        listMatrix = [self FEMElementUtils_makeListMatrixInModel:model solution:solution mesh:mesh reorder:perm sizeOfReorder:permSize localNodes:k equation:equation dgSolver:&dg globalBubbles:&gb nodalDofsOnly:nodalDofsOnly projectorDofs:projectorDofs];
        n = [optimizeBW optimizeBandwidthInListMatrix:listMatrix permutation:perm sizeOfPerm:permSize invInitialReorder:invInitialReorder localNodes:k optimize:optimizeBandwidth useOptimized:useOptimized equation:equation];
    } else {
        listMatrix = [self FEMElementUtils_makeListMatrixInModel:model solution:solution mesh:mesh reorder:perm sizeOfReorder:permSize localNodes:k equation:nil dgSolver:&dg globalBubbles:&gb nodalDofsOnly:nodalDofsOnly projectorDofs:projectorDofs];
        n = [optimizeBW optimizeBandwidthInListMatrix:listMatrix permutation:perm sizeOfPerm:permSize invInitialReorder:invInitialReorder localNodes:k optimize:optimizeBandwidth useOptimized:useOptimized equation:@"[empty field]"];
    }
    
    // Create and initialize the matrix
    switch (matrixFormat) {
        case MATRIX_CRS:
            crsMatrix = [[FEMMatrixCRS alloc] init];
            matrix = [crsMatrix createMatrixWithNumberOfRows:dofs*k totalNonZeros:model.totalMatrixElements rowNonZeros:modelContainers->rowNonZeros degreesFreedom:dofs reorder:perm sizeOfReorder:permSize allocateValues:YES];
            matrix.format = matrixFormat;
            [self FEMElementUtils_initializeMatrix:matrix size:k list:listMatrix reorder:perm invInitialReorder:invInitialReorder dofs:dofs];
            break;
            
        case MATRIX_BAND:
            bandMatrix = [[FEMMatrixBand alloc] init];
            matrix = [bandMatrix createMatrixWithNumberOfRows:dofs*k subBand:dofs*n symmetric:NO allocateValues:YES];
            break;
            
        case MATRIX_SBAND:
            bandMatrix = [[FEMMatrixBand alloc] init];
            matrix = [bandMatrix createMatrixWithNumberOfRows:dofs*k subBand:dofs*n symmetric:YES allocateValues:YES];
            break;
    }
    
    FEMListMatrix *list = [[FEMListMatrix alloc] init];
    [list freeMatrix:listMatrix size:k];
    
    matrix.dgMatrix = dg;
    matrix.subband = dofs * n;
    matrix.complexMatrix = NO;
    matrix.format = matrixFormat;
    
    if (solution.solutionInfo[@"constraint dofs"] != nil) {
        n = [solution.solutionInfo[@"constraint dofs"] intValue];
    } n = 0;
    
    if (n > 0) {
        constraint = [[FEMMatrix alloc] init];
        matrix.constraint = constraint;
        matContainers = matrix.constraint.getContainers;
        
        matContainers->Rows = intvec(0, (n+1)-1);
        matContainers->Diag = intvec(0, n-1);
        matContainers->RHS = doublevec(0, n-1);
        
        for (i=0; i<n; i++) {
            str1 = [NSMutableString stringWithString:@"constraint dof "];
            str2 = [NSMutableString stringWithFormat:@"%d",i+1];
            [str1 appendString:str2];
            [str1 appendString:@" value"];
            matContainers->RHS[i] = [listUtils listGetConstReal:model inArray:solution.valuesList forVariable:str1 info:&found minValue:NULL maxValue:NULL];
        }
        
        crsMatrix = [[FEMMatrixCRS alloc] init];
        
        cols = 0;
        matContainers->Rows[0] = 0;
        for (i=0; i<n; i++) {
            str1 = [NSMutableString stringWithString:@"constraint dof "];
            str2 = [NSMutableString stringWithFormat:@"%d",i+1];
            [str1 appendString:str2];
            [str1 appendString:@" body"];
            if ([listUtils listGetIntegerArray:model inArray:solution.valuesList forVariable:str1 buffer:&ivals] == YES) {
                memset( invInitialReorder, -1, k*sizeof(int) );
                elements = solution.mesh.getElements;
                for (j=0; j<solution.mesh.numberOfBulkElements; j++) {
                    all = YES;
                    for (l=0; l<ivals.m; l++) {
                        if (ivals.ivector[l] != elements[j].BodyID) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) continue;
                    
                    all = YES;
                    for (l=0; l<elements[j].Type.NumberOfNodes; l++) {
                        if (perm[elements[j].NodeIndexes[l]] >= 0) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) {
                        for (l=0; l<elements[j].Type.NumberOfNodes; l++) {
                            invInitialReorder[elements[j].NodeIndexes[l]] = 0;
                        }
                    }
                }
                cols = cols + dofs * count(invInitialReorder, '=', 0, k);
            }
            
            str1 = [NSMutableString stringWithString:@"constraint dof "];
            str2 = [NSMutableString stringWithFormat:@"%d",i+1];
            [str1 appendString:str2];
            [str1 appendString:@" bc"];
            if ([listUtils listGetIntegerArray:model inArray:solution.valuesList forVariable:str1 buffer:&ivals] == YES) {
                memset( invInitialReorder, -1, k*sizeof(int) );
                elements = solution.mesh.getElements;
                for (j=solution.mesh.numberOfBulkElements; j<solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements; j++) {
                    all = YES;
                    for (l=0; l<ivals.m; l++) {
                        if (ivals.ivector[l] != elements[j].BoundaryInfo->Constraint) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) continue;
                    
                    all = YES;
                    for (l=0; l<elements[j].Type.NumberOfNodes; l++) {
                        if (perm[elements[j].NodeIndexes[l]] >= 0) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) {
                        for (l=0; l<elements[j].Type.NumberOfNodes; l++) {
                            invInitialReorder[elements[j].NodeIndexes[l]] = 0;
                        }
                    }
                }
                cols = cols + dofs * count(invInitialReorder, '=', 0, k);
            }
            
            matContainers->Rows[i+1] = matContainers->Rows[i]+cols;
        }
        
        matContainers->Cols = intvec(0, cols-1);
        matContainers->Values = doublevec(0, cols-1);
        matContainers->sizeCols = cols;
        matContainers->sizeValues = cols;
        memset( matContainers->Cols, -1, cols*sizeof(int) );
        memset( matContainers->Values, 0.0, cols*sizeof(double) );
        
        for (i=0; i<n; i++) {
            str1 = [NSMutableString stringWithString:@"constraint dof "];
            str2 = [NSMutableString stringWithFormat:@"%d",i+1];
            [str1 appendString:str2];
            [str1 appendString:@" body"];
            if ([listUtils listGetIntegerArray:model inArray:solution.valuesList forVariable:str1 buffer:&ivals] == YES) {
                elements = solution.mesh.getElements;
                for (j=0; j<solution.mesh.numberOfBulkElements; j++) {
                    all = YES;
                    for (l=0; l<ivals.m; l++) {
                        if (ivals.ivector[l] != elements[j].BodyID) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) continue;
                    
                    all = YES;
                    for (l=0; l<elements[j].Type.NumberOfNodes; l++) {
                        if (perm[elements[j].NodeIndexes[l]] >= 0) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) {
                        for (p=0; p<elements[j].Type.NumberOfNodes; p++) {
                            l = perm[elements[j].NodeIndexes[p]];
                            for (m=0; m<dofs; m++) {
                                k1 = dofs * l + m;
                                [crsMatrix makeMatrixIndex:matrix.constraint row:i col:k1];
                            }
                        }
                    }
                }
            }
            
            str1 = [NSMutableString stringWithString:@"constraint dof "];
            str2 = [NSMutableString stringWithFormat:@"%d",i+1];
            [str1 appendString:str2];
            [str1 appendString:@" bc"];
            if ([listUtils listGetIntegerArray:model inArray:solution.valuesList forVariable:str1 buffer:&ivals] == YES) {
                elements = solution.mesh.getElements;
                for (j=solution.mesh.numberOfBulkElements; j<solution.mesh.numberOfBulkElements+solution.mesh.numberOfBoundaryElements; j++) {
                    all = YES;
                    for (l=0; l<ivals.m; l++) {
                        if (ivals.ivector[l] != elements[j].BoundaryInfo->Constraint) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) continue;
                    
                    all = YES;
                    for (l=0; l<elements[j].Type.NumberOfNodes; l++) {
                        if (perm[elements[j].NodeIndexes[l]] >= 0) {
                            continue;
                        } else {
                            all = NO;
                            break;
                        }
                    }
                    if (all == YES) {
                        for (p=0; p<elements[j].Type.NumberOfNodes; p++) {
                            l = perm[elements[j].NodeIndexes[p]];
                            for (m=0; m<dofs; m++) {
                                k1 = dofs * l + m;
                                [crsMatrix makeMatrixIndex:matrix.constraint row:i col:k1];
                            }
                        }
                    }
                }
            }
        }
        [crsMatrix sortMatrix:matrix.constraint alsoValues:NULL];
    }
    
    free_ivector(modelContainers->rowNonZeros, 0, modelContainers->sizeRowNonZeros-1);
    modelContainers->rowNonZeros = NULL;
    free_ivector(invInitialReorder, 0, k-1);
    if (ivals.ivector != NULL) {
        free_ivector(ivals.ivector, 0, ivals.m-1);
        ivals.vector = NULL;
    }
    
    return matrix;
}

/*************************************************************************************************
    Given the normal, return the tangent directions. The first tangent direction will always be on 
    the xy-plane if also the normal is in the xy-plane.
*************************************************************************************************/
-(void)tangentDirectionsForNormal:(double * _Nonnull)normal tangent1:(double * _Nonnull)tangent1 tangent2:(double * _Nonnull)tangent2 {
    
    int i;
    double n1, n2, n3, sum;
    
    n1 = fabs(normal[0]);
    n2 = fabs(normal[1]);
    n3 = fabs(normal[2]);
    
    if (n1 <= n3 && n2 <= n3) {
        tangent1[0] = 0.0;
        tangent1[1] = -normal[2];
        tangent1[2] = normal[1];
    } else {
        tangent1[0] = -normal[1];
        tangent1[1] = normal[0];
        tangent1[2] = 0.0;
    }
    
    vDSP_svesqD(tangent1, 1, &sum, 3);
    for (i=0; i<3; i++) {
        tangent1[i] = tangent1[i] / sqrt(sum);
    }
    
    tangent2[0] = normal[1]*tangent1[2] - normal[2]*tangent1[1];
    tangent2[1] = normal[2]*tangent1[0] - normal[0]*tangent1[2];
    tangent2[2] = normal[0]*tangent1[1] - normal[1]*tangent1[0];
    
    vDSP_svesqD(tangent2, 1, &sum, 3);
    for (i=0; i<3; i++) {
        tangent2[i] = tangent2[i] / sqrt(sum);
    }
}

-(double)elementArea:(Element_t * _Nonnull)element numberOfNodes:(int)n mesh:(FEMMesh * _Nonnull)mesh nodel:(FEMModel * _Nonnull)model {
    
    int i, t;
    double a, detJ, sqrtMetric, nx[n], ny[n], nz[n], u, v, w, x, y, z;
    BOOL stat;
    Nodes_t *meshNodes, nodes;
    FEMNumericIntegration *integration;
    FEMCoordinateSystems *coordinateSystem;
    GaussIntegrationPoints *IP = NULL;
    
    nodes.x = nx;
    nodes.y = ny;
    nodes.z = nz;

    meshNodes = mesh.getNodes;
    for (i=0; i<element->Type.NumberOfNodes; i++) {
        nodes.x[i] = meshNodes->x[element->NodeIndexes[i]];
        nodes.y[i] = meshNodes->y[element->NodeIndexes[i]];
        nodes.z[i] = meshNodes->z[element->NodeIndexes[i]];
    }
    
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) fatal("FEMElementUtils:elementArea", "Allocation error in FEMNumericIntegration.");
    IP = GaussQuadrature(element, NULL, NULL);
    
    // Start integrating
    coordinateSystem = [[FEMCoordinateSystems alloc] init];
    a = 0.0;
    for (t=0; t<IP->n; t++) {
        // Integration stuff
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];
        
        // Basis function values & derivatives at the integration point
        stat = [integration setBasisForElement:element elementNodes:&nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:&nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        
        // Coordinate system dependent info
        if (model.coordinates != cartesian) {
            x = 0.0;
            y = 0.0;
            z = 0.0;
            for (i=0; i<n; i++) {
                x = x + (nodes.x[i]*integration.basis[i]);
                y = y + (nodes.y[i]*integration.basis[i]);
                z = z + (nodes.z[i]*integration.basis[i]);
            }
            sqrtMetric = [coordinateSystem coordinateSquareRootMetricModel:model coordX:x coordY:y coordZ:z];
            a = a + sqrtMetric * detJ * IP->s[t];
        } else {
            a = a + detJ * IP->s[t];
        }
    }
    
    [integration deallocation:mesh];
    return a;
}

@end
