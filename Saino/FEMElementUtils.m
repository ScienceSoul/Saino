//
//  FEMElementUtils.m
//  Saino
//
//  Created by Seddik hakime on 28/12/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMElementUtils.h"
#import "FEMListUtilities.h"
#import "FEMUtilities.h"
#import "FEMListMatrix.h"
#import "FEMBoundaryCondition.h"
#import "FEMBandwidthOptimize.h"
#import "FEMMatrixCRS.h"
#import "FEMMatrixBand.h"
#import "FEMNumericIntegration.h"
#import "FEMCoordinateSystems.h"
#import "GaussIntegration.h"
#import "Utils.h"

@interface FEMElementUtils ()
-(void)FEMElementUtils_makeListMatrixInModel:(FEMModel *)model solution:(FEMSolution *)solution mesh:(FEMMesh *)mesh list:(ListMatrix_t *)list reorder:(int *)reorder sizeOfReorder:(int)sizeOfReorder localNodes:(int)localNodes equation:(NSString *)equation dgSolver:(BOOL *)dgSolver globalBubbles:(BOOL *)globalBubbles;
-(void)FEMElementUtils_initializeMatrix:(FEMMatrix *)matrix size:(int)n list:(ListMatrix_t *)list reorder:(int *)reorder invInitialReorder:(int *)invInitialReorder dofs:(int)dofs;
@end

@implementation FEMElementUtils

#pragma mark Private methods

/************************************************************************************************************************
    Create a list a matrix given the mesh, the active domains and the element type related to the solver. The list
    matrix is flexible since it can account for any entries. Also constraints and periodic BCs may give rise to 
    entries in the list matrix topology.
************************************************************************************************************************/
-(void)FEMElementUtils_makeListMatrixInModel:(FEMModel *)model solution:(FEMSolution *)solution mesh:(FEMMesh *)mesh list:(ListMatrix_t *)list reorder:(int *)reorder sizeOfReorder:(int)sizeOfReorder localNodes:(int)localNodes equation:(NSString *)equation dgSolver:(BOOL *)dgSolver globalBubbles:(BOOL *)globalBubbles {
    
    int t, i, j, k, l, m, k1, k2, n, edofs, fdofs, bdofs, this;
    int indexSize, numberOfFactors;
    int invPerm[localNodes], *indexes;
    BOOL foundDG, gb, found, radiation;
    FEMMatrix *projector;
    FEMListUtilities *listUtilities;
    FEMListMatrix *listMatrix;
    FEMBoundaryCondition *boundaryConditions;
    ListMatrixEntry_t *cList, *lptr;
    Element_t *elements, *element, *edges, *faces;
    matrixArraysContainer *matContainers = NULL;
    modelArraysContainer *modelContainers = NULL;
    
    gb = NO;
    if (globalBubbles != NULL) gb = *globalBubbles;
    
    listUtilities = [[FEMListUtilities alloc] init];
    
    listMatrix = [[FEMListMatrix alloc] init];
    list = [listMatrix allocateMatrix:localNodes];
    
    bdofs = mesh.maxBdofs;
    edofs = mesh.maxEdgeDofs;
    fdofs = mesh.maxFaceDofs;
    
    indexSize = 128;
    indexes = intvec(0, indexSize-1);
    
    elements = mesh.getElements;
    faces = mesh.getFaces;
    edges = mesh.getEdges;
    
    foundDG = NO;
    if (*dgSolver == YES) {
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
                for (j = 0; j<n; j++) {
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
            
            if (gb == YES) n = n + element->BDOFs;
            
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
            
            if (gb == YES && element->BubbleIndexes != NULL) {
                for (i=0; i<element->BDOFs; i++) {
                    indexes[n] = fdofs * mesh.numberOfFaces + mesh.numberOfNodes + edofs * mesh.numberOfEdges + element->BubbleIndexes[i];
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
        
        if (elements[i].Type.ElementCode == 102) {
            k2 = reorder[elements[i].NodeIndexes[1]];
            if (k2 >= 0) lptr = [listMatrix getMatrixIndexInListMatrix:list atIndex:k2 andIndex:k2];
        }
    }
    
    for (this=0; this<model.numberOfBoundaryConditions; this++) {
        boundaryConditions = (model.boundaryConditions)[this];
        projector = boundaryConditions.pMatrix;
        if (projector == nil) continue;
        
        if ([listUtilities listGetLogical:model inArray:boundaryConditions.valuesList forVariable:@"periodic bc explicit" info:&found] == YES) continue;
        
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
}

/***************************************************************************************************************
    Initialize a CRS format matrix to the effect that it will be ready to accept values when the method
    glueLocalMatrix (defined in the FEMMatrixCRS class) is called (build up the index tables of 
    a CRS format matrix).
****************************************************************************************************************/
-(void)FEMElementUtils_initializeMatrix:(FEMMatrix *)matrix size:(int)n list:(ListMatrix_t *)list reorder:(int *)reorder invInitialReorder:(int *)invInitialReorder dofs:(int)dofs {
    
    int i, j, k, l, m, k1, k2;
    ListMatrixEntry_t *cList;
    FEMMatrixCRS *crsMatrix;
    
    crsMatrix = [[FEMMatrixCRS alloc] init];
    
    for (i=0; i<n; i++) {
        cList = list[i].Head;
        j = reorder[invInitialReorder[i]];
        while (cList != NULL) {
            k = reorder[invInitialReorder[cList->Index]];
            for (l=0; l<dofs; l++) {
                for (m=0; m<dofs; m++) {
                    k1 = dofs * j + l;
                    k2 = dofs * k * m;
                    [crsMatrix makeMatrixIndex:matrix atIndex:k1 andIndex:k2];
                }
            }
            cList = cList->Next;
        }
    }
    
    if (matrix.format == MATRIX_CRS) [crsMatrix sortInMatrix:matrix alsoValues:NULL];
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

-(FEMMatrix *)createMatrixInModel:(FEMModel *)model forSolution:(FEMSolution *)solution mesh:(FEMMesh *)mesh dofs:(int)dofs permutation:(int *)perm sizeOfPermutation:(int)permSize matrixFormat:(int)matrixFormat optimizeBandwidth:(BOOL)optimizeBandwidth equationName:(NSString *)equation discontinuousGalerkinSolution:(BOOL *)dgSolution globalBubbles:(BOOL *)gbBubbles {
    
    int i, j, k, l, m, n, p, k1, edofs, fdofs, bdofs, cols;
    int *invInitialReorder;
    BOOL dg, gb, useOptimized, found, all;
    NSMutableString *str1, *str2;
    Element_t *elements, *edges, *faces;
    modelArraysContainer *modelContainers = NULL;
    ListMatrix_t *listMatrix;
    FEMMatrix *matrix, *constraint;
    FEMListUtilities *listUtils;
    FEMUtilities *utils;
    FEMBandwidthOptimize *optimizeBW;
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    FEMListMatrix *list;
    matrixArraysContainer *matContainers;
    listBuffer ivals = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    listUtils = [[FEMListUtilities alloc] init];
    utils = [[FEMUtilities alloc] init];
    list = [[FEMListMatrix alloc] init];
    
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
        [self FEMElementUtils_makeListMatrixInModel:model solution:solution mesh:mesh list:listMatrix reorder:perm sizeOfReorder:permSize localNodes:k equation:equation dgSolver:&dg globalBubbles:&gb];
        n = [optimizeBW optimizeBandwidthInListMatrix:listMatrix permutation:perm sizeOfPerm:permSize invInitialReorder:invInitialReorder localNodes:k optimize:optimizeBandwidth useOptimized:useOptimized equation:equation];
    } else {
        [self FEMElementUtils_makeListMatrixInModel:model solution:solution mesh:mesh list:listMatrix reorder:perm sizeOfReorder:permSize localNodes:k equation:nil dgSolver:&dg globalBubbles:&gb];
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
    
    [list freeMatrix:listMatrix size:k];
    
    matrix.dgMatrix = dg;
    matrix.subband = dofs * n;
    matrix.complexMatrix = NO;
    matrix.format = matrixFormat;
    
    n = [listUtils listGetInteger:model inArray:solution.valuesList forVariable:@"constraint dofs" info:&found minValue:NULL maxValue:NULL];
    if (n > 0) {
        constraint = [[FEMMatrix alloc] init];
        matrix.constraint = constraint;
        matContainers = matrix.constraint.getContainers;
        
        matContainers->Rows = intvec(0, (n+1)-1);
        matContainers->Diag = intvec(0, n-1);
        matContainers->RHS = doublevec(0, n-1);
        
        for (i=0; i<n; i++) {
            str1 = [NSMutableString stringWithString:@"constraint dof "];
            str2 = [NSString stringWithFormat:@"%d",i+1];
            [str1 appendString:str2];
            [str1 appendString:@" value"];
            matContainers->RHS[i] = [listUtils listGetConstReal:model inArray:solution.valuesList forVariable:str1 info:&found minValue:NULL maxValue:NULL];
        }
        
        crsMatrix = [[FEMMatrixCRS alloc] init];
        
        cols = 0;
        matContainers->Rows[0] = 0;
        for (i=0; i<n; i++) {
            str1 = [NSMutableString stringWithString:@"constraint dof "];
            str2 = [NSString stringWithFormat:@"%d",i+1];
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
            str2 = [NSString stringWithFormat:@"%d",i+1];
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
            str2 = [NSString stringWithFormat:@"%d",i+1];
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
                                [crsMatrix makeMatrixIndex:matrix.constraint atIndex:i andIndex:k1];
                            }
                        }
                    }
                }
            }
            
            str1 = [NSMutableString stringWithString:@"constraint dof "];
            str2 = [NSString stringWithFormat:@"%d",i+1];
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
                                [crsMatrix makeMatrixIndex:matrix.constraint atIndex:i andIndex:k1];
                            }
                        }
                    }
                }
            }
        }
        [crsMatrix sortInMatrix:matrix.constraint alsoValues:NULL];
    }
    
    free_ivector(modelContainers->rowNonZeros, 0, modelContainers->sizeRowNonZeros-1);
    modelContainers->rowNonZeros = NULL;
    free_ivector(invInitialReorder, 0, k-1);
    free_ivector(ivals.ivector, 0, ivals.m-1);
    
    return matrix;
}

/*************************************************************************************************
    Given the normal, return the tangent directions. The first tangent direction will always be on 
    the xy-plane if also the normal is in the xy-plane.
*************************************************************************************************/
-(void)tangentDirectionsForNormal:(double *)normal tangent1:(double *)tangent1 tangent2:(double *)tangent2 {
    
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
    
    sum = 0.0;
    for (i=0; i<2; i++) {
        sum = sum + pow(tangent1[i], 2.0);
    }
    for (i=0; i<2; i++) {
        tangent1[i] = tangent1[i] / sqrt(sum);
    }
    
    tangent2[0] = normal[1]*tangent1[2] - normal[2]*tangent1[1];
    tangent2[1] = normal[2]*tangent1[0] - normal[0]*tangent1[2];
    tangent2[2] = normal[0]*tangent1[1] - normal[1]*tangent1[0];
    
    sum = 0.0;
    for (i=0; i<2; i++) {
        sum = sum + pow(tangent2[i], 2.0);
    }
    for (i=0; i<2; i++) {
        tangent2[i] = tangent2[i] / sqrt(sum);
    }
}

-(double)elementArea:(Element_t *)element numberOfNodes:(int)n mesh:(FEMMesh *)mesh nodel:(FEMModel *)model {
    
    int i, t;
    double a, detJ, sqrtMetric, sum, nx[n], ny[n], nz[n], u, v, w, x, y, z;
    BOOL stat;
    Nodes_t *meshNodes, nodes;
    FEMNumericIntegration *integration;
    FEMCoordinateSystems *coordinateSystem;
    GaussIntegrationPoints *IP;
    
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
    if ([integration allocation:mesh] == NO) errorfunct("FEMModel_localMatrix", "Allocation error in FEMNumericIntegration!");
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
            sum = 0.0;
            for (i=0; i<n; i++) {
                sum = sum + (nodes.x[i]*integration.basis[i]);
            }
            x = sum;
            
            sum = 0.0;
            for (i=0; i<n; i++) {
                sum = sum + (nodes.y[i]*integration.basis[i]);
            }
            y = sum;
            
            sum = 0.0;
            for (i=0; i<n; i++) {
                sum = sum + (nodes.z[i]*integration.basis[i]);
            }
            z = sum;
            
            sqrtMetric = [coordinateSystem coordinateSquareRootMetricModel:model coordX:x coordY:y coordZ:z];
            a = a + sqrtMetric * detJ * IP->s[t];
        } else {
            a = a + detJ * IP->s[t];
        }
    }
    
    return a;
}

@end
