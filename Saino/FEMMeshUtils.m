//
//  FEMMeshUtils.m
//  Saino
//
//  Created by Seddik hakime on 28/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMeshUtils.h"
#import "Utils.h"
#import "FEMElementUtils.h"
#import "FEMListUtilities.h"
#import "FEMBoundaryCondition.h"
#import "FEMPElementMaps.h"
#import "FEMElementDescription.h"
#import "FEMUtilities.h"

static double AEPS = 10.0 * DBL_EPSILON;

@interface FEMMeshUtils ()
-(void)FEMMeshUtils_assignConstraints:(FEMMesh *)mesh;
-(void)FEMMeshUtils_fixFaceEdges:(FEMMesh *)mesh;
-(Element_t *)FEMMeshUtils_getEntityForElement:(Element_t *)element edge:(int)number inMesh:(FEMMesh *)mesh;
@end

@implementation FEMMeshUtils

#pragma mark Private methods

-(void)FEMMeshUtils_assignConstraints:(FEMMesh *)mesh {
    
    int i, j, k, l, n, nfound, count;
    int *faceIndexes;
    Element_t *elements, *element, *boundary, *faces, *face;
    
    elements = mesh.getElements;
    
    for (i=0; i<mesh.numberOfBoundaryElements; i++) {
            
        boundary = &elements[mesh.numberOfBulkElements];
        
        element = boundary->BoundaryInfo->Left;
        if (element == NULL) element = boundary->BoundaryInfo->Right;
        if (element == NULL) continue;
        
        switch (boundary->Type.dimension) {
            case 1:
                faces = mesh.getEdges;
                faceIndexes = element->EdgeIndexes;
                count = element->sizeEdgeIndexes;
                break;
            case 2:
                faces = mesh.getFaces;
                faceIndexes = element->FaceIndexes;
                count = element->sizeFaceIndexes;
                break;
            default:
                faces = NULL;
                faceIndexes = NULL;
                break;
        }
        
        if (faces == NULL || faceIndexes == NULL) continue;

        for (j=0; j<count; j++) {
            face = &faces[faceIndexes[j]];
            
            n = boundary->Type.NumberOfNodes;
            if (n == 0) continue;
            nfound = 0;
            for (k=0; k<n; k++) {
                for (l=0; l<n; l++) {
                    if (boundary->NodeIndexes[k] == face->NodeIndexes[l]) nfound++;
                }
            }
            if (nfound == n) {
                face->BoundaryInfo = boundary->BoundaryInfo;
                break;
            }
        }
    }
}

-(void)FEMMeshUtils_fixFaceEdges:(FEMMesh *)mesh {
    
    int i, j, k, l, n, swap;
    int edgeind[4], i1[2], i2[2];
    Element_t *edges, *faces;
    BOOL condition;
    
    faces = mesh.getFaces;
    edges = mesh.getEdges;
    
    for (i=0; i<mesh.numberOfFaces; i++) {
        n = faces[i].Type.NumberOfEdges;
        for (l=0; l<n; l++) {
            edgeind[l] = faces[i].EdgeIndexes[l];
        }
        for (j=0; j<n; j++) {
            for (l=0; l<2; l++) {
                i1[l] = edges[edgeind[j]].NodeIndexes[l];
            }
            if (i1[0] > i1[1]) {
                swap = i1[0];
                i1[0] = i1[1];
                i1[1] = swap;
            }
            for (k=0; k<n; k++) {
                i2[0] = k;
                i2[1] = k+1;
                if (i2[1] > n-1) i2[1] = 0;
                for (l=0; l<2; l++) {
                    i2[l] = faces[i].NodeIndexes[i2[l]];
                }
                if (i2[0] > i2[1]) {
                    swap = i2[0];
                    i2[0] = i2[1];
                    i2[1] = swap;
                }
                condition = YES;
                for (l=0; l<2; l++) {
                    if (i1[l] != i2[l]) {
                        condition = NO;
                        break;
                    }
                }
                if (condition == YES) {
                    faces[i].EdgeIndexes[k] = edgeind[k];
                    break;
                }
            }
        }
    }
    faces = NULL;
    edges = NULL;
}

-(Element_t *)FEMMeshUtils_getEntityForElement:(Element_t *)element edge:(int)number inMesh:(FEMMesh *)mesh {
    
    Element_t *edges, *faces, *entity;
    
    edges = mesh.getEdges;
    faces = mesh.getFaces;
    
    entity = NULL;
    // Switch by element dimension
    switch (element->Type.dimension) {
        case 2:
            entity = &edges[element->EdgeIndexes[number]];
            break;
        case 3:
            entity = &faces[element->FaceIndexes[number]];
            break;
        default:
            NSLog(@"FEMMeshUtils_getEntityForElement: Unsupported dimension.");
            return entity;
            break;
    }
    return entity;
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

/************************************************************************

    Find 2D mesh edges.
 
************************************************************************/
-(void)findEdges2DInMesh:(FEMMesh *)mesh {
    
    int i, j, k, n, numbOfEdges;
    int node1, node2, edge, swap, degree;
    HashTable_t *hashTable;
    HashEntry_t *hashPtr, *hashPtr1;
    Element_t *elements, *element, *edges;
    ElementType_t *elmType;
    FEMElementDescription *elmDescription;
    BOOL found;
    
    elmDescription = [[FEMElementDescription alloc] init];
    
    elements = mesh.getElements;
    edges = mesh.getEdges;
    edges = (Element_t*)malloc(sizeof(Element_t) * (4*mesh.numberOfBulkElements) );
    
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        
        switch (element->Type.ElementCode / 100) {
            case 3:
                n = 3;
                break;
            case 4:
                n = 4;
                break;
        }
        
        if (element->EdgeIndexes == NULL) element->EdgeIndexes = intvec(0, n-1);
        memset( element->EdgeIndexes, 0, (n*sizeof(element->EdgeIndexes)) );
    }
    
    hashTable = (HashTable_t*)malloc(sizeof(HashTable_t) * (mesh.numberOfNodes) );
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashTable[i].head = NULL;
    }
    
    // Loop over elements
    numbOfEdges = 0;
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        
        switch (element->Type.ElementCode / 100) {
            case 3:
                n = 3;
                break;
            case 4:
                n = 4;
                break;
        }
        
        // Loop over edge of every element
        for (k=0; k<n; k++) {
            // We use min(node1, node2) as the hash table key
            node1 = element->NodeIndexes[k];
            if (k < n-1) {
                node2 = element->NodeIndexes[k+1];
            } else {
                node2 = element->NodeIndexes[0];
            }
            
            if (node2 < node1) {
                swap = node1;
                node1 = node2;
                node2 = swap;
            }
            
            // Look for the edge for the hash table
            hashPtr = hashTable[node1].head;
            found = NO;
            while (hashPtr != NULL) {
                if (hashPtr->node1 == node2) {
                    found = YES;
                    edge = hashPtr->edge;
                    break;
                }
                hashPtr = hashPtr->next;
            }
            
            // Existing edge, update structures
            if (found == YES) {
                element->EdgeIndexes[k] = edge;
                edges[edge].BoundaryInfo->Right = element;
            } else {
                
                //Edge not there yet, create it
                edge = numbOfEdges;
                degree = element->Type.BasisFunctionDegree;
                
                edges[edge].ElementIndex = edge+1;
                edges[edge].NodeIndexes = intvec(0, (degree+1)-1);
                edges[edge].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(edges[edge].BoundaryInfo);
                elmType = [elmDescription getElementType:(201+degree) inMesh:mesh stabilization:NO];
                edges[edge].Type = *elmType;
                
                edges[edge].NodeIndexes[0] = element->NodeIndexes[k];
                if (k < n-1) {
                    edges[edge].NodeIndexes[1] = element->NodeIndexes[k+1];
                } else {
                    edges[edge].NodeIndexes[1] = element->NodeIndexes[0];
                }
                
                for (j=1; j<degree; j++) {
                    edges[edge].NodeIndexes[j+1] = element->NodeIndexes[k+n+j-1];
                }
                
                // Create P element definitions if necessary
                if (element->Pdefs != NULL) {
                    edges[edge].Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
                    edges[edge].Pdefs->p = 0;
                } else {
                    edges[edge].Pdefs = NULL;
                }
                
                edges[edge].NDOFs = 0;
                if (element->NDOFs != 0) edges[edge].NDOFs = edges[edge].Type.NumberOfNodes;
                edges[edge].BDOFs = 0;
                edges[edge].DGDOFs = 0;
                edges[edge].EdgeIndexes = NULL;
                edges[edge].FaceIndexes = NULL;
                
                element->EdgeIndexes[k] = edge;
                
                edges[edge].BoundaryInfo->Left = element;
                edges[edge].BoundaryInfo->Right = NULL;
                
                // Update the hash table
                hashPtr = (HashEntry_t *)malloc(sizeof(HashEntry_t));
                hashPtr->edge = edge;
                hashPtr->node1 = node2;
                hashPtr->next = hashTable[node1].head;
                hashTable[node1].head = hashPtr;
                
                numbOfEdges++;
            }
        }
    }
    
    mesh.numberOfEdges = numbOfEdges;
    
    // Delete the hash table
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashPtr = hashTable[i].head;
        while (hashPtr != NULL) {
            hashPtr1 = hashPtr->next;
            free(hashPtr);
            hashPtr = hashPtr1;
        }
    }
    free(hashTable);
    
    [elmDescription deallocation];
}

/************************************************************************
 
 Find 3D mesh edges.
 
************************************************************************/
-(void)findEdges3DInMesh:(FEMMesh *)mesh {
    
    int i, ii, jj, k, n, numbOfEdges;
    int n1, n2;
    int node1, node2, edge, degree;
    int **edgeMap, **faceEdgeMap;
    int **tetraEdgeMap, **tetraFaceMap, **brickEdgeMap, **wedgeEdgeMap, **pyramidEdgeMap;
    int **tetraFaceEdgeMap, **brickFaceEdgeMap, **wedgeFaceEdgeMap, **pyramidFaceEdgeMap;
    HashTable_t *hashTable;
    HashEntry_t *hashPtr, *hashPtr1;
    Element_t *elements, *element, *edges, *faces;
    ElementType_t *elmType;
    FEMPElementMaps *elementMaps;
    FEMElementDescription *elmDescription;
    BOOL found;

    tetraEdgeMap = intmatrix(0, 5, 0, 2);
    tetraFaceMap = intmatrix(0, 3, 0, 5);
    brickEdgeMap = intmatrix(0, 5, 0, 8);
    wedgeEdgeMap = intmatrix(0, 4, 0, 3);
    pyramidEdgeMap = intmatrix(0, 4, 0, 7);
    
    tetraFaceEdgeMap = intmatrix(0, 3, 0, 2);
    brickFaceEdgeMap = intmatrix(0, 7, 0, 3);
    wedgeFaceEdgeMap = intmatrix(0, 5, 0, 3);
    pyramidFaceEdgeMap = intmatrix(0, 4, 0, 3);
    
    tetraFaceMap[0][0] = 0; tetraFaceMap[0][1] = 1; tetraFaceMap[0][2] = 2; tetraFaceMap[0][3] = 4; tetraFaceMap[0][4] = 5; tetraFaceMap[0][5] = 6;
    tetraFaceMap[1][0] = 0; tetraFaceMap[1][1] = 1; tetraFaceMap[1][2] = 3; tetraFaceMap[1][3] = 4; tetraFaceMap[1][4] = 8; tetraFaceMap[1][5] = 7;
    tetraFaceMap[2][0] = 1; tetraFaceMap[2][1] = 2; tetraFaceMap[2][2] = 3; tetraFaceMap[2][3] = 5; tetraFaceMap[2][4] = 9; tetraFaceMap[2][5] = 8;
    tetraFaceMap[3][0] = 2; tetraFaceMap[3][1] = 0; tetraFaceMap[3][2] = 3; tetraFaceMap[3][3] = 6; tetraFaceMap[3][4] = 7; tetraFaceMap[3][5] = 9;
    
    tetraFaceEdgeMap[0][0] = 0; tetraFaceEdgeMap[0][1] = 1; tetraFaceEdgeMap[0][2] = 2;
    tetraFaceEdgeMap[1][0] = 0; tetraFaceEdgeMap[1][1] = 4; tetraFaceEdgeMap[1][2] = 3;
    tetraFaceEdgeMap[2][0] = 1; tetraFaceEdgeMap[2][1] = 5; tetraFaceEdgeMap[2][2] = 4;
    tetraFaceEdgeMap[3][0] = 2; tetraFaceEdgeMap[3][1] = 3; tetraFaceEdgeMap[3][2] = 5;
    
    tetraEdgeMap[0][0] = 0; tetraEdgeMap[0][1] = 1; tetraEdgeMap[0][2] = 4;
    tetraEdgeMap[1][0] = 1; tetraEdgeMap[1][1] = 2; tetraEdgeMap[1][2] = 5;
    tetraEdgeMap[2][0] = 2; tetraEdgeMap[2][1] = 0; tetraEdgeMap[2][2] = 6;
    tetraEdgeMap[3][0] = 0; tetraEdgeMap[3][1] = 3; tetraEdgeMap[3][2] = 7;
    tetraEdgeMap[4][0] = 1; tetraEdgeMap[4][1] = 3; tetraEdgeMap[4][2] = 8;
    tetraEdgeMap[5][0] = 2; tetraEdgeMap[5][1] = 3; tetraEdgeMap[5][2] = 9;
    
    pyramidEdgeMap[0][0] = 0; pyramidEdgeMap[0][1] = 1; pyramidEdgeMap[0][2] = 0;
    pyramidEdgeMap[1][0] = 1; pyramidEdgeMap[1][1] = 2; pyramidEdgeMap[1][2] = 0;
    pyramidEdgeMap[2][0] = 2; pyramidEdgeMap[2][1] = 3; pyramidEdgeMap[2][2] = 0;
    pyramidEdgeMap[3][0] = 3; pyramidEdgeMap[3][1] = 0; pyramidEdgeMap[3][2] = 0;
    pyramidEdgeMap[4][0] = 0; pyramidEdgeMap[4][1] = 4; pyramidEdgeMap[4][2] = 0;
    pyramidEdgeMap[5][0] = 1; pyramidEdgeMap[5][1] = 4; pyramidEdgeMap[5][2] = 0;
    pyramidEdgeMap[6][0] = 2; pyramidEdgeMap[6][1] = 4; pyramidEdgeMap[6][2] = 0;
    pyramidEdgeMap[7][0] = 3; pyramidEdgeMap[7][1] = 4; pyramidEdgeMap[7][2] = 0;
    
    pyramidFaceEdgeMap[0][0] = 0; pyramidFaceEdgeMap[0][1] = 1; pyramidFaceEdgeMap[0][2] = 2; pyramidFaceEdgeMap[0][3] = 3;
    pyramidFaceEdgeMap[1][0] = 0; pyramidFaceEdgeMap[1][1] = 5; pyramidFaceEdgeMap[1][2] = 4; pyramidFaceEdgeMap[1][3] = -1;
    pyramidFaceEdgeMap[2][0] = 1; pyramidFaceEdgeMap[2][1] = 6; pyramidFaceEdgeMap[2][2] = 5; pyramidFaceEdgeMap[2][3] = -1;
    pyramidFaceEdgeMap[3][0] = 2; pyramidFaceEdgeMap[3][1] = 7; pyramidFaceEdgeMap[3][2] = 6; pyramidFaceEdgeMap[3][3] = -1;
    pyramidFaceEdgeMap[4][0] = 3; pyramidFaceEdgeMap[4][1] = 4; pyramidFaceEdgeMap[4][2] = 7; pyramidFaceEdgeMap[4][3] = -1;
    
    wedgeEdgeMap[0][0] = 0; wedgeEdgeMap[0][1] = 1; wedgeEdgeMap[0][2] = 0;
    wedgeEdgeMap[1][0] = 1; wedgeEdgeMap[1][1] = 2; wedgeEdgeMap[1][2] = 0;
    wedgeEdgeMap[2][0] = 0; wedgeEdgeMap[2][1] = 2; wedgeEdgeMap[2][2] = 0;
    wedgeEdgeMap[3][0] = 3; wedgeEdgeMap[3][1] = 4; wedgeEdgeMap[3][2] = 0;
    wedgeEdgeMap[4][0] = 4; wedgeEdgeMap[4][1] = 5; wedgeEdgeMap[4][2] = 0;
    wedgeEdgeMap[5][0] = 5; wedgeEdgeMap[5][1] = 3; wedgeEdgeMap[5][2] = 0;
    wedgeEdgeMap[6][0] = 0; wedgeEdgeMap[6][1] = 3; wedgeEdgeMap[6][2] = 0;
    wedgeEdgeMap[7][0] = 1; wedgeEdgeMap[7][1] = 4; wedgeEdgeMap[7][2] = 0;
    wedgeEdgeMap[8][0] = 2; wedgeEdgeMap[8][1] = 5; wedgeEdgeMap[8][2] = 0;
    
    wedgeFaceEdgeMap[0][0] = 0; wedgeFaceEdgeMap[0][1] = 1; wedgeFaceEdgeMap[0][2] = 2; wedgeFaceEdgeMap[0][3] = -1;
    wedgeFaceEdgeMap[1][0] = 3; wedgeFaceEdgeMap[1][1] = 4; wedgeFaceEdgeMap[1][2] = 5; wedgeFaceEdgeMap[1][3] = -1;
    wedgeFaceEdgeMap[2][0] = 0; wedgeFaceEdgeMap[2][1] = 7; wedgeFaceEdgeMap[2][2] = 3; wedgeFaceEdgeMap[2][3] = 6;
    wedgeFaceEdgeMap[3][0] = 1; wedgeFaceEdgeMap[3][1] = 8; wedgeFaceEdgeMap[3][2] = 4; wedgeFaceEdgeMap[3][3] = 7;
    wedgeFaceEdgeMap[4][0] = 2; wedgeFaceEdgeMap[4][1] = 6; wedgeFaceEdgeMap[4][2] = 5; wedgeFaceEdgeMap[4][3] = 8;
    
    brickEdgeMap[0][0] = 0; brickEdgeMap[0][1] = 1; brickEdgeMap[0][2] = 8;
    brickEdgeMap[1][0] = 1; brickEdgeMap[1][1] = 2; brickEdgeMap[1][2] = 9;
    brickEdgeMap[2][0] = 3; brickEdgeMap[2][1] = 2; brickEdgeMap[2][2] = 10;
    brickEdgeMap[3][0] = 0; brickEdgeMap[3][1] = 3; brickEdgeMap[3][2] = 11;
    brickEdgeMap[4][0] = 4; brickEdgeMap[4][1] = 5; brickEdgeMap[4][2] = 12;
    brickEdgeMap[5][0] = 5; brickEdgeMap[5][1] = 6; brickEdgeMap[5][2] = 13;
    brickEdgeMap[6][0] = 7; brickEdgeMap[6][1] = 6; brickEdgeMap[6][2] = 14;
    brickEdgeMap[7][0] = 4; brickEdgeMap[7][1] = 7; brickEdgeMap[7][2] = 15;
    brickEdgeMap[8][0] = 0; brickEdgeMap[8][1] = 4; brickEdgeMap[8][2] = 16;
    brickEdgeMap[9][0] = 1; brickEdgeMap[9][1] = 5; brickEdgeMap[9][2] = 17;
    brickEdgeMap[10][0] = 2; brickEdgeMap[10][1] = 6; brickEdgeMap[10][2] = 18;
    brickEdgeMap[11][0] = 3; brickEdgeMap[11][1] = 7; brickEdgeMap[11][2] = 19;
    
    brickFaceEdgeMap[0][0] = 0; brickFaceEdgeMap[0][1] = 1; brickFaceEdgeMap[0][2] = 2; brickFaceEdgeMap[0][3] = 3;
    brickFaceEdgeMap[1][0] = 4; brickFaceEdgeMap[1][1] = 5; brickFaceEdgeMap[1][2] = 6; brickFaceEdgeMap[1][3] = 7;
    brickFaceEdgeMap[2][0] = 0; brickFaceEdgeMap[2][1] = 9; brickFaceEdgeMap[2][2] = 4; brickFaceEdgeMap[2][3] = 8;
    brickFaceEdgeMap[3][0] = 1; brickFaceEdgeMap[3][1] = 10; brickFaceEdgeMap[3][2] = 5; brickFaceEdgeMap[3][3] = 9;
    brickFaceEdgeMap[4][0] = 2; brickFaceEdgeMap[4][1] = 11; brickFaceEdgeMap[4][2] = 6; brickFaceEdgeMap[4][3] = 10;
    brickFaceEdgeMap[5][0] = 3; brickFaceEdgeMap[5][1] = 8; brickFaceEdgeMap[5][2] = 7; brickFaceEdgeMap[5][3] = 11;
    
    elmDescription = [[FEMElementDescription alloc] init];
    elementMaps = [[FEMPElementMaps alloc] init];
    
    elements = mesh.getElements;
    edges = mesh.getEdges;
    edges = (Element_t*)malloc(sizeof(Element_t) * (12*mesh.numberOfBulkElements) );
    faces = mesh.getFaces;
    
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        if (element->EdgeIndexes == NULL) element->EdgeIndexes = intvec(0, 11);
        memset( element->EdgeIndexes, 0, (12*sizeof(element->EdgeIndexes)) );
    }
    
    hashTable = (HashTable_t*)malloc(sizeof(HashTable_t) * (mesh.numberOfNodes) );
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashTable[i].head = NULL;
    }

    // Loop over elements
    numbOfEdges = 0;
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        
        // For p elements, mappings are different
        if (element->Pdefs != NULL) {
            [elementMaps getEdgeMapForElement:element edgeMap:edgeMap];
            [elementMaps getFaceEdgeMapForElement:element faceEdgeMap:faceEdgeMap];
            n = element->Type.NumberOfEdges;
        } else {
            switch (element->Type.ElementCode / 100) {
                case 5:
                    n = 6;
                    edgeMap = tetraEdgeMap;
                    faceEdgeMap = tetraFaceEdgeMap;
                    break;
                case 6:
                    n = 8;
                    edgeMap = pyramidEdgeMap;
                    faceEdgeMap = pyramidFaceEdgeMap;
                    break;
                case 7:
                    n = 9;
                    edgeMap = wedgeEdgeMap;
                    faceEdgeMap = wedgeFaceEdgeMap;
                    break;
                case 8:
                    n = 12;
                    edgeMap = brickEdgeMap;
                    faceEdgeMap = brickFaceEdgeMap;
                    break;
                default:
                    printf("Element type %d not implemented.\n", element->Type.ElementCode);
                    errorfunct("findEdges3DInMesh", "Program terminating now.");
                    break;
            }
        }
        
        // Loop over every edge of every element
        for (k=0; k<n; k++) {
         
            // Use min(node1,node2) as key to hash table
            n1 = element->NodeIndexes[edgeMap[k][0]];
            n2 = element->NodeIndexes[edgeMap[k][1]];
            if (n1 < n2) {
                node1 = n1;
                node2 = n2;
            } else {
                node1 = n2;
                node2 = n1;
            }
            
            // Look the edge from the hash table
            hashPtr = hashTable[node1].head;
            found = NO;
            while (hashPtr != NULL) {
                if (hashPtr->node1 == node2) {
                    found = YES;
                    edge = hashPtr->edge;
                    break;
                }
                hashPtr = hashPtr->next;
            }
            
            // Existing edge, update structures
            if (found) {
                element->EdgeIndexes[k] = edge;
                
                // Mark as an edge of pyramid square face
                if ([elementMaps isPPyramid:element] == YES && k < 4) {
                    edges[edge].Pdefs->PyramidQuadEdge = true;
                }
                
                if (faces != NULL) {
                    for (ii=0; ii<element->Type.NumberOfFaces; ii++) {
                        if (faces[element->FaceIndexes[ii]].EdgeIndexes == NULL) {
                            faces[element->FaceIndexes[ii]].EdgeIndexes = intvec(0, faces[element->FaceIndexes[ii]].Type.NumberOfEdges-1);
                            memset( faces[element->FaceIndexes[ii]].EdgeIndexes, 0, (faces[element->FaceIndexes[ii]].Type.NumberOfEdges*sizeof(faces[element->FaceIndexes[ii]].EdgeIndexes)) );
                        }
                        for (jj=0; jj<faces[element->FaceIndexes[ii]].Type.NumberOfEdges; jj++) {
                            if (faceEdgeMap[ii][jj] == k) {
                                faces[element->FaceIndexes[ii]].EdgeIndexes[jj] = edge;
                                if (edges[edge].BoundaryInfo->Left == NULL) {
                                    edges[edge].BoundaryInfo->Left = &faces[element->FaceIndexes[ii]];
                                } else {
                                    edges[edge].BoundaryInfo->Right = &faces[element->FaceIndexes[ii]];
                                }
                            }
                        }
                    }
                }
            } else {
                
                //Edge not yet there, create it
                edge = numbOfEdges;
                edges[edge].ElementIndex = edge+1;
                degree = element->Type.BasisFunctionDegree;
                
                // Edge is always a line segment with deg+1 nodes
                elmType = [elmDescription getElementType:(201+degree) inMesh:mesh stabilization:NO];
                edges[edge].Type = *elmType;
                
                edges[edge].NDOFs = 0;
                if (element->NDOFs != 0) edges[edge].NDOFs = edges[edge].Type.NumberOfNodes;
                edges[edge].BDOFs = 0;
                edges[edge].DGDOFs = 0;
                edges[edge].EdgeIndexes = NULL;
                edges[edge].FaceIndexes = NULL;
                
                edges[edge].NodeIndexes = intvec(0, (degree+1)-1);
                for (n2=0; n2<degree+1; n2++) {
                    edges[edge].NodeIndexes[n2] = element->NodeIndexes[edgeMap[k][n2]];
                }
                
                element->EdgeIndexes[k] = edge;
                edges[edge].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(edges[edge].BoundaryInfo);
                edges[edge].BoundaryInfo->Left = NULL;
                edges[edge].BoundaryInfo->Right = NULL;
                
                // Allocate p element definitions
                if (element->Pdefs != NULL) {
                    edges[edge].Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
                    edges[edge].Pdefs->p = 0;
                    edges[edge].Pdefs->PyramidQuadEdge = false;
                    // Here mark edge as edge of pyramid if needed (or set as not)
                    if ([elementMaps isPPyramid:element] == YES && k < 4) {
                        edges[edge].Pdefs->PyramidQuadEdge = true;
                    }
                } else {
                    edges[edge].Pdefs = NULL;
                }
                
                if (faces != NULL) {
                    for (ii=0; ii<element->Type.NumberOfFaces; ii++) {
                        if (faces[element->FaceIndexes[ii]].EdgeIndexes == NULL) {
                            faces[element->FaceIndexes[ii]].EdgeIndexes = intvec(0, faces[element->FaceIndexes[ii]].Type.NumberOfEdges-1);
                            memset(faces[element->FaceIndexes[ii]].EdgeIndexes, 0, (faces[element->FaceIndexes[ii]].Type.NumberOfEdges*sizeof(faces[element->FaceIndexes[ii]].EdgeIndexes)) );
                        }
                        for (jj=0; jj<faces[element->FaceIndexes[ii]].Type.NumberOfEdges; jj++) {
                            if (faceEdgeMap[ii][jj] == k) {
                                faces[element->FaceIndexes[ii]].EdgeIndexes[jj] = edge;
                                if (edges[edge].BoundaryInfo->Left == NULL) {
                                    edges[edge].BoundaryInfo->Left = &faces[element->FaceIndexes[ii]];
                                } else {
                                    edges[edge].BoundaryInfo->Right = &faces[element->FaceIndexes[ii]];
                                }
                            }
                        }
                    }
                }
                
                // Update the hash table
                hashPtr = (HashEntry_t *)malloc(sizeof(HashEntry_t));
                hashPtr->edge = edge;
                hashPtr->node1 = node2;
                hashPtr->next = hashTable[node1].head;
                hashTable[node1].head = hashPtr;

                
                numbOfEdges++;
            }
        }
    }
    
    mesh.numberOfEdges = numbOfEdges;
    
    // Delete the hash table
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashPtr = hashTable[i].head;
        while (hashPtr != NULL) {
            hashPtr1 = hashPtr->next;
            free(hashPtr);
            hashPtr = hashPtr1;
        }
    }
    free(hashTable);
    
    if (faces != NULL) [self FEMMeshUtils_fixFaceEdges:mesh];
    
    free_imatrix(tetraEdgeMap, 0, 5, 0, 2);
    free_imatrix(tetraFaceMap, 0, 3, 0, 5);
    free_imatrix(brickEdgeMap, 0, 5, 0, 8);
    free_imatrix(wedgeEdgeMap, 0, 4, 0, 3);
    free_imatrix(pyramidEdgeMap, 0, 4, 0, 7);
    
    free_imatrix(tetraFaceEdgeMap, 0, 3, 0, 2);
    free_imatrix(brickFaceEdgeMap, 0, 7, 0, 3);
    free_imatrix(wedgeFaceEdgeMap, 0, 5, 0, 3);
    free_imatrix(pyramidFaceEdgeMap, 0, 4, 0, 3);

    [elmDescription deallocation];
    [elementMaps deallocation];
}

/************************************************************************
 
 Find 3D mesh faces.
 
 ************************************************************************/
-(void)findFaces3DInMesh:(FEMMesh *)mesh {
    
    int i, j, k, n, numbOfFaces;
    int n1, n2;
    int node1, node2, node3, face, degree;
    int **faceMap;
    int **tetraFaceMap, **brickFaceMap, **wedgeFaceNap, **pyramidFaceMap;
    int *nf;
    HashTable_t *hashTable;
    HashEntry_t *hashPtr, *hashPtr1;
    Element_t *elements, *element, *faces;
    ElementType_t *elmType;
    FEMPElementMaps *elementMaps;
    FEMElementDescription *elmDescription;
    BOOL found;
    
    elementMaps = [[FEMPElementMaps alloc] init];
    elmDescription = [[FEMElementDescription alloc] init];
    
    tetraFaceMap = intmatrix(0, 3, 0, 5);
    brickFaceMap = intmatrix(0, 5, 0, 8);
    wedgeFaceNap = intmatrix(0, 4, 0, 3);
    pyramidFaceMap = intmatrix(0, 4, 0, 7);
    nf = intvec(0, 3);
    
    tetraFaceMap[0][0] = 0; tetraFaceMap[0][1] = 1; tetraFaceMap[0][2] = 2; tetraFaceMap[0][3] = 4; tetraFaceMap[0][4] = 5; tetraFaceMap[0][5] = 6;
    tetraFaceMap[1][0] = 0; tetraFaceMap[1][1] = 1; tetraFaceMap[1][2] = 3; tetraFaceMap[1][3] = 4; tetraFaceMap[1][4] = 8; tetraFaceMap[1][5] = 7;
    tetraFaceMap[2][0] = 1; tetraFaceMap[2][1] = 2; tetraFaceMap[2][2] = 3; tetraFaceMap[2][3] = 5; tetraFaceMap[2][4] = 9; tetraFaceMap[2][5] = 8;
    tetraFaceMap[3][0] = 2; tetraFaceMap[3][1] = 0; tetraFaceMap[3][2] = 3; tetraFaceMap[3][3] = 6; tetraFaceMap[3][4] = 7; tetraFaceMap[3][5] = 9;
    
    wedgeFaceNap[0][0] = 0; wedgeFaceNap[0][1] = 1; wedgeFaceNap[0][2] = 2; wedgeFaceNap[0][3] = -1;
    wedgeFaceNap[1][0] = 3; wedgeFaceNap[1][1] = 4; wedgeFaceNap[1][2] = 5; wedgeFaceNap[1][3] = -1;
    wedgeFaceNap[2][0] = 0; wedgeFaceNap[2][1] = 1; wedgeFaceNap[2][2] = 4; wedgeFaceNap[2][3] = 3;
    wedgeFaceNap[3][0] = 2; wedgeFaceNap[3][1] = 1; wedgeFaceNap[3][2] = 4; wedgeFaceNap[3][3] = 5;
    wedgeFaceNap[4][0] = 2; wedgeFaceNap[4][1] = 0; wedgeFaceNap[4][2] = 3; wedgeFaceNap[4][3] = 5;
    
    pyramidFaceMap[0][0] = 0; pyramidFaceMap[0][1] = 1; pyramidFaceMap[0][2] = 2; pyramidFaceMap[0][3] = 3; pyramidFaceMap[0][4] = 5;
    pyramidFaceMap[0][5] = 6; pyramidFaceMap[0][6] = 7; pyramidFaceMap[0][7] = 8;
    pyramidFaceMap[1][0] = 0; pyramidFaceMap[1][1] = 1; pyramidFaceMap[1][2] = 4; pyramidFaceMap[1][3] = 5; pyramidFaceMap[1][4] = 10;
    pyramidFaceMap[1][5] = 9; pyramidFaceMap[1][6] = -1; pyramidFaceMap[1][7] = -1;
    pyramidFaceMap[2][0] = 1; pyramidFaceMap[2][1] = 2; pyramidFaceMap[2][2] = 4; pyramidFaceMap[2][3] = 6; pyramidFaceMap[2][4] = 11;
    pyramidFaceMap[2][5] = 10; pyramidFaceMap[2][6] = -1; pyramidFaceMap[2][7] = -1;
    pyramidFaceMap[3][0] = 2; pyramidFaceMap[3][1] = 3; pyramidFaceMap[3][2] = 4; pyramidFaceMap[3][3] = 7; pyramidFaceMap[3][4] = 12;
    pyramidFaceMap[3][5] = 11; pyramidFaceMap[3][6] = -1; pyramidFaceMap[3][7] = -1;
    pyramidFaceMap[4][0] = 3; pyramidFaceMap[4][1] = 0; pyramidFaceMap[4][2] = 4; pyramidFaceMap[4][3] = 8; pyramidFaceMap[4][4] = 9;
    pyramidFaceMap[4][5] = 12; pyramidFaceMap[4][6] = -1; pyramidFaceMap[4][7] = -1;
    
    brickFaceMap[0][0] = 0; brickFaceMap[0][1] = 1; brickFaceMap[0][2] = 2; brickFaceMap[0][3] = 3; brickFaceMap[0][4] = 8; brickFaceMap[0][5] = 9;
    brickFaceMap[0][6] = 10; brickFaceMap[0][7] = 11; brickFaceMap[0][8] = 24;
    brickFaceMap[1][0] = 4; brickFaceMap[1][1] = 5; brickFaceMap[1][2] = 6; brickFaceMap[1][3] = 7; brickFaceMap[1][4] = 16; brickFaceMap[1][5] = 17;
    brickFaceMap[1][6] = 18; brickFaceMap[1][7] = 19; brickFaceMap[1][8] = 25;
    brickFaceMap[2][0] = 0; brickFaceMap[2][1] = 1; brickFaceMap[2][2] = 5; brickFaceMap[2][3] = 4; brickFaceMap[2][4] = 8; brickFaceMap[2][5] = 13;
    brickFaceMap[2][6] = 16; brickFaceMap[2][7] = 12; brickFaceMap[2][8] = 20;
    brickFaceMap[3][0] = 1; brickFaceMap[3][1] = 2; brickFaceMap[3][2] = 6; brickFaceMap[3][3] = 5; brickFaceMap[3][4] = 9; brickFaceMap[3][5] = 14;
    brickFaceMap[3][6] = 16; brickFaceMap[3][7] = 13; brickFaceMap[3][8] = 21;
    brickFaceMap[4][0] = 2; brickFaceMap[4][1] = 3; brickFaceMap[4][2] = 7; brickFaceMap[4][3] = 6; brickFaceMap[4][4] = 10; brickFaceMap[4][5] = 15;
    brickFaceMap[4][6] = 18; brickFaceMap[4][7] = 14; brickFaceMap[4][8] = 22;
    brickFaceMap[5][0] = 3; brickFaceMap[5][1] = 0; brickFaceMap[5][2] = 4; brickFaceMap[5][3] = 7; brickFaceMap[5][4] = 11; brickFaceMap[5][5] = 12;
    brickFaceMap[5][6] = 19; brickFaceMap[5][7] = 15; brickFaceMap[5][8] = 23;
    
    elements = mesh.getElements;
    faces = mesh.getFaces;
    faces = (Element_t*)malloc(sizeof(Element_t) * (6*mesh.numberOfBulkElements) );
    
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        if (element->FaceIndexes == NULL) element->FaceIndexes = intvec(0, 5);
        memset( element->FaceIndexes, 0, (6*sizeof(element->FaceIndexes)) );
    }
    
    hashTable = (HashTable_t*)malloc(sizeof(HashTable_t) * (mesh.numberOfNodes) );
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashTable[i].head = NULL;
    }
    
    // Loop over elements
    numbOfFaces = 0;
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        
        // For p elements, mappings are different
        if (element->Pdefs != NULL) {
            [elementMaps getFaceMapForElement:element faceMap:faceMap];
            n = element->Type.NumberOfFaces;
        } else {
            switch (element->Type.ElementCode / 100) {
                case 5:
                    n = 4;
                    faceMap = tetraFaceMap;
                    break;
                case 6:
                    n = 5;
                    faceMap = pyramidFaceMap;
                    break;
                case 7:
                    n = 5;
                    faceMap = wedgeFaceNap;
                    break;
                case 8:
                    n = 6;
                    faceMap = brickFaceMap;
                    break;
                default:
                    continue;
                    break;
            }
        }
        
        // Loop over every face of every element
        for (k=0; k<n; k++) {
            
            // We use min(node1,node2,node3) as the hash table key
            switch (element->Type.ElementCode / 100) {
                case 5:
                    
                    // Tetras:
                    // =======
                    for (j=0; j<3; j++) {
                        nf[j] = element->NodeIndexes[faceMap[k][j]];
                    }
                    sort(3, nf-1);
                    break;
                    
                case 6:
                    
                    // Pyramids:
                    // ========
                    if (k == 0) {
                        for (j=0; j<4; j++) {
                            nf[j] = element->NodeIndexes[faceMap[k][j]];
                        }
                        sort(4, nf-1);
                    } else {
                        for (j=0; j<3; j++) {
                            nf[j] = element->NodeIndexes[faceMap[k][j]];
                        }
                        sort(3, nf-1);
                    }
                    break;
                    
                case 7:
                    
                    // Wedges:
                    // ======
                    if (k <= 1) {
                        for (j=0; j<3; j++) {
                            nf[j] = element->NodeIndexes[faceMap[k][j]];
                        }
                        sort(3, nf-1);
                    } else {
                        for (j=0; j<4; j++) {
                            nf[j] = element->NodeIndexes[faceMap[k][j]];
                        }
                        sort(4, nf-1);
                    }
                    break;
                    
                case 8:
                    
                    // Bricks:
                    // ======
                    for (j=0; j<4; j++) {
                        nf[j] = element->NodeIndexes[faceMap[k][j]];
                    }
                    sort(4, nf-1);
                    break;
                    
                default:
                    printf("Element type %d not implemented.\n", element->Type.ElementCode);
                    errorfunct("findFaces3DForMesh", "Program terminating now.");
                    break;
            }
            
            node1 = nf[0];
            node2 = nf[1];
            node3 = nf[2];
            
            // Look the face from the hash table
            hashPtr = hashTable[node1].head;
            found = NO;
            while (hashPtr != NULL) {
                if (hashPtr->node1 == node2 && hashPtr->node2 == node3) {
                    found = YES;
                    face = hashPtr->face;
                    break;
                }
                hashPtr = hashPtr->next;
            }
            
            // Existing face, update structures
            if (found == YES) {
                element->FaceIndexes[k] = face;
                faces[face].BoundaryInfo->Right = element;
            } else {
                
                // Face not yet there, create it
                face = numbOfFaces;
                faces[face].ElementIndex = face+1;
                
                degree = element->Type.BasisFunctionDegree;
                switch (element->Type.ElementCode / 100) {
                    case 5:
                        // Tetras:
                        // ======
                        switch (degree) {
                            case 1:
                                n1 = 3;
                                break;
                            case 2:
                                n1 = 6;
                                break;
                            case 3:
                                n1 = 10;
                                break;
                        }
                        n1 = 3; // TODO: The switch stastement above is then useless?
                        elmType = [elmDescription getElementType:(300+n1) inMesh:mesh stabilization:NO];
                        faces[face].Type = *elmType;
                        break;
                        
                    case 6:
                        // Pyramid (only 605 supported)
                        // ============================
                        if (k == 0) {
                            n1 = 4 ;
                            elmType = [elmDescription getElementType:(400+n1) inMesh:mesh stabilization:NO];
                            faces[face].Type = *elmType;
                        } else {
                            n1 = 3;
                            elmType = [elmDescription getElementType:(300+n1) inMesh:mesh stabilization:NO];
                            faces[face].Type = *elmType;;
                        }
                        break;
                        
                    case 7:
                        // Wedges (only 706 supported)
                        // ===========================
                        if (k <= 1) {
                            n1 = 3;
                            elmType = [elmDescription getElementType:303 inMesh:mesh stabilization:NO];
                            face[faces].Type = *elmType;
                        } else {
                            n1 = 4;
                            elmType = [elmDescription getElementType:404 inMesh:mesh stabilization:NO];
                            faces[face].Type = *elmType;
                        }
                        break;
                        
                    case 8:
                        // Bricks
                        // ======
                        switch (element->Type.NumberOfNodes) {
                            case 8:
                                n1 = 4;
                                break;
                            case 20:
                                n1 = 8;
                                break;
                            case 27:
                                n1 = 9;
                                break;
                        }
                        elmType = [elmDescription getElementType:(400+n1) inMesh:mesh stabilization:NO];
                        faces[face].Type = *elmType;
                        
                    default:
                        printf("Element type %d not implemented.\n", element->Type.ElementCode);
                        errorfunct("findFaces3DForMesh", "Program terminating now.");
                        break;
                }
                
                // Create P element definitions if necessary
                if (element->Pdefs != NULL) {
                    faces[face].Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
                    faces[face].Pdefs->p = 0;
                } else {
                    faces[face].Pdefs = NULL;
                }
                
                faces[face].NDOFs = 0;
                if (element->NDOFs != 0) faces[face].NDOFs = faces[face].Type.NumberOfNodes;
                faces[face].BDOFs = 0;
                faces[face].DGDOFs = 0;
                faces[face].EdgeIndexes = NULL;
                faces[face].FaceIndexes = NULL;
                
                faces[face].NodeIndexes = intvec(0, n1-1);
                for (n2=0; n2<n1; n2++) {
                    faces[face].NodeIndexes[n2] = element->NodeIndexes[faceMap[k][n2]];
                }
                
                element->FaceIndexes[k] = face;
                
                faces[face].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(faces[face].BoundaryInfo);
                faces[face].BoundaryInfo->Left = element;
                faces[face].BoundaryInfo->Right = NULL;
                
                // Update the hash table
                hashPtr = (HashEntry_t *)malloc(sizeof(HashEntry_t));
                hashPtr->face = face;
                hashPtr->node1 = node2;
                hashPtr->node2 = node3;
                hashPtr->next = hashTable[node1].head;
                hashTable[node1].head = hashPtr;
                
                numbOfFaces++;
            }
        }
    }
    
    mesh.numberOfFaces = numbOfFaces;
    
    // Delete the hash table
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashPtr = hashTable[i].head;
        while (hashPtr != NULL) {
            hashPtr1 = hashPtr->next;
            free(hashPtr);
            hashPtr = hashPtr1;
        }
    }
    free(hashTable);
    
    free_imatrix(tetraFaceMap, 0, 3, 0, 5);
    free_imatrix(brickFaceMap, 0, 5, 0, 8);
    free_imatrix(wedgeFaceNap, 0, 4, 0, 3);
    free_imatrix(pyramidFaceMap, 0, 4, 0, 7);
    free_ivector(nf, 0, 3);
    
    [elementMaps deallocation];
    [elmDescription deallocation];
}

/************************************************************************
 
    Generate element edge (face in 3D) tables for given mesh.
    Only for triangles and tetras. If mesh already has edges, do nothing

************************************************************************/
-(void)findEdgesForMesh:(FEMMesh *)mesh findEdges:(BOOL *)present {
    
    Element_t *edges, *faces;
    BOOL findEdges3D;
    
    if (present != NULL) {
        findEdges3D = *present;
    } else {
        findEdges3D = YES;
    }
    
    edges = mesh.getEdges;
    faces = mesh.getFaces;
    
    switch (mesh.dimension) {
        case 2:
            if (edges == NULL) [self findEdges2DInMesh:mesh];
            break;
        case 3:
            if (faces == NULL) [self findFaces3DInMesh:mesh];
            if (findEdges3D == YES) {
                if (edges == NULL) [self findEdges3DInMesh:mesh];
            }
            break;
    }
    
    [self FEMMeshUtils_assignConstraints:mesh];
    
    edges = NULL;
    faces = NULL;
}

/*************************************************************************
    
    Assign local number of edge to given boundary element. Also copies all
    p element attributes from element edge to boundary edge

*************************************************************************/
-(void)assignLocalNumberToEdgeElement:(Element_t *)edge fromElement:(Element_t *)element inMesh:(FEMMesh *)mesh {
    
    int i, j, n, edgeNumber, numbEdges, bMap[4];
    Element_t *entity;
    FEMPElementMaps *pMaps;
    
    pMaps = [[FEMPElementMaps alloc] init];
    
    numbEdges = 0;
    switch (element->Type.dimension) {
        case 2:
            numbEdges = element->Type.NumberOfEdges;
            break;
        case 3:
            numbEdges = element->Type.NumberOfFaces;
            break;
        default:
            NSLog(@"assignLocalNumberToEdgeElement: Unsupported dimension.");
            return;
            break;
    }
    
    // For each edge or face in element, try to find local number
    for (edgeNumber=0; edgeNumber<numbEdges; edgeNumber++) {
        // If edges have been created, stop search. Actually, this should not happen
        if (element->EdgeIndexes == NULL) {
             NSLog(@"assignLocalNumberToEdgeElement: edges have not been creates. Returning now!");
            return;
        }
        
        entity = [self FEMMeshUtils_getEntityForElement:element edge:edgeNumber inMesh:mesh];
        
        // Edge element not found. This should not be possible, unless there is
        // an error in the mesh read in process
        if (entity == NULL) {
            NSLog(@"assignLocalNumberToEdgeElement: edge element not found.");
            return;
        }
        
        n = 0;
        // For each element node
        for (i=0; i<entity->Type.NumberOfNodes; i++) {
            // For each node in edge element
            for (j=0; j<edge->Type.NumberOfNodes; j++) {
                // If entity and edge element node match, increment counter
                if (entity->NodeIndexes[i] == edge->NodeIndexes[j]) n++;
            }
        }
        
        // If all nodes are on boundary, edge was found
        if (n == edge->Type.NumberOfNodes) {
            edge->Pdefs->LocalNumber = edgeNumber;
            
            // Change ordering of global nodes to match that of element
            [pMaps getBoundaryMapForElement:element localNumber:edgeNumber resultMap:bMap];
            for (i=0; i<n; i++) {
                edge->NodeIndexes[i] = element->NodeIndexes[bMap[i]];
            }
            
            // Copy attributes of edge element to boundary element
            // Misc attributes
            edge->Pdefs->isEdge = entity->Pdefs->isEdge;
            
            // Gauss points
            edge->Pdefs->GaussPoints = entity->Pdefs->GaussPoints;
            
            // Element p (and boundary bubble dofs)
            edge->BDOFs = entity->BDOFs;
            edge->Pdefs->p = entity->Pdefs->p;
            
            // If this boundary has edges copy edge indexes
            if (entity->EdgeIndexes != NULL) {
                // Allocate element edges to element
                n = entity->Type.NumberOfEdges;
                [pMaps getFaceEdgeMapForElement:element index:edgeNumber resultMap:bMap];
             
                if (edge->EdgeIndexes != NULL) {
                    free_ivector(edge->EdgeIndexes, 0, edge->sizeEdgeIndexes-1);
                }
                edge->EdgeIndexes = intvec(0, n-1);
                edge->sizeEdgeIndexes = n;
                // Copy edges from edge to boundary edge
                for (i=0; i<n; i++) {
                    edge->EdgeIndexes[i] = element->EdgeIndexes[bMap[i]];
                }
            }
            
            // Edge fields copied and local edge found so return
            return;
        }
    }
    
    // If we are here, local number not found
    NSLog(@"assignLocalNumberToEdgeElement: Unable to find local edge.");
}

-(FEMMatrix *)periodicProjectorInModel:(FEMModel *)model forMesh:(FEMMesh *)mesh boundary:(int)this target:(int)trgt {
    
    int i, j, k, n, n1, n2, k1, k2, constraint, dim, numbBothActive, count;
    int *perm1, *perm2, *invPerm1, *invPerm2;
    double dot1min, dot2min, sum, alpha, s1, s2;
    double normals[3], normals1[3], normals2[3], x1Min[3], x1Max[3], x2Min[3], x2Max[3], x2rMin[3], x2rMax[3], x[4], identity[4][4], rotMatrix[4][4],
           trsMatrix[4][4], sclMatrix[4][4], trfMatrix[4][4], angles[3], normals2r[3], scl[3];
    listBuffer pArray = { NULL, NULL, NULL, NULL, 0, 0, 0};
    FEMMesh *bmesh1, *bmesh2;
    FEMMatrix *projector;
    FEMBoundaryCondition *boundaryConditionAtId;
    FEMElementDescription *elmDescription;
    FEMListUtilities *listUtil;
    FEMUtilities *utilities;
    Element_t *elements, *mesh1elements, *mesh2elements;
    Nodes_t *nodes, *elementNodes, *mesh1nodes, *mesh2nodes;
    matrixArraysContainer *projectorContainers = NULL;
    BOOL thisActive, targetActive, constantNormals, check, gotRotate, useQuadrantTree, found;
    
    dim = model.dimension;
    
    if (this < 0 || trgt < 0) return nil;
    
    elmDescription = [[FEMElementDescription alloc] init];
    listUtil = [[FEMListUtilities alloc] init];
    utilities = [[FEMUtilities alloc] init];
    
    elements = mesh.getElements;
    nodes = mesh.getNodes;
    
    // Search elements in this boundary and its periodic counterpart:
    n1 = 0;
    n2 = 0;
    for (i=mesh.numberOfBulkElements; i<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; i++) {
        constraint = elements[i].BoundaryInfo->Constraint;
        boundaryConditionAtId = (model.boundaryConditions)[this];
        if (boundaryConditionAtId.tag == constraint) n1++;
        
        boundaryConditionAtId = (model.boundaryConditions)[trgt];
        if (boundaryConditionAtId.tag == constraint) n2++;
    }
    
    if (n1 <= 0 || n2 <= 0) return nil;
    NSLog(@"FEMMesh_periodicProjector: starting to build...");
    
    // Initialize mesh structures for boundaries, this is
    // for getting the mesh projector
    bmesh1 = [[FEMMesh alloc] init];
    bmesh2 = [[FEMMesh alloc] init];
    
    mesh1elements = bmesh1.getElements;
    mesh2elements = bmesh2.getElements;
    
    mesh1elements = (Element_t*) malloc( sizeof(Element_t) * n1 );
    initElements(mesh1elements, n1);
    mesh2elements = (Element_t*) malloc( sizeof(Element_t) * n2 );
    initElements(mesh2elements, n2);
    
    perm1 = intvec(0, mesh.numberOfNodes-1);
    perm2 = intvec(0, mesh.numberOfNodes-1);
    
    n = mesh.maxElementNodes;
    
    elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(elementNodes);
    elementNodes->x = doublevec(0, n-1);
    elementNodes->y = doublevec(0, n-1);
    elementNodes->z = doublevec(0, n-1);
    
    // Fill in the mesh element structures with the boundary elements
    n1 = 0;
    n2 = 0;
    memset( perm1, 0, (mesh.numberOfNodes*sizeof(perm1)) );
    memset( perm2, 0, (mesh.numberOfNodes*sizeof(perm2)) );
    bmesh1.maxElementNodes = 0;
    bmesh2.maxElementNodes = 0;
    dot1min = HUGE_VAL;
    dot2min = HUGE_VAL;
    numbBothActive = 0;
    
    for (i=mesh.numberOfBulkElements; i<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; i++) {
        
        if (dim > 1 && elements[i].Type.ElementCode < 200) continue;
        
        constraint = elements[i].BoundaryInfo->Constraint;
        
        boundaryConditionAtId = (model.boundaryConditions)[this];
        thisActive = (boundaryConditionAtId.tag == constraint) ? YES : NO;
        
        boundaryConditionAtId = (model.boundaryConditions)[trgt];
        targetActive = (boundaryConditionAtId.tag == constraint) ? YES : NO;
        
        if (thisActive == YES || targetActive == YES) {
            n = elements[i].Type.NumberOfNodes;
            for (j=0; j<n; j++) {
                elementNodes->x[j] = nodes->x[elements[i].NodeIndexes[j]];
                elementNodes->y[j] = nodes->y[elements[i].NodeIndexes[j]];
                elementNodes->z[j] = nodes->z[elements[i].NodeIndexes[j]];
            }
            
            // Angle smaller than 180 is always chosen
            check = NO;
            [elmDescription normalVectorForBDElement:&elements[i] boundaryNodes:elementNodes mesh:mesh paraU:NULL paraV:NULL check:&check normals:normals];
            
            if (thisActive == YES && targetActive == YES) {
                numbBothActive++;
            }
        }
        
        if (thisActive == YES) {
            if (n1 == 0) {
                memcpy(normals1, normals, sizeof(normals1));
            } else {
                sum = 0.0;
                for (j=0; j<3; j++) {
                    sum = sum + normals[j]*normals1[j];
                }
                dot1min = min(dot1min, sum);
            }
            
            bmesh1.maxElementNodes = max(bmesh1.maxElementNodes, n);
            mesh1elements[n1] = elements[i];
            mesh1elements[n1].NodeIndexes = intvec(0, n-1);
            mesh1elements[n1].sizeNodeIndexes = n;
            mesh1elements[n1].NodeIndexes = elements[i].NodeIndexes;
            mesh1elements[n1].EdgeIndexes = NULL;
            mesh1elements[n1].FaceIndexes = NULL;
            for (j=0; j<n; j++) {
                perm1[elements[i].NodeIndexes[j]] = 1;
            }
            n1++;
        }
        
        if (targetActive == YES) {
            if (n2 == 0) {
                memcpy(normals2, normals, sizeof(normals2));
            } else {
                sum = 0.0;
                for (j=0; j<3; j++) {
                    sum = sum + normals[j]*normals2[j];
                }
                dot2min = min(dot2min, sum);
            }
            
            bmesh2.maxElementNodes = max(bmesh2.maxElementNodes, n);
            mesh2elements[n2] = elements[i];
            mesh2elements[n2].NodeIndexes = intvec(0, n-1);
            mesh2elements[n2].sizeNodeIndexes = n;
            mesh2elements[n2].NodeIndexes = elements[i].NodeIndexes;
            mesh2elements[n2].EdgeIndexes = NULL;
            mesh2elements[n2].FaceIndexes = NULL;
            for (j=0; j<n; j++) {
                perm2[elements[i].NodeIndexes[j]] = 1;
            }
            n2++;
        }
    }
    
    if (numbBothActive > 0) {
        NSLog(@"periodicProjectorInModel: Nodes belonging to both master and target: %d\n", numbBothActive);
    }
    
    constantNormals = ( (1.0 - dot1min < 1.0e-6) && (1.0 - dot2min < 1.0e-6) ) ? YES : NO;
    if (constantNormals == YES) {
        // The full angle between the two normals
        sum = 0.0;
        for (i = 0; i<3; i++) {
            sum = sum + normals1[i] * normals2[i];
        }
        alpha = acos(sum) * 180 / pi;
        NSLog(@"periodicProjectorInModel: Suggested angle between two normals in degrees: %e\n", alpha);
    }
    
    bmesh1.numberOfBulkElements = n1;
    bmesh2.numberOfBulkElements = n2;
    
    // Fill in the mesh node structures with the boundary nodes
    count = 0;
    for (i=0; i<mesh.numberOfNodes; i++) {
        if (perm2[i] == 1) count++;
    }
    bmesh2.numberOfNodes = count;
    
    count = 0;
    for (i=0; i<mesh.numberOfNodes; i++) {
        if (perm1[i] == 1) count++;
    }
    bmesh1.numberOfNodes = count;
    
    if (bmesh1.numberOfNodes == 0 || bmesh2.numberOfNodes == 0) {
        NSLog(@"periodicProjectorInModel: No active nodes on periodic boundary!");
        
        for (i=0; i<bmesh1.numberOfBulkElements; i++) {
            free_ivector(mesh1elements[i].NodeIndexes, 0, mesh1elements[i].sizeNodeIndexes-1);
        }
        for (i=0; i<bmesh2.numberOfBulkElements; i++) {
            free_ivector(mesh2elements[i].NodeIndexes, 0, mesh1elements[i].sizeNodeIndexes-1);
        }
        
        free(mesh1elements);
        free(mesh2elements);
        free_ivector(perm1, 0, mesh.numberOfNodes-1);
        free_ivector(perm2, 0, mesh.numberOfNodes-1);
        perm1 = NULL; perm2 = NULL;
        bmesh1 = nil;
        bmesh2 = nil;
    } else {
        NSLog(@"periodicProjectorInModel: Number of periodic nodes: %d %d\n", bmesh1.numberOfNodes, bmesh2.numberOfNodes);
    }
    
    mesh1nodes = bmesh1.getNodes;
    mesh2nodes = bmesh2.getNodes;
    
    mesh1nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(mesh1nodes);
    mesh1nodes->x = doublevec(0, bmesh1.numberOfNodes-1);
    mesh1nodes->y = doublevec(0, bmesh1.numberOfNodes-1);
    mesh1nodes->z = doublevec(0, bmesh1.numberOfNodes-1);
    
    mesh2nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(mesh2nodes);
    mesh2nodes->x = doublevec(0, bmesh2.numberOfNodes-1);
    mesh2nodes->y = doublevec(0, bmesh2.numberOfNodes-1);
    mesh2nodes->z = doublevec(0, bmesh2.numberOfNodes-1);
    
    invPerm1 = intvec(0, bmesh1.numberOfNodes-1);
    invPerm2 = intvec(0, bmesh2.numberOfNodes-1);
    
    for (i=0; i<3; i++) {
        x1Min[i] = HUGE_VAL;
        x1Max[i] = HUGE_VAL;
        x2Min[i] = HUGE_VAL;
        x2Max[i] = HUGE_VAL;
        x2rMin[i] = HUGE_VAL;
        x2rMax[i] = HUGE_VAL;
    }
    
    // Initialize the mapping matrices
    memset( x, 0.0, sizeof(x) );
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            identity[i][j] = 0.0;
        }
    }
    for (i=0; i<4; i++) {
        identity[i][i] = 1.0;
    }
    for (i=0; i<4; i++) {
        for (j=0; j<4; j++) {
            trsMatrix[i][j] = identity[i][j];
            rotMatrix[i][j] = identity[i][j];
            sclMatrix[i][j] = identity[i][j];
        }
    }
    
    // Rotations:
    // These are called first since they are not accounted for in the
    // automatic scaling and translation.
    boundaryConditionAtId = (model.boundaryConditions)[this];
    gotRotate = [listUtil listGetConstRealArray:model inArray:boundaryConditionAtId.valuesList forVariable:@"periodic bc rotate" buffer:&pArray];
    if (gotRotate == YES) {
        for (i=0; i<3; i++) {
            angles[i] = pArray.matrix[i][0];
        }
    } else {
        memset( angles, 0.0, sizeof(angles) );
        if ([listUtil listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"periodic bc rotate automatic" info:&found] == YES) {
            
            if (constantNormals == NO) {
                errorfunct("periodicProjectorInModel", "Normals are not constant, can not test for rotation!");
            } else if (alpha > DBL_EPSILON) {
                // Rotation should be performed
                for (i=0; i<3; i++) {
                    if (fabs(normals1[i] - normals2[i]) < DBL_EPSILON) {
                        gotRotate = YES;
                        NSLog(@"periodicProjectorInModel: Rotation around axis %d in degrees %f.\n", i, alpha);
                        angles[i] = alpha;
                        break;
                    }
                }
                if (gotRotate == NO) {
                    errorfunct("periodicProjectorInModel", "could not define rotation axis, improve algorithm.");
                }
            }
        }
    }
    
    if (gotRotate == YES) {
        NSLog(@"periodicProjectorInModel: rotating target with: %f %f %f\n", angles[0], angles[1], angles[2]);
        

        
        for (i=0; i<pArray.m; i++) {
            alpha = angles[i] * pi / 180;
            if (fabs(alpha) < DBL_MIN) continue;
            
            memcpy(trfMatrix, identity, sizeof(identity));
            switch (i) {
                case 0:
                    trfMatrix[1][1] = cos(alpha);
                    trfMatrix[1][2] = -sin(alpha);
                    trfMatrix[2][1] = sin(alpha);
                    trfMatrix[2][2] = cos(alpha);
                    break;
                case 1:
                    trfMatrix[0][0] = cos(alpha);
                    trfMatrix[0][2] = -sin(alpha);
                    trfMatrix[2][0] = sin(alpha);
                    trfMatrix[2][2] = cos(alpha);
                    break;
                case 2:
                    trfMatrix[0][0] = cos(alpha);
                    trfMatrix[0][1] = -sin(alpha);
                    trfMatrix[1][0] = sin(alpha);
                    trfMatrix[1][1] = cos(alpha);
                    break;
            }
            
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, 4, 4, 1.0, (double *)rotMatrix, 4, (double *)trfMatrix, 4, 0.0, (double *)rotMatrix, 4);
        }
        
        // Check the rotated normal vector of the first target element.
        // This is maily done for debugging purposes
        
        if (constantNormals == YES) {
            for (i=mesh.numberOfBulkElements; i<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; i++) {
                if (dim > 1 && elements[i].Type.ElementCode < 200) continue;
                
                constraint = elements[i].BoundaryInfo->Constraint;
                
                boundaryConditionAtId = (model.boundaryConditions)[trgt];
                if (boundaryConditionAtId.tag != constraint) continue;
                
                n = elements[i].Type.NumberOfNodes;
                
                for (j=0; j<n; j++) {
                    elementNodes->x[j] = nodes->x[elements[i].NodeIndexes[j]];
                    elementNodes->y[j] = nodes->y[elements[i].NodeIndexes[j]];
                    elementNodes->z[j] = nodes->z[elements[i].NodeIndexes[j]];
                }
                
                for (j=0; j<n; j++) {
                    k = elements[i].NodeIndexes[j];
                    x[0] = nodes->x[k];
                    x[1] = nodes->y[k];
                    x[2] = nodes->z[k];
                    x[3] = 1.0;
                    
                    cblas_dgemv(CblasRowMajor, CblasNoTrans, 4, 4, 1.0, (double *)rotMatrix, 4, x, 1, 0.0, x, 1);
                    
                    elementNodes->x[j] = x[0];
                    elementNodes->y[j] = x[1];
                    elementNodes->z[j] = x[2];
                    
                }
                
                check = NO;
                [elmDescription normalVectorForBDElement:&elements[i] boundaryNodes:elementNodes mesh:mesh paraU:NULL paraV:NULL check:&check normals:normals2r];
                break;
            }
            
            NSLog(@"periodicProjectorInModel: master normal: %f %f %f\n", normals1[0], normals1[1], normals1[2]);
            NSLog(@"periodicProjectorInModel: initial target normal: %f %f %f\n", normals2[0], normals2[1], normals2[2]);
            NSLog(@"periodicProjectorInModel: rotated target normal: %f %f %f\n", normals2r[0], normals2r[1], normals2r[2]);
            sum = 0.0;
            for (i=0; i<3; i++) {
                sum = sum + normals1[i] * normals2[i];
            }
            alpha = acos(sum) * 180.0 / pi;
            NSLog(@"periodicProjectorInModel: angle between master and initial target: %e\n", alpha);
            
            sum = 0.0;
            for (i=0; i<3; i++) {
                sum = sum + normals1[i] * normals2r[i];
            }
            alpha = acos(sum) * 180.0 / pi;
            NSLog(@"periodicProjectorInModel: angle between master and rotated target: %e\n", alpha);
            
            if (fabs(alpha) > 1.0e-2 && fabs(alpha - 180.0) > 1.0e-2) {
                NSLog(@"periodicProjectorInModel: rotation may be incorrect!");
            }
        }
    }
    
    // Create the master and target meshes that only include the active elements
    k1 = 0; k2 = 0;
    for (i=0; i<mesh.numberOfNodes; i++) {
        if (perm1[i] > 0) {
            perm1[i] = k1;
            invPerm1[k1] = i;
            
            mesh1nodes->x[k1] = nodes->x[i];
            mesh1nodes->y[k1] = nodes->y[i];
            mesh1nodes->z[k1] = nodes->z[i];
            
            x[0] = nodes->x[i];
            x[1] = nodes->y[i];
            x[2] = nodes->z[i];
            
            for (j=0; j<3; j++) {
                x1Min[j] = min(x1Min[j], x[j]);
                x1Max[j] = max(x1Max[j], x[j]);
            }
            k1++;
        }
        
        if (perm2[i] > 0) {
            perm2[i] = k2;
            invPerm1[k2] = i;
            
            x[0] = nodes->x[i];
            x[1] = nodes->y[i];
            x[2] = nodes->z[i];
            
            for (j=0; j<3; j++) {
                x2Min[j] = min(x2Min[j], x[j]);
                x2Max[j] = max(x2Max[j], x[j]);
            }
            
            if (gotRotate == YES) {
                x[3] = 1.0;
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 4, 4, 1.0, (double *)rotMatrix, 4, x, 1, 0.0, x, 1);
                for (j=0; j<3; j++) {
                    x2rMin[j] = min(x2rMin[j], x[j]);
                    x2rMax[j] = max(x2rMax[j], x[j]);
                }
            }
            
            mesh2nodes->x[k2] = x[0];
            mesh2nodes->y[k2] = x[1];
            mesh2nodes->z[k2] = x[2];
            k2++;
        }
    }
    
    NSLog(@"periodicProjectorInModel: minimum values for this periodic BC: %f %f %f\n", x1Min[0], x1Min[1], x1Min[2]);
    NSLog(@"periodicProjectorInModel: minimum values for target periodic BC: %f %f %f\n", x2Min[0], x2Min[1], x2Min[2]);
    if (gotRotate == YES) {
        NSLog(@"periodicProjectorInModel: minimum values for rotated target: %f %f %f\n", x2rMin[0], x2rMin[1], x2rMin[2]);
    }
    NSLog(@"periodicProjectorInModel: maximum values for this periodic BC: %f %f %f\n", x1Max[0], x1Max[1], x1Max[2]);
    NSLog(@"periodicProjectorInModel: maximum values for target periodic BC: %f %f %f\n", x2Max[0], x2Max[1], x2Max[2]);
    if (gotRotate == YES) {
        NSLog(@"periodicProjectorInModel: maximum values for rotated target: %f %f %f\n", x2rMax[0], x2rMax[1], x2rMax[2]);
    }
    
    NSLog(@"periodicProjectorInModel: bounding box for this periodic BC: %f %f %f\n", x1Max[0]-x1Min[0], x1Max[1]-x1Min[1], x1Max[2]-x1Min[2]);
    NSLog(@"periodicProjectorInModel: bounding box for target periodic BC: %f %f %f\n", x2Max[0]-x2Min[0], x2Max[1]-x2Min[1], x2Max[2]-x2Min[2]);
    if (gotRotate == YES) {
        NSLog(@"periodicProjectorInModel: bounding box for rotated target: %f %f %f\n", x2rMax[0]-x2rMin[0], x2rMax[1]-x2rMin[1], x2rMax[2]-x2rMin[2]);
    }
    
    if (gotRotate == NO) {
        memcpy(x2rMin, x2Min, sizeof(x2Min));
        memcpy(x2rMax, x2Max, sizeof(x2Max));
    }
    
    // Scales
    if (pArray.matrix != NULL) {
        free_dmatrix(pArray.matrix, 0, pArray.m-1, 0, pArray.n-1);
        pArray.matrix = NULL;
    }
    boundaryConditionAtId = (model.boundaryConditions)[this];
    found = [listUtil listGetConstRealArray:model inArray:boundaryConditionAtId.valuesList forVariable:@"periodic bc scale" buffer:&pArray];

    if (found == YES) {
        for (i=0; i<pArray.m; i++) {
            sclMatrix[i][i] = pArray.matrix[i][0];
        }
    } else {
        // Define scaling from the bounding boxes
        if (NO) {
            // This makes scaling component wise, currently disabled
            for (i=0; i<3; i++) {
                scl[i] = x2rMax[i] -  x2rMin[i];
                
                if (fabs(scl[i]) > AEPS) {
                    scl[i] = ( x1Max[i] - x1Min[i] ) / scl[i];
                } else {
                    scl[i] = 1.0;
                }
            }
        } else {
            // This assumes isotropic scaling since the previous was prone
            // to errors when the surfaces were almost but not quite aligned
            // with some coordinate direction.
            sum = 0.0;
            for (i=0; i<3; i++) {
                sum = sum + pow((x1Max[i]-x1Min[i]), 2.0);
            }
            s1 = sum;
            
            sum = 0.0;
            for (i=0; i<3; i++) {
                sum = sum + pow((x2rMax[i]-x2rMin[i]), 2.0);
            }
            s2 = sum;
            if (s2 > DBL_EPSILON) {
                for (i=0; i<3; i++) {
                    scl[i] = sqrt(s1/s2);
                }
            } else {
                for (i=0; i<3; i++) {
                    scl[i] = 1.0;
                }
            }
        }
        
        NSLog(@"periodicProjectorInModel: scaling with: %f %f %f\n", scl[0], scl[1], scl[2]);
        for (i=0; i<3; i++) {
            sclMatrix[i][i] = scl[i];
        }
    }
    
    // Transform the target boundary on top of this boundary
    // Translations:
    if (pArray.matrix != NULL) {
        free_dmatrix(pArray.matrix, 0, pArray.m-1, 0, pArray.n-1);
        pArray.matrix = NULL;
    }
    boundaryConditionAtId = (model.boundaryConditions)[this];
    found = [listUtil listGetConstRealArray:model inArray:boundaryConditionAtId.valuesList forVariable:@"periodic bc translate" buffer:&pArray];
    // Define translations so that the lower left corner is the same
    if (found == NO) {
        for (i=0; i<3; i++) {
            trsMatrix[3][i] = x1Min[i] - sclMatrix[i][i] * x2rMin[i];
        }
        NSLog(@"periodicProjectorInModel: translation: %f %f %f\n", trsMatrix[3][0], trsMatrix[3][1], trsMatrix[3][2]);
    } else {
        for (i=0; i<pArray.m; i++) {
            trsMatrix[3][i] = pArray.matrix[i][0];
        }
    }
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, 4, 4, 1.0, (double *)sclMatrix, 4, (double *)trsMatrix, 4, 0.0, (double *)trfMatrix, 4);
    
    // If whole transf. matrix given, it will override all the other settings
    if (pArray.matrix != NULL) {
        free_dmatrix(pArray.matrix, 0, pArray.m-1, 0, pArray.n-1);
        pArray.matrix = NULL;
    }
    boundaryConditionAtId = (model.boundaryConditions)[this];
    found = [listUtil listGetConstRealArray:model inArray:boundaryConditionAtId.valuesList forVariable:@"periodic bc matrix" buffer:&pArray];
    if (found == YES) {
        for (i=0; i<pArray.m; i++) {
            for (j=0; j<pArray.n; j++) {
                trfMatrix[i][j] = pArray.matrix[j][i];
                
            }
        }
    }
    
    if (pArray.matrix != NULL) {
        free_dmatrix(pArray.matrix, 0, pArray.m-1, 0, pArray.n-1);
        pArray.matrix = NULL;
    }
    
    // Re-number the element node pointers to use only boundary nodes
    for (i=0; i<n1; i++) {
        for (j=0; j<mesh1elements[i].sizeNodeIndexes; j++) {
            mesh1elements[i].NodeIndexes[j] = perm1[mesh1elements[i].NodeIndexes[j]];
        }
    }
    
    for (i=0; i<n2; i++) {
        for (j=0; j<mesh2elements[i].sizeNodeIndexes; j++) {
            mesh2elements[i].NodeIndexes[j] = perm2[mesh2elements[i].NodeIndexes[j]];
        }
    }
    
    // Now transform the coordinates
    for (i=0; i<bmesh2.numberOfNodes; i++) {
        x[0] = mesh2nodes->x[i];
        x[1] = mesh2nodes->y[i];
        x[2] = mesh2nodes->z[i];
        x[3] = 1.0;
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 4, 4, 1.0, (double *)trfMatrix, 4, x, 1, 0.0, x, 1);
        mesh2nodes->x[i] = x[0] / x[3];
        mesh2nodes->y[i] = x[1] / x[3];
        mesh2nodes->z[i] = x[2] / x[3];
    }
    
    // Get the mesh projector which now contains weights between the boundary nodes
    useQuadrantTree = [listUtil listGetLogical:model inArray:model.simulation.valuesList forVariable:@"use quadrant tree" info:&found];
    if (found == NO) useQuadrantTree = YES;
    
    projector = [utilities meshProjector:bmesh2 secondmesh:bmesh1 model:model useQuadrantTree:&useQuadrantTree transpose:NULL];
    projectorContainers = projector.getContainers;
    projectorContainers->InvPerm = invPerm1;
    for (i=0; i<projector.numberOfRows; i++) {
        for (j=projectorContainers->Rows[i]; j<=projectorContainers->Rows[i+1]-1; j++) {
            k = projectorContainers->Cols[j];
            if (k >= 0) projectorContainers->Cols[j] = invPerm2[k];
        }
    }
    
    // Deallocate mesh structures
    free_dvector(elementNodes->x, 0, mesh.maxElementNodes-1);
    free_dvector(elementNodes->y, 0, mesh.maxElementNodes-1);
    free_dvector(elementNodes->z, 0, mesh.maxElementNodes-1);
    free(elementNodes);
    
    for (i=0; i<bmesh1.numberOfBulkElements; i++) {
        free_ivector(mesh1elements[i].NodeIndexes, 0, mesh1elements[i].sizeNodeIndexes-1);
    }
    [bmesh1 deallocateQuadrantTree];
    
    for (i=0; i<bmesh2.numberOfBulkElements; i++) {
        free_ivector(mesh2elements[i].NodeIndexes, 0, mesh2elements[i].sizeNodeIndexes-1);
    }
    [bmesh2 deallocateQuadrantTree];
    
    free_dvector(mesh1nodes->x, 0, bmesh1.numberOfNodes-1);
    free_dvector(mesh1nodes->y, 0, bmesh1.numberOfNodes-1);
    free_dvector(mesh1nodes->z, 0, bmesh1.numberOfNodes-1);
    free(mesh1nodes);
    
    free_dvector(mesh2nodes->x, 0, bmesh2.numberOfNodes-1);
    free_dvector(mesh2nodes->y, 0, bmesh2.numberOfNodes-1);
    free_dvector(mesh2nodes->z, 0, bmesh2.numberOfNodes-1);
    free(mesh2nodes);
    
    if (perm1 != NULL) {
        free_ivector(perm1, 0, mesh.numberOfNodes-1);
    }
    if (perm2 != NULL) {
        free_ivector(perm2, 0, mesh.numberOfNodes-1);
    }
    invPerm1 = NULL;
    free_ivector(invPerm2, 0, bmesh2.numberOfNodes-1);
    bmesh1 = nil;
    bmesh2 = nil;
    
    [elmDescription deallocation];
    return projector;
}

/********************************************************************************************
 
    Split a mesh equally to smaller pieces by performing a uniform split. Also know as 
    mesh mulyiplication. A 2D element splits into 4 elemenrs of same form and a #d element 
    into 8 elements.
    Currently works only for linear elements
 
********************************************************************************************/
-(FEMMesh *)splitMeshEqual:(FEMMesh *)mesh model:(FEMModel *)model nodal:(double *)h sizeNodal:(int *)sizeNodal {
    
    int i, j, k, l, n, newElCnt, nodeCnt, edgeCnt, node, parentID, diag, minVal, maxVal;
    int faceNumber, edge1, edge2, edge3, edge4, node12, node23, node34, node41, node31;
    int n1, n2, n3, *eoldNodes, *faceNodes, *edgeNodes;        // Only linear so far
    int **child;
    double *nodals, *xh = NULL;
    double dxyz[3][3], dist[3], r, s, t, h1, h2, sum;
    BOOL found;
    Element_t *elements, *newElements, *faces, *edges, *eptr, *eParent;
    Nodes_t *nodes, *newNodes;
    FEMMesh *newMesh;
    FEMBoundaryCondition *boundaryConditionAtId;
    FEMElementDescription *elementDescription;
    FEMListUtilities *listUtil;
    
    if (mesh == nil) return nil;
    
    newMesh = [[FEMMesh alloc] init];
    
    elementDescription = [[FEMElementDescription alloc] init];
    
    [self findEdgesForMesh:mesh findEdges:NULL];
    
    NSLog(@"splitMeshEqual: ********** Old mesh **********\n");
    NSLog(@"Nodes: %d\n", mesh.numberOfNodes);
    NSLog(@"splitMeshEqual: bulk elements: %d\n", mesh.numberOfBulkElements);
    NSLog(@"splitMeshEqual: boundary elements: %d\n", mesh.numberOfBoundaryElements);
    NSLog(@"splitMeshEqual: edges: %d\n", mesh.numberOfEdges);
    NSLog(@"splitMeshEqual: faces: %d\n", mesh.numberOfFaces);
    
    // Update nodal coordinates
    nodeCnt = mesh.numberOfNodes + mesh.numberOfEdges;
    
    elements = mesh.getElements;
    faces = mesh.getFaces;
    edges = mesh.getEdges;
    nodes = mesh.getNodes;
    
    // For bricks, count faces
    for (i=0; i<mesh.numberOfFaces; i++) {
        if (faces[i].Type.NumberOfNodes == 4) nodeCnt++;
    }
    
    // For quads and bricks, count centerpoints
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        switch (elements[i].Type.ElementCode / 100) {
            case 4:
            case 8:
                nodeCnt++;
                break;
        }
    }
    
    // New mesh node coordinate arrays
    newNodes = newMesh.getNodes;
    newNodes->x = doublevec(0, nodeCnt-1);
    newNodes->y = doublevec(0, nodeCnt-1);
    newNodes->z = doublevec(0, nodeCnt-1);
    
    // New mesh includes old mesh nodes
    for (i=0; i<mesh.numberOfNodes; i++) {
        newNodes->x[i] = nodes->x[i];
        newNodes->y[i] = nodes->y[i];
        newNodes->z[i] = nodes->z[i];
    }
    
    if (h != NULL) {
        xh = doublevec(0, nodeCnt-1);
        for (i=0; i<*sizeNodal; i++) {
            xh[i] = h[i];
        }
    }
    
    // Add edge centers
    for (i=0; i<mesh.numberOfEdges; i++) {
        k = edges[i].Type.NumberOfNodes;
        j = i + mesh.numberOfNodes;
        if (h != NULL) {
            h1 = h[edges[i].NodeIndexes[0]];
            h2 = h[edges[i].NodeIndexes[1]];
            r = 1.0 / (1.0 + h1/h2);
            newNodes->x[j] = r*nodes->x[edges[i].NodeIndexes[0]] + (1.0-r)*nodes->x[edges[i].NodeIndexes[1]];
            newNodes->y[j] = r*nodes->y[edges[i].NodeIndexes[0]] + (1.0-r)*nodes->y[edges[i].NodeIndexes[1]];
            newNodes->z[j] = r*nodes->z[edges[i].NodeIndexes[0]] + (1.0-r)*nodes->z[edges[i].NodeIndexes[1]];
            xh[j] = r*h1+(1.0-r)*h2;
        } else {
            sum = 0.0;
            for (l=0; l<edges[i].Type.NumberOfNodes; l++) {
                sum = sum + nodes->x[edges[i].NodeIndexes[l]];
            }
            newNodes->x[j] = sum / k;
            
            sum = 0.0;
            for (l=0; l<edges[i].Type.NumberOfNodes; l++) {
                sum = sum + nodes->y[edges[i].NodeIndexes[l]];
            }
            newNodes->y[j] = sum / k;
            
            sum = 0.0;
            for (l=0; l<edges[i].Type.NumberOfNodes; l++) {
                sum = sum + nodes->z[edges[i].NodeIndexes[l]];
            }
            newNodes->z[j] = sum / k;
        }
    }
    
    nodals = doublevec(0, mesh.maxElementNodes-1);
    memset( nodals, 0.0, (mesh.maxElementNodes*sizeof(nodals)) );
    
    // Add face centers for bricks
    for (i=0; i<mesh.numberOfFaces; i++) {
        k = faces[i].Type.NumberOfNodes;
        if (k == 4) {
            j = i + mesh.numberOfNodes + mesh.numberOfEdges;
            if (h != NULL) {
                n = mesh.numberOfNodes;
                h1 = xh[n+faces[i].EdgeIndexes[1]];
                h2 = xh[n+faces[i].EdgeIndexes[3]];
                r = 2.0/(1.0 + h1/h2) - 1.0 ;
                h1 = xh[n+faces[i].EdgeIndexes[2]];
                h2 = xh[n+faces[i].EdgeIndexes[0]];
                s = 2.0/(1.0 + h1/h2) - 1.0;
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    nodals[l] = nodes->x[faces[i].NodeIndexes[l]];
                }
                newNodes->x[j] = [elementDescription interpolate2DInElement:&faces[i] nodalValues:nodals evaluatedAt:r andAt:s];
                
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    nodals[l] = nodes->y[faces[i].NodeIndexes[l]];
                }
                newNodes->y[j] = [elementDescription interpolate2DInElement:&faces[i] nodalValues:nodals evaluatedAt:r andAt:s];
                
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    nodals[l] = nodes->z[faces[i].NodeIndexes[l]];
                }
                newNodes->z[j] = [elementDescription interpolate2DInElement:&faces[i] nodalValues:nodals evaluatedAt:r andAt:s];
                
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    nodals[l] = h[faces[i].NodeIndexes[l]];
                }
                xh[j] = [elementDescription interpolate2DInElement:&faces[i] nodalValues:nodals evaluatedAt:r andAt:s];
            } else {
                sum = 0.0;
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    sum = sum + nodes->x[faces[i].NodeIndexes[l]];
                }
                newNodes->x[j] = sum / k;
                
                sum = 0.0;
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    sum = sum + nodes->y[faces[i].NodeIndexes[l]];
                }
                newNodes->y[j] = sum / k;
                
                sum = 0.0;
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    sum = sum + nodes->z[faces[i].NodeIndexes[l]];
                }
                newNodes->z[j] = sum / k;
            }
        }
    }
    
    memset( nodals, 0.0, (mesh.maxElementNodes*sizeof(nodals)) );
    
    // Add centerpoint for quads & bricks
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        k = elements[i].Type.NumberOfNodes;
        switch (elements[i].Type.ElementCode / 100) {
            case 4:
                j++;
                if (h != NULL) {
                    n = mesh.numberOfNodes;
                    h1 = xh[n+elements[i].EdgeIndexes[1]];
                    h2 = xh[n+elements[i].EdgeIndexes[3]];
                    r = 2.0/(1.0 + h1/h2) - 1.0;
                    h1 = xh[n+elements[i].EdgeIndexes[2]];
                    h2 = xh[n+elements[i].EdgeIndexes[0]];
                    s = 2.0/(1.0 + h1/h2) - 1.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->x[elements[i].NodeIndexes[l]];
                    }
                    newNodes->x[j] = [elementDescription interpolate2DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s];
                    
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->y[elements[i].NodeIndexes[l]];
                    }
                    newNodes->y[j] = [elementDescription interpolate2DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s];
                    
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->z[elements[i].NodeIndexes[l]];
                    }
                    newNodes->z[j] = [elementDescription interpolate2DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s];
                } else {
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->x[elements[i].NodeIndexes[l]];
                    }
                    newNodes->x[j] = sum / k;
                    
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->y[elements[i].NodeIndexes[l]];
                    }
                    newNodes->y[j] = sum / k;
                    
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->z[elements[i].NodeIndexes[l]];
                    }
                    newNodes->z[j] = sum / k;
                }
                break;
                
            case 8:
                j++;
                if (h != NULL) {
                    n = mesh.numberOfNodes + mesh.numberOfEdges;
                    h1 = xh[n+elements[i].FaceIndexes[3]];
                    h2 = xh[n+elements[i].FaceIndexes[5]];
                    r = 2.0/(1.0 + h1/h2) - 1.0;
                    
                    h1 = xh[n+elements[i].FaceIndexes[4]];
                    h2 = xh[n+elements[i].FaceIndexes[2]];
                    s = 2.0/(1.0 + h1/h2) - 1.0;
                    
                    h1 = xh[n+elements[i].FaceIndexes[1]];
                    h2 = xh[n+elements[i].FaceIndexes[0]];
                    t = 2.0/(1.0 + h1/h2) - 1.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->x[elements[i].NodeIndexes[l]];
                    }
                    newNodes->x[j] = [elementDescription interpolate3DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s andAt:t];
                    
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->y[elements[i].NodeIndexes[l]];
                    }
                    newNodes->y[j] = [elementDescription interpolate3DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s andAt:t];
                    
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->z[elements[i].NodeIndexes[l]];
                    }
                    newNodes->z[j] = [elementDescription interpolate3DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s andAt:t];
                } else {
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->x[elements[i].NodeIndexes[l]];
                    }
                    newNodes->x[j] = sum / k;
                    
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->y[elements[i].NodeIndexes[l]];
                    }
                    newNodes->y[j] = sum / k;
                    
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->z[elements[i].NodeIndexes[l]];
                    }
                    newNodes->z[j] = sum / k;
                }
                break;
        }
    }
    
    if (xh != NULL) {
        free_dvector(xh, 0, nodeCnt-1);
    }
    free_dvector(nodals, 0, mesh.maxElementNodes-1);
    [elementDescription deallocation];
    
    // Update new mesh node count
    newMesh.numberOfEdges = 0;
    newMesh.numberOfFaces = 0;
    newMesh.maxBdofs = mesh.maxBdofs;
    newMesh.maxEdgeDofs = mesh.maxEdgeDofs;
    newMesh.maxFaceDofs = mesh.maxFaceDofs;
    newMesh.maxElementDofs = mesh.maxElementDofs;
    newMesh.dimension = mesh.dimension;
    
    newMesh.numberOfNodes = nodeCnt;
    newNodes->numberOfNodes = nodeCnt;
    
    // Update bulk elements
    
    // First count new elements
    newElCnt = 0;
    for (i=0; i<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; i++) {
        
        switch (elements[i].Type.ElementCode / 100) {
            // Each element will be divided into 2^dim new elements
            case 2:
                newElCnt = newElCnt + 2;  // Lines
                break;
            case 3:
                newElCnt = newElCnt + 4;  // Trias
                break;
            case 4:
                newElCnt = newElCnt + 4;  // Quads
                break;
            case 5:
                newElCnt = newElCnt + 8;  // Tetras
                break;
            case 8:
                newElCnt = newElCnt + 8;  // Hexas
                break;
        }
    }
    newElements = newMesh.getElements;
    newElements = (Element_t*) malloc( sizeof(Element_t) * newElCnt);
    initElements(newElements, newElCnt);
    for (i=0; i<newElCnt; i++) {
        newElements[i].EdgeIndexes = NULL;
        newElements[i].FaceIndexes = NULL;
    }
    
    child = intmatrix(0, mesh.numberOfBulkElements-1, 0, 7);
    newElCnt = 0;
    nodeCnt = mesh.numberOfNodes;
    edgeCnt = mesh.numberOfEdges;
    
    // Index to old quad/hexa centerpoint node in the new nodal arrays
    node = nodeCnt + mesh.numberOfEdges + mesh.numberOfFaces;

    // Now update all new mesh elements
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        
        switch (elements[i].Type.ElementCode) {
            case 303:
                // Split triangle to four triangles from
                // edge centerpoints
                
                // 1st new element
                child[i][0] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                
                // 2nd new element
                child[i][1] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElCnt++;
                
                // 3rd new element
                child[i][2] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                
                // 4th new element
                child[i][3] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                break;
                
            case 404:
                // Index to old quad centerpoint node in the
                // new mesh nodal arrays
                // node = node not node = node + 1 as it's the case in Elmer
                
                // Spit quad to four new quads from edge
                // centerpoints and centerpoint of the element
                
                // 1st new element
                child[i][0] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = node;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElCnt++;
                
                // 2nd new element
                child[i][1] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = node;
                newElCnt++;
                
                // 3rd new element
                child[i][2] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = node;
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                
                // 4th new element
                child[i][3] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = node;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].NodeIndexes[3];
                newElCnt++;
                break;
                
            case 504:
                // Split tetra to 8 new elements from
                // corners and edge centerpoints
                
                // 1st new element
                child[i][0] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElCnt++;
                
                // 2nd new element
                child[i][1] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[4] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElCnt++;
                
                // 3rd new element
                child[i][2] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[5] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                
                // 4th new element
                child[i][3] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[3];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[5] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[4] + nodeCnt;
                newElCnt++;
                
                // Then the annoying part; we will have to split the remaining octahedron into four elements. This can be done in three ways
                // of which only one preserves the minimum angle condition (Delaunay splitting)
                dxyz[0][0] = newNodes->x[elements[i].EdgeIndexes[3] + nodeCnt] - newNodes->x[elements[i].EdgeIndexes[1] + nodeCnt];
                dxyz[1][0] = newNodes->y[elements[i].EdgeIndexes[3] + nodeCnt] - newNodes->y[elements[i].EdgeIndexes[1] + nodeCnt];
                dxyz[2][0] = newNodes->z[elements[i].EdgeIndexes[3] + nodeCnt] - newNodes->z[elements[i].EdgeIndexes[1] + nodeCnt];
                
                dxyz[0][1] = newNodes->x[elements[i].EdgeIndexes[4] + nodeCnt] - newNodes->x[elements[i].EdgeIndexes[2] + nodeCnt];
                dxyz[1][1] = newNodes->y[elements[i].EdgeIndexes[4] + nodeCnt] - newNodes->y[elements[i].EdgeIndexes[2] + nodeCnt];
                dxyz[2][1] = newNodes->z[elements[i].EdgeIndexes[4] + nodeCnt] - newNodes->z[elements[i].EdgeIndexes[2] + nodeCnt];
                
                dxyz[0][2] = newNodes->x[elements[i].EdgeIndexes[5] + nodeCnt] - newNodes->x[elements[i].EdgeIndexes[0] + nodeCnt];
                dxyz[1][2] = newNodes->y[elements[i].EdgeIndexes[5] + nodeCnt] - newNodes->y[elements[i].EdgeIndexes[0] + nodeCnt];
                dxyz[2][2] = newNodes->z[elements[i].EdgeIndexes[5] + nodeCnt] - newNodes->z[elements[i].EdgeIndexes[0] + nodeCnt];
                
                dist[0] = sqrt( pow(dxyz[0][0], 2.0) + pow(dxyz[1][0], 2.0) + pow(dxyz[2][0], 2.0) );
                dist[1] = sqrt( pow(dxyz[0][1], 2.0) + pow(dxyz[1][1], 2.0) + pow(dxyz[2][1], 2.0) );
                dist[2] = sqrt( pow(dxyz[0][2], 2.0) + pow(dxyz[1][2], 2.0) + pow(dxyz[2][2], 2.0) );
                
                diag = 1; // The default diaginal for splitting is between edges 2-4
                if (dist[1] < dist[0] && dist[1] < dist[2]) diag = 2; // Edges 3-5
                if (dist[2] < dist[0] && dist[2] < dist[1]) diag = 3; // Edges 1-6
                
                switch (diag) {
                    case 1:
                        // 5th new element
                        child[i][4] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElCnt++;
                        
                        // 6th new element
                        child[i][5] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElCnt++;
                        
                        // 7th new element
                        child[i][6] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElCnt++;
                        
                        // 8th new element
                        child[i][7] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElCnt++;
                        break;
                        
                    case 2:
                        // 5th new element
                        child[i][4] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElCnt++;
                        
                        // 6th new element
                        child[i][5] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElCnt++;
                        
                        // 7th new element
                        child[i][6] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElCnt++;
                        
                        // 8th new element
                        child[i][7] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElCnt++;
                        break;
                        
                    case 3:
                        // 5th new element
                        child[i][4] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElCnt++;
                        
                        // 6th new element
                        child[i][5] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElCnt++;
                        
                        // 7th new element
                        child[i][6] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElCnt++;
                        
                        // 8th new element
                        child[i][7] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElCnt++;
                        break;
                }
                
            case 808:
                // Index to old quad centerpoints node in the new mesh nodal arrays
                // node = node not node = node + 1 as it's the case in Elmer
                
                // Split brick to 8 new bricks from edge centerpoints
                // and centerpoint of the element
                
                // 1st new element
                child[i][0] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].FaceIndexes[0] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].EdgeIndexes[8] + nodeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].FaceIndexes[2] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[6] = node;
                newElements[newElCnt].NodeIndexes[7] = elements[i].FaceIndexes[5] + nodeCnt + edgeCnt;
                newElCnt++;
                
                // 2nd new element
                child[i][1] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].FaceIndexes[0] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].FaceIndexes[2] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].EdgeIndexes[9] + nodeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].FaceIndexes[3] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[7] = node;
                newElCnt++;
                
                // 3rd new element
                child[i][2] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].FaceIndexes[0] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].NodeIndexes[3];
                newElements[newElCnt].NodeIndexes[4] = elements[i].FaceIndexes[5] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[5] = node;
                newElements[newElCnt].NodeIndexes[6] = elements[i].FaceIndexes[4] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].EdgeIndexes[11] + nodeCnt;
                newElCnt++;
                
                // 4th new element
                child[i][3] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].FaceIndexes[0] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElements[newElCnt].NodeIndexes[4] = node;
                newElements[newElCnt].NodeIndexes[5] = elements[i].FaceIndexes[3] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].EdgeIndexes[10] + nodeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].FaceIndexes[4] + nodeCnt + edgeCnt;
                newElCnt++;
                
                // 5th new element
                child[i][4] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[8] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].FaceIndexes[2] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[2] = node;
                newElements[newElCnt].NodeIndexes[3] = elements[i].FaceIndexes[5] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].NodeIndexes[4];
                newElements[newElCnt].NodeIndexes[5] = elements[i].EdgeIndexes[4] + nodeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].FaceIndexes[1] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].EdgeIndexes[7] + nodeCnt;
                newElCnt++;
                
                // 6th new element
                child[i][5] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].FaceIndexes[2] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[9] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].FaceIndexes[3] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[3] = node;
                newElements[newElCnt].NodeIndexes[4] = elements[i].EdgeIndexes[4] + nodeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].NodeIndexes[5];
                newElements[newElCnt].NodeIndexes[6] = elements[i].EdgeIndexes[5] + nodeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].FaceIndexes[1] + nodeCnt + edgeCnt;
                newElCnt++;
                
                // 7th new element
                child[i][6] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].FaceIndexes[5] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[1] = node;
                newElements[newElCnt].NodeIndexes[2] = elements[i].FaceIndexes[4] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[11] + nodeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].EdgeIndexes[7] + nodeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].FaceIndexes[1] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].EdgeIndexes[6] + nodeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].NodeIndexes[7];
                newElCnt++;
                
                // 8th new element
                child[i][7] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = node;
                newElements[newElCnt].NodeIndexes[1] = elements[i].FaceIndexes[3] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[10] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].FaceIndexes[4] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].FaceIndexes[1] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].EdgeIndexes[5] + nodeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].NodeIndexes[6];
                newElements[newElCnt].NodeIndexes[7] = elements[i].EdgeIndexes[6] + nodeCnt;
                newElCnt++;
                break;
                
            default:
                printf("Element type %d not supported by the multigrid solver.\n", elements[i].Type.ElementCode);
                errorfunct("splitMeshEqual", "Program terminating now...");
        }
    }
    
    
    // Update new mesh element counts
    newMesh.numberOfBulkElements = newElCnt;
    
    // Update boundary elements
    // Note:  internal boundaries not taken care of....!!!
    eoldNodes = intvec(0, 3);
    faceNodes = intvec(0, 3);
    edgeNodes = intvec(0, 1);
    for (i=0; i<mesh.numberOfBoundaryElements; i++) {
                
        // Get parent of the boundary element
        eParent = elements[i+mesh.numberOfBulkElements].BoundaryInfo->Left;
        if (eParent == NULL) eParent = elements[i+mesh.numberOfBulkElements].BoundaryInfo->Right;
        if (eParent == NULL) continue;
        
        parentID = eParent->ElementIndex-1;
        
        switch (elements[i+mesh.numberOfBulkElements].Type.ElementCode / 100) {
            case 2:
                // Line segments
                
                // Which edge of the parent element are we?
                for (edge1=0; edge1<eParent->sizeEdgeIndexes; edge1++) {
                    if ((elements[i+mesh.numberOfBulkElements].NodeIndexes[0] == edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0]
                         && elements[i+mesh.numberOfBulkElements].NodeIndexes[1] == edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]) ||
                        (elements[i+mesh.numberOfBulkElements].NodeIndexes[1] == edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0] &&
                         elements[i+mesh.numberOfBulkElements].NodeIndexes[0] == edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1])) break;
                }
                
                // Index of the old edge centerpoint in the
                // new mesh nodal arrays
                node = eParent->EdgeIndexes[edge1] + mesh.numberOfNodes;
                
                // 1st new element
                newElements[newElCnt] = elements[i+mesh.numberOfBulkElements];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 1);
                newElements[newElCnt].sizeNodeIndexes = 2;
                newElements[newElCnt].NodeIndexes[0] = elements[i+mesh.numberOfBulkElements].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = node;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfBulkElements].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the
                // children of the old mesh parent element
                for (j=0; j<4; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    found = NO;
                    for (k=0; k<n-1; k++) {
                        if ((newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[k] &&
                            newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[k+1]) ||
                            (newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[k] &&
                             newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[k+1])) {
                                found = YES;
                                break;
                        }
                    }
                    if ((newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[n-1] &&
                         newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[0]) ||
                        (newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[n-1] &&
                         newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[0])) {
                            found = YES;
                    }
                    if (found == YES) break;
                }
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 2nd new element
                newElements[newElCnt] = elements[i+mesh.numberOfBulkElements];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 1);
                newElements[newElCnt].sizeNodeIndexes = 2;
                newElements[newElCnt].NodeIndexes[0] = node;
                newElements[newElCnt].NodeIndexes[1] = elements[i+mesh.numberOfBulkElements].NodeIndexes[1];
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfBulkElements].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the
                // children of the old mesh parent element
                for (j=0; j<4; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    found = NO;
                    for (k=0; k<n-1; k++) {
                        if ((newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[k] &&
                             newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[k+1]) ||
                            (newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[k] &&
                             newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[k+1])) {
                                found = YES;
                                break;
                        }
                    }
                    if ((newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[n-1] &&
                         newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[0]) ||
                        (newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[n-1] &&
                         newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[0])) {
                            found = YES;
                    }
                    if (found == YES) break;
                }
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                break;
                
            case 3:
                // Trias
                
                // On which face of the parent element are we?
                for (j=0; j<3; j++) {
                    eoldNodes[j] = elements[i+mesh.numberOfBulkElements].NodeIndexes[j];
                }
                sort(3, eoldNodes-1);
                
                for (faceNumber=0; faceNumber<eParent->sizeFaceIndexes; faceNumber++) {
                    for (j=0; j<3; j++) {
                        faceNodes[j] = faces[eParent->FaceIndexes[faceNumber]].NodeIndexes[j];
                    }
                    sort(3, faceNodes-1);
                    
                    if (eoldNodes[0] == faceNodes[0] && eoldNodes[1] == faceNodes[1] && eoldNodes[2] == faceNodes[2]) break;
                }
                
                // Then what are the edges on this face?
                
                // First edge
                eoldNodes[0] = min(elements[i+mesh.numberOfBulkElements].NodeIndexes[0], elements[i+mesh.numberOfBulkElements].NodeIndexes[1]);
                eoldNodes[1] = max(elements[i+mesh.numberOfBulkElements].NodeIndexes[0], elements[i+mesh.numberOfBulkElements].NodeIndexes[1]);
                for (edge1=0; edge1<eParent->sizeEdgeIndexes; edge1++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Second edge
                eoldNodes[0] = min(elements[i+mesh.numberOfBulkElements].NodeIndexes[1], elements[i+mesh.numberOfBulkElements].NodeIndexes[2]);
                eoldNodes[1] = max(elements[i+mesh.numberOfBulkElements].NodeIndexes[1], elements[i+mesh.numberOfBulkElements].NodeIndexes[2]);
                for (edge2=0; edge2<eParent->sizeEdgeIndexes; edge2++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge2]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge2]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge2]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge2]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Third edge
                eoldNodes[0] = min(elements[i+mesh.numberOfBulkElements].NodeIndexes[2], elements[i+mesh.numberOfBulkElements].NodeIndexes[0]);
                eoldNodes[1] = max(elements[i+mesh.numberOfBulkElements].NodeIndexes[2], elements[i+mesh.numberOfBulkElements].NodeIndexes[0]);
                for (edge3=0; edge3<eParent->sizeEdgeIndexes; edge3++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge3]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge3]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge3]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge3]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Index of the old face and edge centerpoints
                // in the new mesh nodal arrays
                node12 = eParent->EdgeIndexes[edge1] + mesh.numberOfNodes;
                node23 = eParent->EdgeIndexes[edge2] + mesh.numberOfNodes;
                node31 = eParent->EdgeIndexes[edge3] + mesh.numberOfNodes;
                
                // 1st new element
                newElements[newElCnt] = elements[i+mesh.numberOfBulkElements];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i+mesh.numberOfBulkElements].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = node12;
                newElements[newElCnt].NodeIndexes[2] = node31;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfBulkElements].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0; // Count matches
                    for (n1=0; n1<3; n1++) {
                        for (n2=0; n2<4; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) errorfunct("splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 2nd new element
                newElements[newElCnt] = elements[i+mesh.numberOfBulkElements];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = node12;
                newElements[newElCnt].NodeIndexes[1] = elements[i+mesh.numberOfBulkElements].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = node23;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfBulkElements].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<3; n1++) {
                        for (n2=0; n2<4; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) errorfunct("splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 3rd new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = node12;
                newElements[newElCnt].NodeIndexes[1] = node23;
                newElements[newElCnt].NodeIndexes[2] = node31;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the children
                // of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<3; n1++) {
                        for (n2=0; n2<4; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) errorfunct("splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 4th new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = node31;
                newElements[newElCnt].NodeIndexes[1] = node23;
                newElements[newElCnt].NodeIndexes[2] = elements[i+mesh.numberOfNodes].NodeIndexes[2];
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among
                // the children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<3; n1++) {
                        for (n2=0; n2<4; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) errorfunct("splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                break;
                
            case 4:
                // Quads
                
                // On which face of the parent element are we?
                for (j=0; j<4; j++) {
                    eoldNodes[j] = elements[i+mesh.numberOfNodes].NodeIndexes[j];
                }
                sort(4, eoldNodes-1);
                
                for (faceNumber=0; faceNumber<eParent->sizeFaceIndexes; faceNumber++) {
                    for (j=0; j<4; j++) {
                        faceNodes[j] = faces[eParent->FaceIndexes[faceNumber]].NodeIndexes[j];
                    }
                    sort(4, faceNodes-1);
                    
                    if (eoldNodes[0] == faceNodes[0] && eoldNodes[1] == faceNodes[1] && eoldNodes[2] == faceNodes[2]) break;
                }
                
                // Then what are the edges on this face?
                
                // First edge
                eoldNodes[0] = min(elements[i+mesh.numberOfNodes].NodeIndexes[0], elements[i+mesh.numberOfNodes].NodeIndexes[1]);
                eoldNodes[1] = max(elements[i+mesh.numberOfNodes].NodeIndexes[0], elements[i+mesh.numberOfNodes].NodeIndexes[1]);
                for (edge1=0; edge1<eParent->sizeEdgeIndexes; edge1++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Second edge
                eoldNodes[0] = min(elements[i+mesh.numberOfNodes].NodeIndexes[1], elements[i+mesh.numberOfNodes].NodeIndexes[2]);
                eoldNodes[1] = max(elements[i+mesh.numberOfNodes].NodeIndexes[1], elements[i+mesh.numberOfNodes].NodeIndexes[2]);
                for (edge2=0; edge2<eParent->sizeEdgeIndexes; edge2++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge2]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge2]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge2]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge2]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Third edge
                eoldNodes[0] = min(elements[i+mesh.numberOfNodes].NodeIndexes[2], elements[i+mesh.numberOfNodes].NodeIndexes[3]);
                eoldNodes[1] = max(elements[i+mesh.numberOfNodes].NodeIndexes[2], elements[i+mesh.numberOfNodes].NodeIndexes[3]);
                for (edge3=0; edge3<eParent->sizeEdgeIndexes; edge3++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge3]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge3]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge3]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge3]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Fourth edge
                eoldNodes[0] = min(elements[i+mesh.numberOfNodes].NodeIndexes[3], elements[i+mesh.numberOfNodes].NodeIndexes[0]);
                eoldNodes[1] = max(elements[i+mesh.numberOfNodes].NodeIndexes[3], elements[i+mesh.numberOfNodes].NodeIndexes[0]);
                for (edge4=0; edge4<eParent->sizeEdgeIndexes; edge4++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge4]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge4]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge4]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge4]].NodeIndexes[1]);
                    
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Index of the old face and edge centerpoints
                // in the new mesh nodal arrays
                node = eParent->FaceIndexes[faceNumber] + // Face mid-point
                        mesh.numberOfNodes + mesh.numberOfEdges;
                node12 = eParent->EdgeIndexes[edge1] + mesh.numberOfNodes;
                node23 = eParent->EdgeIndexes[edge2] + mesh.numberOfNodes;
                node34 = eParent->EdgeIndexes[edge3] + mesh.numberOfNodes;
                node41 = eParent->EdgeIndexes[edge4] + mesh.numberOfNodes;
                
                // 1st new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i+mesh.numberOfNodes].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = node12;
                newElements[newElCnt].NodeIndexes[2] = node;
                newElements[newElCnt].NodeIndexes[3] = node41;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mech parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<4; n1++) {
                        for (n2=0; n2<8; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) errorfunct("splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 2nd new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = node12;
                newElements[newElCnt].NodeIndexes[1] = elements[i+mesh.numberOfNodes].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = node23;
                newElements[newElCnt].NodeIndexes[3] = node;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mech parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<4; n1++) {
                        for (n2=0; n2<8; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) errorfunct("splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 3rd new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = node41;
                newElements[newElCnt].NodeIndexes[1] = node;
                newElements[newElCnt].NodeIndexes[2] = node34;
                newElements[newElCnt].NodeIndexes[3] = elements[i+mesh.numberOfNodes].NodeIndexes[3];
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mech parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<4; n1++) {
                        for (n2=0; n2<8; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2)  break;
                }
                if (n3 < 3) errorfunct("splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 4th new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = node;
                newElements[newElCnt].NodeIndexes[1] = node23;
                newElements[newElCnt].NodeIndexes[2] = elements[i+mesh.numberOfNodes].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[3] = node34;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mech parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0; // Count matches
                    for (n1=0; n1<4; n1++) {
                        for (n2=0; n2<8; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) errorfunct("splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                break;
        }
        
    }
    free_ivector(eoldNodes, 0, 3);
    free_ivector(faceNodes, 0, 3);
    free_ivector(edgeNodes, 0, 1);
    
    // Update new mesh boundary element counts
    newMesh.numberOfBoundaryElements = newElCnt - newMesh.numberOfBulkElements;
    newMesh.maxElementDofs = mesh.maxElementDofs;
    newMesh.maxElementNodes = mesh.maxElementNodes;
    
    for (i=0; i<newMesh.numberOfBulkElements+newMesh.numberOfBoundaryElements; i++) {
        newElements[i].EdgeIndexes = NULL;
        newElements[i].FaceIndexes = NULL;
    }
    
    NSLog(@"splitMeshEqual: ********** New mesh **********\n");
    NSLog(@"Nodes: %d\n", newMesh.numberOfNodes);
    NSLog(@"splitMeshEqual: bulk elements: %d\n", newMesh.numberOfBulkElements);
    NSLog(@"splitMeshEqual: boundary elements: %d\n", newMesh.numberOfBoundaryElements);
    
    // TODO: add support for parallel run
    
    // If periodic BC given, compute boundary mesh projector
    listUtil = [[FEMListUtilities alloc] init];
    for (i=0; i<model.numberOfBoundaryConditions; i++) {
        boundaryConditionAtId = (model.boundaryConditions)[i];
        minVal = 1;
        maxVal = model.numberOfBoundaryConditions;
        k = [listUtil listGetInteger:model inArray:boundaryConditionAtId.valuesList forVariable:@"periodic bc" info:&found minValue:&minVal maxValue:&maxVal];
        boundaryConditionAtId.pMatrix = [self periodicProjectorInModel:model forMesh:mesh boundary:i target:k-1];
    }
    
    free_imatrix(child, 0, mesh.numberOfBulkElements-1, 0, 7);
    
    [mesh deallocateMeshEdgeTables];
    [mesh deallocateMeshFaceTables];
    
    return newMesh;
}

-(void)SetStabilizationParametersInMesh:(FEMMesh *)mesh model:(FEMModel *)model {
    
    int i, j, n;
    Element_t *elements;
    Nodes_t *nodes, *meshNodes;
    FEMElementDescription *elmDescription;
    
    for (FEMSolution *solution in model.solutions) {
        if (mesh == solution.mesh) {
            mesh.stabilize = (mesh.stabilize == YES || [(solution.solutionInfo)[@"stabilize"] boolValue] == YES) ? YES : NO;
            mesh.stabilize = (mesh.stabilize == YES || [(solution.solutionInfo)[@"stabilization method"] isEqualToString:@"vms"] == YES) ? YES : NO;
            mesh.stabilize = (mesh.stabilize == YES || [(solution.solutionInfo)[@"stabilization method"] isEqualToString:@"stabilized"] == YES) ? YES : NO;
        }
    }
    
    nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    nodes->x = doublevec(0, mesh.maxElementNodes-1);
    nodes->y = doublevec(0, mesh.maxElementNodes-1);
    nodes->z = doublevec(0, mesh.maxElementNodes-1);
    
    elements = mesh.getElements;
    meshNodes = mesh.getNodes;
    
    elmDescription = [[FEMElementDescription alloc] init];
    
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        n = elements[i].Type.NumberOfNodes;
        for (j=0; i<n; j++) {
            nodes->x[j] = meshNodes->x[elements[i].NodeIndexes[j]];
            nodes->y[j] = meshNodes->y[elements[i].NodeIndexes[j]];
            nodes->z[j] = meshNodes->z[elements[i].NodeIndexes[j]];
        }
        if (mesh.stabilize == YES) {
            [elmDescription computeStabilizationParameterInElement:&elements[i] nodes:nodes mesh:mesh numberOfNodes:n mk:elements[i].StabilizationMK hk:&elements[i].hK];
        } else {
            [elmDescription elementDiameter:&elements[i] nodes:nodes];
        }
    }
    
    free_dvector(nodes->x, 0, mesh.maxElementNodes-1);
    free_dvector(nodes->y, 0, mesh.maxElementNodes-1);
    free_dvector(nodes->z, 0, mesh.maxElementNodes-1);
    [elmDescription deallocation];
}

-(void)setCurrentMesh:(FEMMesh *)mesh inModel:(FEMModel *)model {
    
    Element_t *modelElements, *meshElements;
    Nodes_t *modelNodes, *meshNodes;
    
    model.variables = mesh.variables;
    model.mesh = mesh;
    
    modelElements = model.getElements;
    modelNodes = model.getNodes;
    meshElements = mesh.getElements;
    meshNodes = mesh.getNodes;
    
    modelNodes = meshNodes;
    model.numberOfNodes = mesh.numberOfNodes;
    modelNodes->numberOfNodes = meshNodes->numberOfNodes;
    
    modelElements = meshElements;
    model.maxElementNodes = mesh.maxElementNodes;
    model.numberOfBulkElements = mesh.numberOfBulkElements;
    model.numberOfBoundaryElements = mesh.numberOfBoundaryElements;
}

-(void)updateMesh:(FEMMesh *)mesh inSolution:(FEMSolution *)solution model:(FEMModel *)model {
    
    int i, j, k, n, dofs, n1, n2;
    int *permutation;
    double *work;
    BOOL onlySearch, found, optimizeBandwidth;
    NSString *str;
    FEMUtilities *utilities;
    FEMElementUtils *elementUtils;
    FEMMatrix *matrix;
    variableArraysContainer *saveVarContainers = NULL, *varContainers = NULL;
    matrixArraysContainer *matContainers = NULL, *solutionMatContainers = NULL;
    
    saveVarContainers = solution.variable.getContainers;
    dofs = solution.variable.dofs;
    
    solution.mesh = mesh;
    [self setCurrentMesh:mesh inModel:model];
    
    // Create matrix and variable structures for
    // current equation on the mesh
    utilities = [[FEMUtilities alloc] init];
    onlySearch = NO;
    solution.variable = [utilities getVariableFrom:mesh.variables model:model name:solution.variable.name onlySearch:&onlySearch maskName:NULL info:&found];
    
    varContainers = solution.variable.getContainers;
    permutation = intvec(0, varContainers->sizePerm-1);
    
    if ((solution.solutionInfo)[@"optimize bandwidth"] != nil) {
        optimizeBandwidth = [(solution.solutionInfo)[@"optimize bandwidth"] boolValue];
    } else optimizeBandwidth = YES;
    
    elementUtils = [[FEMElementUtils alloc] init];
    str = (solution.solutionInfo)[@"equation"];
    matrix = [elementUtils createMatrixInModel:model forSolution:solution mesh:mesh dofs:dofs permutation:permutation sizeOfPermutation:varContainers->sizePerm matrixFormat:MATRIX_CRS optimizeBandwidth:optimizeBandwidth equationName:str discontinuousGalerkinSolution:NULL globalBubbles:NULL];
    
    if ((solution.solutionInfo)[@"linear system symmetric"] != nil) {
        matrix.symmetric = [(solution.solutionInfo)[@"linear system symmetric"] boolValue];
    }
    
    if ((solution.solutionInfo)[@"lumped mass matrix"] != nil) {
        matrix.lumped = [(solution.solutionInfo)[@"lumped mass matrix"] boolValue];
    }
    
    work = doublevec(0, varContainers->sizeValues-1);
    memcpy(work, varContainers->Values, varContainers->sizeValues);
    for (k=1; k<=dofs; k++) {
        for (i=0; i<varContainers->sizePerm; i++) {
            if (permutation[i] >= 0) {
                varContainers->Values[dofs*(permutation[i]+1)-k] = work[dofs*(varContainers->Perm[i]+1)-k];
            }
        }
    }
    
    if (varContainers->PrevValues != NULL) {
        for (j=0; j<varContainers->size2PrevValues; j++) {
            for (i=0; i<varContainers->size1PrevValues; i++) {
                work[i] = varContainers->PrevValues[i][j];
            }
             for (k=1; k<=dofs; k++) {
                 for (i=0; i<varContainers->sizePerm; i++) {
                     if (permutation[i] >= 0) {
                         varContainers->PrevValues[dofs*(permutation[i]+1)-k][j] = work[dofs*(varContainers->Perm[i]+1)-k];
                     }
                 }
             }
        }
    }
    free_dvector(work, 0, varContainers->sizeValues-1);
    
    memcpy(varContainers->Perm, permutation, varContainers->sizePerm);
    solution.variable.solution = solution;
    free_ivector(permutation, 0, varContainers->sizePerm-1);
    
    matContainers = matrix.getContainers;
    matContainers->RHS = doublevec(0, matrix.numberOfRows-1);
    matContainers->sizeRHS = matrix.numberOfRows;
    
    solutionMatContainers = solution.matrix.getContainers;
    
    if (saveVarContainers->EigenValues != NULL) {
        n = saveVarContainers->sizeEigenValues;
        
        if (n > 0) {
            solution.nOfEigenValues = n;
            varContainers->EigenValues = cdoublevec(0, n-1);
            varContainers->sizeEigenValues = n;
            varContainers->EigenVectors = cdoublematrix(0, n-1, 0, varContainers->sizeValues-1);
            varContainers->size1EigenVectors = n;
            varContainers->size2EigenVectors = varContainers->sizeValues;
            
            for (i=0; i<n; i++) {
                varContainers->EigenValues[i] = 0.0;
                for (j=0; j<varContainers->sizeValues; j++) {
                    varContainers->EigenVectors[i][j] = 0.0;
                }
            }
            matContainers->MassValues = doublevec(0, matContainers->sizeValues-1);
            matContainers->sizeMassValues = matContainers->sizeValues;
            memset( matContainers->MassValues, 0.0, (matContainers->sizeMassValues*sizeof(matContainers->MassValues)) );
        }
    } else if (solutionMatContainers->Force != NULL) {
        n1 = matrix.numberOfRows;
        n2 = solutionMatContainers->size2Force;
        matContainers->Force = doublematrix(0, n1-1, 0, n2-1);
        matContainers->size1force = n1;
        matContainers->size2Force = n2;
        for (i=0; i<n1; i++) {
            for (j=0; j<n2; j++) {
                matContainers->Force[i][j] = 0.0;
            }
        }
    }
    
    solution.matrix = matrix;
    solution.mesh.changed = YES;
}

@end

