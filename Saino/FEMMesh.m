//
//  FEMMesh.m
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import "FEMMesh.h"

#include "memory.h"

@implementation FEMMesh

-(void)Simple2DMesh:Borders:(double*)borders withSize:(int*) intervals elemetCode:(int) elementID {
    
    int i, j, k, l;
    int m, n, o, p;
    double minx, maxx, miny, maxy, incri, incrj;
    double meshsize1, meshsize2;
    int domainShape;
    int **side1, **side2, **side3, **side4;
    BOOL ct;
    
    minx = borders[0];
    maxx = borders[1];
    miny = borders[2];
    maxy = borders[3];
    
    // Rectangle or square?
    if ( (maxx-minx) > 0 && (maxy-miny) > 0 && (maxx-minx) == (maxy-miny) ) {
        domainShape = 1;                                                               // valid square
    } else  if ((maxx-minx) > 0 && (maxy-miny) > 0 && (maxx-minx) != (maxy-miny)) {
        domainShape = 2;                                                              // valid rectangle
    }
    else {
        domainShape = -1;                                                             // Something wrong with borders input
    }
    
    [self setDimension:2];
    
    // Build nodes
    switch (domainShape) {
        case 1:
            
            meshsize1 = (maxx-minx) / intervals[0];
            meshsize2 = meshsize1;
            intervals[1] = intervals[0];
            
            // Allocate memory for the GlobalNodes
            [self AllocateNodes:0 :((intervals[0]+1)*(intervals[0]+1))-1];
            
            break;
            
        case 2:
            
            if (intervals[1] <= 0) {
                errorfunct("Simple2DMesh", "No intervals given for second direction discretization in rectangle domain.");
            }
            
            meshsize1 = (maxx-minx) / intervals[0];
            meshsize2 = (maxy-miny) / intervals[1];
            
            // Allocate memory for the GlobalNodes
            [self AllocateNodes:0 :((intervals[0]+1)*(intervals[1]+1))-1];
            
            break;
            
        default:
            errorfunct("Simple2DMesh", "Failure in method Simple2DMesh. Cant initialize mesh.");
            break;
    }
    
    l = 1;
    k = 0;
    incri = minx;
    incrj = miny;
    
    for (i=0; i<=intervals[0]; i++) {
        for (j=0; j<=intervals[1]; j++) {
            
            [self setGlabalNodes_ID:k :l];
            [self setGlabalNodes_p:k :-1];
            [self setGlobalNodes_x:k :incri];
            [self setGlobalNodes_y:k :incrj];
            [self setGlobalNodes_z:k :0.0];
            incrj = incrj + meshsize2;
            l++;
            k++;
            
        }
        incrj = miny;
        incri = incri + meshsize1;
    }
            
             
    // Build elements and remember boundary elements
    switch (domainShape) {
        case 1:
            
            // Allocate memory for the Elements
            //Elements = ElementVec(0, (intervals[0]*intervals[0]));
            //Elements->NodeIndexes = intvec(0, 3);
            
            side1 = intmatrix(0, intervals[0], 0, 2);
            side2 = intmatrix(0, intervals[0], 0, 2);
            side3 = intmatrix(0, intervals[0], 0, 2);
            side4 = intmatrix(0, intervals[0], 0, 2);
            
            break;
            
        case 2:
            
            //Elements = ElementVec(0, (intervals[0]*intervals[1]));
            //Elements->NodeIndexes = intvec(0, 3);
            
            side1 = intmatrix(0, intervals[0], 0, 2);
            side2 = intmatrix(0, intervals[1], 0, 2);
            side3 = intmatrix(0, intervals[0], 0, 2);
            side4 = intmatrix(0, intervals[1], 0, 2);
            
            break;
            
    }

    k = 0;
    l = 1;
    
    m = 0;
    n = 0;
    o = 0;
    p = 0;
    
    for (i=0; i<((intervals[0]+1)*(intervals[1]+1))-(intervals[1]+2); i++) {
        
        // Top boundary elements
        for (j=1;j<(intervals[1]+1);j++) {
            if ( i == 0+(j*(intervals[1]+1)-1) ) {
                side3[o][0] = [self GlobalNodes_ID:0+((j+1)*(intervals[1]+1)-1)];
                side3[o][1] = [self GlobalNodes_ID:i];
                side3[o][2] = k-1;
                o++;
                ct = YES;
                break;
            }
        }
        
        if (ct == YES) {
            ct = NO;
            continue;
        }
        
        [self setElements_ElementIndex:k :l];
        [self setElements_BodyID:k :1];
        [self setElements_Type_ElementCode:k :404];
        [self setElements_NodeIndexes:k :0 :[self GlobalNodes_ID:i]];
        [self setElements_NodeIndexes:k :1 :[self GlobalNodes_ID:i+(intervals[1]+1)]];
        [self setElements_NodeIndexes:k :2 :[self GlobalNodes_ID:i+((intervals[1]+1)+1)]];
        [self setElements_NodeIndexes:k :3 :[self GlobalNodes_ID:i+1]];
        
        // Left boundary elements
        if ( i < (intervals[1]+1) ) {
            side4[p][0] = [self GlobalNodes_ID:i+1];
            side4[p][1] = [self GlobalNodes_ID:i];
            side4[p][2] = k;
            p++;
        }
        
        // Bottom boundary elements
        for ( j=0;j<(intervals[1]+1);j++ ) {
            if ( i == 0+j*(intervals[1]+1) ) {
                side1[m][0] = [self GlobalNodes_ID:i];
                side1[m][1] = [self GlobalNodes_ID:i+(intervals[1]+1)];
                side1[m][2] = k;
                m++;
                break;
            }
        }
        
        // Right boundary elements
        if ( i >= ((intervals[0]+1)*(intervals[1]+1))-(2*(intervals[1]+1)) ) {
            side2[n][0] = [self GlobalNodes_ID:i+(intervals[1]+1)];
            side2[n][1] = [self GlobalNodes_ID:i+((intervals[1]+1)+1)];
            side2[n][2] = k;
            n++;
            break;
        }
        
        k++;
        l++;
    }
            
    // Build boundary elements
    
    switch (domainShape) {
        case 1:
            
            // Allocate memory for the boundary elements
            //BDElements = BDElementVec(0, (intervals[0]*4)-1);
            //BDElements->NodeIndexes = intvec(0, 1);
            
            break;
            
        case 2:
            
            //BDElements = BDElementVec(0, ((intervals[0]*2)+(intervals[1]*2))-1);
            //BDElements->NodeIndexes = intvec(0, 1);
            
            break;
            
    }
    
    k = 0;
    l = 1;
    
    // Boundaries are numbered anti-clockwise
    // First boundary -> bottom side
    for (i=0; i<m; i++) {
        
        [self setBDElements_ID:k :l];
        [self setBDElements_boundary:k :1];
        [self setBDElements_p1:k :side1[m][2]];
        [self setBDElements_p2:k :0];
        [self setDElements_code:k :202];
        [self setBDElements_NodeIndexes:k :0 :side1[m][0]];
        [self setBDElements_NodeIndexes:k :1 :side1[m][1]];
        k++;
        l++;
    }
    
    // Second boundary -> right side
    for (i=0; i<n; i++) {
        
        [self setBDElements_ID:k :l];
        [self setBDElements_boundary:k :2];
        [self setBDElements_p1:k :side2[n][2]];
        [self setBDElements_p2:k :0];
        [self setDElements_code:k :202];
        [self setBDElements_NodeIndexes:k :0 :side2[n][0]];
        [self setBDElements_NodeIndexes:k :1 :side2[n][1]];
        k++;
        l++;
    }
    
    // Third boundary -> top side
    for (i=0; i<o; i++) {
        
        [self setBDElements_ID:k :l];
        [self setBDElements_boundary:k :3];
        [self setBDElements_p1:k :side3[o][2]];
        [self setBDElements_p2:k :0];
        [self setDElements_code:k :202];
        [self setBDElements_NodeIndexes:k :0 :side3[o][0]];
        [self setBDElements_NodeIndexes:k :1 :side3[o][1]];
        k++;
        l++;
    }
    
    // Fourth boundary -> left side
    for (i=0; i<p; i++) {
        
        [self setBDElements_ID:k :l];
        [self setBDElements_boundary:k :4];
        [self setBDElements_p1:k :side4[p][2]];
        [self setBDElements_p2:k :0];
        [self setDElements_code:k :202];
        [self setBDElements_NodeIndexes:k :0 :side4[p][0]];
        [self setBDElements_NodeIndexes:k :1 :side4[p][1]];
        k++;
        l++;
    }
    
    switch (domainShape) {
        case 1:

            free_imatrix(side1, 0, intervals[0], 0, 2);
            free_imatrix(side2, 0, intervals[0], 0, 2);
            free_imatrix(side3, 0, intervals[0], 0, 2);
            free_imatrix(side4, 0, intervals[0], 0, 2);
            
            break;
            
        case 2:
            
            free_imatrix(side1, 0, intervals[0], 0, 2);
            free_imatrix(side2, 0, intervals[1], 0, 2);
            free_imatrix(side3, 0, intervals[0], 0, 2);
            free_imatrix(side4, 0, intervals[1], 0, 2);

            break;
    }
    
}   

-(int)dimension {
    
    return dimension;
}

-(int)numberOfNodes {
    
    return numberOfNodes;
}

-(int)numberOfBulkElements {
    
    return numberOfBulkElements;
}

-(int)numberOfEdges {
    
    return numberOfEdges;
}

-(int)numberOfFaces {
    
    return numberOfFaces;
}

-(int)numberOfBoundaryElements {
    
    return numberOfBoundaryElements;
}

-(int)maxElementNodes {
    
    return maxElementNodes;
}

-(int)maxElementDofs {
    
    return maxElementDofs;
}

-(int)maxFaceDofs {
    
    return maxFaceDofs;
}

-(int)maxEdgeDofs {
    
    return maxEdgeDofs;
}

-(int)maxBdofs {
    
    return maxBdofs;
}

-(void)setDimension:(int)n {
    
    dimension = n;
}

-(void)setNumberOfNodes:(int)n {
    
    numberOfNodes = n;
}

-(void)setNumberOfBulkElements:(int)n {
    
    numberOfBulkElements = n;
}

-(void)setNumberOfEdges:(int)n {
    
    numberOfEdges = n;
}

-(void)setNumberOfFaces:(int)n {
    
    numberOfFaces = n;
}

-(void)setNumberOfBoundaryElements:(int)n {
    
    numberOfBoundaryElements = n;
}

-(void)setMaxElementNodes:(int)n {
    
    maxElementNodes = n;
}

-(void)setMaxElementDofs:(int)n {
    
    maxElementDofs = n;
}

-(void)setMaxFaceDofs:(int)n {
    
    maxFaceDofs = n;
}

-(void)setMaxEdgeDofs:(int)n {
    
    maxEdgeDofs = n;
}

-(void)setMaxBdofs:(int)n {
    
    maxBdofs = n;
}

#pragma mark Nodes_t and Element_t allocations

-(BOOL)AllocateNodes:(int)nl :(int)nh {
    
    GlobalNodes = nodesvec(nl,nh-1);
    if (GlobalNodes == NULL) {
        return NO;
    }
    else {
        return YES;
    }
}

#pragma mark Nodes fields getters/setters

-(int)GlabalNodes_ID:(int)i {
    
    return GlobalNodes[i].ID;
}

-(int)GlabalNodes_p:(int)i {
    
    return GlobalNodes[i].p;
}

-(double)GlobalNodes_x:(int)i {
    
    return GlobalNodes[i].x;
}

-(double)GlobalNodes_y:(int)i {
    
    return GlobalNodes[i].y;
}

-(double)GlobalNodes_z:(int)i {
    
    return GlobalNodes[i].z;
}

-(void)setGlabalNodes_ID:(int)i :(int)n {
    
    GlobalNodes[i].ID = n;
}

-(void)setGlabalNodes_p:(int)i :(int)n {
    
    GlobalNodes[i].p = n;
}

-(void)setGlobalNodes_x:(int)i :(double)n {
    
    GlobalNodes[i].x = n;
}

-(void)setGlobalNodes_y:(int)i :(double)n {
    
    GlobalNodes[i].y = n;
}

-(void)setGlobalNodes_z:(int)i :(double)n {
    
    GlobalNodes[i].z = n;
}

#pragma mark Elements getters/setters

-(Element_t *)returnElementAtIndex:(int)i {
    
    return &Elements[i];
}

-(int)Elements_NodeIndexes:(int)i :(int)j {
    
    return Elements[i].NodeIndexes[j];
}

-(int)Elements_BodyID:(int)i {
    
    return Elements[i].BodyID;
}

-(int)Elements_ElementIndex:(int)i {
    
    return Elements[i].ElementIndex;
}

-(int)Elements_NDOFs:(int)i {
    
    return Elements[i].NDOFs;
    
}

-(double)Elements_StabilizationMK:(int)i {
    
    return Elements[i].StabilizationMK;
}

-(double)Elements_hk:(int)i {
    
    return Elements[i].hK;
}

-(int)Elements_Type_ElementCode:(int)i {
    
    return Elements[i].Type.ElementCode;
}

-(int)Elements_Type_BasisFunctionDegree:(int)i {
    
    return Elements[i].Type.BasisFunctionDegree;
}

-(int)Elements_Type_NumberOfNodes:(int)i {
    
    return Elements[i].Type.NumberOfNodes;
}

-(int)Elements_Type_NumberOfEdges:(int)i {
    
    return Elements[i].Type.NumberOfEdges;
}

-(int)Elements_Type_NumberOfFaces:(int)i {
    
    return Elements[i].Type.NumberOfFaces;
}

-(int)Elements_Type_dimension:(int)i {
    
    return Elements[i].Type.dimension;
}

-(int)Elements_Type_GaussPoints:(int)i {
    
    return Elements[i].Type.GaussPoints;
}

-(int)Elements_Type_GaussPoints2:(int)i {
    
    return Elements[i].Type.GaussPoints2;
}

-(double)Elements_Type_StabilizationMK:(int)i {
    
    return Elements[i].Type.StabilizationMK;
}

-(double)Elements_Type_NodeU:(int)i :(int)j {
    
    return Elements[i].Type.NodeU[j];
}

-(double)Elements_Type_NodeV: (int)i: (int)j {
    
    return Elements[i].Type.NodeV[j];
}

-(double)Elements_Type_NodeW: (int)i: (int)j {
    
    return Elements[i].Type.NodeW[j];
}

-(int)Elements_Type_BasisFunctions_n:(int)i: (int)j {
    
    return Elements[i].Type.BasisFunctions[j].n;
}

-(int)Elements_Type_BasisFunctions_p:(int)i: (int)j: (int)k {
    
    return Elements[i].Type.BasisFunctions[j].p[k];
}

-(int)Elements_Type_BasisFunctions_q:(int)i: (int)j: (int)k {
    
    return Elements[i].Type.BasisFunctions[j].q[k];
}

-(int)Elements_Type_BasisFunctions_r:(int)i: (int)j: (int)k {
    
    return Elements[i].Type.BasisFunctions[j].r[k];
}

-(double)Elements_Type_BasisFunctions_coeff:(int)i: (int)j: (int)k {
    
    return Elements[i].Type.BasisFunctions[j].coeff[k];
}

-(void)setElements_NodeIndexes:(int)i :(int)j :(int)n {
    
    Elements[i].NodeIndexes[j] = n;
}

-(void)setElements_BodyID:(int)i :(int)n {
    
    Elements[i].BodyID = n;
}

-(void)setElements_ElementIndex:(int)i :(int)n {
    
    Elements[i].ElementIndex = n;
}

-(void)setElements_NDOFs:(int)i :(int)n {
    
    Elements[i].NDOFs = n;
}

-(void)setElements_StabilizationMK:(int)i :(double)n {
    
    Elements[i].StabilizationMK = n;
}

-(void)setElements_hK:(int)i :(double)n {
    
    Elements[i].hK = n;
}

-(void)setElements_Type_ElementCode:(int)i :(int)n {
    
    Elements[i].Type.ElementCode = n;
}

-(void)setElements_Type_BasisFunctionDegree:(int)i :(int)n {
    
    Elements[i].Type.BasisFunctionDegree = n;
}

-(void)setElements_Type_NumberOfNodes:(int)i :(int)n {
    
    Elements[i].Type.NumberOfNodes = n;
}

-(void)setElements_Type_NumberOfEdges:(int)i :(int)n {
    
    Elements[i].Type.NumberOfEdges = n;
}

-(void)setElements_Type_NumberOfFAces:(int)i :(int)n {
    
    Elements[i].Type.NumberOfFaces = n;
}

-(void)setElements_Type_dimension: (int)i: (int)n {
    
    Elements[i].Type.dimension = n;
}

-(void)setElements_Type_GaussPoints:(int)i :(int)n {
    
    Elements[i].Type.GaussPoints = n;
}

-(void)setElements_Type_GaussPoints2:(int)i :(int)n {
    
    Elements[i].Type.GaussPoints2 = n;
}

-(void)setElements_Type_StabilizationMK:(int)i :(double)n {
    
    Elements[i].Type.StabilizationMK = n;
}

-(void)setElements_Type_NodeU:(int)i :(int)j :(double)n {
    
    Elements[i].Type.NodeU[j] = n;
}

-(void)setElements_Type_NodeV:(int)i :(int)j :(double)n {
    
    Elements[i].Type.NodeV[j] = n;
}

-(void)setElements_Type_NodeW:(int)i :(int)j :(double)n {
    
    Elements[i].Type.NodeW[j] = n;
}

-(void)setElements_Type_BasisFunctions_n:(int)i :(int)j :(int)n {
    
    Elements[i].Type.BasisFunctions[j].n = n;
}

-(void)setElements_Type_BasisFunctions_p:(int)i :(int)j :(int)k :(int)n {
    
    Elements[i].Type.BasisFunctions[j].p[k] = n;
}

-(void)setElements_Type_BasisFunctions_q: (int)i: (int)j: (int)k: (int)n {
    
    Elements[i].Type.BasisFunctions[j].q[k] = n;
}

-(void)setElements_Type_BasisFunctions_r: (int)i: (int)j: (int)k: (int)n {
    
    Elements[i].Type.BasisFunctions[j].r[k] = n;
}

-(void)setElements_Type_BasisFunctions_coeff: (int)i: (int)j: (int)k: (double)n {
    
    Elements[i].Type.BasisFunctions[j].coeff[k] = n;
}

#pragma mark Edges getters/setters

-(Element_t *)returnEdgeAtIndex:(int)i {
    
    return &Edges[i];
}

-(int)Edges_NodeIndexes:(int)i :(int)j {
    
    return Edges[i].NodeIndexes[j];
}

-(int)Edges_Type_NumberOfNodes:(int)i {
    
    return Edges[i].Type.NumberOfNodes;
}

-(int)Edges_BDOFs:(int)i {
    
    return Edges[i].BDOFs;
    
}

-(void)setEdges_NodeIndexes:(int)i :(int)j :(int)n {
    
    Edges[i].NodeIndexes[j] = n;
}

-(void)setEdges_Type_NumberOfNodes:(int)i :(int)n {
    
    Edges[i].Type.NumberOfNodes = n;
}

-(void)setEdges_BDOFs: (int)i: (int)n {
    
    Edges[i].BDOFs = n;
}

#pragma mark Faces getters/setters

-(Element_t *)returnFaceAtIndex:(int)i {
    
    return &Faces[i];
}

-(int)Faces_NodeIndexes:(int)i :(int)j {
    
    return Faces[i].NodeIndexes[j];
}

-(int)Faces_Type_NumberOfNodes:(int)i {
    
    return Faces[i].Type.NumberOfNodes;
}

-(int)Faces_BDOFs:(int)i {
    
    return Faces[i].BDOFs;
    
}

-(void)setFaces_NodeIndexes:(int)i :(int)j :(int)n {
    
    Faces[i].NodeIndexes[j] = n;
}

-(void)setFaces_Type_NumberOfNodes:(int)i :(int)n {
    
    Faces[i].Type.NumberOfNodes = n;
}

-(void)setFaces_BDOFs: (int)i: (int)n {
    
    Faces[i].BDOFs = n;
}

#pragma mark BDElement_t fields getters/setters

-(void)setBDElements_ID:(int)i :(int)n {
    
    BDElements[i].ID = n;
}

-(void)setBDElements_boundary:(int)i :(int)n {
    
    BDElements[i].boundary = n;
}

-(void)setBDElements_p1:(int)i :(int)n {
    
    BDElements[i].p1 = n;
}

-(void)setBDElements_p2:(int)i :(int)n {
    
    BDElements[i].p2 = n;
}

-(void)setDElements_code:(int)i :(int)n {
    
    BDElements[i].code = n;
}

-(void)setBDElements_NodeIndexes:(int)i :(int)j :(int)n {
    
    BDElements[i].NodeIndexes[j] = n;
}

#pragma mark variables instance methods

-(Variable_t *)returnPointerToVariables {
    
    return variables;
}

#pragma mark Test associativity

-(BOOL)isAssociatedEdges {
    
    if (Edges != NULL) {
        return YES;
    } else {
        return NO;
    }
    
}

-(BOOL)isAssociatedFaces {
    
    if (Faces != NULL) {
        return YES;
    } else {
        return NO;
    }

}

#pragma mark Sizes

-(int)sizeOfGlobalNodes {
    
    return sizeOfGlobalNodes;
}

-(void)setSizeOfGlobalNodes:(int)n {
    
    sizeOfGlobalNodes = n;
}



@end
