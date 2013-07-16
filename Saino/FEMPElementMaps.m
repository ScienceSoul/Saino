//
//  FEMPElementMaps.m
//  Saino
//
//  Created by Seddik hakime on 15/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMPElementMaps.h"
#import "Utils.h"

@implementation FEMPElementMaps {
    
    int **_quadEdgeMap;
    int **_triangleEdgeMap;
    int **_tetraEdgeMap1;
    int **_tetraFaceMap1;
    int **_tetraFaceEdgeMap1;
    int **_tetraEdgeMap2;
    int **_tetraFaceMap2;
    int **_tetraFaceEdgeMap2;
    int **_brickEdgeMap;
    int **_brickFaceMap;
    int **_brickFaceEdgeMap;
    int **_wedgeEdgeMap;
    int **_wedgeFaceMap;
    int **_wedgeFaceEdgeMap;
    int **_pyramidEdgeMap;
    int **_pyramidFaceMap;
    int **_pyramidFaceEdgeMap;
}

- (id)init
{    
    self = [super init];
    if (self) {
        _quadEdgeMap = intmatrix(0, 3, 0, 1);
        _triangleEdgeMap = intmatrix(0, 2, 0, 1);
        _tetraEdgeMap1 = intmatrix(0, 5, 0, 1);
        _tetraFaceMap1 = intmatrix(0, 3, 0, 2);
        _tetraFaceEdgeMap1 = intmatrix(0, 3, 0, 2);
        _tetraEdgeMap2 = intmatrix(0, 5, 0, 1);
        _tetraFaceMap2 = intmatrix(0, 3, 0, 2);
        _tetraFaceEdgeMap2 = intmatrix(0, 3, 0, 2);
        _brickEdgeMap = intmatrix(0, 11, 0, 1);
        _brickFaceMap = intmatrix(0, 5, 0, 3);
        _brickFaceEdgeMap = intmatrix(0, 5, 0, 3);
        _wedgeEdgeMap = intmatrix(0, 8, 0, 1);
        _wedgeFaceMap = intmatrix(0, 4, 0, 3);
        _wedgeFaceEdgeMap = intmatrix(0, 4, 0, 3);
        _pyramidEdgeMap = intmatrix(0, 7, 0, 1);
        _pyramidFaceMap = intmatrix(0, 4, 0, 3);
        _pyramidFaceEdgeMap = intmatrix(0, 4, 0, 3);
        
        NSLog(@"FEMPElementMaps:FEMPElementMaps: initalizing mappings for elements...\n");
        
        // Quad edge mappings
        _quadEdgeMap[0][0] = 0; _quadEdgeMap[0][1] = 1;
        _quadEdgeMap[1][0] = 1; _quadEdgeMap[1][1] = 2;
        _quadEdgeMap[2][0] = 3; _quadEdgeMap[2][1] = 2;
        _quadEdgeMap[3][0] = 0; _quadEdgeMap[3][1] = 3;
        
        // Triangle edge mappings
        _triangleEdgeMap[0][0] = 0; _triangleEdgeMap[0][1] = 1;
        _triangleEdgeMap[1][0] = 1; _triangleEdgeMap[1][1] = 2;
        _triangleEdgeMap[2][0] = 2; _triangleEdgeMap[2][1] = 0;
        
        // Brick edge mappings
        _brickEdgeMap[0][0] = 0; _brickEdgeMap[0][1] = 1;
        _brickEdgeMap[1][0] = 1; _brickEdgeMap[1][1] = 2;
        _brickEdgeMap[2][0] = 3; _brickEdgeMap[2][1] = 2;
        _brickEdgeMap[3][0] = 0; _brickEdgeMap[3][1] = 3;
        _brickEdgeMap[4][0] = 4; _brickEdgeMap[4][1] = 5;
        _brickEdgeMap[5][0] = 5; _brickEdgeMap[5][1] = 6;
        _brickEdgeMap[6][0] = 7; _brickEdgeMap[6][1] = 6;
        _brickEdgeMap[7][0] = 4; _brickEdgeMap[7][1] = 7;
        _brickEdgeMap[8][0] = 0; _brickEdgeMap[8][1] = 4;
        _brickEdgeMap[9][0] = 1; _brickEdgeMap[9][1] = 5;
        _brickEdgeMap[10][0] = 2; _brickEdgeMap[10][1] = 6;
        _brickEdgeMap[11][0] = 3; _brickEdgeMap[11][1] = 7;
        
        // Brick face mappings
        _brickFaceMap[0][0] = 0; _brickFaceMap[0][1] = 1; _brickFaceMap[0][2] = 2; _brickFaceMap[0][3] = 3;  // xi, eta
        _brickFaceMap[1][0] = 4; _brickFaceMap[1][1] = 5; _brickFaceMap[1][2] = 6; _brickFaceMap[1][3] = 7;  // xi, eta
        _brickFaceMap[2][0] = 0; _brickFaceMap[2][1] = 1; _brickFaceMap[2][2] = 5; _brickFaceMap[2][3] = 4;  // xi, eta
        _brickFaceMap[3][0] = 1; _brickFaceMap[3][1] = 2; _brickFaceMap[3][2] = 6; _brickFaceMap[3][3] = 5;  // eta, zeta
        _brickFaceMap[4][0] = 3; _brickFaceMap[4][1] = 2; _brickFaceMap[4][2] = 6; _brickFaceMap[4][3] = 7;
        _brickFaceMap[5][0] = 0; _brickFaceMap[5][1] = 3; _brickFaceMap[5][2] = 7; _brickFaceMap[5][3] = 4;
        
        _brickFaceEdgeMap[0][0] = 0; _brickFaceEdgeMap[0][1] = 1; _brickFaceEdgeMap[0][2] = 2; _brickFaceEdgeMap[0][3] = 3;
        _brickFaceEdgeMap[1][0] = 4; _brickFaceEdgeMap[1][1] = 5; _brickFaceEdgeMap[1][2] = 6; _brickFaceEdgeMap[1][3] = 7;
        _brickFaceEdgeMap[2][0] = 0; _brickFaceEdgeMap[2][1] = 9; _brickFaceEdgeMap[2][2] = 4; _brickFaceEdgeMap[2][3] = 8;
        _brickFaceEdgeMap[3][0] = 1; _brickFaceEdgeMap[3][1] = 10; _brickFaceEdgeMap[3][2] = 5; _brickFaceEdgeMap[3][3] = 9;
        _brickFaceEdgeMap[4][0] = 2; _brickFaceEdgeMap[4][1] = 10; _brickFaceEdgeMap[4][2] = 6; _brickFaceEdgeMap[4][3] = 11;
        _brickFaceEdgeMap[5][0] = 3; _brickFaceEdgeMap[5][1] = 11; _brickFaceEdgeMap[5][2] = 7; _brickFaceEdgeMap[5][3] = 8;
        
        // Tetra edge mappings (not need for enforcing parity)
        // Type 1
        _tetraEdgeMap1[0][0] = 0; _tetraEdgeMap1[0][1] = 1;
        _tetraEdgeMap1[1][0] = 1; _tetraEdgeMap1[1][1] = 2;
        _tetraEdgeMap1[2][0] = 0; _tetraEdgeMap1[2][1] = 2;
        _tetraEdgeMap1[3][0] = 0; _tetraEdgeMap1[3][1] = 3;
        _tetraEdgeMap1[4][0] = 1; _tetraEdgeMap1[4][1] = 3;
        _tetraEdgeMap1[5][0] = 2; _tetraEdgeMap1[5][1] = 3;
        // Type 2
        _tetraEdgeMap2[0][0] = 0; _tetraEdgeMap2[0][1] = 1;
        _tetraEdgeMap2[1][0] = 2; _tetraEdgeMap2[1][1] = 1;
        _tetraEdgeMap2[2][0] = 0; _tetraEdgeMap2[2][1] = 2;
        _tetraEdgeMap2[3][0] = 0; _tetraEdgeMap2[3][1] = 3;
        _tetraEdgeMap2[4][0] = 1; _tetraEdgeMap2[4][1] = 3;
        _tetraEdgeMap2[5][0] = 2; _tetraEdgeMap2[5][1] = 3;
        
        // Tetra face mappings (not need for enforcing parity)
        // Type 1
        _tetraFaceMap1[0][0] = 0; _tetraEdgeMap2[0][1] = 1; _tetraEdgeMap2[0][2] = 2;
        _tetraEdgeMap2[1][0] = 0; _tetraEdgeMap2[1][1] = 1; _tetraEdgeMap2[1][2] = 3;
        _tetraEdgeMap2[2][0] = 1; _tetraEdgeMap2[2][1] = 2; _tetraEdgeMap2[2][2] = 3;
        _tetraEdgeMap2[3][0] = 0; _tetraEdgeMap2[3][1] = 2; _tetraEdgeMap2[3][2] = 3;
        // Type 2
        _tetraFaceMap2[0][0] = 0; _tetraFaceMap2[0][1] = 2; _tetraFaceMap2[0][2] = 1;
        _tetraFaceMap2[1][0] = 0; _tetraFaceMap2[1][1] = 1; _tetraFaceMap2[1][2] = 3;
        _tetraFaceMap2[2][0] = 2; _tetraFaceMap2[2][1] = 1; _tetraFaceMap2[2][2] = 3;
        _tetraFaceMap2[3][0] = 0; _tetraFaceMap2[3][1] = 2; _tetraFaceMap2[3][2] = 3;
        
        // Type 1
        _tetraFaceEdgeMap1[0][0] = 0; _tetraFaceEdgeMap1[0][1] = 1; _tetraFaceEdgeMap1[0][2] = 2;
        _tetraFaceEdgeMap1[1][0] = 0; _tetraFaceEdgeMap1[1][1] = 4; _tetraFaceEdgeMap1[1][2] = 3;
        _tetraFaceEdgeMap1[2][0] = 1; _tetraFaceEdgeMap1[2][1] = 5; _tetraFaceEdgeMap1[2][2] = 4;
        _tetraFaceEdgeMap1[3][0] = 2; _tetraFaceEdgeMap1[3][1] = 5; _tetraFaceEdgeMap1[3][2] = 3;
        // Type 2
        _tetraFaceEdgeMap2[0][0] = 2; _tetraFaceEdgeMap2[0][1] = 1; _tetraFaceEdgeMap2[0][2] = 0;
        _tetraFaceEdgeMap2[1][0] = 0; _tetraFaceEdgeMap2[1][1] = 4; _tetraFaceEdgeMap2[1][2] = 3;
        _tetraFaceEdgeMap2[2][0] = 1; _tetraFaceEdgeMap2[2][1] = 4; _tetraFaceEdgeMap2[2][2] = 5;
        _tetraFaceEdgeMap2[3][0] = 2; _tetraFaceEdgeMap2[3][1] = 5; _tetraFaceEdgeMap2[3][3] = 3;
        
        // Wedge edge mappings
        _wedgeEdgeMap[0][0] = 0; _wedgeEdgeMap[0][1] = 1;
        _wedgeEdgeMap[1][0] = 1; _wedgeEdgeMap[1][1] = 2;
        _wedgeEdgeMap[2][0] = 2; _wedgeEdgeMap[2][1] = 0;
        _wedgeEdgeMap[3][0] = 3; _wedgeEdgeMap[3][1] = 4;
        _wedgeEdgeMap[4][0] = 4; _wedgeEdgeMap[4][1] = 5;
        _wedgeEdgeMap[5][0] = 5; _wedgeEdgeMap[5][1] = 3;
        _wedgeEdgeMap[6][0] = 0; _wedgeEdgeMap[6][1] = 3;
        _wedgeEdgeMap[7][0] = 1; _wedgeEdgeMap[7][1] = 4;
        _wedgeEdgeMap[8][0] = 2; _wedgeEdgeMap[8][1] = 5;
        
        // Wedge face mappings
        _wedgeFaceMap[0][0] = 0; _wedgeFaceMap[0][1] = 1; _wedgeFaceMap[0][2] = 2; _wedgeFaceMap[0][3] = -1; // -1 -> not used
        _wedgeFaceMap[1][0] = 3; _wedgeFaceMap[1][1] = 4; _wedgeFaceMap[1][2] = 5; _wedgeFaceMap[1][3] = -1;
        _wedgeFaceMap[2][0] = 0; _wedgeFaceMap[2][1] = 1; _wedgeFaceMap[2][2] = 4; _wedgeFaceMap[2][3] = 3;
        _wedgeFaceMap[3][0] = 1; _wedgeFaceMap[3][1] = 2; _wedgeFaceMap[3][2] = 5; _wedgeFaceMap[3][3] = 4;
        _wedgeFaceMap[4][0] = 2; _wedgeFaceMap[4][1] = 0; _wedgeFaceMap[4][2] = 3; _wedgeFaceMap[4][3] = 5;
        
        _wedgeFaceEdgeMap[0][0] = 0; _wedgeFaceEdgeMap[0][1] = 1; _wedgeFaceEdgeMap[0][2] = 2; _wedgeFaceEdgeMap[0][3] = -1;
        _wedgeFaceEdgeMap[1][0] = 3; _wedgeFaceEdgeMap[1][1] = 4; _wedgeFaceEdgeMap[1][2] = 5; _wedgeFaceEdgeMap[1][3] = -1;
        _wedgeFaceEdgeMap[2][0] = 0; _wedgeFaceEdgeMap[2][1] = 7; _wedgeFaceEdgeMap[2][2] = 3; _wedgeFaceEdgeMap[2][3] = 6;
        _wedgeFaceEdgeMap[3][0] = 1; _wedgeFaceEdgeMap[3][1] = 8; _wedgeFaceEdgeMap[3][2] = 4; _wedgeFaceEdgeMap[3][3] = 7;
        _wedgeFaceEdgeMap[4][0] = 2; _wedgeFaceEdgeMap[4][1] = 6; _wedgeFaceEdgeMap[4][2] = 5; _wedgeFaceEdgeMap[4][3] = 8;
        
        // Pyramid edge mappings
        _pyramidEdgeMap[0][0] = 0; _pyramidEdgeMap[0][1] = 1;
        _pyramidEdgeMap[1][0] = 1; _pyramidEdgeMap[1][1] = 2;
        _pyramidEdgeMap[2][0] = 3; _pyramidEdgeMap[2][1] = 2;
        _pyramidEdgeMap[3][0] = 0; _pyramidEdgeMap[3][1] = 3;
        _pyramidEdgeMap[4][0] = 0; _pyramidEdgeMap[4][1] = 4;
        _pyramidEdgeMap[5][0] = 1; _pyramidEdgeMap[5][1] = 4;
        _pyramidEdgeMap[6][0] = 2; _pyramidEdgeMap[6][1] = 4;
        _pyramidEdgeMap[7][0] = 3; _pyramidEdgeMap[7][1] = 4;
        
        // Pyramid face mappings
        _pyramidFaceMap[0][0] = 0; _pyramidFaceMap[0][1] = 1; _pyramidFaceMap[0][2] = 2; _pyramidFaceMap[0][3] = 3;
        _pyramidFaceMap[1][0] = 0; _pyramidFaceMap[1][1] = 1; _pyramidFaceMap[1][2] = 4; _pyramidFaceMap[1][3] = -1;
        _pyramidFaceMap[2][0] = 1; _pyramidFaceMap[2][1] = 2; _pyramidFaceMap[2][2] = 4; _pyramidFaceMap[2][3] = -1;
        _pyramidFaceMap[3][0] = 2; _pyramidFaceMap[3][1] = 3; _pyramidFaceMap[3][2] = 4; _pyramidFaceMap[3][3] = -1;
        _pyramidFaceMap[4][0] = 3; _pyramidFaceMap[4][1] = 0; _pyramidFaceMap[4][2] = 4; _pyramidFaceMap[4][4] = -1;
        
        _pyramidFaceEdgeMap[0][0] = 0; _pyramidFaceEdgeMap[0][1] = 1; _pyramidFaceEdgeMap[0][2] = 2; _pyramidFaceEdgeMap[0][3] = 3;
        _pyramidFaceEdgeMap[1][0] = 0; _pyramidFaceEdgeMap[1][1] = 5; _pyramidFaceEdgeMap[1][2] = 4; _pyramidFaceEdgeMap[1][3] = -1;
        _pyramidFaceEdgeMap[2][0] = 1; _pyramidFaceEdgeMap[2][1] = 6; _pyramidFaceEdgeMap[2][2] = 5; _pyramidFaceEdgeMap[2][3] = -1;
        _pyramidFaceEdgeMap[3][0] = 2; _pyramidFaceEdgeMap[3][1] = 7; _pyramidFaceEdgeMap[3][2] = 6; _pyramidFaceEdgeMap[3][3] = -1;
        _pyramidFaceEdgeMap[4][0] = 3; _pyramidFaceEdgeMap[4][1] = 4; _pyramidFaceEdgeMap[4][2] = 7; _pyramidFaceEdgeMap[4][3] = -1;
    }
    
    return self;
}

-(void)deallocation {
    free_imatrix(_quadEdgeMap, 0, 3, 0, 1);
    free_imatrix(_triangleEdgeMap, 0, 2, 0, 1);
    free_imatrix(_tetraEdgeMap1, 0, 5, 0, 1);
    free_imatrix(_tetraFaceMap1, 0, 3, 0, 2);
    free_imatrix(_tetraFaceEdgeMap1, 0, 3, 0, 2);
    free_imatrix(_tetraEdgeMap2, 0, 5, 0, 1);
    free_imatrix(_tetraFaceMap2, 0, 3, 0, 2);
    free_imatrix(_tetraFaceEdgeMap2, 0, 3, 0, 2);
    free_imatrix(_brickEdgeMap, 0, 11, 0, 1);
    free_imatrix(_brickFaceMap, 0, 5, 0, 3);
    free_imatrix(_brickFaceEdgeMap, 0, 5, 0, 3);
    free_imatrix(_wedgeEdgeMap, 0, 8, 0, 1);
    free_imatrix(_wedgeFaceMap, 0, 4, 0, 3);
    free_imatrix(_wedgeFaceEdgeMap, 0, 4, 0, 3);
    free_imatrix(_pyramidEdgeMap, 0, 7, 0, 1);
    free_imatrix(_pyramidFaceMap, 0, 4, 0, 3);
    free_imatrix(_pyramidFaceEdgeMap, 0, 4, 0, 3);
}


/******************************************************************
 
    Based on element bubble polynomial degree p, return degrees of
    freedom for given element bubbles.

******************************************************************/
-(int)bubbleDofsForElement:(Element_t *)element degree:(int)p {
    
    int bubbleDofs;
    
    // Not defined for non p elements
    if (element->Pdefs == NULL) {
        return bubbleDofs = 0;
    }
    
    // Select by element type
    bubbleDofs = 0;
    switch (element->Type.ElementCode / 100) {
        case 2: // Line
            if (p >= 2) bubbleDofs = p-1;
            break;
        case 3: // Triangle
            if (p >= 3) bubbleDofs = (p-1)*(p-2)/2;
            break;
        case 4: // Quad
            if (p >= 4) bubbleDofs = (p-2)*(p-3)/2;
            break;
        case 5: // Tetrahedron
            if (p >= 4) bubbleDofs = (p-1)*(p-2)*(p-3)/6;
            break;
        case 6: // Pyramid
            if (p >= 4) bubbleDofs = (p-1)*(p-2)*(p-3)/6;
            break;
        case 7: // Wedge
            if (p >= 5) bubbleDofs = (p-2)*(p-3)*(p-4)/6;
            break;
        case 8: // Brick
            if (p >= 6) bubbleDofs = (p-3)*(p-4)*(p-5)/6;
            break;
        default:
            NSLog(@"FEMPElementMaps:bubbleDofsForElement: unsupported p element type.\n");
            bubbleDofs = p;
            break;
    }
    
    bubbleDofs = max(0, bubbleDofs);
    return bubbleDofs;
}

-(void)getTriangleEdgeMap:(int *)edge index:(int)i {
    
    int j;
    
    for (j=0; j<2; j++) {
        edge[j] = _triangleEdgeMap[i][j];
    }
}

-(void)getQuadEdgeMap:(int *)edge index:(int)i {
    
    int j;
    
    for (j=0; j<2; j++) {
        edge[j] = _quadEdgeMap[i][j];
    }
}

-(void)getBrickEdgeMap:(int *)edge index:(int)i {
    
    int j;
    
    for (j=0; j<2; j++) {
        edge[j] = _brickEdgeMap[i][j];
    }
}

-(void)getBrickFaceMap:(int *)face index:(int)i {
    
    int j;
    
    for (j=0; j<4; j++) {
        face[j] = _brickFaceMap[i][j];
    }
}

-(int)getBrickFaceEdgeMap:(int)face localNode:(int)node {
    
    int localEdge;
    
    localEdge = _brickFaceEdgeMap[face][node];
    
    if (localEdge < 0) {
        NSLog(@"FEMPElementMaps:getBrickFaceEdgeMap: unknown combination node for (face,node): %d %d.\n", face, node);
        errorfunct("FEMPElementMaps:getBrickFaceEdgeMap", "Terminating now!!!");
    }    
    return localEdge;
}


-(void)getTetraEdgeMap:(int *)edge index:(int)i type:(int *)type {
    
    int j, t;
    
    // If not given, use default (1)
    t = 1;
    if (type != NULL) t = *type;
    
    //Select face type by tetra type
    switch (t) {
        case 1:
            for (j=0; j<2; j++) {
                edge[j] = _tetraEdgeMap1[i][j];
            }
            break;
        case 2:
            for (j=0; j<2; j++) {
                edge[j] = _tetraEdgeMap2[i][j];
            }
            break;
        default:
            errorfunct("FEMPElementMaps:getTetraEdgeMap", "Unknown tetra type.");
            break;
    }
}

-(void)getTetraFaceMap:(int *)face index:(int)i type:(int *)type {
    
    int j, t;
    
    // If not given, use default (1)
    t = 1;
    if (type != NULL) t = *type;
    
    //Select face type by tetra type
    switch (t) {
        case 1:
            for (j=0; j<3; j++) {
                face[j] = _tetraFaceMap1[i][j];
            }
            break;
        case 2:
            for (j=0; j<3; j++) {
                face[j] = _tetraFaceMap2[i][j];
            }
            break;
        default:
            errorfunct("FEMPElementMaps:getTetraFaceMap", "Unknown tetra type.");
            break;
    }
}

-(void)getWedgeEdgeMap:(int *)edge index:(int)i {
    
    int j;
    
    for (j=0; j<2; j++) {
        edge[j] = _wedgeEdgeMap[i][j];
    }
}

-(void)getWedgeFaceMap:(int *)face index:(int)i {
    
    int j;
    
    for (j=0; j<4; j++) {
        face[j] = _wedgeFaceMap[i][j];
    }
}

-(void)getPyramidEdgeMap:(int *)edge index:(int)i {
    
    int j;
    
    for (j=0; j<2; j++) {
        edge[j] = _pyramidEdgeMap[i][j];
    }
}

-(void)getPyramidFaceMap:(int *)face index:(int)i {
    
    int j;
    
    for (j=0; j<4; j++) {
        face[j] = _pyramidFaceMap[i][j];
    }
}

-(BOOL)isPElement:(Element_t *)element {
    
    if (element->Pdefs != NULL) {
        return YES;
    } else return NO;
}

/***************************************************************
 
 This method gets the reference p element nodes since they are 
 not yet defined in the FEMElementDescription class
 
***************************************************************/
-(void)getRefPElementNodesForElement:(Element_t *)element nodeU:(double *)u nodeV:(double *)v nodeW:(double *)w {
    
    int n;
    
    if ([self isPElement:element] == NO) {
        NSLog(@"FEMPElementMaps:getRefPElementNodesForElement: element given not a p element!\n");
        return;
    }
    
    // Reverse space for element nodes
    n = element->Type.NumberOfNodes;
    
    // Select by element type
    switch (element->Type.ElementCode / 100) {
        case 2: // Line
            u[0] = -1.0;
            u[1] = 1.0;
            break;
        case 3: // Triangle
            u[0] = -1.0; u[1] = 1.0; u[2] = 0.0;
            v[0] = 0.0; v[1] = 0.0; v[2] = sqrt(3.0);
            break;
        case 4: // Quad
            u[0] = -1.0; u[1] = 1.0; u[2] = 1.0; u[3] = -1.0;
            v[0] = -1.0; v[1] = -1.0; v[2] = 1.0; v[3] = 1.0;
            break;
        case 5: // Tetrahedron
            u[0] = -1.0; u[1] = 1.0; u[2] = 0.0; u[3] = 0.0;
            v[0] = 0.0; v[1] = 0.0; v[2] = sqrt(3.0); v[3] = 1.0/sqrt(3.0);
            w[0] = 0.0; w[1] = 0.0; w[2] = 0.0; w[3] = 2.0*sqrt(2.0/3.0);
            break;
        case 6: // Pyramid
            u[0] = -1.0; u[1] = 1.0; u[2] = 1.0; u[3] = -1.0; u[4] = 0.0;
            v[0] = -1.0; v[1] = -1.0; v[2] = 1.0; v[3] = 1.0; v[4] = 0.0;
            w[0] = 0.0; w[1] = 0.0; w[2] = 0.0; w[3] = 0.0; w[4] = sqrt(2.0);
            break;
        case 7: // Wedge
            u[0] = -1.0; u[1] = 1.0; u[2] = 0.0; u[3] = -1.0; u[4] = 1.0; u[5] = 0.0;
            v[0] = 0.0; v[1] = 0.0; v[2] = sqrt(3.0); v[3] = 0.0; v[4] = 0.0; v[5] = sqrt(3.0);
            w[0] = -1.0; w[1] = -1.0; w[2] = -1.0; w[3] = 1.0; w[4] = 1.0; w[5] = 1.0;
            break;
        case 8: // Brick
            u[0] = -1.0; u[1] = 1.0; u[2] = 1.0; u[3] = -1.0; u[4] = -1.0; u[5] = 1.0; u[6] = 1.0; u[7] = -1.0;
            v[0] = -1.0; v[1] = -1.0; v[2] = 1.0; v[3] = 1.0; v[4] = -1.0; v[5] = -1.0; v[6] = 1.0; v[7] = 1.0;
            w[0] = -1.0; w[1] = -1.0; w[2] = -1.0; w[3] = -1.0; w[4] = 1.0; w[5] = 1.0; w[6] = 1.0; w[7] = 1.0;
            break;
        default:
            NSLog(@"FEMPElementMaps:getRefPElementNodesForElement: unknown element type.\n");
            break;
    }
}

/**************************************************************************************
 
    Get mappings for given elemen to element faces and their nodes. Given the element,
    this method returns a map containing nodes (endpoints) of element face.
 
    Arguments:
        Element_t *element -> Element to get map for
        int **map          -> Map containing local numbers of local faces

***************************************************************************************/
-(void)getFaceMapForElement:(Element_t *)element faceMap:(int **)map {
    
    // Not defined for non p elements
    if (element->Pdefs == NULL) {
        NSLog(@"FEMPElementMaps:getFaceMapForElement: element is not a p element.");
        map = NULL;
        return;
    }
    
    switch (element->Type.ElementCode / 100) {
        case 5:
            switch (element->Pdefs->TetraType) {
                case 1:
                    map = _tetraFaceMap1;
                    break;
                case 2:
                    map = _tetraFaceMap2;
                    break;
                default:
                    errorfunct("FEMPElementMaps:getFaceMapForElement", "Unknown tetra type for p element.");
                    break;
            }
            break;
        case 6:
            map = _pyramidFaceMap;
            break;
        case 7:
            map = _wedgeFaceMap;
            break;
        case 8:
            map = _brickFaceMap;
            break;
        default:
            errorfunct("FEMPElementMaps:getFaceMapForElement", "Unsupported element type.");
            break;
    }
}

/**************************************************************************************
 
 Get mappings for given element to element edges and their nodes. Given the element, 
 this method returns a map containing nodes (endpoints) of elements edges.
 
 Arguments:
 Element_t *element -> Element to get map for
 int **map          -> Map containing local numbers of local faces
 
 ***************************************************************************************/
-(void)getEdgeMapForElement:(Element_t *)element edgeMap:(int **)map {
    
    // Not defined for non p elements
    if (element->Pdefs == NULL) {
        NSLog(@"FEMPElementMaps:getEdgeMapForElement: element is not a p element.");
        map = NULL;
        return;
    }
    
    switch (element->Type.ElementCode / 100) {
        case 3:
            map = _triangleEdgeMap;
            break;
        case 4:
            map = _quadEdgeMap;
            break;
        case 5:
            switch (element->Pdefs->TetraType) {
                case 1:
                    map = _tetraEdgeMap1;
                    break;
                case 2:
                    map = _tetraEdgeMap2;
                    break;
                default:
                     errorfunct("FEMPElementMaps:getEdgeMapForElement", "Unknown tetra type for p element.");
                    break;
            }
            break;
        case 6:
            map = _pyramidEdgeMap;
            break;
        case 7:
            map = _wedgeEdgeMap;
            break;
        case 8:
            map =_brickEdgeMap;
            break;
        default:
            errorfunct("FEMPElementMaps:getEdgeMapForElement", "Unsupported element type.");
            break;
    }
}

/**************************************************************************************
 
 Get mappings for given element to element faces and their edges. Given the element,
 this method returns a map containing local edge number of elements faces.
 
 Arguments:
 Element_t *element -> Element to get map for
 int **map          -> Map containing local numbers of local faces
 
 ***************************************************************************************/
-(void)getFaceEdgeMapForElement:(Element_t *)element faceEdgeMap:(int **)map {
    
    // Not defined for non p elements
    if (element->Pdefs == NULL) {
        NSLog(@"FEMPElementMaps:getFaceEdgeMapForElement: element is not a p element.");
        map = NULL;
        return;
    }
    
    switch (element->Type.ElementCode / 100) {
        case 5:
            switch (element->Pdefs->TetraType) {
                case 1:
                    map = _tetraFaceEdgeMap1;
                    break;
                case 2:
                    map = _tetraFaceEdgeMap2;
                    break;
                default:
                    errorfunct("FEMPElementMaps:getFaceEdgeMapForElement", "Unknown tetra type for p element.");
                    break;
            }
            break;
        case 6:
            map = _pyramidFaceEdgeMap;
            break;
        case 7:
            map = _wedgeFaceEdgeMap;
            break;
        case 8:
            map =_brickFaceEdgeMap;
            break;
        default:
            errorfunct("FEMPElementMaps:getFaceEdgeMapForElement", "Unsupported element type.");
            break;
    }
}

/**********************************************************************
 
    Method checks if given element is a p element pyramid

**********************************************************************/
-(BOOL)isPPyramid:(Element_t *)element {
    
    return (element->Type.ElementCode/100 == 6 && [self isPElement:element] == YES) ? YES : NO;
}

/**********************************************************************
 
    Mapping from element local edge or face number to nodes contained
    in that edge or face.
 
**********************************************************************/
-(void)getBoundaryMapForElement:(Element_t *)element localNumber:(int)i resultMap:(int *)map {
    
    memset( map, 0, 4*sizeof(int) );
    
    // Method not defined for non p elements
    if (element->Pdefs == NULL) {
        NSLog(@"FEMPElementMaps:getBoundaryMapForElement: element not p element.");
        return;
    }
    
    switch (element->Type.ElementCode / 100) {
        case 3:
            [self getTriangleEdgeMap:map index:i];
            break;
        case 4:
            [self getQuadEdgeMap:map index:i];
            break;
        case 5:
            [self getTetraFaceMap:map index:i type:&element->Pdefs->TetraType];
            break;
        case 6:
            [self getPyramidFaceMap:map index:i];
            break;
        case 7:
            [self getWedgeFaceMap:map index:i];
            break;
        case 8:
            [self getBrickFaceMap:map index:i];
            break;
        default:
            errorfunct("FEMPElementMaps:getBoundaryMapForElement", "Unsupported element type.");
            break;
    }
}

/**********************************************************************
 
    Mapping from element local face to local edges contained in face.
    Given element and local face number, this routine returns numbers
    of local edges on face.
 
**********************************************************************/
-(void)getFaceEdgeMapForElement:(Element_t *)element index:(int)i resultMap:(int *)map {
    
    int j;
    
    memset( map, 0, 4*sizeof(int) );
    
    // Method not defined for non p elements
    if (element->Pdefs == NULL) {
        NSLog(@"FEMPElementMaps:getFaceEdgeMapForElement: element not p element.");
        return;
    }
    
    switch (element->Type.ElementCode / 100) {
        case 5:
            switch (element->Pdefs->TetraType) {
                case 1:
                    for (j=0; j<3; j++) {
                        map[j] = _tetraFaceEdgeMap1[i][j];
                    }
                    break;
                case 2:
                    for (j=0; j<3; j++) {
                        map[j] = _tetraFaceEdgeMap2[i][j];
                    }
                    break;
                default:
                    errorfunct("FEMPElementMaps:getFaceEdgeMapForElement", "Unknown tetra type.");
                    break;
            }
            break;
        case 6:
            for (j=0; j<4; j++) {
                map[j] = _pyramidFaceEdgeMap[i][j];
            }
            break;
        case 7:
            for (j=0; j<4; j++) {
                map[j] = _wedgeFaceEdgeMap[i][j];
            }
            break;
        case 8:
            for (j=0; j<4; j++) {
                map[j] = _brickFaceEdgeMap[i][j];
            }
            break;
        default:
            errorfunct("FEMPElementMaps:getFaceEdgeMapForElement", "Unsupported element type.");
            break;
    }
}

/**********************************************************************
 
    Based on element face polynomial degree p, return degrees of
    freedom for given face.
 
    Element_t *element  ->  element to get face dofs
    int p               ->  face polynomial degree p
    int number          ->  local number of face for element 
                            (important for wedges and pyramids)
 
    Return number of face dofs for element.
 
**********************************************************************/
-(int)getFaceDofsForElement:(Element_t *)element polyDegree:(int)p faceNumber:(int)number {
    
    int faceDofs;
    
    // This function is not defined for non p elements
    if (element->Pdefs == NULL) {
        return faceDofs = 0;
    }
    
    faceDofs = 0;
    switch (element->Type.ElementCode / 100) {
        case 5: // Tetrahedron
            if (p >=3) faceDofs = (p-1)*(p-2)/2;
            break;
        case 6: // Pyramid
            switch (number) {
                case 0:
                    if (p >= 4) faceDofs = (p-2)*(p-3)/2;
                    break;
                case 1:
                case 4:
                    if (p >= 3) faceDofs = (p-1)*(p-2)/2;
                    break;
            }
            break;
        case 7: // Wedge
            switch (number) {
                case 0:
                case 1:
                    if (p >= 3) faceDofs = (p-1)*(p-2)/2;
                    break;
                case 2:
                case 4:
                    if (p >= 4) faceDofs = (p-2)*(p-3)/2;
                    break;
            }
            break;
        case 8: // Brick
            if (p >= 4) faceDofs = (p-2)*(p-3)/2;
            break;
        default:
            NSLog(@"FEMPElementMaps:getFaceDofsForElement: unsupported p element type.");
            faceDofs = p;
            break;
    }
    
    return faceDofs = max(0, faceDofs);
}

-(int)getNumberOfGaussPointsForFace:(Element_t *)face inMesh:(FEMMesh *)mesh {
    
    int i, edgep, maxp;
    Element_t *edges;
    
    edges = mesh.getEdges;
    
    // Get max p of edges contained in face
    edgep = 0;
    for (i=0; i<face->Type.NumberOfEdges; i++) {
        edgep = max(edgep, edges[face->EdgeIndexes[i]].Pdefs->p);
    }
    
    // If no face dofs use max edge dofs as number of points
    if (face->BDOFs <= 0) {
        return pow((edgep+1), 2);
    }
    
    edges = NULL;
    
    // Get number of Gauss points
    maxp = max(edgep, face->Pdefs->p) + 1;
    return pow(maxp, 2);
}

-(int)getEdgePForElement:(Element_t *)element inMesh:(FEMMesh *)mesh {
    
    int i, edgep;
    Element_t *edges;
    
    if (element->Pdefs == NULL) {
        NSLog(@"FEMPElementMaps:getEdgePForElement: element not p element");
        return edgep = 0;
    }
    
    edges = mesh.getEdges;
    
    // Get max p of edges of element if any
    edgep = 0;
    if (element->EdgeIndexes != NULL) {
        for (i=0; element->Type.NumberOfEdges; i++) {
            // If edge has no dofs, it effectively has degree of 0
            if (edges[element->EdgeIndexes[i]].BDOFs <= 0) continue;
            edgep = max(edgep, edges[element->EdgeIndexes[i]].Pdefs->p);
        }
    }
    
    edges = NULL;
    return edgep;
}

-(int)getFacePForElement:(Element_t *)element inMesh:(FEMMesh *)mesh {
    
    int i, facep;
    Element_t *faces;
    
    if (element->Pdefs == NULL) {
        NSLog(@"getEdgePForElement: element not p element");
        return facep = 0;
    }
    
    faces = mesh.getFaces;
    
    facep = 0;
    if (element->FaceIndexes != NULL) {
        for (i=0; i<element->Type.NumberOfFaces; i++) {
            // If face has no dofs, it effectively has degree of 0
            if (faces[element->FaceIndexes[i]].BDOFs <= 0) continue;
            facep = max(facep, faces[element->FaceIndexes[i]].Pdefs->p);
        }
    }
    
    faces= NULL;
    return facep;
}

/****************************************************************************
  
    Based on element bubble polynomial degree p, return degrees of freedom
    for given elements bubbles.
 
    ELement_t *element   ->  element to get bubble dofs to
    int p                ->  element polynomial degree p
 
    Returns number of bubble dofs for element

*****************************************************************************/
-(int)getBubbleDofsForElement:(Element_t *)element polyDegree:(int)p {
    
    int bubbleDofs;
    
    if (element->Pdefs == NULL) {
        NSLog(@"FEMPElementMaps:getBubbleDofsForElement: element not p element");
        return bubbleDofs = 0;
    }
    
    bubbleDofs = 0;
    switch (element->Type.ElementCode / 100) {
        case 2: // Line
            if (p >= 2) bubbleDofs = p - 1;
            break;
        case 3: // Triangle
            if (p >= 3) bubbleDofs = (p-1)*(p-2)/2;
            break;
        case 4: // Quad
            if (p >= 4) bubbleDofs = (p-2)*(p-3)/2;
            break;
        case 5: // Tetrahedron
            if (p >= 4) bubbleDofs = (p-1)*(p-2)*(p-3)/6;
            break;
        case 6: // Pyramid
            if (p >= 4) bubbleDofs = (p-1)*(p-2)*(p-3)/6;
            break;
        case 7: // Wedge
            if (p >= 5) bubbleDofs = (p-2)*(p-3)*(p-4)/6;
            break;
        case 8: // Brick
            if (p >= 6) bubbleDofs = (p-3)*(p-4)*(p-5)/6;
            break;
        default:
            NSLog(@"getBubbleDofsForElement: unsupported p element type");
            bubbleDofs = p;
            break;
    }
    
    return bubbleDofs = max(0,bubbleDofs);
}

-(int)getNumberOfGaussPointsForElement:(Element_t *)element inMesh:(FEMMesh *)mesh {
    
    int edgep, facep, bubblep, nb, maxp;
    
    if (element->Pdefs == NULL) {
        NSLog(@"FEMPElementMaps:getNumberOfGaussPointsForElement: element not p element.");
        return 0;
    }
    
    // Max p of edges
    edgep = 0;
    if (element->Type.dimension == 2 || element->Type.dimension == 3) {
        edgep = [self getEdgePForElement:element inMesh:mesh];
    }
    
    // Max p of faces
    facep = 0;
    if (element->Type.dimension == 3) {
        facep = [self getFacePForElement:element inMesh:mesh];
    }
    
    // Element bubble p
    bubblep = 0;
    if (element->BDOFs > 0) {
        bubblep = element->Pdefs->p;
        
        switch (element->Type.ElementCode / 100) {
            case 3:
                nb = max([self getBubbleDofsForElement:element polyDegree:bubblep], element->BDOFs);
                bubblep = (int)round( (3.0+sqrt(1.0+8.0*nb)) / 2.0 );
                break;
            case 4:
                nb = max([self getBubbleDofsForElement:element polyDegree:bubblep], element->BDOFs);
                bubblep = (int)round( (5.0+sqrt(1.0+8.0*nb)) / 2.0 );
                break;
            case 5:
                nb = max([self getBubbleDofsForElement:element polyDegree:bubblep], element->BDOFs);
                bubblep = (int)round(1.0/3.0*pow(81.0*nb+3.0*sqrt(-3.0+729.0*pow(nb, 2.0)), 1.0/3.0) + 1.0 / pow(81.0*nb+3.0*sqrt(-3.0+729.0*pow(nb, 2.0)), 1.0/3.0) + 2.0);
                break;
            case 6:
                nb = max([self getBubbleDofsForElement:element polyDegree:bubblep], element->BDOFs);
                bubblep = (int)round(1.0/3.0*pow(81.0*nb+3.0*sqrt(-3.0+729.0*pow(nb, 2.0)), 1.0/3.0) + 1.0 / pow(81.0*nb+3.0*sqrt(-3.0+729.0*pow(nb, 2.0)), 1.0/3.0) + 2.0);
                break;
            case 7:
                nb = max([self getBubbleDofsForElement:element polyDegree:bubblep], element->BDOFs);
                bubblep = (int)round(1.0/3.0*pow(81.0*nb+3.0*sqrt(-3.0+729.0*pow(nb, 2.0)), 1.0/3.0) + 1.0 / pow(81.0*nb+3.0*sqrt(-3.0+729.0*pow(nb, 2.0)), 1.0/3.0) + 3.0);
                break;
            case 8:
                nb = max([self getBubbleDofsForElement:element polyDegree:bubblep], element->BDOFs);
                bubblep = (int)round(1.0/3.0*pow(81.0*nb+3.0*sqrt(-3.0+729.0*pow(nb, 2.0)), 1.0/3.0) + 1.0 / pow(81.0*nb+3.0*sqrt(-3.0+729.0*pow(nb, 2.0)), 1.0/3.0) + 4.0);
                break;
        }
    }
    
    // Get number of Gauss points. Number needed is 2*max(p) = 2*(max(p)+1)/2
    maxp = max(1, edgep, facep, bubblep) + 1;
    return pow(maxp, element->Type.dimension);
}

-(int)getEdgeDofsForElement:(Element_t *)element polyDegree:(int)p {
    
    if (element->Pdefs == NULL) {
        return 0;
    }
    
    return max(0, p-1);
}

@end
