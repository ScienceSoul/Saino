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

@synthesize dimension = _dimension;
@synthesize numberOfNodes = _numberOfNodes;
@synthesize numberOfBulkElements = _numberOfBulkElements;
@synthesize numberOfEdges = _numberOfEdges;
@synthesize numberOfFaces = _numberOfFaces;
@synthesize numberOfBoundaryElements = _numberOfBoundaryElements;
@synthesize maxElementNodes = _maxElementNodes;
@synthesize maxElementDofs = _maxElementDofs;
@synthesize maxEdgeDofs = _maxEdgeDofs;
@synthesize maxFaceDofs = _maxFaceDofs;
@synthesize maxBdofs = _maxBdofs;
@synthesize sizeOfGlobalNodes = _sizeOfGlobalNodes;

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
            
            _globalNodes[k].ID = l;
            _globalNodes[k].p = -1;
            _globalNodes[k].x = incri;
            _globalNodes[k].y = incrj;
            _globalNodes[k].z = 0.0;
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
                side3[o][0] = _globalNodes[0+((j+1)*(intervals[1]+1)-1)].ID;
                side3[o][1] = _globalNodes[i].ID;
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
        
        _elements[k].ElementIndex = l;
        _elements[k].BodyID = 1;
        _elements[k].Type.ElementCode = 404;
        _elements[k].NodeIndexes[0] =  _globalNodes[i].ID;
        _elements[k].NodeIndexes[1] = _globalNodes[i+(intervals[1]+1)].ID;
        _elements[k].NodeIndexes[2] = _globalNodes[i+((intervals[1]+1)+1)].ID;
        _elements[k].NodeIndexes[3] = _globalNodes[i+1].ID;
        
        // Left boundary elements
        if ( i < (intervals[1]+1) ) {
            side4[p][0] = _globalNodes[i+1].ID;
            side4[p][1] = _globalNodes[i].ID;
            side4[p][2] = k;
            p++;
        }
        
        // Bottom boundary elements
        for ( j=0;j<(intervals[1]+1);j++ ) {
            if ( i == 0+j*(intervals[1]+1) ) {
                side1[m][0] = _globalNodes[i].ID;
                side1[m][1] = _globalNodes[i+(intervals[1]+1)].ID;
                side1[m][2] = k;
                m++;
                break;
            }
        }
        
        // Right boundary elements
        if ( i >= ((intervals[0]+1)*(intervals[1]+1))-(2*(intervals[1]+1)) ) {
            side2[n][0] = _globalNodes[i+(intervals[1]+1)].ID;
            side2[n][1] = _globalNodes[i+((intervals[1]+1)+1)].ID;
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
        
        //TODO: 
    }
    
    // Second boundary -> right side
    for (i=0; i<n; i++) {
        
        //TODO:
    }
    
    // Third boundary -> top side
    for (i=0; i<o; i++) {
        
        //TODO:
    }
    
    // Fourth boundary -> left side
    for (i=0; i<p; i++) {
        
        //TODO:
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

#pragma mark Nodes_t and Element_t allocations

-(BOOL)AllocateNodes:(int)nl :(int)nh {
    
    _globalNodes = nodesvec(nl,nh-1);
    if (_globalNodes == NULL) {
        return NO;
    }
    else {
        return YES;
    }
}

#pragma mark Nodes getter

-(Nodes_t *)getNodes {
    
    return _globalNodes;
}

#pragma mark Elements getter

-(Element_t *)getElements {
    
    return _elements;
}

#pragma mark Edges getter

-(Element_t *)getEdges {
    
    return _edges;
}

#pragma mark Faces getter

-(Element_t *)getFaces {
    
    return _faces;
}

#pragma mark variables instance methods

-(Variable_t *)getVariables {
    
    return _variables;
}

#pragma mark Test associativity

-(BOOL)isAssociatedEdges {
    
    if (_edges != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedFaces {
    
    if (_faces != NULL) {
        return YES;
    } else {
        return NO;
    }
}

@end
