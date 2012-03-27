//
//  FEMMesh.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"

@interface FEMMesh : NSObject {
  
    int dimension;
    int numberOfNodes, numberOfBulkElements, numberOfEdges, numberOfFaces, 
        numberOfBoundaryElements, maxElementNodes, maxElementDofs, maxEdgeDofs, 
        maxFaceDofs, maxBdofs;
    
    Element_t *Elements, *Edges, *Faces;
    BDElement_t *BDElements;
    Nodes_t *GlobalNodes;
    Variable_t *variables;
    
}

-(void)Simple2DMesh:Borders:(double*)borders withSize:(int*) meshSize elemetCode:(int) elementID;

-(int)dimension;
-(int)numberOfNodes;
-(int)numberOfEdges;
-(int)numberOfFaces;
-(int)maxElementNodes;
-(int)maxFaceDofs;
-(int)maxEdgeDofs;
-(int)maxBdofs;

-(void)setDimension: (int) n;
-(void)setNumberOfNodes:(int)n;
-(void)setNumberOfEdges:(int)n;
-(void)setNumberOfFaces:(int)n;
-(void)setMaxElementNodes:(int)n;
-(void)setMaxFaceDofs: (int)n;
-(void)setMaxEdgeDofs: (int)n;
-(void)setMaxBdofs: (int)n;


//Allocations methods
-(BOOL)AllocateNodes: (int)nl: (int)nh;
// To do's: Add diallocation methods........

// Getters and setters for node structure fields
-(int) GlobalNodes_ID: (int)i;
-(int) GlobalNodes_p: (int)i;
-(double) GlobalNodes_x: (int)i;
-(double) GlobalNodes_y: (int)i;
-(double) GlobalNodes_z: (int)i;

-(void)setGlabalNodes_ID: (int)i: (int)n;
-(void)setGlabalNodes_p: (int)i: (int)n;
-(void)setGlobalNodes_x: (int)i: (double)n;
-(void)setGlobalNodes_y: (int)i: (double)n;
-(void)setGlobalNodes_z: (int)i: (double)n;

// Getters and setters for element structure fields
-(int)Elements_NodeIndexes:(int)i: (int)j;
-(int)Elements_BodyID:(int)i;
-(int)Elements_ElementIndex:(int)i;
-(int)Elements_NDOFs:(int)i;
-(double)Elements_StabilizationMK:(int)i;
-(double)Elements_hk:(int)i;
-(int)Elements_Type_ElementCode: (int)i;
-(int)Elements_Type_BasisFunctionDegree: (int)i;
-(int)Elements_Type_NumberOfNodes: (int)i;
-(int)Elements_Type_NumberOfEdges: (int)i;
-(int)Elements_Type_NumberOfFaces: (int)i;
-(int)Elements_Type_dimension: (int)i;
-(int)Elements_Type_GaussPoints: (int)i;
-(int)Elements_Type_GaussPoints2: (int)i;
-(double)Elements_Type_StabilizationMK: (int)i;
-(double)Elements_Type_NodeU: (int)i: (int)j;
-(double)Elements_Type_NodeV: (int)i: (int)j;
-(double)Elements_Type_NodeW: (int)i: (int)j;
-(int)Elements_Type_BasisFunctions_n:(int)i: (int)j;
-(int)Elements_Type_BasisFunctions_p:(int)i: (int)j: (int)k;
-(int)Elements_Type_BasisFunctions_q:(int)i: (int)j: (int)k;
-(int)Elements_Type_BasisFunctions_r:(int)i: (int)j: (int)k;
-(double)Elements_Type_BasisFunctions_coeff:(int)i: (int)j: (int)k;

-(void)setElements_NodeIndexes: (int)i: (int)j: (int)n;
-(void)setElements_BodyID: (int)i: (int)n;
-(void)setElements_ElementIndex: (int)i: (int)n;
-(void)setElements_NDOFs: (int)i: (int)n;
-(void)setElements_StabilizationMK: (int)i: (double)n;
-(void)setElements_hK: (int)i: (double)n;
-(void)setElements_Type_ElementCode: (int)i: (int)n;
-(void)setElements_Type_BasisFunctionDegree: (int)i: (int)n;
-(void)setElements_Type_NumberOfNodes: (int)i: (int)n;
-(void)setElements_Type_NumberOfEdges: (int)i: (int)n;
-(void)setElements_Type_NumberOfFaces: (int)i: (int)n;
-(void)setElements_Type_dimension: (int)i: (int)n;
-(void)setElements_Type_GaussPoints: (int)i: (int)n;
-(void)setElements_Type_GaussPoints2: (int)i: (int)n;
-(void)setElements_Type_StabilizationMK: (int)i: (double)n;
-(void)setElements_Type_NodeU: (int)i: (int)j: (double)n;
-(void)setElements_Type_NodeV: (int)i: (int)j: (double)n;
-(void)setElements_Type_NodeW: (int)i: (int)j: (double)n;
-(void)setElements_Type_BasisFunctions_n: (int)i: (int)j: (int)n;
-(void)setElements_Type_BasisFunctions_p: (int)i: (int)j: (int)k: (int)n;
-(void)setElements_Type_BasisFunctions_q: (int)i: (int)j: (int)k: (int)n;
-(void)setElements_Type_BasisFunctions_r: (int)i: (int)j: (int)k: (int)n;
-(void)setElements_Type_BasisFunctions_coeff: (int)i: (int)j: (int)k: (double)n;


// Getters and setters for edges structure fields
-(int)Edges_NodeIndexes:(int)i: (int)j;
-(int)Edges_Type_NumberOfNodes: (int)i;
-(int)Edges_BDOFs:(int)i;

-(void)setEdges_NodeIndexes: (int)i: (int)j: (int)n;
-(void)setEdges_Type_NumberOfNodes: (int)i: (int)n;
-(void)setEdges_BDOFs: (int)i: (int)n;


// Getters and setters for faces structure fields
-(int)Faces_NodeIndexes:(int)i: (int)j;
-(int)Faces_Type_NumberOfNodes: (int)i;
-(int)Faces_BDOFs:(int)i;

-(void)setFaces_NodeIndexes: (int)i: (int)j: (int)n;
-(void)setFaces_Type_NumberOfNodes: (int)i: (int)n;
-(void)setFaces_BDOFs: (int)i: (int)n;


// Getters and setters for boundary element structure fields
-(void)setBDElements_ID: (int)i: (int)n;
-(void)setBDElements_boundary: (int)i: (int)n;
-(void)setBDElements_p1: (int)i: (int)n;
-(void)setBDElements_p2: (int)i: (int)n;
-(void)setDElements_code:(int)i: (int)n;
-(void)setBDElements_NodeIndexes: (int)i: (int)j: (int)n;

// Getters for Variables
-(Variable_t *)returnPointerToVariables;




@end
