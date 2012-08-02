//
//  FEMMesh.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMVariable.h"
#import "Constructors.h"

@interface FEMMesh : NSObject {
  
    int _dimension;
    int _numberOfNodes; 
    int _numberOfBulkElements; 
    int _numberOfEdges; 
    int _numberOfFaces;
    int _numberOfBoundaryElements; 
    int _maxElementNodes; 
    int _maxElementDofs; 
    int _maxEdgeDofs;
    int _maxFaceDofs; 
    int _maxBdofs;
    int _sizeOfGlobalNodes;
    
    Element_t *_elements, *_edges, *_faces;
    Nodes_t *_globalNodes;
    
    NSMutableDictionary *_variables;                // Mutable dictionary holding FEMVariable classes
}

@property(nonatomic, assign) int dimension;
@property(nonatomic, assign) int numberOfNodes;
@property(nonatomic, assign) int numberOfBulkElements;
@property(nonatomic, assign) int numberOfEdges;
@property(nonatomic, assign) int numberOfFaces;
@property(nonatomic, assign) int numberOfBoundaryElements;
@property(nonatomic, assign) int maxElementNodes;
@property(nonatomic, assign) int maxElementDofs;
@property(nonatomic, assign) int maxEdgeDofs;
@property(nonatomic, assign) int maxFaceDofs;
@property(nonatomic, assign) int maxBdofs;
@property(nonatomic, assign) int sizeOfGlobalNodes;
@property(nonatomic, strong) NSMutableDictionary *variables;

//Allocations methods
-(BOOL)AllocateNodes: (int)nl: (int)nh;
// TODO: Add diallocation methods........

// Nodes getter
-(Nodes_t *)getNodes;

// Elements getter
-(Element_t *)getElements;

// Edges getter 
-(Element_t *)getEdges;

// Faces getter
-(Element_t *)getFaces;

// Test associativity
-(BOOL)isAssociatedEdges;
-(BOOL)isAssociatedFaces;

-(void)Simple2DMesh:Borders:(double*)borders withSize:(int*) meshSize elemetCode:(int) elementID;







@end
