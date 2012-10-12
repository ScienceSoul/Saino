//
//  FEMMesh.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"
#import "FEMVariable.h"
#import "SIOMeshIO.h"
#import "Constructors.h"

#import "Utils.h"

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
    BOOL _savesDone;
    BOOL _outputActive;
    BOOL _adaptiveMesh;
    BOOL _changed;
    BOOL _stabilize;
    
    Element_t *_elements, *_edges, *_faces;
    Nodes_t *_globalNodes;
    Quadrant_t *_rootQuadrant;
    
    NSMutableArray *_variables;            // Mutable array of FEMVariable classes
    NSMutableArray *_projectors;           // Mutable array of FEMProjector classes
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
@property(nonatomic, assign, getter = isSavesDone) BOOL savesDone;
@property(nonatomic, assign, getter = isOutputActive) BOOL outputActive;
@property(nonatomic, assign, getter = isAdaptiveMesh) BOOL adaptiveMesh;
@property(nonatomic, assign, getter = isChanged) BOOL changed;
@property(nonatomic, assign, getter = isStabilize) BOOL stabilize;
@property(nonatomic, strong) NSMutableArray *variables;
@property(nonatomic, strong) NSMutableArray *projectors;

-(void)loadMeshForModel:(FEMModel *)model meshDirectory:(NSString *)dir meshName:(NSString *)name boundariesOnly:(BOOL)bd numberOfPartitions:(int *)numParts partitionID:(int *)partID definitions:(int *)defDofs;

//Allocations methods
// TODO: Add diallocation methods........
-(BOOL)AllocateNodes: (int)nl: (int)nh;
-(void)allocatePDefinitionsForElement:(Element_t *)element;

// Nodes getter
-(Nodes_t *)getNodes;

// Elements getter
-(Element_t *)getElements;

// Edges getter 
-(Element_t *)getEdges;

// Faces getter
-(Element_t *)getFaces;

// Quadrant getter
-(Quadrant_t *)getQuadrant;

// Test associativity
-(BOOL)isAssociatedEdges;
-(BOOL)isAssociatedFaces;

// Deallocation
-(void)deallocateQuadrantTree;

-(void)Simple2DMesh:Borders:(double*)borders withSize:(int*) meshSize elemetCode:(int) elementID;







@end
