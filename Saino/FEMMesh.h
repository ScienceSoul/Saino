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

@interface FEMMesh : NSObject {
  
    int _dimension;
    int _numberOfNodes;
    int _numberOfElements;
    int _numberOfBulkElements; 
    int _numberOfEdges; 
    int _numberOfFaces;
    int _numberOfBoundaryElements;
    int _numberOfViewFactors;
    int _maxElementNodes; 
    int _maxElementDofs; 
    int _maxEdgeDofs;
    int _maxFaceDofs; 
    int _maxBdofs;
    int _numberOfPassiveBCs;
    int _savesDone;
    int _numberOfColors;
    BOOL _outputActive;
    BOOL _adaptiveMesh;
    BOOL _changed;
    BOOL _stabilize;
    
    Element_t *_elements;
    Element_t *_edges;
    Element_t *_faces;
    Nodes_t *_globalNodes;
    Quadrant_t *_rootQuadrant;
    Factors_t *_viewFactors;
    
    NSMutableString *_name;
    NSMutableArray *_variables;              // Mutable array of FEMVariable classes
    NSMutableArray *_projectors;             // Mutable array of FEMProjector classes
    NSMutableArray *_next;                   // Mutable array of FEMMesh classes
    NSMutableArray *_colors;                 // Mutable array of mutable arrays for each color. Each array stores the number of elements for that color,
                                             // The color index and the corresponding RGB values for that color as: [number of elements, color index, R, G, B]
    FEMMesh *_parent;
    FEMMesh *_child;
    
    int *_colorMapping;
}

@property(nonatomic, assign) int dimension;
@property(nonatomic, assign) int numberOfNodes;
@property(nonatomic, assign) int numberOfElements;
@property(nonatomic, assign) int numberOfBulkElements;
@property(nonatomic, assign) int numberOfEdges;
@property(nonatomic, assign) int numberOfFaces;
@property(nonatomic, assign) int numberOfBoundaryElements;
@property(nonatomic, assign) int numberOfViewFactors;
@property(nonatomic, assign) int maxElementNodes;
@property(nonatomic, assign) int maxElementDofs;
@property(nonatomic, assign) int maxEdgeDofs;
@property(nonatomic, assign) int maxFaceDofs;
@property(nonatomic, assign) int maxBdofs;
@property(nonatomic, assign) int numberOfPassiveBCs;
@property(nonatomic, assign) int savesDone;
@property(nonatomic, assign) int numberOfColors;
@property(nonatomic, assign, getter = isOutputActive) BOOL outputActive;
@property(nonatomic, assign, getter = isAdaptiveMesh) BOOL adaptiveMesh;
@property(nonatomic, assign, getter = isChanged) BOOL changed;
@property(nonatomic, assign, getter = isStabilize) BOOL stabilize;
@property(nonatomic, strong) NSMutableString *name;
@property(nonatomic, strong) NSMutableArray *variables;
@property(nonatomic, strong) NSMutableArray *projectors;
@property(nonatomic, strong) NSMutableArray *next;
@property(nonatomic, strong) NSMutableArray *colors;
@property(nonatomic, strong) FEMMesh *parent;
@property(nonatomic, strong) FEMMesh *child;

-(void)loadMeshForModel:(FEMModel *)model meshDirectory:(NSString *)dir meshName:(NSString *)name boundariesOnly:(BOOL)bd numberOfPartitions:(int *)numParts partitionID:(int *)partID definitions:(int *)defDofs;

//Allocations methods
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

// Color mapping getter
-(int *)getColorMapping;

// Test associativity
-(BOOL)isAssociatedEdges;
-(BOOL)isAssociatedFaces;

// Deallocation
-(void)deallocation;
-(void)deallocationMeshVariables;
-(void)deallocateMeshEdgeTables;
-(void)deallocateMeshFaceTables;
-(void)deallocateMeshViewFactorTables;
-(void)deallocateQuadrantTree;

-(void)Simple2DMeshBorders:(double*)borders withSize:(int*) meshSize elemetCode:(int) elementID;







@end
