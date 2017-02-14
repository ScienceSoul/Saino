//===----------------------------------------------------------------------===//
//  FEMMesh.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright (c) 2011 ScienceSoul. All rights reserved.
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

#import <Foundation/Foundation.h>
#import "FEMModel.h"
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
    BOOL _adaptiveMesh;
    BOOL _changed;
    BOOL _discontinuousMesh;
    BOOL _stabilize;
    BOOL _outputActive;
    
    Element_t * __nullable _elements;
    Element_t * __nullable _edges;
    Element_t * __nullable _faces;
    Nodes_t * __nullable _globalNodes;
    Quadrant_t * __nullable _rootQuadrant;
    Factors_t * __nullable _viewFactors;
    
    NSMutableString *_name;
    NSMutableArray *_variables;
    NSMutableArray *_projectors;               // Mutable array of FEMProjector classes
    NSMutableArray *_next;
    NSMutableArray *_colors;                   // Mutable arrays for each color. Each array stores the number of elements for that color,
                                               // the color index and the corresponding RGB values for that color as: [number of elements, color index, R, G, B]
    FEMMesh *_parent;
    FEMMesh *_child;
    
    int * __nullable _colorMapping;
    int * __nullable _elementNodeIndexesStore;
    int * __nullable _discontinousPerm;
    int * __nullable _invPerm;
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
@property(nonatomic, assign, getter = isDiscontinuousMesh) BOOL discontinuousMesh;
@property(nonatomic, strong, nullable) NSMutableString *name;
@property(nonatomic, strong, nullable) NSMutableArray <FEMVariable *> *variables;
@property(nonatomic, strong, nonnull) NSMutableArray *projectors;
@property(nonatomic, strong, nonnull) NSMutableArray <FEMMesh *> *next;
@property(nonatomic, strong, nonnull) NSMutableArray <NSMutableArray *> *colors;
@property(nonatomic, strong, nullable) FEMMesh *parent;
@property(nonatomic, strong, nullable) FEMMesh *child;

-(void)loadMeshForModel:(FEMModel * __nonnull)model meshDirectory:(NSString * __nonnull)dir meshName:(NSString * __nonnull)name boundariesOnly:(BOOL)bd numberOfPartitions:(int * __nullable)numParts partitionID:(int * __nullable)partID definitions:(int * __nullable)defDofs;

// Allocations methods
-(void)allocatePDefinitionsForElement:(Element_t * __nonnull)element;

// Assign nodes
-(void)assignNodes:(Nodes_t * __nonnull)nodes;

// Assign elements
-(void)assignElements:(Element_t * __nonnull)elements;

// Assign faces
-(void)assignFaces:(Element_t * __nonnull)faces;

// Assign edges
-(void)assignEdges:(Element_t * __nonnull)edges;

// Assign quadrants
-(void)assignQuadrant:(Quadrant_t * __nonnull)quadrant;

// Assign view factors
-(void)assignViewFactors:(Factors_t * __nonnull)factors;

// Assign inverse permutation
-(void)assignInvPerm:(int * __nonnull)perm;

// Assign color mapping
-(void)assignColorMapping:(int * __nonnull)colorMap;

// Assign indexes store
-(void)assignElementNodeIndexesStore:(int * __nonnull)elementNodeIndexesStore;

// Nodes getter
-(Nodes_t * __nullable)getNodes;

// Elements getter
-(Element_t * __nullable)getElements;

// Edges getter 
-(Element_t * __nullable)getEdges;

// Faces getter
-(Element_t * __nullable)getFaces;

// Quadrant getter
-(Quadrant_t * __nullable)getQuadrant;

// View Factors getter
-(Factors_t * __nullable)getViewFactors;

// Color mapping getter
-(int * __nullable)getColorMapping;

// Element permutation store getter
-(int * __nullable)getElementNodeIndexesStore;

// Discontinous permutation getter
-(int * __nullable)getDiscontinousPerm;

// Inverse permutation getter
-(int * __nullable)getInvPerm;

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

@end
