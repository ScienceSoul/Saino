//
//  FEMModel.h
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"
#import "FEMSimulation.h"
#import "FEMConstants.h"
#import "FileReader.h"

@interface FEMModel : NSObject {
    
    int _dimension;
    int _numberOfNodes;
    int _numberOfBulkElements;
    int _numberOfBoundaryElements;
    int _numberOfBodies;
    int _numberOfBodyForces;
    int _numberOfBoundaryConditions;
    int _numberOfBoundaries;
    int _numberOfInitialConditions;
    int _numberOfSolutions;
    int _numberOfEquations;
    int _numberOfMaterials;
    int _coordinates;
    int _totalMatrixElements;
    int _maxElementNodes;
    
    id __weak _mesh;
    id __weak _solution;
    NSMutableString *_meshDir;
    NSMutableString *_meshName;
    NSMutableString *_outputPath;
    NSArray *_boundaryID;                           // Array of NSNumbers for boundaries ID
    NSArray *_solutions;                            // Array of all solutions in the model
    NSMutableArray *_meshes;                        // Array of meshes used by the model
    NSArray *_bodies;                               // Array of dictionaries
    NSArray *_bodyForces;                           // Array of FEMBodyForce objects
    NSArray *_boundaryConditions;                   // Array of FEMBoundaryCondition objects
    NSArray *_boundaries;                           // Array of FEMBoundary objects
    NSArray *_equations;                            // Array of FEMEquation objects
    NSArray *_initialConditions;                    // Array of FEMInitialConditions objects
    NSArray *_materials;                            // Array of FEMMaterial objects
    NSMutableArray *_variables;                     // Mutable array of FEMVariable classes
    FEMSimulation *_simulation;
    FEMConstants *_constants;
    FileReader *_mdf;
    
    Element_t * __nullable _elements;
    Element_t * __nullable _currentElement;
    Nodes_t * __nullable _nodes;
    modelArraysContainer * __nonnull _containers;
}

@property(nonatomic, assign) int dimension;
@property(nonatomic, assign) int numberOfNodes;
@property(nonatomic, assign) int numberOfBulkElements;
@property(nonatomic, assign) int numberOfBoundaryElements;
@property(nonatomic, assign) int numberOfBodies;
@property(nonatomic, assign) int numberOfBodyForces;
@property(nonatomic, assign) int numberOfBoundaryConditions;
@property(nonatomic, assign) int numberOfBoundaries;
@property(nonatomic, assign) int numberOfInitialConditions;
@property(nonatomic, assign) int numberOfSolutions;
@property(nonatomic, assign) int numberOfEquations;
@property(nonatomic, assign) int numberOfMaterials;
@property(nonatomic, assign) int coordinates;
@property(nonatomic, assign) int totalMatrixElements;
@property(nonatomic, assign) int maxElementNodes;
@property(nonatomic, weak, nullable) id mesh;
@property(nonatomic, weak, nullable) id solution;
@property(nonatomic, strong, nonnull) NSMutableString *meshDir;
@property(nonatomic, strong, nonnull) NSMutableString *meshName;
@property(nonatomic, strong, nonnull) NSMutableString *outputPath;
@property(nonatomic, strong, nonnull) NSArray *boundaryID;
@property(nonatomic, strong, nullable) NSArray *solutions;
@property(nonatomic, strong, nonnull) NSMutableArray *meshes;
@property(nonatomic, strong, nullable) NSArray *bodies;
@property(nonatomic, strong, nullable) NSArray *bodyForces;
@property(nonatomic, strong, nullable) NSArray *boundaryConditions;
@property(nonatomic, strong, nullable) NSArray *boundaries;
@property(nonatomic, strong, nullable) NSArray *equations;
@property(nonatomic, strong, nullable) NSArray *initialConditions;
@property(nonatomic, strong, nullable) NSArray *materials;
@property(nonatomic, strong, nonnull) NSMutableArray *variables;
@property(nonatomic, strong, nonnull) FEMSimulation *simulation;
@property(nonatomic, strong, nonnull) FEMConstants *constants;
@property(nonatomic, strong, nullable) FileReader *mdf;

-(void)deallocation;

// The dummy arguments are placeholders for future arguments related to MPI support
-(void)loadModelName:(NSString * __nonnull)name boundariesOnly:(BOOL)bd dummy:(int * __nullable)d1 dummy:(int * __nullable)d2;

// Elements getter
-(Element_t * __nullable)getElements;
-(Element_t * __nullable)getCurrentElement;

// Nodes getter
-(Nodes_t * __nullable)getNodes;

// Elements and nodes setters
-(void)SetElements:(Element_t * __nonnull)elements;
-(void)setNodes:(Nodes_t * __nonnull)nodes;

-(modelArraysContainer * __nonnull)getContainers;

@end
