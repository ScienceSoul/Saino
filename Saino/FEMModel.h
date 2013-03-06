//
//  FEMModel.h
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"
#import "FEMVariable.h"
#import "FEMSimulation.h"

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
    
    id _mesh;
    id _solution;
    NSMutableString *_outputPath;
    NSArray *_boundaryID;                           // Array of NSNumbers for boundaries ID
    NSArray *_solutions;                            // Array of all solutions in the model
    NSMutableArray *_meshes;                        // Array of meshes used by the model
    NSArray *_bodies;                               // Array of dictionaries
    NSArray *_bodyForces;                           // Array of FEMBodyForce objects
    NSArray *_boundaryConditions;                   // Array of FEMBoundaryCondition objects
    NSArray *_boundaries;                           // Array of FEMBoundary objects
    NSArray *_equations;                            // Array of FEMEquation objects
    NSMutableArray *_variables;                     // Mutable array of FEMVariable classes
    FEMSimulation *_simulation;                   
    
    Element_t *_elements;
    Element_t *_currentElement;
    Nodes_t *_nodes;
    modelArraysContainer *_containers;
    
    /* Initialize bodies like this:
     NSArray* names = [NSArray arrayWithObjects:
     [NSDictionary dictionaryWithObjectsAndKeys:
     @"Joe",@"firstname",
     @"Bloggs",@"surname",
     nil],
     [NSDictionary dictionaryWithObjectsAndKeys:
     @"Simon",@"firstname",
     @"Templar",@"surname",
     nil],
     [NSDictionary dictionaryWithObjectsAndKeys:
     @"Amelia",@"firstname",
     @"Pond",@"surname",
     nil],
     nil];
     */
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
@property(nonatomic, strong) id mesh;
@property(nonatomic, strong) id solution;
@property(nonatomic, strong) NSMutableString *outputPath;
@property(nonatomic, strong) NSArray *boundaryID;
@property(nonatomic, strong) NSArray *solutions;
@property(nonatomic, strong) NSMutableArray *meshes;
@property(nonatomic, strong) NSArray *bodies;
@property(nonatomic, strong) NSArray *bodyForces;
@property(nonatomic, strong) NSArray *boundaryConditions;
@property(nonatomic, strong) NSArray *boundaries;
@property(nonatomic, strong) NSArray *equations;
@property(nonatomic, strong) FEMSimulation *simulation;
@property(nonatomic, strong) NSMutableArray *variables;

-(void)deallocation;

// The dummy arguments are placeholders for future arguments related to MPI support
-(void)loadModelName:(NSString *)name boundariesOnly:(BOOL)bd dummy:(int *)d1 dummy:(int *)d2;

// Elements getter
-(Element_t *)getElements;
-(Element_t *)getCurrentElement;

// Nodes getter
-(Nodes_t *)getNodes;

-(modelArraysContainer *)getContainers;

@end
