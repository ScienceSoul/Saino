//
//  FEMModel.h
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMVariable.h"
#import "FEMSimulation.h"
#import "Constructors.h"

@interface FEMModel : NSObject {
    
    int _dimension;
    int _numberOfNodes;
    int _numberOfBulkElements;
    int _numberOfBoundaryElements;
    int _numberOfBodies;
    int _numberOfBodyForces;
    int _numberOfBoundaryConditions;
    int _numberOfBoundaries;
    int _numberOfICs;
    int _numberOfSolutions;
    int _numberOfEquations;
    
    id _mesh;
    NSArray *_boundaryID;                           // Array of NSNumbers for boundaries ID
    NSArray *_solutions;                            // Array of all solutions in the model
    NSArray *_meshes;                               // Array of meshes used by the model
    NSArray *_bodies;                               // Array of dictionaries
    NSArray *_bodyForces;                           // Array of FEMBodyForce objects
    NSArray *_boundaryConditions;                   // Array of FEMBoundaryCondition objects
    NSArray *_boundaries;                           // Array of FEMBoundary objects
    NSArray *_equations;                            // Array of FEMEquation objects
    FEMSimulation *_simulation;                   
    NSMutableArray *_variables;                     // Mutable array of FEMVariable classes
    
    Element_t *_currentElement;
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
@property(nonatomic, assign) int numberOfICs;
@property(nonatomic, assign) int numberOfSolutions;
@property(nonatomic, assign) int numberOfEquations;
@property(nonatomic, strong) id mesh;
@property(nonatomic, strong) NSArray *boundaryID;
@property(nonatomic, strong) NSArray *solutions;
@property(nonatomic, strong) NSArray *meshes;
@property(nonatomic, strong) NSArray *bodies;
@property(nonatomic, strong) NSArray *bodyForces;
@property(nonatomic, strong) NSArray *boundaryConditions;
@property(nonatomic, strong) NSArray *boundaries;
@property(nonatomic, strong) NSArray *equations;
@property(nonatomic, strong) FEMSimulation *simulation;
@property(nonatomic, strong) NSMutableArray *variables;

-(void)deallocation;

// Elements getter
-(Element_t *)getCurrentElement;

-(modelArraysContainer *)getContainers;

@end
