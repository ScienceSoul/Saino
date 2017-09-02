//===----------------------------------------------------------------------===//
//  FEMModel.h
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
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
#import "Constructors.h"
#import "FEMVariable.h"
#import "FEMSimulation.h"
#import "FEMConstants.h"
#import "FEMBodyForce.h"
#import "FEMBoundaryCondition.h"
#import "FEMBoundary.h"
#import "FEMEquation.h"
#import "FEMInitialConditions.h"
#import "FEMMaterial.h"
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
    NSArray *_boundaryID;                        // Boundaries ID
    NSArray *_solutions;                         // Array of all solutions in the model
    NSMutableArray *_meshes;                     // Array of meshes used by the model
    NSArray *_bodies;
    NSArray *_bodyForces;
    NSArray *_boundaryConditions;
    NSArray *_boundaries;
    NSArray *_equations;
    NSArray *_initialConditions;
    NSArray *_materials;
    NSMutableArray *_variables;
    FEMSimulation *_simulation;
    FEMConstants *_constants;
    FileReader *_mdf;
    
    Element_t * _Nullable _elements;
    Element_t * _Nullable _currentElement;
    Nodes_t * _Nullable _nodes;
    modelArraysContainer * _Nonnull _containers;
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
@property(nonatomic, strong, nonnull) NSArray <NSNumber *> *boundaryID;
@property(nonatomic, strong, nullable) NSArray *solutions;
@property(nonatomic, strong, nonnull) NSMutableArray *meshes;
@property(nonatomic, strong, nullable) NSArray <NSDictionary *> *bodies;
@property(nonatomic, strong, nullable) NSArray <FEMBodyForce *> *bodyForces;
@property(nonatomic, strong, nullable) NSArray <FEMBoundaryCondition *> *boundaryConditions;
@property(nonatomic, strong, nullable) NSArray <FEMBoundary *> *boundaries;
@property(nonatomic, strong, nullable) NSArray <FEMEquation *> *equations;
@property(nonatomic, strong, nullable) NSArray <FEMInitialConditions *> *initialConditions;
@property(nonatomic, strong, nullable) NSArray <FEMMaterial *> *materials;
@property(nonatomic, strong, nonnull) NSMutableArray <FEMVariable *> *variables;
@property(nonatomic, strong, nonnull) FEMSimulation *simulation;
@property(nonatomic, strong, nonnull) FEMConstants *constants;
@property(nonatomic, strong, nullable) FileReader *mdf;

-(void)deallocation;

// The dummy arguments are placeholders for future arguments related to MPI support
-(void)loadModelName:(NSString * _Nonnull)name boundariesOnly:(BOOL)bd dummy:(int * _Nullable)d1 dummy:(int * _Nullable)d2;

// Elements getter
-(Element_t * _Nullable)getElements;
-(Element_t * _Nullable)getCurrentElement;

// Nodes getter
-(Nodes_t * _Nullable)getNodes;

// Elements and nodes setters
-(void)SetElements:(Element_t * _Nonnull)elements;
-(void)setNodes:(Nodes_t * _Nonnull)nodes;

-(modelArraysContainer * _Nonnull)getContainers;

@end
