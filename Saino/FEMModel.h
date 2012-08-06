//
//  FEMModel.h
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMVariable.h"
#import "Constructors.h"

@interface FEMModel : NSObject {
    
    int _dimension;
    int _numberOfNodes;
    int _numberOfBulkElements;
    int _numberOfBoundaryElements;
    int _numberOfBodyForces;
    int _numberOfBoundaries;
    
    NSArray *_bodies;                             // Array of dictionaries
    NSArray *_bodyForces;                         // Array of FEMBodyForce objects
    NSArray *_boundaries;                         // Array of FEMBoundaryCondition objects
    NSArray *_simulations;                        // Array of FEMSimulation objects
    NSMutableDictionary *_variables;              // Mutable dictionary holding FEMVariable classes
    
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
@property(nonatomic, assign) int numberOfBodyForces;
@property(nonatomic, assign) int numberOfBoundaries;
@property(nonatomic, strong) NSArray *bodies;
@property(nonatomic, strong) NSArray *bodyForces;
@property(nonatomic, strong) NSArray *boundaries;
@property(nonatomic, strong) NSArray *simulations;
@property(nonatomic, strong) NSMutableDictionary *variables;

@end
