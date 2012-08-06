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
    
    int dimension;
    int numberOfNodes, numberOfBulkElements, numberOfBoundaryElements;
    int numberOfBodyForces;
    int numberOfBoundaryConditions;
    
    
    NSArray *bodies;                            // Array of dictionaries
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
    
    NSArray *_bodyForces;                        // Array of FEMBodyForce objects
    NSArray *_boundaries;                // Array of FEMBoundaryCondition objects
    NSArray *_simulations;                        // Array of FEMSimulation objects
    NSMutableDictionary *_variables;             // Mutable dictionary holding FEMVariable classes
}

@property(nonatomic, strong) NSMutableDictionary *variables;
@property(nonatomic, strong) NSArray *bodyForces;
@property(nonatomic, strong) NSArray *boundaries;
@property(nonatomic, strong) NSArray *simulations;

-(int)dimension;
-(int)numberOfNodes;
-(int)numberOfBulkElements;
-(int)numberOfBoundaryElements;
-(int)numberOfBodyForces;
-(int)numberOfBoundaryConditions;
-(NSArray *)bodies;

-(void)setDimension:(int)n;
-(void)setNumberOfNodes:(int)n;
-(void)setNumberOfBulkElements:(int)n;
-(void)setNumberOfBoundaryElements:(int)n;
-(void)setNumberOfBodyForces:(int)n;
-(void)setNumberOfBoundaryConditions:(int)n;

@end
