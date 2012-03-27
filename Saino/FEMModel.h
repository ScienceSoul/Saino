//
//  FEMModel.h
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"

@interface FEMModel : NSObject {
    
    int dimension;
    int numberOfBodyForces;
    int numberOfBoundaryConditions;
    
    
    NSArray *bodies;     // Array of dictionaries
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
    
    NSArray *bodyForces; // Array of FEMBodyForce objects
    NSArray *boundaryConditions; // Array of FEMBoundaryCondition objects
    
    Variable_t *variables;
}

-(int)dimension;
-(int)numberOfBodyForces;
-(int)numberOfBoundaryConditions;
-(NSArray *)bodiesValuesForKey:(NSString *)key;
-(NSArray *)returnBodyForces;
-(NSArray *)returnBoundaryConditions;

-(void)setDimension:(int)n;
-(void)setNumberOfBodyForces:(int)n;
-(void)setNumberOfBoundaryConditions:(int)n;

-(Variable_t *)returnPointerToVariables;

@end
