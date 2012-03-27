//
//  FEMModel.m
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMModel.h"

@implementation FEMModel

-(int)dimension {
    
    return dimension;
}

-(int)numberOfBodyForces {
    
    return numberOfBodyForces;
}

-(int)numberOfBoundaryConditions {
    
    return numberOfBoundaryConditions;
}

-(NSArray *)bodiesValuesForKey:(NSString *)key {
    
    return [bodies valueForKey:key];
}

-(NSArray *)returnBodyForces {
    
    return bodyForces;
}

-(NSArray *)returnBoundaryConditions {
    
    return boundaryConditions;
}

-(void)setDimension:(int)n {
    
    dimension = n;
}

-(void)setNumberOfBodyForces:(int)n {
    
    numberOfBodyForces = n;
}

-(void)setNumberOfBoundaryConditions:(int)n {
    
    numberOfBoundaryConditions = n;
}

-(Variable_t *)returnPointerToVariables {
    
    return variables;
}

@end
