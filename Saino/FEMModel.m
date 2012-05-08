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

-(int)numberOfNodes {
    
    return numberOfNodes;
}

-(int)numberOfBulkElements {
    
    return numberOfBulkElements;
}

-(int)numberOfBoundaryElements {
    
    return numberOfBoundaryElements;
}

-(int)numberOfBodyForces {
    
    return numberOfBodyForces;
}

-(int)numberOfBoundaryConditions {
    
    return numberOfBoundaryConditions;
}

-(NSArray *)bodies {
    
    return bodies;
}

-(NSArray *)bodyForces {
    
    return bodyForces;
}

-(NSArray *)boundaryConditions {
    
    return boundaryConditions;
}

-(NSArray *)simulations {
    
    return simulations;
}

-(void)setDimension:(int)n {
    
    dimension = n;
}

-(void)setNumberOfNodes:(int)n {
    
    numberOfNodes = n;
}

-(void)setNumberOfBulkElements:(int)n {
    
    numberOfBulkElements = n;
}

-(void)setNumberOfBoundaryElements:(int)n {
    
    numberOfBoundaryElements = n;
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
