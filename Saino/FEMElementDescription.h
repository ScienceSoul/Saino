//
//  FEMElementDescription.h
//  Saino
//
//  Created by Hakime Seddik on 05/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <Accelerate/Accelerate.h>

#import "FEMSolution.h"
#import "FEMNumericIntegration.h"
#import "Constructors.h"
#import "memory.h"

@interface FEMElementDescription : NSObject {
@private
    ElementType_t *elementTypeList;
    
    int **point;
    int **line;
    int **triangle;
    int **quad;
    int **tetra;
    int **prism;
    int **wedge;
    int **brick;
    
    BOOL initialized[8];
    
}

-(int **)getEdgeMap:(int)elementFamily;
-(double)elementDiameter:(Element_t *)element: (Nodes_t *)nodes;
-(void)computeStabilizationParameter:(Element_t *)element: (Nodes_t *)nodes: (FEMMesh *)mesh: (int)n: (double)mk: (double *)hk;
-(ElementType_t *)getElementType:(int)code inMesh:(FEMMesh *)mesh stabilization:(BOOL *)computeStab;

@end
