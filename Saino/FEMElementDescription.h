//
//  FEMElementDescription.h
//  Saino
//
//  Created by Hakime Seddik on 05/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <Accelerate/Accelerate.h>

#import "FEMElementsDefinition.h"
#import "FEMSolution.h"
#import "FEMNumericIntegration.h"

#import "Constructors.h"
#import "memory.h"

@interface FEMElementDescription : NSObject {
@private
    
    int maxDeg, maxDeg3, maxDeg2;
    
    int MAX_ELEMEMT_NODES;
        
    //------------------------------
    // Number of elements definition
    //------------------------------
    int numberOfElementDefs;
    
    //---------------------------------------
    // List of supported elements definition
    //---------------------------------------
    FEMElementsDefinition *_1nodePoint;
    FEMElementsDefinition *_2nodePeriodic;
    FEMElementsDefinition *_2nodeLine;
    FEMElementsDefinition *_3nodeLine;
    FEMElementsDefinition *_4nodeLine;
    FEMElementsDefinition *_3nodeTriangle;
    FEMElementsDefinition *_6nodeTriangle;
    FEMElementsDefinition *_10nodeTriangle;
    FEMElementsDefinition *_4nodeQuadrilateral;
    FEMElementsDefinition *_8nodeQuadrilateral;
    FEMElementsDefinition *_9nodeQuadrilateral;
    FEMElementsDefinition *_12nodeQuadrilateral;
    FEMElementsDefinition *_16nodeQuadrilateral;
    FEMElementsDefinition *_4nodeTetrahedron;
    FEMElementsDefinition *_10nodeTetrahedron;
    FEMElementsDefinition *_5nodePyramid;
    FEMElementsDefinition *_13nodePyramid;
    FEMElementsDefinition *_6nodeWedge;
    FEMElementsDefinition *_8nodeOctahedron;
    FEMElementsDefinition *_20nodeOctahedron;
    FEMElementsDefinition *_27nodeOctahedron;
    
    NSArray *listOfDefinitions;
    
    //-----------------------------------------------------
    // Linked list strurture which contains the elements 
    // definition and their allocations. The elements 
    // definition are provided by the class itself.
    //-----------------------------------------------------
    ElementType_t *elementTypeList;
    BOOL isTypeListInitialized;
    
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

-(void)deallocation;
-(void)addDescriptionOfElement:(ElementType_t)element withBasisTerms:(int *)terms;
-(void)initElementDescriptions;
-(int **)getEdgeMap:(int)elementFamily;
-(double)elementDiameter:(Element_t *)element: (Nodes_t *)nodes;
-(void)computeStabilizationParameter:(Element_t *)element: (Nodes_t *)nodes: (FEMMesh *)mesh: (int)n: (double)mk: (double *)hk;
-(ElementType_t *)getElementType:(int)code inMesh:(FEMMesh *)mesh stabilization:(BOOL *)computeStab;

@end
