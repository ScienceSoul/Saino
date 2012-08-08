//
//  FEMElementDescription.m
//  Saino
//
//  Created by Hakime Seddik on 05/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMElementDescription.h"
#import "FEMLinearAlgebra.h"

#import "Numerics.h"
#import "Utils.h"
#import "GaussIntegration.h"

static double AEPS = 10.0 * DBL_EPSILON;

@interface FEMElementDescription ()

-(void)FEMElementDescription_initElementsDefinition;
-(void)FEMElementDescription_compute1DPBasis:(double **)basis: (int)n;

@end

@implementation FEMElementDescription {
    
    int _maxDeg, _maxDeg3, _maxDeg2;
    int _maxElementNodes;
    
    //------------------------------
    // Number of elements definition
    //------------------------------
    int _numberOfElementDefs;
    
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
    
    NSArray *_listOfDefinitions;
    
    //-----------------------------------------------------
    // Linked list strurture which contains the elements
    // definition and their allocations. The elements
    // definition are provided by the class itself.
    //-----------------------------------------------------
    ElementType_t *_elementTypeList;
    BOOL _isTypeListInitialized;
    
    int **_point;
    int **_line;
    int **_triangle;
    int **_quad;
    int **_tetra;
    int **_prism;
    int **_wedge;
    int **_brick;
    
    BOOL _initialized[8];
}

#pragma mark Private methods...

-(void)FEMElementDescription_initElementsDefinition {
    
    //---------------------
    // 1 node point element
    //---------------------
    _1nodePoint = [[FEMElementsDefinition alloc] init];
    _1nodePoint.dimension = [NSNumber numberWithInt:1];
    _1nodePoint.topology = @"point";
    _1nodePoint.code = [NSNumber numberWithInt:101];
    _1nodePoint.nodes = [NSNumber numberWithInt:1];
    _1nodePoint.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], nil];
    _1nodePoint.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], nil];
    _1nodePoint.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:0], nil];
    _1nodePoint.stabilization = [NSNumber numberWithDouble:0.0];
    
    //-------------------------
    // 2 nodes periodic element
    //-------------------------
    _2nodePeriodic = [[FEMElementsDefinition alloc] init];
    _2nodePeriodic.dimension = [NSNumber numberWithInt:1];
    _2nodePeriodic.topology = @"point";
    _2nodePeriodic.code = [NSNumber numberWithInt:102];
    _2nodePeriodic.nodes = [NSNumber numberWithInt:2];
    _2nodePeriodic.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil];
    _2nodePeriodic.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], nil];
    _2nodePeriodic.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:0], [NSNumber numberWithInt:0], nil];
    _2nodePeriodic.stabilization = [NSNumber numberWithDouble:0.0];
    
    //---------------------
    // 2 nodes line element
    //---------------------
    _2nodeLine = [[FEMElementsDefinition alloc] init];
    _2nodeLine.dimension = [NSNumber numberWithInt:1];
    _2nodeLine.topology = @"line";
    _2nodeLine.code = [NSNumber numberWithInt:202];
    _2nodeLine.nodes = [NSNumber numberWithInt:2];
    //          1    2 
    _2nodeLine.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], nil];
    //          1  2   3   4
    //          1  u  u^2 u^3
    _2nodeLine.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], nil];
    _2nodeLine.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:2], [NSNumber numberWithInt:4], nil];
    _2nodeLine.stabilization = [NSNumber numberWithDouble:0.33333333333333333333333];
    
    //---------------------
    // 3 nodes line element
    //---------------------
    _3nodeLine = [[FEMElementsDefinition alloc] init];
    _3nodeLine.dimension = [NSNumber numberWithInt:1];
    _3nodeLine.topology = @"line";
    _3nodeLine.code = [NSNumber numberWithInt:203];
    _3nodeLine.nodes = [NSNumber numberWithInt:3];
    //          1     2     3
    _3nodeLine.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil];
    //          1  2   3   4
    //          1  u  u^2 u^3
    _3nodeLine.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], nil];
    _3nodeLine.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:3], [NSNumber numberWithInt:5], nil];
    _3nodeLine.stabilization = [NSNumber numberWithDouble:0.1666666666666666666666];
    
    //---------------------
    // 4 nodes line element
    //---------------------
    _4nodeLine = [[FEMElementsDefinition alloc] init];
    _4nodeLine.dimension = [NSNumber numberWithInt:1];
    _4nodeLine.topology = @"line";
    _4nodeLine.code = [NSNumber numberWithInt:204];
    _4nodeLine.nodes = [NSNumber numberWithInt:4];
    //          1     2     3     
    _4nodeLine.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-0.333333333333333333],
                        [NSNumber numberWithDouble:0.3333333333333333333], nil];
    //          1  2   3   4
    //          1  u  u^2 u^3
    _4nodeLine.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], [NSNumber numberWithInt:4], nil];
    _4nodeLine.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:4], [NSNumber numberWithInt:6], nil];
    _4nodeLine.stabilization = [NSNumber numberWithDouble:0.03333333333333333333333];
    
    //-------------------------
    // 3 nodes triangle element
    //-------------------------
    _3nodeTriangle = [[FEMElementsDefinition alloc] init];
    _3nodeTriangle.dimension = [NSNumber numberWithInt:2];
    _3nodeTriangle.topology = @"triangle";
    _3nodeTriangle.code = [NSNumber numberWithInt:303];
    _3nodeTriangle.nodes = [NSNumber numberWithInt:3];
    //          1     2     3
    _3nodeTriangle.nodeU = [NSArray  arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil];
    _3nodeTriangle.nodeV = [NSArray  arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _3nodeTriangle.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:5], nil];
    _3nodeTriangle.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:3], [NSNumber numberWithInt:20], nil];
    _3nodeTriangle.stabilization = [NSNumber numberWithDouble:0.333333333333333333333];
    
    //-------------------------
    // 6 nodes triangle element
    //-------------------------
    _6nodeTriangle = [[FEMElementsDefinition alloc] init];
    _6nodeTriangle.dimension = [NSNumber numberWithInt:2];
    _6nodeTriangle.topology = @"triangle";
    _6nodeTriangle.code = [NSNumber numberWithInt:306];
    _6nodeTriangle.nodes = [NSNumber numberWithInt:6];
    //          1     2     3     4     5     6
    _6nodeTriangle.nodeU = [NSArray  arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], 
                            [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.0], nil];
    _6nodeTriangle.nodeV = [NSArray  arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], 
                            [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.5], nil];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _6nodeTriangle.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], 
                            [NSNumber numberWithInt:5], [NSNumber numberWithInt:6], [NSNumber numberWithInt:9], nil];
    _6nodeTriangle.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:6], [NSNumber numberWithInt:17], nil];
    _6nodeTriangle.stabilization = [NSNumber numberWithDouble:0.041666666666666666];
    
    //--------------------------
    // 10 nodes triangle element
    //--------------------------
    _10nodeTriangle = [[FEMElementsDefinition alloc] init];
    _10nodeTriangle.dimension = [NSNumber numberWithInt:2];
    _10nodeTriangle.topology = @"triangle";
    _10nodeTriangle.code = [NSNumber numberWithInt:310];
    _10nodeTriangle.nodes = [NSNumber numberWithInt:10];
    //          1     2     3     ......
    _10nodeTriangle.nodeU = [NSArray  arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], 
                             [NSNumber numberWithDouble:0.333333333333333333], [NSNumber numberWithDouble:0.666666666666666667], 
                             [NSNumber numberWithDouble:0.666666666666666667], [NSNumber numberWithDouble:0.333333333333333333], 
                             [NSNumber numberWithDouble:0.000000000000000000], [NSNumber numberWithDouble:0.000000000000000000], 
                             [NSNumber numberWithDouble:0.333333333333333333], nil];
    _10nodeTriangle.nodeV = [NSArray  arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], 
                             [NSNumber numberWithDouble:0.000000000000000000], [NSNumber numberWithDouble:0.000000000000000000], 
                             [NSNumber numberWithDouble:0.333333333333333333], [NSNumber numberWithDouble:0.666666666666666667], 
                             [NSNumber numberWithDouble:0.666666666666666667], [NSNumber numberWithDouble:0.333333333333333333], 
                             [NSNumber numberWithDouble:0.333333333333333333], nil];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _10nodeTriangle.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], 
                             [NSNumber numberWithInt:4], [NSNumber numberWithInt:5], [NSNumber numberWithInt:6], [NSNumber numberWithInt:7], 
                             [NSNumber numberWithInt:9], [NSNumber numberWithInt:10], [NSNumber numberWithInt:13], nil];
    _10nodeTriangle.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:25], [NSNumber numberWithInt:64], nil];
    _10nodeTriangle.stabilization = [NSNumber numberWithDouble:0.01341555597798937329];
       
    //------------------------------
    // 4 nodes quadrilateral element
    //------------------------------
    _4nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _4nodeQuadrilateral.dimension = [NSNumber numberWithInt:2];
    _4nodeQuadrilateral.topology = @"quadrilateral";
    _4nodeQuadrilateral.code = [NSNumber numberWithInt:404];
    _4nodeQuadrilateral.nodes = [NSNumber numberWithInt:4];
    //          1     2     3     4
    _4nodeQuadrilateral.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0],
                                 [NSNumber numberWithDouble:-1.0], nil];
    _4nodeQuadrilateral.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0],
                                 [NSNumber numberWithDouble:1.0], nil];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _4nodeQuadrilateral.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:5],
                                 [NSNumber numberWithInt:6], nil];
    _4nodeQuadrilateral.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:4], [NSNumber numberWithInt:16], nil];
    _4nodeQuadrilateral.stabilization = [NSNumber numberWithDouble:0.333333333333333333333];
    
    //------------------------------
    // 8 nodes quadrilateral element
    //------------------------------
    _8nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _8nodeQuadrilateral.dimension = [NSNumber numberWithInt:2];
    _8nodeQuadrilateral.topology = @"quadrilateral";
    _8nodeQuadrilateral.code = [NSNumber numberWithInt:408];
    _8nodeQuadrilateral.nodes = [NSNumber numberWithInt:8];
    //          1     2     3     4     5     6      7      8
    _8nodeQuadrilateral.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                                 [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0],
                                 [NSNumber numberWithDouble:-1.0], nil];
    _8nodeQuadrilateral.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                                 [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0],
                                 [NSNumber numberWithDouble:0.0], nil];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _8nodeQuadrilateral.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], [NSNumber numberWithInt:5],
                                 [NSNumber numberWithInt:6], [NSNumber numberWithInt:7], [NSNumber numberWithInt:9], [NSNumber numberWithInt:10], nil];
    _8nodeQuadrilateral.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:9], [NSNumber numberWithInt:25], nil];
    _8nodeQuadrilateral.stabilization = [NSNumber numberWithDouble:0.0833333333333333333333];
    
    //------------------------------
    // 9 nodes quadrilateral element
    //------------------------------
    _9nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _9nodeQuadrilateral.dimension = [NSNumber numberWithInt:2];
    _9nodeQuadrilateral.topology = @"quadrilateral";
    _9nodeQuadrilateral.code = [NSNumber numberWithInt:409];
    _9nodeQuadrilateral.nodes = [NSNumber numberWithInt:9];
    //          1     2     3     4     5     6      7      8      9
    _9nodeQuadrilateral.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                                 [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0],
                                 [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], nil];
    _9nodeQuadrilateral.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                                 [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0],
                                 [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], nil];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _9nodeQuadrilateral.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], [NSNumber numberWithInt:5],
                                 [NSNumber numberWithInt:6], [NSNumber numberWithInt:7], [NSNumber numberWithInt:9], [NSNumber numberWithInt:10], 
                                 [NSNumber numberWithInt:11], nil];
    _9nodeQuadrilateral.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:9], [NSNumber numberWithInt:25], nil];
    _9nodeQuadrilateral.stabilization = [NSNumber numberWithDouble:0.0833333333333333333333];
    
    //-------------------------------
    // 12 nodes quadrilateral element
    //-------------------------------
    _12nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _12nodeQuadrilateral.dimension = [NSNumber numberWithInt:2];
    _12nodeQuadrilateral.topology = @"quadrilateral";
    _12nodeQuadrilateral.code = [NSNumber numberWithInt:412];
    _12nodeQuadrilateral.nodes = [NSNumber numberWithInt:12];
    //          1     2     3     ......
    _12nodeQuadrilateral.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                                  [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-0.3333333333333333333], [NSNumber numberWithDouble:0.3333333333333333333],
                                  [NSNumber numberWithDouble:1.0000000000000000000], [NSNumber numberWithDouble:1.0000000000000000000], 
                                  [NSNumber numberWithDouble:0.3333333333333333333], [NSNumber numberWithDouble:-0.3333333333333333333], 
                                  [NSNumber numberWithDouble:-1.0000000000000000000], [NSNumber numberWithDouble:-1.0000000000000000000], nil];
    _12nodeQuadrilateral.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                                  [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0000000000000000000], [NSNumber numberWithDouble:-1.0000000000000000000],
                                  [NSNumber numberWithDouble:-0.3333333333333333333], [NSNumber numberWithDouble:0.3333333333333333333], 
                                  [NSNumber numberWithDouble:1.0000000000000000000], [NSNumber numberWithDouble:1.0000000000000000000], 
                                  [NSNumber numberWithDouble:0.3333333333333333333], [NSNumber numberWithDouble:-0.3333333333333333333], nil];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _12nodeQuadrilateral.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], [NSNumber numberWithInt:4],
                                  [NSNumber numberWithInt:5], [NSNumber numberWithInt:6], [NSNumber numberWithInt:7], [NSNumber numberWithInt:8], 
                                  [NSNumber numberWithInt:9], [NSNumber numberWithInt:10], [NSNumber numberWithInt:13], [NSNumber numberWithInt:14], nil];
    _12nodeQuadrilateral.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:16], [NSNumber numberWithInt:36], nil];
    _12nodeQuadrilateral.stabilization = [NSNumber numberWithDouble:0.0];
    
    //-------------------------------
    // 16 nodes quadrilateral element
    //-------------------------------
    _16nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _16nodeQuadrilateral.dimension = [NSNumber numberWithInt:2];
    _16nodeQuadrilateral.topology = @"quadrilateral";
    _16nodeQuadrilateral.code = [NSNumber numberWithInt:416];
    _16nodeQuadrilateral.nodes = [NSNumber numberWithInt:16];
    //          1     2     3     ......
    _16nodeQuadrilateral.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                                  [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-0.3333333333333333333], [NSNumber numberWithDouble:0.3333333333333333333],
                                  [NSNumber numberWithDouble:1.0000000000000000000], [NSNumber numberWithDouble:1.0000000000000000000], 
                                  [NSNumber numberWithDouble:0.3333333333333333333], [NSNumber numberWithDouble:-0.3333333333333333333], 
                                  [NSNumber numberWithDouble:-1.0000000000000000000], [NSNumber numberWithDouble:-1.0000000000000000000], 
                                  [NSNumber numberWithDouble:-0.3333333333333333333], [NSNumber numberWithDouble:0.3333333333333333333], 
                                  [NSNumber numberWithDouble:0.3333333333333333333], [NSNumber numberWithDouble:-0.3333333333333333333], nil];
    _16nodeQuadrilateral.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                                  [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0000000000000000000], [NSNumber numberWithDouble:-1.0000000000000000000],
                                  [NSNumber numberWithDouble:-0.3333333333333333333], [NSNumber numberWithDouble:0.3333333333333333333], 
                                  [NSNumber numberWithDouble:1.0000000000000000000], [NSNumber numberWithDouble:1.0000000000000000000], 
                                  [NSNumber numberWithDouble:0.3333333333333333333], [NSNumber numberWithDouble:-0.3333333333333333333], 
                                  [NSNumber numberWithDouble:-0.3333333333333333333], [NSNumber numberWithDouble:-0.3333333333333333333], 
                                  [NSNumber numberWithDouble:0.3333333333333333333], [NSNumber numberWithDouble:0.3333333333333333333], nil];
    _16nodeQuadrilateral.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], [NSNumber numberWithInt:4],
                                  [NSNumber numberWithInt:5], [NSNumber numberWithInt:6], [NSNumber numberWithInt:7], [NSNumber numberWithInt:8], 
                                  [NSNumber numberWithInt:9], [NSNumber numberWithInt:10], [NSNumber numberWithInt:11], [NSNumber numberWithInt:12], 
                                  [NSNumber numberWithInt:13], [NSNumber numberWithInt:14], [NSNumber numberWithInt:15], [NSNumber numberWithInt:16], nil];
    _16nodeQuadrilateral.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:16], [NSNumber numberWithInt:36], nil];
    _16nodeQuadrilateral.stabilization = [NSNumber numberWithDouble:0.01766875890919188522];
    
    //----------------------------
    // 4 nodes tetrahedron element
    //----------------------------
    _4nodeTetrahedron = [[FEMElementsDefinition alloc] init];
    _4nodeTetrahedron.dimension = [NSNumber numberWithInt:3];
    _4nodeTetrahedron.topology = @"tetrahedron";
    _4nodeTetrahedron.code = [NSNumber numberWithInt:504];
    _4nodeTetrahedron.nodes = [NSNumber numberWithInt:4];
    //          1     2     3     4
    _4nodeTetrahedron.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], 
                               [NSNumber numberWithDouble:0.0], nil];
    _4nodeTetrahedron.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:0.0], nil];
    _4nodeTetrahedron.nodeW = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], 
                               [NSNumber numberWithDouble:1.0], nil];
    // 1       2     3      4     5       6        7       8         9
    // 1       u    u^2    u^3    v       uv      u^2v   u^3v       v^2 
    //
    // 10      11   12      13    14       15      16       17        18 
    // uv^2 u^2v^2 u^3v^2   v^3  uv^3     u^2v^3 u^3v^3     w         uw
    //
    // 19     20     21     22    23      24       25      26        27
    // u^2w  u^3w    vw    uvw   u^2vw   u^3vw    v^2w    uv^2w     u^2v^2w
    //
    //   28    
    // u^3v^2w v^3w  uv^3w u^2v^3w u^3v^3w ...
    _4nodeTetrahedron.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:5], [NSNumber numberWithInt:17], nil];
    _4nodeTetrahedron.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:4], [NSNumber numberWithInt:11], nil];
    _4nodeTetrahedron.stabilization = [NSNumber numberWithDouble:0.333333333333333333];
    
    //-----------------------------
    // 10 nodes tetrahedron element
    //-----------------------------
    _10nodeTetrahedron = [[FEMElementsDefinition alloc] init];
    _10nodeTetrahedron.dimension = [NSNumber numberWithInt:3];
    _10nodeTetrahedron.topology = @"tetrahedron";
    _10nodeTetrahedron.code = [NSNumber numberWithInt:510];
    _10nodeTetrahedron.nodes = [NSNumber numberWithInt:10];
    //          1     2     3     4     5     6      7      8      9     10
    _10nodeTetrahedron.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], 
                                [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.0], 
                                [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.0], nil];
    _10nodeTetrahedron.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], 
                                [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.5], 
                                [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.5], nil];
    _10nodeTetrahedron.nodeW = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], 
                                [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], 
                                [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.5], nil];
    // 1       2     3      4     5       6        7       8         9
    // 1       u    u^2    u^3    v       uv      u^2v   u^3v       v^2 
    //
    // 10      11   12      13    14       15      16       17        18 
    // uv^2 u^2v^2 u^3v^2   v^3  uv^3     u^2v^3 u^3v^3     w         uw
    //
    // 19     20     21     22    23      24       25      26        27
    // u^2w  u^3w    vw    uvw   u^2vw   u^3vw    v^2w    uv^2w     u^2v^2w
    //
    //   28    
    // u^3v^2w v^3w  uv^3w u^2v^3w u^3v^3w ...
    _10nodeTetrahedron.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], [NSNumber numberWithInt:5], 
                                [NSNumber numberWithInt:6], [NSNumber numberWithInt:9], [NSNumber numberWithInt:17], [NSNumber numberWithInt:18], 
                                [NSNumber numberWithInt:21], [NSNumber numberWithInt:33], nil];
    _10nodeTetrahedron.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:27], [NSNumber numberWithInt:64], nil];
    _10nodeTetrahedron.stabilization = [NSNumber numberWithDouble:0.0416666666666666666];
    

    
    //------------------------
    // 5 nodes pyramid element
    //------------------------
    _5nodePyramid = [[FEMElementsDefinition alloc] init];
    _5nodePyramid.dimension = [NSNumber numberWithInt:3];
    _5nodePyramid.topology = @"pyramid";
    _5nodePyramid.code = [NSNumber numberWithInt:605];
    _5nodePyramid.nodes = [NSNumber numberWithInt:5];
    //          1     2     3     4     5
    _5nodePyramid.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                           [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], nil];
    _5nodePyramid.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                           [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil];
    _5nodePyramid.nodeW = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], 
                           [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil];
    // 1       2     3      4     5       6        7       8         9
    // 1       u    u^2    u^3    v       uv      u^2v   u^3v       v^2 
    //
    // 10      11   12      13    14       15      16       17        18 
    // uv^2 u^2v^2 u^3v^2   v^3  uv^3     u^2v^3 u^3v^3     w         uw
    //
    // 19     20     21     22    23      24       25      26        27
    // u^2w  u^3w    vw    uvw   u^2vw   u^3vw    v^2w    uv^2w     u^2v^2w
    //
    //   28    
    // u^3v^2w v^3w  uv^3w u^2v^3w u^3v^3w ...
    _5nodePyramid.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:5], [NSNumber numberWithInt:6], 
                           [NSNumber numberWithInt:17], nil];
    _5nodePyramid.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:8], [NSNumber numberWithInt:64], nil];
    _5nodePyramid.stabilization = [NSNumber numberWithDouble:0.333333333333333333];
    
    //-------------------------
    // 13 nodes pyramid element
    //-------------------------
    _13nodePyramid = [[FEMElementsDefinition alloc] init];
    _13nodePyramid.dimension = [NSNumber numberWithInt:3];
    _13nodePyramid.topology = @"pyramid";
    _13nodePyramid.code = [NSNumber numberWithInt:613];
    _13nodePyramid.nodes = [NSNumber numberWithInt:13];
    //          1     2     3     4     5     6      7      8      9     10    11    12    13
    _13nodePyramid.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                            [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], 
                            [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-0.5], [NSNumber numberWithDouble:0.5], 
                            [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:-0.5], nil];
    _13nodePyramid.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                            [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], 
                            [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:-0.5], [NSNumber numberWithDouble:-0.5], 
                            [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.5], nil];
    _13nodePyramid.nodeW = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], 
                            [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], 
                            [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.5], 
                            [NSNumber numberWithDouble:0.5], [NSNumber numberWithDouble:0.5], nil];
    // 1       2     3      4     5       6        7       8         9
    // 1       u    u^2    u^3    v       uv      u^2v   u^3v       v^2 
    //
    // 10      11   12      13    14       15      16       17        18 
    // uv^2 u^2v^2 u^3v^2   v^3  uv^3     u^2v^3 u^3v^3     w         uw
    //
    // 19     20     21     22    23      24       25      26        27
    // u^2w  u^3w    vw    uvw   u^2vw   u^3vw    v^2w    uv^2w     u^2v^2w
    //
    //   28    
    // u^3v^2w v^3w  uv^3w u^2v^3w u^3v^3w ...
    _13nodePyramid.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], [NSNumber numberWithInt:5], 
                            [NSNumber numberWithInt:6], [NSNumber numberWithInt:7], [NSNumber numberWithInt:9], [NSNumber numberWithInt:10], [NSNumber numberWithInt:17],
                            [NSNumber numberWithInt:18], [NSNumber numberWithInt:21], [NSNumber numberWithInt:22], [NSNumber numberWithInt:33], nil];
    _13nodePyramid.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:27], [NSNumber numberWithInt:125], nil];
    _13nodePyramid.stabilization = [NSNumber numberWithDouble:0.333333333333333333];
    
    //----------------------
    // 6 nodes wedge element
    //----------------------
    _6nodeWedge = [[FEMElementsDefinition alloc] init];
    _6nodeWedge.dimension = [NSNumber numberWithInt:3];
    _6nodeWedge.topology = @"wedge";
    _6nodeWedge.code = [NSNumber numberWithInt:706];
    _6nodeWedge.nodes = [NSNumber numberWithInt:6];
    //          1     2     3     4     5     6
    _6nodeWedge.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], 
                         [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil];
    _6nodeWedge.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], 
                         [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], nil];
    _6nodeWedge.nodeW = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], 
                         [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], nil];
    // 1       2     3      4     5       6        7       8         9
    // 1       u    u^2    u^3    v       uv      u^2v   u^3v       v^2 
    //
    // 10      11   12      13    14       15      16       17        18 
    // uv^2 u^2v^2 u^3v^2   v^3  uv^3     u^2v^3 u^3v^3     w         uw
    //
    // 19     20     21     22    23      24       25      26        27
    // u^2w  u^3w    vw    uvw   u^2vw   u^3vw    v^2w    uv^2w     u^2v^2w
    //
    //   28    
    // u^3v^2w v^3w  uv^3w u^2v^3w u^3v^3w ...
    _6nodeWedge.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:5], [NSNumber numberWithInt:17], 
                         [NSNumber numberWithInt:18], [NSNumber numberWithInt:21], nil];
    _6nodeWedge.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:8], [NSNumber numberWithInt:64], nil];
    _6nodeWedge.stabilization = [NSNumber numberWithDouble:0.333333333333333333];
    
    //---------------------------
    // 8 nodes octahedron element
    //---------------------------
    _8nodeOctahedron = [[FEMElementsDefinition alloc] init];
    _8nodeOctahedron.dimension = [NSNumber numberWithInt:3];
    _8nodeOctahedron.topology = @"brick";
    _8nodeOctahedron.code = [NSNumber numberWithInt:808];
    _8nodeOctahedron.nodes = [NSNumber numberWithInt:8];
    //          1     2     3     4     5     6     7     8
    _8nodeOctahedron.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                              [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                              [NSNumber numberWithDouble:-1.0], nil];
    _8nodeOctahedron.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                              [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                              [NSNumber numberWithDouble:1.0], nil];
    _8nodeOctahedron.nodeW = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], 
                              [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                              [NSNumber numberWithDouble:1.0], nil];
    // 1       2     3      4     5       6        7       8         9
    // 1       u    u^2    u^3    v       uv      u^2v   u^3v       v^2 
    //
    // 10      11   12      13    14       15      16       17        18 
    // uv^2 u^2v^2 u^3v^2   v^3  uv^3     u^2v^3 u^3v^3     w         uw
    //
    // 19     20     21     22    23      24       25      26        27
    // u^2w  u^3w    vw    uvw   u^2vw   u^3vw    v^2w    uv^2w     u^2v^2w
    //
    //   28    
    // u^3v^2w v^3w  uv^3w u^2v^3w u^3v^3w  w^2 uw^2 u^2w^2 u^3w^3 
    //
    //
    // vw^2 uvw^2 u^2vw^2 u^3vw^2 v^2w^2 uv^2w^2 u^2v^2w^2 u^3v^2w^2
    _8nodeOctahedron.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:5], [NSNumber numberWithInt:6],
                              [NSNumber numberWithInt:17], [NSNumber numberWithInt:18], [NSNumber numberWithInt:21], [NSNumber numberWithInt:22], nil];
    _8nodeOctahedron.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:8], [NSNumber numberWithInt:64], nil];
    _8nodeOctahedron.stabilization = [NSNumber numberWithDouble:0.166666666666666666667];
    
    //----------------------------
    // 20 nodes octahedron element
    //----------------------------
    _20nodeOctahedron = [[FEMElementsDefinition alloc] init];
    _20nodeOctahedron.dimension = [NSNumber numberWithInt:3];
    _20nodeOctahedron.topology = @"brick";
    _20nodeOctahedron.code = [NSNumber numberWithInt:820];
    _20nodeOctahedron.nodes = [NSNumber numberWithInt:20];
    //          1     2     3     4     5     6      7      8      9     10    11    12    13    14    15   16   17    18    19    20
    _20nodeOctahedron.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0],
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], 
                               [NSNumber numberWithDouble:-1.0], nil];
    _20nodeOctahedron.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:0.0], nil];
    _20nodeOctahedron.nodeW = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], 
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0],
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0],
                               [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:1.0], nil];
    // 1       2     3      4     5       6        7       8         9
    // 1       u    u^2    u^3    v       uv      u^2v   u^3v       v^2 
    //
    // 10      11   12      13    14       15      16       17        18 
    // uv^2 u^2v^2 u^3v^2   v^3  uv^3     u^2v^3 u^3v^3     w         uw
    //
    // 19     20     21     22    23      24       25      26        27
    // u^2w  u^3w    vw    uvw   u^2vw   u^3vw    v^2w    uv^2w     u^2v^2w
    //
    //   28    
    // u^3v^2w v^3w  uv^3w u^2v^3w u^3v^3w  w^2 uw^2 u^2w^2 u^3w^3 
    //
    //
    // vw^2 uvw^2 u^2vw^2 u^3vw^2 v^2w^2 uv^2w^2 u^2v^2w^2 u^3v^2w^2
    _20nodeOctahedron.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], [NSNumber numberWithInt:5],
                               [NSNumber numberWithInt:6], [NSNumber numberWithInt:7], [NSNumber numberWithInt:9], [NSNumber numberWithInt:10], [NSNumber numberWithInt:17],
                               [NSNumber numberWithInt:18], [NSNumber numberWithInt:19], [NSNumber numberWithInt:21], [NSNumber numberWithInt:22], 
                               [NSNumber numberWithInt:23], [NSNumber numberWithInt:25], [NSNumber numberWithInt:26], [NSNumber numberWithInt:33], 
                               [NSNumber numberWithInt:34], [NSNumber numberWithInt:37], [NSNumber numberWithInt:38], nil];
    _20nodeOctahedron.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:27], [NSNumber numberWithInt:125], nil];
    _20nodeOctahedron.stabilization = [NSNumber numberWithDouble:0.08148148148148];
    
    //----------------------------
    // 27 nodes octahedron element
    //----------------------------
    _27nodeOctahedron = [[FEMElementsDefinition alloc] init];
    _27nodeOctahedron.dimension = [NSNumber numberWithInt:3];
    _27nodeOctahedron.topology = @"brick";
    _27nodeOctahedron.code = [NSNumber numberWithInt:827];
    _27nodeOctahedron.nodes = [NSNumber numberWithInt:27];
    //          1     2     3     4     5     6      7      8      9     10    11    12    13    14    15   16   17    18    19    20    21    22    23    24    25
    //          26    27
    _27nodeOctahedron.nodeU = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0],
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], 
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], 
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], nil];
    _27nodeOctahedron.nodeV = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0],
                               [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], nil];
    _27nodeOctahedron.nodeW = [NSArray arrayWithObjects:[NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], 
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:-1.0], 
                               [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0],
                               [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:1.0], 
                               [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:0.0],
                               [NSNumber numberWithDouble:0.0], [NSNumber numberWithDouble:-1.0], [NSNumber numberWithDouble:1.0], [NSNumber numberWithDouble:0.0], nil];
    // 1       2     3      4     5       6        7       8         9
    // 1       u    u^2    u^3    v       uv      u^2v   u^3v       v^2 
    //
    // 10      11   12      13    14       15      16       17        18 
    // uv^2 u^2v^2 u^3v^2   v^3  uv^3     u^2v^3 u^3v^3     w         uw
    //
    // 19     20     21     22    23      24       25      26        27
    // u^2w  u^3w    vw    uvw   u^2vw   u^3vw    v^2w    uv^2w     u^2v^2w
    //
    //   28    
    // u^3v^2w v^3w  uv^3w u^2v^3w u^3v^3w  w^2 uw^2 u^2w^2 u^3w^3 
    //
    //
    // vw^2 uvw^2 u^2vw^2 u^3vw^2 v^2w^2 uv^2w^2 u^2v^2w^2 u^3v^2w^2
    _27nodeOctahedron.basis = [NSArray arrayWithObjects:[NSNumber numberWithInt:1], [NSNumber numberWithInt:2], [NSNumber numberWithInt:3], [NSNumber numberWithInt:5], 
                               [NSNumber numberWithInt:6], [NSNumber numberWithInt:7], [NSNumber numberWithInt:9], [NSNumber numberWithInt:10], 
                               [NSNumber numberWithInt:11], [NSNumber numberWithInt:17], [NSNumber numberWithInt:18], [NSNumber numberWithInt:19], 
                               [NSNumber numberWithInt:21], [NSNumber numberWithInt:22], [NSNumber numberWithInt:23], [NSNumber numberWithInt:25], 
                               [NSNumber numberWithInt:26], [NSNumber numberWithInt:27], [NSNumber numberWithInt:33], [NSNumber numberWithInt:34], 
                               [NSNumber numberWithInt:35], [NSNumber numberWithInt:37], [NSNumber numberWithInt:38], [NSNumber numberWithInt:39], 
                               [NSNumber numberWithInt:41], [NSNumber numberWithInt:42], [NSNumber numberWithInt:43], nil];
    _27nodeOctahedron.gaussPoints = [NSArray arrayWithObjects:[NSNumber numberWithInt:27], [NSNumber numberWithInt:125], nil];
    _27nodeOctahedron.stabilization = [NSNumber numberWithDouble:0.08148148148148];
    
}

-(void)FEMElementDescription_compute1DPBasis:(double **)basis :(int)n {
/***************************************************************************
    Method to compute the 1D P-basis from Legendre polynomials
***************************************************************************/
 
    int i, j, k, l;
    double s, *p, *q, *p0, *p1;
    
    p = doublevec(1, n+1);
    q = doublevec(1, n);
    p0 = doublevec(1, n);
    p1 = doublevec(1, n+1);
    
    if (n <= 1) {
        basis[0][0] = 1.0;
        return;
    }
    
    //--------------------------------------------------------------------------
    // Compute coefficients of nth Legendre polynomial from the recurrence:
    // (i+1)P_{i+1}(x) = (2i+1)*x*P_i(x) - i*P{i-1}(x), P_{0} = 1; P_{1} = x;
    // IMPORTANT: Computed coefficients inaccurate for n > ~15
    //--------------------------------------------------------------------------
    memset( p, 0.0, ((n+1)*sizeof(p)) );
    memset( p0, 0.0, (n*sizeof(p0)) );
    memset( p1, 0.0, ((n+1)*sizeof(p1)) );
    p0[1] = 1;
    p1[1] = 1;
    p1[2] = 0;
    
    basis[0][0] = 0.5;
    basis[0][1] = -0.5;
    
    basis[1][0] = 0.5;
    basis[1][1] = 0.5;
    
    for (k=2; k<=n; k++) {
        if (k > 2) {
            s = sqrt( (2.0*(k-1)-1) / 2.0 );
            for (j=1; j<=k-1; j++) {
                basis[k-1][(k-1)-j+1] = s * p0[j] / (k-j);
                basis[k-1][0] = basis[k-1][0] - s * p0[j]*pow((-1), (j+1)) / (k-j);
            }
        }
        
        i = k - 1;
        for (j=1; j<=i+1; j++) {
            p[j] = (2*i+1) * p[j] / (i+1);
        }
        for (j=3; j<=i+2; j++) {
            for (l=0; l<=i; l++) {
                p[j] = p[j] - i*p0[l] / (i+1);
            }
        }
        for (j=1; j<=i+1; j++) {
            p0[j] = p1[j];
        }
        for (j=1; j<=i+2; j++) {
            p1[j] = p[j];
        }
    }
    
    free_dvector(p, 1, n+1);
    free_dvector(q, 1, n);
    free_dvector(p0, 1, n);
    free_dvector(p1, 1, n+1);
    
}


#pragma mark Public methods...

- (id)init
{
    int i;
    
    self = [super init];
    if (self) {
        // Initialization code here.
        
        _maxElementNodes = 256;
        _maxDeg = 4;
        _maxDeg3 = pow(_maxDeg, 3.0);
        _maxDeg2 = pow(_maxDeg, 2.0);
       
        _point = intmatrix(0, 0, 0, 0);
        _line = intmatrix(0, 0, 0, 1);
        _triangle = intmatrix(0, 2, 0, 1);
        _quad = intmatrix(0, 3, 0, 1);
        _tetra = intmatrix(0, 5, 0, 1);
        _prism = intmatrix(0, 7, 0, 1);
        _wedge = intmatrix(0, 8, 0, 1);
        _brick = intmatrix(0, 11, 0, 1);
        
        for (i=0; i<8; i++) {
            _initialized[i] = NO;
        }
        
        // Construction of definitions
        [self FEMElementDescription_initElementsDefinition];
        
        _numberOfElementDefs = 21;
        
         // Array of objects describing each definition
        _listOfDefinitions = [NSArray arrayWithObjects:_1nodePoint, _2nodePeriodic, _2nodeLine, _3nodeLine, _4nodeLine, _3nodeTriangle, _6nodeTriangle,
                             _10nodeTriangle, _4nodeQuadrilateral, _8nodeQuadrilateral, _9nodeQuadrilateral, _12nodeQuadrilateral, _16nodeQuadrilateral,
                             _4nodeTetrahedron, _10nodeTetrahedron, _5nodePyramid, _13nodePyramid, _6nodeWedge, _8nodeOctahedron, _20nodeOctahedron,
                             _27nodeOctahedron, nil];
        
        _isTypeListInitialized = NO;

    }
    
    return self;
}

-(void)deallocation
{
    free_imatrix(_point, 0, 0, 0, 0);
    free_imatrix(_line, 0, 0, 0, 1);
    free_imatrix(_triangle, 0, 2, 0, 1);
    free_imatrix(_quad, 0, 3, 0, 1);
    free_imatrix(_tetra, 0, 5, 0, 1);
    free_imatrix(_prism, 0, 7, 0, 1);
    free_imatrix(_wedge, 0, 8, 0, 1);
    free_imatrix(_brick, 0, 11, 0, 1);
}

-(void)addDescriptionOfElement:(ElementType_t)element withBasisTerms:(int *)terms {
/***********************************************************************************************
    Add an element description to global list of element types
 
    ElementType_t element  ->  Structure of the element type description
    int *terms             ->  List of terms in the basis function that should be included for
                               this element type. terms(i) is an integer from 1-27 according
                               to the list below
 
***********************************************************************************************/
    
    int i, j, k, n, upow, vpow, wpow;
    double u, v, w;
    double **a;
    FEMLinearAlgebra *algebra;
    ElementType_t *temp;
    
    algebra = [[FEMLinearAlgebra alloc] init];
    
    n = element.NumberOfNodes;
    element.NumberOfEdges = 0;
    element.NumberOfFaces = 0;
    element.BasisFunctionDegree = 0;
    element.BasisFunctions = NULL;
    
    if (element.ElementCode >= 200) {
        
        a = doublematrix(0, n-1, 0, n-1);
        
        //------------------------------
        // 1D line element
        //------------------------------
        if (element.dimension == 1) {
            
            for (i=0; i<n; i++) {
                u = element.NodeU[i];
                for (j=0; j<n; j++) {
                    k = terms[j] -1;
                    upow = k;
                    if (u == 0 && upow == 0) {
                        a[i][j] = 1;
                    } else {
                        a[i][j] = pow(u, upow);
                    }
                    element.BasisFunctionDegree = max(element.BasisFunctionDegree, upow);
                }
            }
        }
        //------------------------------
        // 2D surface element
        //------------------------------
        else if (element.dimension == 2) {
            
            for (i=0; i<n; i++) {
                u = element.NodeU[i];
                v = element.NodeV[i];
                for (j=0; j<n; j++) {
                    k = terms[j] - 1;
                    vpow = k / _maxDeg;
                    upow = k % _maxDeg;
                    
                    if (upow == 0) {
                        a[i][j] = 1;
                    } else {
                        a[i][j] = pow(u, upow);
                    }
                    
                    if (vpow != 0) {
                        a[i][j] = a[i][j] * pow(v, vpow);
                    }
                    
                    element.BasisFunctionDegree = max(element.BasisFunctionDegree, upow);
                    element.BasisFunctionDegree = max(element.BasisFunctionDegree, vpow);
                }
            }
        }
        //------------------------------
        // 3D volume element
        //------------------------------
        else {
            
            for (i=0; i<n; i++) {
                u = element.NodeU[i];
                v = element.NodeV[i];
                w = element.NodeW[i];
                for (j=0; j<n; j++) {
                    k = terms[j] - 1;
                    upow = k % _maxDeg;
                    wpow = k / _maxDeg2;
                    vpow = (k / _maxDeg) % _maxDeg;
                    
                    if (upow == 0) {
                        a[i][j] = 1;
                    } else {
                        a[i][j] = pow(u, upow);
                    }
                    
                    if (vpow != 0) {
                        a[i][j] = a[i][j] * pow(v, vpow);
                    }
                    
                    if (wpow != 0) {
                        a[i][j] = a[i][j] * pow(w, wpow);
                    }
                    
                    element.BasisFunctionDegree = max(element.BasisFunctionDegree, upow);
                    element.BasisFunctionDegree = max(element.BasisFunctionDegree, vpow);
                    element.BasisFunctionDegree = max(element.BasisFunctionDegree, wpow);
                }
            }
        }
        
        // Compute the coefficients of the basis function terms
        [algebra invertMatrix:a ofSize:n];
        
        if (element.ElementCode == 202) {
            element.BasisFunctions = (BasisFunctions_t*) malloc( sizeof(BasisFunctions_t) * 14 );
        } else {
            element.BasisFunctions = (BasisFunctions_t*) malloc( sizeof(BasisFunctions_t) * n );
        }
        
        upow = 0;
        vpow = 0;
        wpow = 0;
        
        for (i=0; i<n; i++) {
            element.BasisFunctions[i].n = n;
            element.BasisFunctions[i].p = intvec(0, n-1);
            element.BasisFunctions[i].q = intvec(0, n-1);
            element.BasisFunctions[i].r = intvec(0, n-1);
            element.BasisFunctions[i].coeff = doublevec(0, n-1);
            
            for (j=0; j<n; j++) {
                k = terms[j] - 1;
                
                switch (element.dimension) {
                    case 1:
                        upow = k;
                        break;
                    case 2:
                        upow = k % _maxDeg;
                        vpow = k / _maxDeg;
                        break;
                    case 3:
                        upow = k % _maxDeg;
                        vpow = (k / _maxDeg) % _maxDeg;
                        wpow = k / _maxDeg2;
                }
                
                element.BasisFunctions[i].p[j] = upow;
                element.BasisFunctions[i].q[j] = vpow;
                element.BasisFunctions[i].r[j] = wpow;
                element.BasisFunctions[i].coeff[j] = a[i][j];
            }
        }
        
        free_dmatrix(a, 0, n-1, 0, n-1);
        
        if (element.ElementCode == 202) {
            a = doublematrix(0, 13, 0, 13);
            for (i=0; i<14; i++) {
                for (j=0; j<14; j++) {
                    a[i][j] = 0.0;
                }
            }
            [self FEMElementDescription_compute1DPBasis:a :14];
            
            for (i=2; i<14; i++) {
                element.BasisFunctions[i].p = intvec(0, i);
                element.BasisFunctions[i].q = intvec(0, i);
                element.BasisFunctions[i].r = intvec(0, i);
                element.BasisFunctions[i].coeff = doublevec(0, i);
                
                k = 0;
                for (j=0; j<=i; j++) {
                    if (a[i][j] != 0.0) {
                        element.BasisFunctions[i].p[k] = j;
                        element.BasisFunctions[i].q[k] = 0;
                        element.BasisFunctions[i].r[k] = 0;
                        element.BasisFunctions[i].coeff[k] = a[i][j];
                        k++;
                    }
                }
                element.BasisFunctions[i].n = k;
            }
            free_dmatrix(a, 0, 13, 0, 13);
        }
        
        switch (element.ElementCode / 100) {
            case 3:
                element.NumberOfEdges = 3;
                break;
            case 4:
                element.NumberOfEdges = 4;
                break;
            case 5:
                element.NumberOfFaces = 4;
                element.NumberOfEdges = 6;
                break;
            case 6:
                element.NumberOfFaces = 5;
                element.NumberOfEdges = 8;
                break;
            case 7:
                element.NumberOfFaces = 5;
                element.NumberOfEdges = 9;
                break;
            case 8:
                element.NumberOfFaces = 6;
                element.NumberOfEdges = 12;
                break;
        }
        
    }
    
    // Finally add the element description to the global list of types
    if (_isTypeListInitialized == NO) {
        _elementTypeList = (ElementType_t *)malloc(sizeof(ElementType_t));
        _elementTypeList = &element;
        _isTypeListInitialized = YES;
        _elementTypeList->NextElementType = NULL;
    } else {
        temp = (ElementType_t *)malloc(sizeof(ElementType_t));
        temp = &element;
        temp->NextElementType = _elementTypeList;
        _elementTypeList = temp;
    }
    
}

-(void)initElementDescriptions {
/***********************************************************************************************

    Read the elements definition provided by the class and add the element types 
    to a global list.
    Users of the class should make sure to call this method in their code.

***********************************************************************************************/
    
    int i, k;
    int *basisTerms;
    ElementType_t element;
    FEMElementsDefinition *defs;
    
    basisTerms = intvec(0, _maxDeg3-1);
    
    // Add the connectivity element types...
    memset( basisTerms, 0.0, (_maxDeg3*sizeof(basisTerms)) );
    element.GaussPoints = 0;
    element.GaussPoints2 = 0;
    element.StabilizationMK = 0;
    element.NodeU = NULL;
    element.NodeV = NULL;
    element.NodeW = NULL;
    for (k=3; k<=64; k++) {
        element.NumberOfNodes = k;
        element.ElementCode = 100 + k;
        [self addDescriptionOfElement:element withBasisTerms:basisTerms];
    }
    
    // ... then the rest of them
    @autoreleasepool {
        
        for (i=0; i<_numberOfElementDefs; i++) {
            
            element.NodeU = NULL;
            element.NodeV = NULL;
            element.NodeW = NULL;
            
            defs = [_listOfDefinitions objectAtIndex:i];
            
            element.dimension = [defs.dimension intValue];
            element.ElementCode = [defs.code intValue];
            element.NumberOfNodes = [defs.nodes intValue];
            if (element.dimension == 1) {
                element.NodeU = doublevec(0, element.NumberOfNodes-1);
                for (k=0; k<element.NumberOfNodes; k++) {
                    element.NodeU[k] = [[defs.nodeU objectAtIndex:k] doubleValue];
                    basisTerms[k] = [[defs.basis objectAtIndex:k] intValue];
                }
            } else if (element.dimension == 2) {
                element.NodeU = doublevec(0, element.NumberOfNodes-1);
                element.NodeV = doublevec(0, element.NumberOfNodes-1);
                for (k=0; k<element.NumberOfNodes; k++) {
                    element.NodeU[k] = [[defs.nodeU objectAtIndex:k] doubleValue];
                    element.NodeV[k] = [[defs.nodeV objectAtIndex:k] doubleValue];
                    basisTerms[k] = [[defs.basis objectAtIndex:k] intValue];
                }
            } else if (element.dimension == 3) {
                element.NodeU = doublevec(0, element.NumberOfNodes-1);
                element.NodeV = doublevec(0, element.NumberOfNodes-1);
                element.NodeW = doublevec(0, element.NumberOfNodes-1);
                for (k=0; k<element.NumberOfNodes; k++) {
                    element.NodeU[k] = [[defs.nodeU objectAtIndex:k] doubleValue];
                    element.NodeV[k] = [[defs.nodeV objectAtIndex:k] doubleValue];
                    element.NodeW[k] = [[defs.nodeW objectAtIndex:k] doubleValue];
                    basisTerms[k] = [[defs.basis objectAtIndex:k] intValue];
                }
            }
            element.GaussPoints = [[defs.gaussPoints objectAtIndex:0] intValue];
            element.GaussPoints2 = [[defs.gaussPoints objectAtIndex:1] intValue];
            element.StabilizationMK = [defs.stabilization doubleValue];
            
            if (element.GaussPoints2 <= 0) element.GaussPoints2 = element.GaussPoints;
            
            // Reassign 0 to stabilization. Don't know why yet!!!
            element.StabilizationMK = 0.0;
            
            if (element.NodeV == NULL) {
                element.NodeV = doublevec(0, element.NumberOfNodes-1);
                memset( element.NodeV, 0.0, (element.NumberOfNodes*sizeof(element.NodeV)) );
            }
            if (element.NodeW == NULL) {
                element.NodeW = doublevec(0, element.NumberOfNodes-1);
                memset( element.NodeW, 0.0, (element.NumberOfNodes*sizeof(element.NodeW)) );
            }
            
            [self addDescriptionOfElement:element withBasisTerms:basisTerms];
            
            free_dvector(element.NodeU, 0, element.NumberOfNodes-1);
            free_dvector(element.NodeV, 0, element.NumberOfNodes-1);
            free_dvector(element.NodeW, 0, element.NumberOfNodes-1);
            
        }
    }
    
}

-(int **)getEdgeMap:(int)elementFamily {
    
    int **edgeMag;
    
    switch (elementFamily) {
        case 1:
            edgeMag = _point;
        case 2:
            edgeMag = _line;
            break;
        case 3:
            edgeMag = _triangle;
            break;
        case 4:
            edgeMag = _quad;
            break;
        case 5:
            edgeMag = _tetra;
            break;
        case 6:
            edgeMag = _prism;
            break;
        case 7:
            edgeMag = _wedge;
            break;
        case 8:
            edgeMag = _brick;
            break;
    }
    
    if (_initialized[elementFamily-1] == NO) {
        _initialized[elementFamily-1] = YES;
        switch (elementFamily) {
            case 1:
                edgeMag[0][0] = 0;
                break;
            case 2:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                break;
            case 3:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 2;
                edgeMag[2][1] = 0;
                break;
            case 4:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 2;
                edgeMag[2][1] = 3;
                
                edgeMag[3][0] = 3;
                edgeMag[3][1] = 0;
                break;
            case 5:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][0] = 2;
                
                edgeMag[2][0] = 2;
                edgeMag[2][1] = 0;
                
                edgeMag[3][0] = 0;
                edgeMag[3][1] = 3;
                
                edgeMag[4][0] = 1;
                edgeMag[4][1] = 3;
                
                edgeMag[5][0] = 2;
                edgeMag[5][1] = 3;
                break;
            case 6:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 3;
                edgeMag[2][1] = 2;
                
                edgeMag[3][0] = 0;
                edgeMag[3][1] = 3;
                
                edgeMag[4][0] = 0;
                edgeMag[4][1] = 4;
                
                edgeMag[5][0] = 1;
                edgeMag[5][1] = 4;
                
                edgeMag[6][0] = 2;
                edgeMag[6][1] = 4;
                
                edgeMag[7][0] = 3;
                edgeMag[7][1] = 4;
                break;
            case 7:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 2;
                edgeMag[2][1] = 0;
                
                edgeMag[3][0] = 3;
                edgeMag[3][1] = 4;
                
                edgeMag[4][0] = 4;
                edgeMag[4][1] = 5;
                
                edgeMag[5][0] = 5;
                edgeMag[5][1] = 3;
                
                edgeMag[6][0] = 0;
                edgeMag[6][1] = 3;
                
                edgeMag[7][0] = 1;
                edgeMag[7][1] = 4;
                
                edgeMag[8][0] = 2;
                edgeMag[8][1] = 5;
                break;
            case 8:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 3;
                edgeMag[2][1] = 2;
                
                edgeMag[3][0] = 0;
                edgeMag[3][1] = 3;
                
                edgeMag[4][0] = 4;
                edgeMag[4][1] = 5;
                
                edgeMag[5][0] = 5;
                edgeMag[5][1] = 6;
                
                edgeMag[6][0] = 7;
                edgeMag[6][1] = 6;
                
                edgeMag[7][0] = 4;
                edgeMag[7][1] = 7;
                
                edgeMag[8][0] = 0;
                edgeMag[8][1] = 4;
                
                edgeMag[9][0] = 1;
                edgeMag[9][1] = 5;
                
                edgeMag[10][0] = 2;
                edgeMag[10][1] = 6;
                
                edgeMag[11][0] = 3;
                edgeMag[11][1] = 7;
                break;
        }
    }
    
    return edgeMag;
    
}

-(double)elementDiameter:(Element_t *)element: (Nodes_t *)nodes {
/*****************************************************************************************
    Retrive element diameter parameter for stabilization
 
    Element_t *element  ->  element
    Nodes_t *nodes      ->  nodal coodinate arrays of the element
 
*****************************************************************************************/
    
    int i, j, k, family;
    double x0, y0, z0, hk, a, s, cx, cy, cz;
    double j11, j12, j13, j21, j22, j23, g11, g12, g22;
    int **edgeMap, size;
    
    family = element->Type.ElementCode / 100;
    
    switch (family) {
        case 1:
            hk = 0.0;
            break;
        case 3: // Triangular element
            j11 = nodes[1].x - nodes[0].x;
            j12 = nodes[1].y - nodes[0].y;
            j13 = nodes[1].z - nodes[0].z;
            j21 = nodes[2].x - nodes[0].x;
            j22 = nodes[2].y - nodes[0].y;
            j23 = nodes[2].z - nodes[0].z;
            g11 = pow(j11, 2.0) + pow(j12, 2.0) + pow(j13, 2.0);
            g12 = j11*j21 + j12+j22 + j13*j23;
            g22 = pow(j21, 2.0) + pow(j22, 2.0) + pow(j23, 2.0);
            a = sqrt(g11*g22 - pow(g12, 2.0)) / 2.0;
            
            cx = ( nodes[0].x + nodes[1].x + nodes[2].x ) / 3.0;
            cy = ( nodes[0].y + nodes[1].y + nodes[2].y ) / 3.0;
            cz = ( nodes[0].z + nodes[1].z + nodes[2].z ) / 3.0;
            
            s = pow((nodes[0].x-cx), 2.0) + pow((nodes[0].y-cy), 2.0) + pow((nodes[0].z-cz), 2.0);
            s = s + pow((nodes[1].x-cx), 2.0) + pow((nodes[1].y-cy), 2.0) + pow((nodes[1].z-cz), 2.0);
            s = s + pow((nodes[2].x-cx), 2.0) + pow((nodes[2].y-cy), 2.0) + pow((nodes[2].z-cz), 2.0);
            
            hk = 16.0 * a * a / (3.0 * s);
            break;
        case 4: // Quadrilateral element
            cx = pow((nodes[1].x-nodes[0].x), 2.0) + pow((nodes[1].y-nodes[0].y), 2.0) + pow((nodes[1].z-nodes[0].z), 2.0);
            cy = pow((nodes[3].x-nodes[0].x), 2.0) + pow((nodes[3].y-nodes[0].y), 2.0) + pow((nodes[3].z-nodes[0].z), 2.0);
            hk = 2.0 * cx * cy / (cx + cy);
        default:
            edgeMap = [self getEdgeMap:family];
            hk = HUGE_VAL;
            switch (family) {
                case 2:
                    size = 1;
                    break;
                case 3:
                    size = 3;
                    break;
                case 4:
                    size = 4;
                    break;
                case 5:
                    size = 6;
                    break;
                case 6:
                    size = 8;
                    break;
                case 7:
                    size = 9;
                    break;
                case 8:
                    size = 12;
                    break;
            }
            for (i=0; i<size; i++) {
                j = edgeMap[i][0];
                k = edgeMap[i][1];
                x0 = nodes[j].x - nodes[k].x;
                y0 = nodes[j].y - nodes[k].y;
                z0 = nodes[j].z - nodes[k].z;
                hk = min( hk, (pow(x0, 2.0)+pow(y0, 2.0)+pow(z0, 2.0)) );
            }
            break;
    }
    edgeMap = NULL;
    return sqrt(hk);
}

-(void)computeStabilizationParameter:(Element_t *)element: (Nodes_t *)nodes: (FEMMesh *)mesh: (int)n: (double)mk: (double *)hk {
/*******************************************************************************************************
    Compute convection diffusion equation stabilization parameter for each and every element of the 
    model by solving the largest eigenvalue of
        
        Lu = \lambda Gu,
        L = (\nabla^2 u,\nabla^2 w), G = (\nabla u,\nabla w)
 
*******************************************************************************************************/
 
    int i, j, k, p, q, t, dim;
    double *eigr, s, *ddp, *ddq, ***dNodalBasisdx;
    double u, v, w, **l, **g, *work, sum1, sum2;
    double *buffer1, *buffer2;
    GaussIntegrationPoints integCompound;
    FEMNumericIntegration *numericIntegration;
    BOOL stat;
    
    // Used for lapack routine
    int itype, order, lda, ldb, lwork, info;
    char *jobz, *uplo;
    
    eigr = doublevec(0, n-1);
    ddp = doublevec(0, 2);
    ddq = doublevec(0, 2);
    dNodalBasisdx = d3tensor(0, n-1, 0, n-1, 0, 2);
    l = doublematrix(0, (n-1)-1, 0, (n-1)-1);
    g = doublematrix(0, (n-1)-1, 0, (n-1)-1);
    work = doublevec(0, (16*n)-1);
    
    numericIntegration = [[FEMNumericIntegration alloc] init];
    [numericIntegration allocation:mesh];
    
    if (element->Type.BasisFunctionDegree <= 1) {
        switch (element->Type.ElementCode) {
            case 202:
            case 303:
            case 404:
            case 504:
            case 605:
            case 706:
                mk = 1.0 / 3.0;
                break;
            case 808:
                mk = 1.0 / 6.0;
        }
        if (hk != NULL) *hk = [self elementDiameter:element :nodes];
        return;
    }
    
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            for (k=0; k<3; k++) {
                dNodalBasisdx[i][j][k] = 0.0;
            }
        }
    }
    for (p=0; p<n; p++) {
        u = element->Type.NodeU[p];
        v = element->Type.NodeV[p];
        w = element->Type.NodeW[p];
        stat = [numericIntegration setBasisFirstDerivativeForElement:element 
                                                        elementNodes:nodes 
                                                              inMesh:mesh 
                                                firstEvaluationPoint:u 
                                               secondEvaluationPoint:v 
                                                thirdEvaluationPoint:w 
                                                         withBubbles:NO 
                                                         basisDegree:NULL];
        for (i=0; i<n; i++) {
            for (j=0; j<3; j++) {
                dNodalBasisdx[i][p][j] = [numericIntegration basisFirstDerivative:i :j];
            }
        }
    }
    
    dim = [mesh dimension];
    integCompound = GaussQuadrature(element);
    
    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1; j++) {
            l[i][j] = 0.0;
            g[i][j] = 0.0;
        }
    }
    for (t =0; t<integCompound.n; t++) {
        u = integCompound.u[t];
        v = integCompound.v[t];
        w = integCompound.w[t];
        
        stat = [numericIntegration setMetricDeterminantForElement:element 
                                                     elementNodes:nodes inMesh:mesh 
                                             firstEvaluationPoint:u 
                                            secondEvaluationPoint:v 
                                             thirdEvaluationPoint:w];
        
        s = [numericIntegration metricDeterminant] * integCompound.s[t];
        
        for (p=1; p<n; p++) {
            for (q=1; q<n; q++) {
               memset( ddp, 0.0, (3*sizeof(ddp)) );
               memset( ddq, 0.0, (3*sizeof(ddq)) );
               for (i=0; i<dim; i++) {
                   g[p-1][q-1] = g[p-1][q-1] + s * [numericIntegration basisFirstDerivative:p :i] * [numericIntegration basisFirstDerivative:q :i];
                   sum1 = 0.0;
                   sum2 = 0.0;
                   for (j=0; j<n; j++) {
                       sum1 = sum1 + dNodalBasisdx[p][j][i] * [numericIntegration basisFirstDerivative:j :i];
                       sum2= sum2 + dNodalBasisdx[q][j][i] * [numericIntegration basisFirstDerivative:j :i];
                   }
                   ddp[i] = ddp[i] + sum1;
                   ddq[i] = ddq[i] + sum2;
               }
                sum1 = 0.0;
                sum2 = 0.0;
               for (i=0; i<3; i++) {
                   sum1 = sum1 + ddp[i];
                   sum2 = sum2 + ddq[i];
               }
                l[p-1][q-1] = l[p-1][q-1] + s * sum1 * sum2; 
            }
        }
    }
    stat = YES;
    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1; j++) {
            if (l[i][j] < AEPS) {
                continue;
            } else {
                stat = NO;
                break;
            }
        }
    }
    if (stat == YES) {
        mk = 1.0 / 3.0;
        if (hk != NULL) {
            *hk = [self elementDiameter:element :nodes];
        }
        return;
    }
    
    // TODO: Here we serialize the arrays l[][] and g[][] because we need to do that for using 
    // them in dsygv_(). Later optimize so that we work initially on serialzed arrays so that 
    // we don't need to do a copy.
    
    buffer1 = doublevec(0, ((n-1)*(n-1))-1);
    buffer2 = doublevec(0, ((n-1)*(n-1))-1);
    
    k = 0;
    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1; j++) {
            buffer1[k] = l[i][j];
            k++;
        }
    }
    k=0;
    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1; j++) {
            buffer2[k] = g[i][j];
            k++;
        }
    }
    
    itype = 1;
    jobz = "N";
    uplo = "U";
    order = n-1;
    lda = n-1;
    ldb = n-1;
    lwork = 12*n;
    dsygv_(&itype, jobz, uplo, &order, buffer1, &lda, buffer2, &ldb, eigr, work, &lwork, &info);
    if (info < 0 || info > 0) {
        warnfunct("computeStabilizationParameter", "Error in lapack routine dsygv. Error code:");
        printf("%d\n", info);
        errorfunct("computeStabilizationParameter", "Program terminating now...");
    }
    mk = eigr[n-2];
    
    free_dvector(buffer1, 0, ((n-1)*(n-1))-1);
    free_dvector(buffer2, 0, ((n-1)*(n-1))-1);
    
    if (mk < 10.0 * AEPS) {
        mk = 1.0 / 3.0;
        if (hk != NULL) {
            *hk = [self elementDiameter:element :nodes];
        }
        return;
    }
    
    if (hk != NULL) {
        *hk = sqrt( 2.0 / (mk * element->Type.StabilizationMK) );
        mk = min( 1.0/3.0, element->Type.StabilizationMK );
    } else {
        switch (element->Type.ElementCode / 100) {
            case 2:
            case 4:
            case 8:
                mk = 4.0 * mk;
                break;
        }
        mk = min( 1.0/3.0, 2.0/mk );
    }

    free_dvector(eigr, 0, n-1);
    free_dvector(ddp, 0, 2);
    free_dvector(ddq, 0, 2);
    free_d3tensor(dNodalBasisdx, 0, n-1, 0, n-1, 0, 2);
    free_dmatrix(l, 0, (n-1)-1, 0, (n-1)-1);
    free_dmatrix(g, 0, (n-1)-1, 0, (n-1)-1);
    free_dvector(work, 0, (16*n)-1);
    
    [numericIntegration deallocation:mesh];
    
}

-(ElementType_t *)getElementType:(int)code inMesh:(FEMMesh *)mesh stabilization:(BOOL *)computeStab {
    
    int i;
    
    Element_t *elm;
    ElementType_t *element;
    Nodes_t *nodes;
    
    element = _elementTypeList;
    
    while (element != NULL) {
        if (code == element->ElementCode) break;
        element = element->NextElementType;
    }
    
    if (element == NULL) {
        errorfunct("getElementType", "Element type code not found:");
        printf("%d\n", code);
        errorfunct("getElementType", "Ignoring element.");
        return NULL;
    }
    
    if (computeStab != NULL) {
        if (*computeStab == NO) return element;
    }
    
    if (element->StabilizationMK == 0.0) {
        elm = (Element_t*) malloc( sizeof(Element_t));
        elm->Type = *element;
        elm->BDOFs = 0;
        elm->DGDOFs = 0;
        elm->Pdefs = NULL;
        elm->DGIndexes = NULL;
        elm->FaceIndexes = NULL;
        elm->BubbleIndexes = NULL;
        nodes = (Nodes_t*)malloc(sizeof(Nodes_t) * element->NumberOfNodes );
        for (i=0; i<element->NumberOfNodes; i++) {
            nodes[i].x = element->NodeU[i];
            nodes[i].y = element->NodeV[i];
            nodes[i].z = element->NodeW[i];
        }
        [self computeStabilizationParameter:elm :nodes :mesh :element->NumberOfNodes :element->StabilizationMK :NULL];
        free(elm);
        free(nodes);
    }
    
    *element = elm->Type;
    return element;
    
}
    

@end
