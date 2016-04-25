//
//  FEMElementDescription.m
//  Saino
//
//  Created by Hakime Seddik on 05/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMElementDescription.h"
#import "FEMLinearAlgebra.h"
#import "FEMElementsDefinition.h"
#import "FEMNumericIntegration.h"
#import "FEMUtilities.h"
#import "GaussIntegration.h"
#import "Utils.h"

@interface FEMElementDescription ()

-(void)FEMElementDescription_initElementsDefinition;
-(void)FEMElementDescription_compute1DPBasis:(double * __nonnull * __nonnull)basis sizeOfBasis:(int)n;
-(double)FEMElementDescription_interpolate1DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evalutationPoint:(double)u;
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
    FEMElementsDefinition * __nonnull _1nodePoint;
    FEMElementsDefinition * __nonnull _2nodePeriodic;
    FEMElementsDefinition * __nonnull _2nodeLine;
    FEMElementsDefinition * __nonnull _3nodeLine;
    FEMElementsDefinition * __nonnull _4nodeLine;
    FEMElementsDefinition * __nonnull _3nodeTriangle;
    FEMElementsDefinition * __nonnull _6nodeTriangle;
    FEMElementsDefinition * __nonnull _10nodeTriangle;
    FEMElementsDefinition * __nonnull _4nodeQuadrilateral;
    FEMElementsDefinition * __nonnull _8nodeQuadrilateral;
    FEMElementsDefinition * __nonnull _9nodeQuadrilateral;
    FEMElementsDefinition * __nonnull _12nodeQuadrilateral;
    FEMElementsDefinition * __nonnull _16nodeQuadrilateral;
    FEMElementsDefinition * __nonnull _4nodeTetrahedron;
    FEMElementsDefinition * __nonnull _10nodeTetrahedron;
    FEMElementsDefinition * __nonnull _5nodePyramid;
    FEMElementsDefinition * __nonnull _13nodePyramid;
    FEMElementsDefinition * __nonnull _6nodeWedge;
    FEMElementsDefinition * __nonnull _8nodeOctahedron;
    FEMElementsDefinition * __nonnull _20nodeOctahedron;
    FEMElementsDefinition * __nonnull _27nodeOctahedron;
    
    NSArray * __nonnull _listOfDefinitions;
    
    //-----------------------------------------------------
    // Linked list strurture which contains the elements
    // definition and their allocations. The elements
    // definition are provided by the class itself.
    //-----------------------------------------------------
    ElementType_t * __nullable _elementTypeList;
    BOOL _isTypeListInitialized;
    
    int * __nonnull * __nonnull _point;
    int * __nonnull * __nonnull _line;
    int * __nonnull * __nonnull _triangle;
    int * __nonnull * __nonnull _quad;
    int * __nonnull * __nonnull _tetra;
    int * __nonnull * __nonnull _prism;
    int * __nonnull * __nonnull _wedge;
    int * __nonnull * __nonnull _brick;
    
    BOOL _initialized[8];
}

#pragma mark Private methods...

-(void)FEMElementDescription_initElementsDefinition {
    
    //---------------------
    // 1 node point element
    //---------------------
    _1nodePoint = [[FEMElementsDefinition alloc] init];
    _1nodePoint.dimension = @1;
    _1nodePoint.topology = @"point";
    _1nodePoint.code = @101;
    _1nodePoint.nodes = @1;
    _1nodePoint.nodeU = @[@0.0];
    _1nodePoint.basis = @[@1];
    _1nodePoint.gaussPoints = @[@1];
    _1nodePoint.stabilization = @0.0;
    
    //-------------------------
    // 2 nodes periodic element
    //-------------------------
    _2nodePeriodic = [[FEMElementsDefinition alloc] init];
    _2nodePeriodic.dimension = @1;
    _2nodePeriodic.topology = @"point";
    _2nodePeriodic.code = @102;
    _2nodePeriodic.nodes = @2;
    _2nodePeriodic.nodeU = @[@0.0, @1.0];
    _2nodePeriodic.basis = @[@1, @2];
    _2nodePeriodic.gaussPoints = @[@0];
    _2nodePeriodic.stabilization = @0.0;
    
    //---------------------
    // 2 nodes line element
    //---------------------
    _2nodeLine = [[FEMElementsDefinition alloc] init];
    _2nodeLine.dimension = @1;
    _2nodeLine.topology = @"line";
    _2nodeLine.code = @202;
    _2nodeLine.nodes = @2;
    //                      1    2 
    _2nodeLine.nodeU = @[@-1.0, @1.0];
    //          1  2   3   4
    //          1  u  u^2 u^3
    _2nodeLine.basis = @[@1, @2];
    _2nodeLine.gaussPoints = @[@2, @4];
    _2nodeLine.stabilization = @0.33333333333333333333333;
    
    //---------------------
    // 3 nodes line element
    //---------------------
    _3nodeLine = [[FEMElementsDefinition alloc] init];
    _3nodeLine.dimension = @1;
    _3nodeLine.topology = @"line";
    _3nodeLine.code = @203;
    _3nodeLine.nodes = @3;
    //                     1     2     3
    _3nodeLine.nodeU = @[@-1.0, @1.0, @0.0];
    //          1  2   3   4
    //          1  u  u^2 u^3
    _3nodeLine.basis = @[@1, @2, @3];
    _3nodeLine.gaussPoints = @[@3, @5];
    _3nodeLine.stabilization = @0.1666666666666666666666;
    
    //---------------------
    // 4 nodes line element
    //---------------------
    _4nodeLine = [[FEMElementsDefinition alloc] init];
    _4nodeLine.dimension = @1;
    _4nodeLine.topology = @"line";
    _4nodeLine.code = @204;
    _4nodeLine.nodes = @4;
    //                      1     2     3                      4    
    _4nodeLine.nodeU = @[@-1.0, @1.0, @-0.333333333333333333, @0.3333333333333333333];
    //          1  2   3   4
    //          1  u  u^2 u^3
    _4nodeLine.basis = @[@1, @2, @3, @4];
    _4nodeLine.gaussPoints = @[@4, @6];
    _4nodeLine.stabilization = @0.03333333333333333333333;
    
    //-------------------------
    // 3 nodes triangle element
    //-------------------------
    _3nodeTriangle = [[FEMElementsDefinition alloc] init];
    _3nodeTriangle.dimension = @2;
    _3nodeTriangle.topology = @"triangle";
    _3nodeTriangle.code = @303;
    _3nodeTriangle.nodes = @3;
    //                        1     2     3
    _3nodeTriangle.nodeU = @[@0.0, @1.0, @0.0];
    _3nodeTriangle.nodeV = @[@0.0, @0.0, @1.0];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _3nodeTriangle.basis = @[@1, @2, @5];
    _3nodeTriangle.gaussPoints = @[@3, @20, @1];
    _3nodeTriangle.stabilization = @0.333333333333333333333;
    
    //-------------------------
    // 6 nodes triangle element
    //-------------------------
    _6nodeTriangle = [[FEMElementsDefinition alloc] init];
    _6nodeTriangle.dimension = @2;
    _6nodeTriangle.topology = @"triangle";
    _6nodeTriangle.code = @306;
    _6nodeTriangle.nodes = @6;
    //                         1     2     3     4     5     6
    _6nodeTriangle.nodeU = @[@0.0, @1.0, @0.0, @0.5, @0.5, @0.0];
    _6nodeTriangle.nodeV = @[@0.0, @0.0, @1.0, @0.0, @0.5, @0.5];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _6nodeTriangle.basis = @[@1, @2, @3, @5, @6, @9];
    _6nodeTriangle.gaussPoints = @[@6, @17];
    _6nodeTriangle.stabilization = @0.041666666666666666;
    
    //--------------------------
    // 10 nodes triangle element
    //--------------------------
    _10nodeTriangle = [[FEMElementsDefinition alloc] init];
    _10nodeTriangle.dimension = @2;
    _10nodeTriangle.topology = @"triangle";
    _10nodeTriangle.code = @310;
    _10nodeTriangle.nodes = @10;
    //                          1     2     3     ......
    _10nodeTriangle.nodeU = @[@0.0, @1.0, @0.0, @0.333333333333333333, @0.666666666666666667,
                             @0.666666666666666667, @0.333333333333333333, @0.000000000000000000, @0.000000000000000000, @0.333333333333333333];
    _10nodeTriangle.nodeV = @[@0.0, @0.0, @1.0, @0.000000000000000000, @0.000000000000000000,
                             @0.333333333333333333, @0.666666666666666667, @0.666666666666666667, @0.333333333333333333, @0.333333333333333333];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _10nodeTriangle.basis = @[@1, @2, @3, @4, @5, @6, @7, @9, @10, @13];
    _10nodeTriangle.gaussPoints = @[@25, @64];
    _10nodeTriangle.stabilization = @0.01341555597798937329;
       
    //------------------------------
    // 4 nodes quadrilateral element
    //------------------------------
    _4nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _4nodeQuadrilateral.dimension = @2;
    _4nodeQuadrilateral.topology = @"quadrilateral";
    _4nodeQuadrilateral.code = @404;
    _4nodeQuadrilateral.nodes = @4;
    //                              1     2     3     4
    _4nodeQuadrilateral.nodeU = @[@-1.0, @1.0, @1.0, @-1.0];
    _4nodeQuadrilateral.nodeV = @[@-1.0, @-1.0, @1.0, @1.0];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _4nodeQuadrilateral.basis = @[@1, @2, @5, @6];
    _4nodeQuadrilateral.gaussPoints = @[@4, @16];
    _4nodeQuadrilateral.stabilization = @0.333333333333333333333;
    
    //------------------------------
    // 8 nodes quadrilateral element
    //------------------------------
    _8nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _8nodeQuadrilateral.dimension = @2;
    _8nodeQuadrilateral.topology = @"quadrilateral";
    _8nodeQuadrilateral.code = @408;
    _8nodeQuadrilateral.nodes = @8;
    //                              1     2     3     4     5     6      7      8
    _8nodeQuadrilateral.nodeU = @[@-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0, @-1.0];
    _8nodeQuadrilateral.nodeV = @[@-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _8nodeQuadrilateral.basis = @[@1, @2, @3, @5, @6, @7, @9, @10];
    _8nodeQuadrilateral.gaussPoints = @[@9, @25];
    _8nodeQuadrilateral.stabilization = @0.0833333333333333333333;
    
    //------------------------------
    // 9 nodes quadrilateral element
    //------------------------------
    _9nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _9nodeQuadrilateral.dimension = @2;
    _9nodeQuadrilateral.topology = @"quadrilateral";
    _9nodeQuadrilateral.code = @409;
    _9nodeQuadrilateral.nodes = @9;
    //                              1     2     3     4     5     6      7      8      9
    _9nodeQuadrilateral.nodeU = @[@-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0, @-1.0, @0.0];
    _9nodeQuadrilateral.nodeV = @[@-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0, @0.0];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _9nodeQuadrilateral.basis = @[@1, @2, @3, @5, @6, @7, @9, @10, @11];
    _9nodeQuadrilateral.gaussPoints = @[@9, @25];
    _9nodeQuadrilateral.stabilization = @0.0833333333333333333333;
    
    //-------------------------------
    // 12 nodes quadrilateral element
    //-------------------------------
    _12nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _12nodeQuadrilateral.dimension = @2;
    _12nodeQuadrilateral.topology = @"quadrilateral";
    _12nodeQuadrilateral.code = @412;
    _12nodeQuadrilateral.nodes = @12;
    //                               1     2     3     ......
    _12nodeQuadrilateral.nodeU = @[@-1.0, @1.0, @1.0, 
                                  @-1.0, @-0.3333333333333333333, @0.3333333333333333333,
                                  @1.0000000000000000000, @1.0000000000000000000, 
                                  @0.3333333333333333333, @-0.3333333333333333333, 
                                  @-1.0000000000000000000, @-1.0000000000000000000];
    _12nodeQuadrilateral.nodeV = @[@-1.0, @-1.0, @1.0, 
                                  @1.0, @-1.0000000000000000000, @-1.0000000000000000000,
                                  @-0.3333333333333333333, @0.3333333333333333333, 
                                  @1.0000000000000000000, @1.0000000000000000000, 
                                  @0.3333333333333333333, @-0.3333333333333333333];
    //       1  2   3   4   5   6   7     8    9     10 
    //       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
    //
    //      11      12     13 14    15     16
    //      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
    _12nodeQuadrilateral.basis = @[@1, @2, @3, @4, @5, @6, @7, @8, @9, @10, @13, @14];
    _12nodeQuadrilateral.gaussPoints = @[@16, @36];
    _12nodeQuadrilateral.stabilization = @0.0;
    
    //-------------------------------
    // 16 nodes quadrilateral element
    //-------------------------------
    _16nodeQuadrilateral = [[FEMElementsDefinition alloc] init];
    _16nodeQuadrilateral.dimension = @2;
    _16nodeQuadrilateral.topology = @"quadrilateral";
    _16nodeQuadrilateral.code = @416;
    _16nodeQuadrilateral.nodes = @16;
    //                               1     2     3     ......
    _16nodeQuadrilateral.nodeU = @[@-1.0, @1.0, @1.0, 
                                  @-1.0, @-0.3333333333333333333, @0.3333333333333333333,
                                  @1.0000000000000000000, @1.0000000000000000000, 
                                  @0.3333333333333333333, @-0.3333333333333333333, 
                                  @-1.0000000000000000000, @-1.0000000000000000000, 
                                  @-0.3333333333333333333, @0.3333333333333333333, 
                                  @0.3333333333333333333, @-0.3333333333333333333];
    _16nodeQuadrilateral.nodeV = @[@-1.0, @-1.0, @1.0, 
                                  @1.0, @-1.0000000000000000000, @-1.0000000000000000000,
                                  @-0.3333333333333333333, @0.3333333333333333333, 
                                  @1.0000000000000000000, @1.0000000000000000000, 
                                  @0.3333333333333333333, @-0.3333333333333333333, 
                                  @-0.3333333333333333333, @-0.3333333333333333333, 
                                  @0.3333333333333333333, @0.3333333333333333333];
    _16nodeQuadrilateral.basis = @[@1, @2, @3, @4, @5, @6, @7, @8, @9, @10, @11, @12, @13, @14, @15, @16];
    _16nodeQuadrilateral.gaussPoints = @[@16, @36];
    _16nodeQuadrilateral.stabilization = @0.01766875890919188522;
    
    //----------------------------
    // 4 nodes tetrahedron element
    //----------------------------
    _4nodeTetrahedron = [[FEMElementsDefinition alloc] init];
    _4nodeTetrahedron.dimension = @3;
    _4nodeTetrahedron.topology = @"tetrahedron";
    _4nodeTetrahedron.code = @504;
    _4nodeTetrahedron.nodes = @4;
    //                           1     2     3     4
    _4nodeTetrahedron.nodeU = @[@0.0, @1.0, @0.0, @0.0];
    _4nodeTetrahedron.nodeV = @[@0.0, @0.0, @1.0, @0.0];
    _4nodeTetrahedron.nodeW = @[@0.0, @0.0, @0.0, @1.0];
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
    _4nodeTetrahedron.basis = @[@1, @2, @5, @17];
    _4nodeTetrahedron.gaussPoints = @[@4, @11, @1];
    _4nodeTetrahedron.stabilization = @0.333333333333333333;
    
    //-----------------------------
    // 10 nodes tetrahedron element
    //-----------------------------
    _10nodeTetrahedron = [[FEMElementsDefinition alloc] init];
    _10nodeTetrahedron.dimension = @3;
    _10nodeTetrahedron.topology = @"tetrahedron";
    _10nodeTetrahedron.code = @510;
    _10nodeTetrahedron.nodes = @10;
    //                           1     2     3     4     5     6      7      8      9     10
    _10nodeTetrahedron.nodeU = @[@0.0, @1.0, @0.0, @0.0, @0.5, @0.5, @0.0, @0.0, @0.5, @0.0];
    _10nodeTetrahedron.nodeV = @[@0.0, @0.0, @1.0, @0.0, @0.0, @0.5, @0.5, @0.0, @0.0, @0.5];
    _10nodeTetrahedron.nodeW = @[@0.0, @0.0, @0.0, @1.0, @0.0, @0.0, @0.0, @0.5, @0.5, @0.5];
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
    _10nodeTetrahedron.basis = @[@1, @2, @3, @5, @6, @9, @17, @18, @21, @33];
    _10nodeTetrahedron.gaussPoints = @[@27, @64];
    _10nodeTetrahedron.stabilization = @0.0416666666666666666;
    

    
    //------------------------
    // 5 nodes pyramid element
    //------------------------
    _5nodePyramid = [[FEMElementsDefinition alloc] init];
    _5nodePyramid.dimension = @3;
    _5nodePyramid.topology = @"pyramid";
    _5nodePyramid.code = @605;
    _5nodePyramid.nodes = @5;
    //                        1     2     3     4     5
    _5nodePyramid.nodeU = @[@-1.0, @1.0, @1.0, @-1.0, @0.0];
    _5nodePyramid.nodeV = @[@-1.0, @-1.0, @1.0, @1.0, @0.0];
    _5nodePyramid.nodeW = @[@0.0, @0.0, @0.0, @0.0, @1.0];
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
    _5nodePyramid.basis = @[@1, @2, @5, @6, @17];
    _5nodePyramid.gaussPoints = @[@8, @64];
    _5nodePyramid.stabilization = @0.333333333333333333;
    
    //-------------------------
    // 13 nodes pyramid element
    //-------------------------
    _13nodePyramid = [[FEMElementsDefinition alloc] init];
    _13nodePyramid.dimension = @3;
    _13nodePyramid.topology = @"pyramid";
    _13nodePyramid.code = @613;
    _13nodePyramid.nodes = @13;
    //                          1     2     3     4     5     6      7      8      9     10    11    12    13
    _13nodePyramid.nodeU = @[@-1.0, @1.0, @1.0, @-1.0, @0.0, @0.0, @1.0, @0.0, @-1.0, @-0.5, @0.5, @0.5, @-0.5];
    _13nodePyramid.nodeV = @[@-1.0, @-1.0, @1.0, @1.0, @0.0, @-1.0, @0.0, @1.0, @0.0, @-0.5, @-0.5, @0.5, @0.5];
    _13nodePyramid.nodeW = @[@0.0, @0.0, @0.0, @0.0, @1.0, @0.0, @0.0, @0.0, @0.0, @0.5, @0.5, @0.5, @0.5];
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
    _13nodePyramid.basis = @[@1, @2, @3, @5, @6, @7, @9, @10, @17, @18, @21, @22, @33];
    _13nodePyramid.gaussPoints = @[@27, @125];
    _13nodePyramid.stabilization = @0.333333333333333333;
    
    //----------------------
    // 6 nodes wedge element
    //----------------------
    _6nodeWedge = [[FEMElementsDefinition alloc] init];
    _6nodeWedge.dimension = @3;
    _6nodeWedge.topology = @"wedge";
    _6nodeWedge.code = @706;
    _6nodeWedge.nodes = @6;
    //                     1     2     3     4     5     6
    _6nodeWedge.nodeU = @[@0.0, @1.0, @0.0, @0.0, @1.0, @0.0];
    _6nodeWedge.nodeV = @[@0.0, @0.0, @1.0, @0.0, @0.0, @1.0];
    _6nodeWedge.nodeW = @[@-1.0, @-1.0, @-1.0, @1.0, @1.0, @1.0];
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
    _6nodeWedge.basis = @[@1, @2, @5, @17, @18, @21];
    _6nodeWedge.gaussPoints = @[@8, @64];
    _6nodeWedge.stabilization = @0.333333333333333333;
    
    //---------------------------
    // 8 nodes octahedron element
    //---------------------------
    _8nodeOctahedron = [[FEMElementsDefinition alloc] init];
    _8nodeOctahedron.dimension = @3;
    _8nodeOctahedron.topology = @"brick";
    _8nodeOctahedron.code = @808;
    _8nodeOctahedron.nodes = @8;
    //                           1     2     3     4     5      6     7     8
    _8nodeOctahedron.nodeU = @[@-1.0, @1.0, @1.0, @-1.0, @-1.0, @1.0, @1.0, @-1.0];
    _8nodeOctahedron.nodeV = @[@-1.0, @-1.0, @1.0, @1.0, @-1.0, @-1.0, @1.0, @1.0];
    _8nodeOctahedron.nodeW = @[@-1.0, @-1.0, @-1.0, @-1.0, @1.0, @1.0, @1.0, @1.0];
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
    _8nodeOctahedron.basis = @[@1, @2, @5, @6, @17, @18, @21, @22];
    _8nodeOctahedron.gaussPoints = @[@8, @64];
    _8nodeOctahedron.stabilization = @0.166666666666666666667;
    
    //----------------------------
    // 20 nodes octahedron element
    //----------------------------
    _20nodeOctahedron = [[FEMElementsDefinition alloc] init];
    _20nodeOctahedron.dimension = @3;
    _20nodeOctahedron.topology = @"brick";
    _20nodeOctahedron.code = @820;
    _20nodeOctahedron.nodes = @20;
    //                            1     2     3     4       5     6      7      8      9     10    11    12    13    14    15     16   17     18    19    20
    _20nodeOctahedron.nodeU = @[@-1.0, @1.0, @1.0,  @-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0, @-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0, @-1.0];
    _20nodeOctahedron.nodeV = @[@-1.0, @-1.0, @1.0, @1.0, @-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0, @-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0];
    _20nodeOctahedron.nodeW = @[@-1.0, @-1.0, @-1.0, @-1.0, @1.0, @1.0, @1.0, @1.0, @-1.0, @-1.0, @-1.0, @-1.0, @0.0, @0.0, @0.0, @0.0, @1.0, @1.0, @1.0, @1.0];
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
    _20nodeOctahedron.basis = @[@1, @2, @3, @5, @6, @7, @9, @10, @17, @18, @19, @21, @22, @23, @25, @26, @33, @34, @37, @38];
    _20nodeOctahedron.gaussPoints = @[@27, @125];
    _20nodeOctahedron.stabilization = @0.08148148148148;
    
    //----------------------------
    // 27 nodes octahedron element
    //----------------------------
    _27nodeOctahedron = [[FEMElementsDefinition alloc] init];
    _27nodeOctahedron.dimension = @3;
    _27nodeOctahedron.topology = @"brick";
    _27nodeOctahedron.code = @827;
    _27nodeOctahedron.nodes = @27;
    //                           1     2     3     4     5     6      7      8      9     10    11    12    13    14    15   16   17    18    19    20    21    22
    //                            23   24    25    26    27
    _27nodeOctahedron.nodeU = @[@-1.0, @1.0, @1.0, @-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0, @-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0, 
                               @-1.0, @0.0, @1.0, @0.0, @-1.0, @0.0, @0.0, @0.0];
    _27nodeOctahedron.nodeV = @[@-1.0, @-1.0, @1.0, @1.0, @-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, @0.0, @-1.0, @-1.0, @1.0, @1.0, @-1.0, @0.0, @1.0, 
                               @0.0, @-1.0, @0.0, @1.0, @0.0, @0.0, @0.0, @0.0];
    _27nodeOctahedron.nodeW = @[@-1.0, @-1.0, @-1.0, @-1.0, @1.0, @1.0, @1.0, @1.0, @-1.0, @-1.0, @-1.0, @-1.0, @0.0, @0.0, @0.0, @0.0, @1.0, @1.0, @1.0, 
                               @1.0, @0.0, @0.0, @0.0, @0.0, @-1.0, @1.0, @0.0];
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
    _27nodeOctahedron.basis = @[@1, @2, @3, @5, @6, @7, @9, @10, @11, @17, @18, @19, @21, @22, @23, @25, @26, @27, @33, @34, @35, @37, @38, @39, @41, @42, @43];
    _27nodeOctahedron.gaussPoints = @[@27, @125];
    _27nodeOctahedron.stabilization = @0.08148148148148;
    
}

/***************************************************************************
 
    Method to compute the 1D P-basis from Legendre polynomials
 
***************************************************************************/
-(void)FEMElementDescription_compute1DPBasis:(double * __nonnull * __nonnull)basis sizeOfBasis:(int)n {
 
    int i, j, k, l;
    double s, p[n+1], p0[n], p1[n+1];
        
    if (n <= 1) {
        basis[0][0] = 1.0;
        return;
    }
    
    //--------------------------------------------------------------------------
    // Compute coefficients of nth Legendre polynomial from the recurrence:
    // (i+1)P_{i+1}(x) = (2i+1)*x*P_i(x) - i*P{i-1}(x), P_{0} = 1; P_{1} = x;
    // IMPORTANT: Computed coefficients inaccurate for n > ~15
    //--------------------------------------------------------------------------
    memset( p, 0.0, sizeof(p) );
    memset( p0, 0.0, sizeof(p0) );
    memset( p1, 0.0, sizeof(p1) );
    p0[0] = 1;
    p1[0] = 1;
    p1[1] = 0;
    
    basis[0][0] = 0.5;
    basis[0][1] = -0.5;
    
    basis[1][0] = 0.5;
    basis[1][1] = 0.5;
    
    for (k=1; k<n; k++) {
        if (k > 1) {
            s = sqrt( (2.0*k-1) / 2.0 );
            for (j=0; j<k; j++) {
                basis[k][k-j] = s * p0[j] / (k-j);
                basis[k][0] = basis[k][0] - s * p0[j]*pow((-1), (j+2)) / (k-j);
            }
        }
        
        i = k;
        for (j=0; j<i+1; j++) {
            p[j] = (2*i+1) * p[j] / (i+1);
        }
        l = 0;
        for (j=2; j<i+2; j++) {
            p[j] = p[j] - i*p0[l] / (i+1);
            l++;
        }
        for (j=0; j<i+1; j++) {
            p0[j] = p1[j];
        }
        for (j=0; j<i+2; j++) {
            p1[j] = p[j];
        }
    }
}

/**********************************************************************************************
 
    Given element structure, return value of a quantity x given at element nodes at local
    cooidinate point u inside the element. Element basis functions are used to compute the
    value. Used for 1d elements
 
    Element_t *element  ->  element structure
    double *x           ->  nodal values of the quantity whose partial derivative is required
    double u            ->  point at which to evaluate the partial derivative
 
    Return y = x(u)
 
**********************************************************************************************/
-(double)FEMElementDescription_interpolate1DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evalutationPoint:(double)u {
    
    int i, n;
    int *p;
    double y, s;
    double *coeff;
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (x[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                s = s + coeff[i] * pow(u, p[i]);
            }
            y = y + s * x[n];
        }
    }
    
    return y;
}


#pragma mark Singleton method

static FEMElementDescription * __nullable sharedElementDescription = nil;
static dispatch_once_t onceToken;

+(id __nonnull)sharedElementDescription {
    
    dispatch_once(&onceToken, ^{
        sharedElementDescription = [[self alloc] init];
    });
    return sharedElementDescription;
}

+(void) selfDestruct {
    sharedElementDescription = nil;
    onceToken = 0;
}

#pragma mark Public methods...

- (id)init
{
    int i;
    
    self = [super init];
    if (self) {
        
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
        
        _elementTypeList = NULL;
        
        // Construction of definitions
        [self FEMElementDescription_initElementsDefinition];
        
        _numberOfElementDefs = 21;
        
         // Array of objects describing each definition
        _listOfDefinitions = @[_1nodePoint, _2nodePeriodic, _2nodeLine, _3nodeLine, _4nodeLine, _3nodeTriangle, _6nodeTriangle,
                             _10nodeTriangle, _4nodeQuadrilateral, _8nodeQuadrilateral, _9nodeQuadrilateral, _12nodeQuadrilateral, _16nodeQuadrilateral,
                             _4nodeTetrahedron, _10nodeTetrahedron, _5nodePyramid, _13nodePyramid, _6nodeWedge, _8nodeOctahedron, _20nodeOctahedron,
                             _27nodeOctahedron];
        
        _isTypeListInitialized = NO;
        
        [self initElementDescriptions];
    }
    
    return self;
}

-(void)deallocation {
    
    int i;
    ElementType_t *ptr;
    
    free_imatrix(_point, 0, 0, 0, 0);
    free_imatrix(_line, 0, 0, 0, 1);
    free_imatrix(_triangle, 0, 2, 0, 1);
    free_imatrix(_quad, 0, 3, 0, 1);
    free_imatrix(_tetra, 0, 5, 0, 1);
    free_imatrix(_prism, 0, 7, 0, 1);
    free_imatrix(_wedge, 0, 8, 0, 1);
    free_imatrix(_brick, 0, 11, 0, 1);
    
    ptr = NULL;
    while (_elementTypeList != NULL) {
        ptr = _elementTypeList->NextElementType;
        if (_elementTypeList->ElementCode >= 200) {
            for (i=0; i<_elementTypeList->NumberOfNodes; i++) {
                free_ivector(_elementTypeList->BasisFunctions[i].p, 0, _elementTypeList->NumberOfNodes-1);
                free_ivector(_elementTypeList->BasisFunctions[i].q, 0, _elementTypeList->NumberOfNodes-1);
                free_ivector(_elementTypeList->BasisFunctions[i].r, 0, _elementTypeList->NumberOfNodes-1);
                free_dvector(_elementTypeList->BasisFunctions[i].coeff, 0, _elementTypeList->NumberOfNodes-1);
            }
            if (_elementTypeList->ElementCode == 202) {
                for (i=2; i<14; i++) {
                    free_ivector(_elementTypeList->BasisFunctions[i].p, 0, i);
                    free_ivector(_elementTypeList->BasisFunctions[i].q, 0, i);
                    free_ivector(_elementTypeList->BasisFunctions[i].r, 0, i);
                    free_dvector(_elementTypeList->BasisFunctions[i].coeff, 0, i);
                }
            }
            free(_elementTypeList->BasisFunctions);
            
            free_dvector(_elementTypeList->NodeU, 0, _elementTypeList->NumberOfNodes-1);
            free_dvector(_elementTypeList->NodeV, 0, _elementTypeList->NumberOfNodes-1);
            free_dvector(_elementTypeList->NodeW, 0, _elementTypeList->NumberOfNodes-1);
        }
        free(_elementTypeList);
        _elementTypeList = ptr;
    }
}

/***********************************************************************************************
    Add an element description to global list of element types
 
    ElementType_t element  ->  Structure of the element type description
    int *terms             ->  List of terms in the basis function that should be included for
                               this element type. terms(i) is an integer from 1-27 according
                               to the list below
 
***********************************************************************************************/
-(void)addDescriptionOfElement:(ElementType_t * __nonnull)element withBasisTerms:(int * __nonnull)terms {
    
    int i, j, k, n, upow, vpow, wpow;
    double u, v, w;
    FEMLinearAlgebra *algebra;
    ElementType_t *temp = NULL;
    
    algebra = [[FEMLinearAlgebra alloc] init];
    
    n = element->NumberOfNodes;
    element->NumberOfEdges = 0;
    element->NumberOfFaces = 0;
    element->BasisFunctionDegree = 0;
    element->BasisFunctions = NULL;
    
    if (element->ElementCode >= 200) {
        
        double **a = doublematrix(0, n-1, 0, n-1);
        
        //------------------------------
        // 1D line element
        //------------------------------
        if (element->dimension == 1) {
            
            for (i=0; i<n; i++) {
                u = element->NodeU[i];
                for (j=0; j<n; j++) {
                    k = terms[j] -1;
                    upow = k;
                    if (u == 0 && upow == 0) {
                        a[i][j] = 1;
                    } else {
                        a[i][j] = pow(u, upow);
                    }
                    element->BasisFunctionDegree = max(element->BasisFunctionDegree, upow);
                }
            }
        }
        //------------------------------
        // 2D surface element
        //------------------------------
        else if (element->dimension == 2) {
            
            for (i=0; i<n; i++) {
                u = element->NodeU[i];
                v = element->NodeV[i];
                for (j=0; j<n; j++) {
                    k = terms[j] - 1;
                    vpow = k / _maxDeg;
                    upow = k % _maxDeg;
                    
                    if (upow == 0) {
                        a[i][j] = 1.0;
                    } else {
                        a[i][j] = pow(u, upow);
                    }
                    
                    if (vpow != 0) {
                        a[i][j] = a[i][j] * pow(v, vpow);
                    }
                    
                    element->BasisFunctionDegree = max(element->BasisFunctionDegree, upow);
                    element->BasisFunctionDegree = max(element->BasisFunctionDegree, vpow);
                }
            }
        }
        //------------------------------
        // 3D volume element
        //------------------------------
        else {
            
            for (i=0; i<n; i++) {
                u = element->NodeU[i];
                v = element->NodeV[i];
                w = element->NodeW[i];
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
                    
                    element->BasisFunctionDegree = max(element->BasisFunctionDegree, upow);
                    element->BasisFunctionDegree = max(element->BasisFunctionDegree, vpow);
                    element->BasisFunctionDegree = max(element->BasisFunctionDegree, wpow);
                }
            }
        }
        
        // Compute the coefficients of the basis function terms
        [algebra invertMatrix:a ofSize:n];
        
        if (element->ElementCode == 202) {
            element->BasisFunctions = (BasisFunctions_t*) malloc( sizeof(BasisFunctions_t) * 14 );
        } else {
            element->BasisFunctions = (BasisFunctions_t*) malloc( sizeof(BasisFunctions_t) * n );
        }
                
        upow = 0;
        vpow = 0;
        wpow = 0;
        
        for (i=0; i<n; i++) {
            element->BasisFunctions[i].n = n;
            element->BasisFunctions[i].p = intvec(0, n-1);
            element->BasisFunctions[i].q = intvec(0, n-1);
            element->BasisFunctions[i].r = intvec(0, n-1);
            element->BasisFunctions[i].coeff = doublevec(0, n-1);
            
            for (j=0; j<n; j++) {
                k = terms[j] - 1;
                
                switch (element->dimension) {
                    case 1:
                        upow = k;
                        break;
                    case 2:
                        vpow = k / _maxDeg;
                        upow = k % _maxDeg;
                        break;
                    case 3:
                        upow = k % _maxDeg;
                        vpow = (k / _maxDeg) % _maxDeg;
                        wpow = k / _maxDeg2;
                }
                
                element->BasisFunctions[i].p[j] = upow;
                element->BasisFunctions[i].q[j] = vpow;
                element->BasisFunctions[i].r[j] = wpow;
                element->BasisFunctions[i].coeff[j] = a[j][i];
            }
        }
        
        free_dmatrix(a, 0, n-1, 0, n-1);
        
        if (element->ElementCode == 202) {
            double **a = doublematrix(0, 13, 0, 13);
            memset( *a, 0.0, (14*14)*sizeof(double) );
            [self FEMElementDescription_compute1DPBasis:a sizeOfBasis:14];
            
            for (i=2; i<14; i++) {
                element->BasisFunctions[i].p = intvec(0, i);
                element->BasisFunctions[i].q = intvec(0, i);
                element->BasisFunctions[i].r = intvec(0, i);
                element->BasisFunctions[i].coeff = doublevec(0, i);
                
                k = 0;
                for (j=0; j<=i; j++) {
                    if (a[i][j] != 0.0) {
                        element->BasisFunctions[i].p[k] = j;
                        element->BasisFunctions[i].q[k] = 0;
                        element->BasisFunctions[i].r[k] = 0;
                        element->BasisFunctions[i].coeff[k] = a[i][j];
                        k++;
                    }
                }
                element->BasisFunctions[i].n = k;
            }
            free_dmatrix(a, 0, 13, 0, 13);
        }
        
        switch (element->ElementCode / 100) {
            case 3:
                element->NumberOfEdges = 3;
                break;
            case 4:
                element->NumberOfEdges = 4;
                break;
            case 5:
                element->NumberOfFaces = 4;
                element->NumberOfEdges = 6;
                break;
            case 6:
                element->NumberOfFaces = 5;
                element->NumberOfEdges = 8;
                break;
            case 7:
                element->NumberOfFaces = 5;
                element->NumberOfEdges = 9;
                break;
            case 8:
                element->NumberOfFaces = 6;
                element->NumberOfEdges = 12;
                break;
        }
        
    }
    
    // Finally add the element description to the global list of types
    if (_isTypeListInitialized == NO) {
        _elementTypeList = element;
        _isTypeListInitialized = YES;
        _elementTypeList->NextElementType = NULL;
    } else {
        temp = element;
        temp->NextElementType = _elementTypeList;
        _elementTypeList = temp;
    }
}

/***********************************************************************************************
 
    Read the elements definition provided by the class and add the element types
    to a global list.
    Users of the class should make sure to call this method in their code.
 
***********************************************************************************************/
-(void)initElementDescriptions {
    
    int i, k;
    int *basisTerms;
    ElementType_t *element = NULL;
    FEMElementsDefinition *defs;
    
    basisTerms = intvec(0, _maxDeg3-1);
    
    // Add the connectivity element types...
    memset( basisTerms, 0, _maxDeg3*sizeof(int) );
    for (k=3; k<=64; k++) {
        element = (ElementType_t *)malloc(sizeof(ElementType_t));
        element->GaussPoints = 0;
        element->GaussPoints2 = 0;
        element->StabilizationMK = 0;
        element->NodeU = NULL;
        element->NodeV = NULL;
        element->NodeW = NULL;
        element->NumberOfNodes = k;
        element->ElementCode = 100 + k;
        [self addDescriptionOfElement:element withBasisTerms:basisTerms];
    }
    
    // ... then the rest of them
        
    for (i=0; i<_numberOfElementDefs; i++) {
        
        element = (ElementType_t *)malloc(sizeof(ElementType_t));
        
        element->NodeU = NULL;
        element->NodeV = NULL;
        element->NodeW = NULL;
        
        defs = _listOfDefinitions[i];
        
        element->dimension = [defs.dimension intValue];
        element->ElementCode = [defs.code intValue];
        element->NumberOfNodes = [defs.nodes intValue];
        if (element->dimension == 1) {
            element->NodeU = doublevec(0, element->NumberOfNodes-1);
            for (k=0; k<element->NumberOfNodes; k++) {
                element->NodeU[k] = [(defs.nodeU)[k] doubleValue];
                basisTerms[k] = [(defs.basis)[k] intValue];
            }
        } else if (element->dimension == 2) {
            element->NodeU = doublevec(0, element->NumberOfNodes-1);
            element->NodeV = doublevec(0, element->NumberOfNodes-1);
            for (k=0; k<element->NumberOfNodes; k++) {
                element->NodeU[k] = [(defs.nodeU)[k] doubleValue];
                element->NodeV[k] = [(defs.nodeV)[k] doubleValue];
                basisTerms[k] = [(defs.basis)[k] intValue];
            }
        } else if (element->dimension == 3) {
            element->NodeU = doublevec(0, element->NumberOfNodes-1);
            element->NodeV = doublevec(0, element->NumberOfNodes-1);
            element->NodeW = doublevec(0, element->NumberOfNodes-1);
            for (k=0; k<element->NumberOfNodes; k++) {
                element->NodeU[k] = [(defs.nodeU)[k] doubleValue];
                element->NodeV[k] = [(defs.nodeV)[k] doubleValue];
                element->NodeW[k] = [(defs.nodeW)[k] doubleValue];
                basisTerms[k] = [(defs.basis)[k] intValue];
            }
        }
        element->GaussPoints = [(defs.gaussPoints)[0] intValue];
        if ([defs.gaussPoints count] > 1) element->GaussPoints2 = [(defs.gaussPoints)[1] intValue];
        if ([defs.gaussPoints count] > 2) element->GaussPoints0 = [(defs.gaussPoints)[2] intValue];
        
        element->StabilizationMK = [defs.stabilization doubleValue];
        
        if (element->GaussPoints2 <= 0) element->GaussPoints2 = element->GaussPoints;
        if (element->GaussPoints0 <= 0) element->GaussPoints0 = element->GaussPoints;
        
        // Reassign 0 to stabilization. Don't know why yet!!!
        element->StabilizationMK = 0.0;
        
        if (element->NodeV == NULL) {
            element->NodeV = doublevec(0, element->NumberOfNodes-1);
            memset( element->NodeV, 0.0, element->NumberOfNodes*sizeof(double) );
        }
        if (element->NodeW == NULL) {
            element->NodeW = doublevec(0, element->NumberOfNodes-1);
            memset( element->NodeW, 0.0, element->NumberOfNodes*sizeof(double) );
        }
        
        [self addDescriptionOfElement:element withBasisTerms:basisTerms];
    }

    free_ivector(basisTerms, 0, _maxDeg3-1);
}

-(int * __nonnull * __nonnull)getEdgeMap:(int)elementFamily {
    
    int **edgeMag = NULL;
    
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

/*****************************************************************************************
 
    Retrive element diameter parameter for stabilization
 
    Element_t *element  ->  element
    Nodes_t *nodes      ->  nodal coodinate arrays of the element
 
*****************************************************************************************/
-(double)elementDiameter:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes {
    
    int i, j, k, family;
    double x0, y0, z0, hk, a, s, cx, cy, cz;
    double j11, j12, j13, j21, j22, j23, g11, g12, g22;
    int **edgeMap = NULL, size=0;
    
    family = element->Type.ElementCode / 100;
    
    switch (family) {
        case 1:
            hk = 0.0;
            break;
        case 3: // Triangular element
            j11 = nodes->x[1] - nodes->x[0];
            j12 = nodes->y[1] - nodes->y[0];
            j13 = nodes->z[1] - nodes->z[0];
            j21 = nodes->x[2] - nodes->x[0];
            j22 = nodes->y[2] - nodes->y[0];
            j23 = nodes->z[2] - nodes->z[0];
            g11 = pow(j11, 2.0) + pow(j12, 2.0) + pow(j13, 2.0);
            g12 = j11*j21 + j12*j22 + j13*j23;
            g22 = pow(j21, 2.0) + pow(j22, 2.0) + pow(j23, 2.0);
            a = sqrt(g11*g22 - pow(g12, 2.0)) / 2.0;
            
            cx = ( nodes->x[0] + nodes->x[1] + nodes->x[2] ) / 3.0;
            cy = ( nodes->y[0] + nodes->y[1] + nodes->y[2] ) / 3.0;
            cz = ( nodes->z[0] + nodes->z[1] + nodes->z[2] ) / 3.0;
            
            s = pow((nodes->x[0]-cx), 2.0) + pow((nodes->y[0]-cy), 2.0) + pow((nodes->z[0]-cz), 2.0);
            s = s + pow((nodes->x[1]-cx), 2.0) + pow((nodes->y[1]-cy), 2.0) + pow((nodes->z[1]-cz), 2.0);
            s = s + pow((nodes->x[2]-cx), 2.0) + pow((nodes->y[2]-cy), 2.0) + pow((nodes->z[2]-cz), 2.0);
            
            hk = 16.0 * a * a / (3.0 * s);
            break;
        case 4: // Quadrilateral element
            cx = pow((nodes->x[1]-nodes->x[0]), 2.0) + pow((nodes->y[1]-nodes->y[0]), 2.0) + pow((nodes->z[1]-nodes->z[0]), 2.0);
            cy = pow((nodes->x[3]-nodes->x[0]), 2.0) + pow((nodes->y[3]-nodes->y[0]), 2.0) + pow((nodes->z[3]-nodes->z[0]), 2.0);
            hk = 2.0 * cx * cy / (cx + cy);
            break;
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
                x0 = nodes->x[j] - nodes->x[k];
                y0 = nodes->y[j] - nodes->y[k];
                z0 = nodes->z[j] - nodes->z[k];
                hk = min( hk, (pow(x0, 2.0)+pow(y0, 2.0)+pow(z0, 2.0)) );
            }
            break;
    }
    return sqrt(hk);
}

/*******************************************************************************************************
 
    Compute convection diffusion equation stabilization parameter for each and every element of the
    model by solving the largest eigenvalue of
 
    Lu = \lambda Gu,
    L = (\nabla^2 u,\nabla^2 w), G = (\nabla u,\nabla w)
 
*******************************************************************************************************/
-(void)computeStabilizationParameterInElement:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes mesh:(FEMMesh * __nonnull)mesh numberOfNodes:(int)n mk:(double * __nonnull)mk hk:(double * __nullable)hk {
 
    int i, j, p, q, t, dim;
    double eigr[n], s, ddp[3], ddq[3], dNodalBasisdx[n][n][3];
    double u, v, w, l[n-1][n-1], g[n-1][n-1], l_transpose[(n-1)*(n-1)], g_transpose[(n-1)*(n-1)], work[16*n], sum1, sum2;
    GaussIntegrationPoints *IP = NULL;
    FEMNumericIntegration *numericIntegration;
    BOOL stat;
    
    // Used for lapack routine
    int itype, order, lda, ldb, lwork, info;
    char *jobz, *uplo;
    
    numericIntegration = [[FEMNumericIntegration alloc] init];
    if ([numericIntegration allocation:mesh] == NO) fatal("FEMElementDescription:computeStabilizationParameter", "Allocation error in FEMNumericIntegration.");
    
    if (element->Type.BasisFunctionDegree <= 1) {
        switch (element->Type.ElementCode) {
            case 202:
            case 303:
            case 404:
            case 504:
            case 605:
            case 706:
                *mk = 1.0 / 3.0;
                break;
            case 808:
                *mk = 1.0 / 6.0;
        }
        if (hk != NULL) *hk = [self elementDiameter:element nodes:nodes];
        return;
    }
    
    memset(**dNodalBasisdx, 0.0, (n*n*3)*sizeof(double) );
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
                dNodalBasisdx[i][p][j] = numericIntegration.basisFirstDerivative[i][j];
            }
        }
    }
    
    dim = [mesh dimension];
    IP = GaussQuadrature(element, NULL, NULL);
    
    memset( *l, 0.0, ((n-1)*(n-1))*sizeof(double) );
    memset( *g, 0.0, ((n-1)*(n-1))*sizeof(double) );
    for (t=0; t<IP->n; t++) {
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];
        
        stat = [numericIntegration setMetricDeterminantForElement:element 
                                                     elementNodes:nodes inMesh:mesh 
                                             firstEvaluationPoint:u 
                                            secondEvaluationPoint:v 
                                             thirdEvaluationPoint:w];
        
        s = [numericIntegration metricDeterminant] * IP->s[t];
        
        for (p=1; p<n; p++) {
            for (q=1; q<n; q++) {
               memset( ddp, 0.0, sizeof(ddp) );
               memset( ddq, 0.0, sizeof(ddq) );
               for (i=0; i<dim; i++) {
                   g[p-1][q-1] = g[p-1][q-1] + s * numericIntegration.basisFirstDerivative[p][i] * numericIntegration.basisFirstDerivative[q][i];
                   sum1 = 0.0;
                   sum2 = 0.0;
                   for (j=0; j<n; j++) {
                       sum1 = sum1 + dNodalBasisdx[p][j][i] * numericIntegration.basisFirstDerivative[j][i];
                       sum2= sum2 + dNodalBasisdx[q][j][i] * numericIntegration.basisFirstDerivative[j][i];
                   }
                   ddp[i] = ddp[i] + sum1;
                   ddq[i] = ddq[i] + sum2;
               }
                vDSP_sveD(ddp, 1, &sum1, 3);
                vDSP_sveD(ddq, 1, &sum2, 3);
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
        *mk = 1.0 / 3.0;
        if (hk != NULL) {
            *hk = [self elementDiameter:element nodes:nodes];
        }
        return;
    }
    
    // We need to transpose our matrices before we pass them to LAPACK since LAPACK needs
    // two-dimensional arrays stored in column-order major.
    // TODO: dsygv assumes symetric matrices, do we really need to transpose?
    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1; j++) {
            l_transpose[j+(n-1)*i] = l[j][i];
            g_transpose[j+(n-1)*i] = g[j][i];
        }
    }
    
    itype = 1;
    jobz = "N";
    uplo = "U";
    order = n-1;
    lda = n-1;
    ldb = n-1;
    lwork = 12*n;
    dsygv_(&itype, jobz, uplo, &order, l_transpose, &lda, g_transpose, &ldb, eigr, work, &lwork, &info);
    if (info < 0 || info > 0) {
        fprintf(stderr, "FEMElementDescription:computeStabilizationParameterInElement: error in lapack routine dsygv. Error code: %d.\n", info);
        fatal("FEMElementDescription:computeStabilizationParameterInElement", "Saino will abort the simulation now...");
    }
    *mk = eigr[n-2];
        
    if (*mk < 10.0 * AEPS) {
        *mk = 1.0 / 3.0;
        if (hk != NULL) {
            *hk = [self elementDiameter:element nodes:nodes];
        }
        return;
    }
    
    if (hk != NULL) {
        *hk = sqrt( 2.0 / (*mk * element->Type.StabilizationMK) );
        *mk = min( 1.0/3.0, element->Type.StabilizationMK );
    } else {
        switch (element->Type.ElementCode / 100) {
            case 2:
            case 4:
            case 8:
                *mk = 4.0 * (*mk);
                break;
        }
        *mk = min( 1.0/3.0, 2.0/(*mk) );
    }
    
    [numericIntegration deallocation:mesh];
}

-(ElementType_t * __nullable)getElementType:(int)code inMesh:(FEMMesh * __nonnull)mesh stabilization:(BOOL * __nullable)computeStab {
        
    Element_t *elm = NULL;
    ElementType_t *element = NULL;
    Nodes_t *nodes = NULL;
    
    element = _elementTypeList;
    
    while (element != NULL) {
        if (code == element->ElementCode) break;
        element = element->NextElementType;
    }
        
    if (element == NULL) {
        fprintf(stdout, "FEMElementDescription:getElementType: element type code not found: %d.\n", code);
        fprintf(stdout, "FEMElementDescription:getElementType: ignoring element.\n");
        return NULL;
    }
        
    if (computeStab != NULL) {
        if (*computeStab == NO) return element;
    }
    
    if (element->StabilizationMK == 0.0) {
        elm = (Element_t*) malloc( sizeof(Element_t));
        initElements(elm, 1);
        elm->Type = *element;
        elm->BDOFs = 0;
        elm->DGDOFs = 0;
        elm->Pdefs = NULL;
        elm->DGIndexes = NULL;
        elm->FaceIndexes = NULL;
        elm->BubbleIndexes = NULL;
        nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
        initNodes(nodes);
        nodes->x = element->NodeU;
        nodes->y = element->NodeV;
        nodes->z = element->NodeW;
        [self computeStabilizationParameterInElement:elm nodes:nodes mesh:mesh numberOfNodes:element->NumberOfNodes mk:&element->StabilizationMK hk:NULL];
        free(elm);
        free(nodes);
    }
    
    return element;
}

/**********************************************************************************************
 
    Given element structure, return value of the first partial derivative with respect
    to local coordinate of a quantity x given at element nodes at local coordinates
    point u inside the element. Element basis functions are used to compute the value.
 
    Element_t *element  ->  element structure
    double u            ->  point at which to evaluate the partial derivative
    double *x           ->  nodal values of the quantity whose partial derivative is required
 
    Return y = @x/@u
 
**********************************************************************************************/
-(double)firstDerivative1DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evalutationPoint:(double)u {
    
    int i, n;
    int *p;
    double y;
    double s;
    double *coeff;
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (x[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                if (p[i] >= 1) {
                    s = s + p[i] * coeff[i] * pow(u, (p[i]-1));
                }
            }
            y = y + s * x[n];
        }
    }
    
    return y;
}

/**********************************************************************************************
 
    Given element structure, return value of the first partial derivative with respect to local
    coordinate u of i quantity x given at element nodes at local coordinate points u, v inside
    the element. Element basis functions are used to compute the value.
 
    Element_t *element  ->  element struncture
    double *x:          ->  nodal values of the quantity whose partial derivative is required
    double u,v          ->  points at which to evaluate the partial derivative
 
    Return first partial derivative in u of quantity x at point u,v y = @x(u,v)/@u
    
**********************************************************************************************/
-(double)firstDerivativeU2DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v {
    
    int i, n;
    int *p, *q;
    double y, s;
    double *coeff;
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (x[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                if (p[i] >= 1) {
                    s = s + p[i] * coeff[i] * pow(u, (p[i]-1)) * pow(v, q[i]);
                }
            }
            y = y + s * x[n];
        }
    }
    
    return y;
}

/**********************************************************************************************
 
    Given element structure, return value of the first partial derivative with respect to local
    coordinate v of i quantity x given at element nodes at local coordinate points u, v inside
    the element. Element basis functions are used to compute the value.
 
    Element_t *element  ->  element struncture
    double *x:          ->  nodal values of the quantity whose partial derivative is required
    double u,v          ->  points at which to evaluate the partial derivative
 
    Return first partial derivative in v of quantity x at point u,v y = @x(u,v)/@v
 
**********************************************************************************************/
-(double)firstDerivativeV2DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v {
    
    int i, n;
    int *p, *q;
    double y, s;
    double *coeff;
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (x[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                if (q[i] >= 1) {
                    s = s + q[i] * coeff[i] * pow(u, p[i]) * pow(v, (q[i]-1));
                }
            }
            y = y + s * x[n];
        }
    }
    
    return y;
}

/**********************************************************************************************
 
    Given element structure, return value of the first partial derivative with respect to
    local coordinate u of a quantity x given at element nodes at local coordinate point
    u, v, w inside the element. Element basis functions are used to compute the value.
 
    Element_t *element  ->  element struncture
    double *x:          ->  nodal values of the quantity whose partial derivative is required
    double u,v, w       ->  points at which to evaluate the partial derivative
 
    Return first partial derivative in u of quantity x at point u,v,w y = @x(u,v,w)/@u
 
**********************************************************************************************/
-(double)firstDerivativeU3DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w {
    
    int i, n;
    int *p, *q, *r;
    double y, s;
    double *coeff;
    
    if (element->Type.ElementCode == 605) {
        if (w == 1) w = 1.0 - 1.0e-12;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * ( -(1.0-v) + v*w * s) / 4.0;
        y = y + x[1] * (  (1.0-v) - v*w * s) / 4.0;
        y = y + x[2] * (  (1.0+v) + v*w * s) / 4.0;
        y = y + x[3] * ( -(1.0+v) - v*w * s) / 4.0;
        return  y;
    } else if (element->Type.ElementCode == 613) {
        if (w == 1) w = 1.0 - 1.0e-12;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * ( -( (1.0-u) * (1.0-v) - w + u*v*w * s ) + (-u-v-1.0) * ( -(1.0-v) + v*w *s ) ) / 4.0;
        
        y = y + x[1] * (  ( (1.0+u) * (1.0-v) - w - u*v*w * s ) + (u-v-1.0) * ( (1.0-v) - v*w *s ) ) / 4.0;
        
        y = y + x[2] * (  ( (1.0+u) * (1.0+v) - w + u*v*w * s ) + (u+v-1.0) * ( (1.0+v) + v*w *s ) ) / 4.0;
        
        y = y + x[3] * ( -( (1.0-u) * (1.0+v) - w - u*v*w * s ) + (-u+v-1.0) * ( -(1.0+v) - v*w *s ) ) / 4.0;
        
        y = y + x[4] * 0.0;
        
        y = y + x[5] * ( (1.0-u-w)*(1.0-v-w) - (1.0+u-w)*(1.0-v-w) ) * s / 2.0;
        y = y + x[6] * ( (1.0+v-w)*(1.0-v-w) ) * s / 2.0;
        y = y + x[7] * ( (1.0-u-w)*(1.0+v-w) - (1.0+u-w)*(1.0+v-w) ) * s  / 2.0;
        y = y + x[8] * ( -(1.0+v-w)*(1.0-v-w) ) * s / 2.0;
        
        y = y - x[9]  * w * (1.0-v-w) * s;
        y = y + x[10] * w * (1.0-v-w) * s;
        y = y + x[11] * w * (1.0+v-w) * s;
        y = y - x[12] * w * (1.0+v-w) * s;
        
        return y;
    }
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (x[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            r = element->Type.BasisFunctions[n].r;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                if (p[i] >= 1) {
                    s = s + p[i] * coeff[i] * pow(u, (p[i]-1.0)) * pow(v, q[i]) * pow(w, r[i]);
                }
            }
            y = y + s * x[n];
        }
    }
    
    return y;
}

/**********************************************************************************************
 
    Given element structure, return value of the first partial derivative with respect to
    local coordinate v of a quantity x given at element nodes at local coordinate point
    u, v, w inside the element. Element basis functions are used to compute the value.
 
    Element_t *element  ->  element struncture
    double *x:          ->  nodal values of the quantity whose partial derivative is required
    double u,v, w       ->  points at which to evaluate the partial derivative
 
    Return first partial derivative in v of quantity x at point u,v,w y = @x(u,v,w)/@v
 
**********************************************************************************************/
-(double)firstDerivativeV3DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w {
    
    int i, n;
    int *p, *q, *r;
    double y, s;
    double *coeff;
    
    if (element->Type.ElementCode == 605) {
        if (w == 1) w = 1.0 - 1.0e-12;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * ( -(1.0-u) + u*w *s ) / 4.0;
        y = y + x[1] * ( -(1.0+u) - u*w *s ) / 4.0;
        y = y + x[2] * (  (1.0+u) + u*w *s ) / 4.0;
        y = y + x[3] * (  (1.0-u) - u*w *s ) / 4.0;
        return y;
    } else if (element->Type.ElementCode == 613) {
        if (w == 1) w = 1.0 - 1.0e-12;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * ( -( (1.0-u) * (1.0-v) - w + u*v*w * s ) + (-u-v-1.0) * ( -(1.0-u) + u*w *s ) ) / 4.0;
        
        y = y + x[1] * ( -( (1.0+u) * (1.0-v) - w - u*v*w * s ) + (u-v-1.0) * ( -(1.0+u) - u*w *s ) ) / 4.0;
        
        y = y + x[2] * (  ( (1.0+u) * (1.0+v) - w + u*v*w * s ) + (u+v-1.0) * ( (1.0+u) + u*w *s ) ) / 4.0;
        
        y = y + x[3] * (  ( (1.0-u) * (1.0+v) - w - u*v*w * s ) + (-u+v-1.0) * ( (1.0-u) - u*w *s ) ) / 4.0;
        
        y = y + x[4] * 0.0;
        
        y = y - x[5] * (1.0+u-w)*(1.0-u-w) * s / 2.0;
        y = y + x[6] * ( (1.0-v-w)*(1.0+u-w) - (1.0+v-w)*(1.0+u-w) ) * s / 2.0;
        y = y + x[7] * (1.0+u-w)*(1.0-u-w) * s / 2.0;
        y = y + x[8] * ( (1.0-v-w)*(1.0-u-w) - (1.0+v-w)*(1.0-u-w) ) * s / 2.0;
        
        y = y - x[9]  * w * (1.0-u-w) * s;
        y = y - x[10] * w * (1.0+u-w) * s;
        y = y + x[11] * w * (1.0+u-w) * s;
        y = y + x[12] * w * (1.0-u-w) * s;
        
        return y;
    }
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (x[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            r = element->Type.BasisFunctions[n].r;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                if (q[i] >= 1) {
                    s = s + q[i] * coeff[i] * pow(u, p[i]) * pow(v, (q[i]-1.0)) * pow(w, r[i]);
                }
            }
            y = y + s * x[n];
        }
    }

    return y;
}

/**********************************************************************************************
 
    Given element structure, return value of the first partial derivative with respect to
    local coordinate v of a quantity x given at element nodes at local coordinate point
    u, v, w inside the element. Element basis functions are used to compute the value.
 
    Element_t *element  ->  element struncture
    double *x:          ->  nodal values of the quantity whose partial derivative is required
    double u,v, w       ->  points at which to evaluate the partial derivative
 
    Return first partial derivative in w of quantity x at point u,v,w y = @x(u,v,w)/@w
 
    Method corresponds to Elmer from git on October 27 2015
 
**********************************************************************************************/
-(double)firstDerivativeW3DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w {
    
    int i, n;
    int *p, *q, *r;
    double y, s;
    double *coeff;
    
    if (element->Type.ElementCode == 605) {
        if (w == 1) w = 1.0 - 1.0e-12;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * ( -1.0 + u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[1] * ( -1.0 - u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[2] * ( -1.0 + u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[3] * ( -1.0 - u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[4];
        return y;
    } else if (element->Type.ElementCode == 613) {
        if (w == 1) w = 1.0 - 1.0e-12;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * (-u-v-1.0) * ( -1.0 + u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[1] * (u-v-1.0)  * ( -1.0 - u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[2] * (u+v-1.0)  * ( -1.0 + u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[3] * (-u+v-1.0) * ( -1.0 - u*v*pow(s, 2.0) ) / 4.0;
        
        y = y + x[4] * (4.0*w-1.0);
        
        y = y + x[5] * ( ( -(1.0-u-w)*(1.0-v-w) - (1.0+u-w)*(1.0-v-w) - (1.0+u-w)*(1.0-u-w) ) * s + (1.0+u-w)*(1.0-u-w)*(1.0-v-w) * pow(s, 2.0) ) / 2.0;
        
        y = y + x[6] * ( ( -(1.0-v-w)*(1.0+u-w) - (1.0+v-w)*(1.0+u-w) - (1.0+v-w)*(1.0-v-w) ) * s + (1.0+v-w)*(1.0-v-w)*(1.0+u-w) * pow(s, 2.0) ) / 2.0;
        
        y = y + x[7] * ( ( -(1.0-u-w)*(1.0+v-w) - (1.0+u-w)*(1.0+v-w) - (1.0+u-w)*(1.0-u-w) ) * s + (1.0+u-w)*(1.0-u-w)*(1.0+v-w) * pow(s, 2.0) ) / 2.0;
        
        y = y + x[8] * ( ( -(1.0-v-w)*(1.0-u-w) - (1.0+v-w)*(1.0-u-w) - (1.0+v-w)*(1.0-v-w) ) * s + (1.0+v-w)*(1.0-v-w)*(1.0-u-w) * pow(s, 2.0) ) / 2.0;
        
        y = y + x[9] * ( ( (1.0-u-w)*(1.0-v-w) - w * (1.0-v-w) - w * (1.0-u-w) ) * s + w * (1.0-u-w) * (1.0-v-w) * pow(s, 2.0) );
        
        y = y + x[10] * ( ( (1.0+u-w)*(1.0-v-w) - w * (1.0-v-w) - w * (1.0+u-w) ) * s + w * (1.0+u-w) * (1.0-v-w) * pow(s, 2.0) );
        
        y = y + x[11] * ( ( (1.0+u-w)*(1.0+v-w) - w * (1.0+v-w) - w * (1.0+u-w) ) * s + w * (1.0+u-w) * (1.0+v-w) * pow(s, 2.0) );
        
        y = y + x[12] * ( ( (1.0-u-w)*(1.0+v-w) - w * (1.0+v-w) - w * (1.0-u-w) ) * s + w * (1.0-u-w) * (1.0+v-w) * pow(s, 2.0) );
        
        return y;
    }
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (x[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            r = element->Type.BasisFunctions[n].r;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                if (r[i] >= 1) {
                    s = s + r[i] * coeff[i] * pow(u, p[i]) * pow(v, q[i]) * pow(w, (r[i]-1.0));
                }
            }
            y = y + s * x[n];
        }
    }

    return y;
}

/**********************************************************************************************

    Given element structure, return value of a quantity x given at element nodes at local
    coordinates point u (1D), (u,v) (2D), (u,v,w) (3D) inside the element. Element basis
    functions are used to compute the value. This wrapper method will call the corresponding
    class private methods according to the element dimension.
 
    Element_t *element  ->  element structure
    double *x           ->  nodal values of the quantity whose partial derivative is required
    double u, v, w      ->  points at which to evaluate the partial derivative
    double *basis       ->  Values of the basis functions at the point (u,v,z) can be given if
                            known, otherwise they will be computed from the definition.
 
    Return y = x(u,v,w)

**********************************************************************************************/
-(double)interpolateInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w withBasis:(double * __nullable)basis {
    
    int n;
    double value=0.0;
    
    if (basis != NULL) {
        // Basis function values given, just sum the result
        n = element->Type.NumberOfNodes;
        return value = cblas_ddot(n, x, 1, basis, 1);
    } else {
        // Otherwise compute from the definition
        switch (element->Type.dimension) {
            case 0:
                value = x[0];
                break;
            case 1:
                value = [self FEMElementDescription_interpolate1DInElement:element nodalValues:x evalutationPoint:u];
                break;
            case 2:
                value = [self interpolate2DInElement:element nodalValues:x evaluatedAt:u andAt:v];
                break;
            case 3:
                value = [self interpolate3DInElement:element nodalValues:x evaluatedAt:u andAt:v andAt:w];
                break;
        }
        return value;
    }
}

/***************************************************************************************
 
    Normals will point from more dense material to less dense material or outwards, if
    no elements of the other side.
 
***************************************************************************************/
-(void)checkNormalDirectionInBDElement:(Element_t * __nonnull)boundary forNormals:(double * __nonnull)normals mesh:(FEMMesh * __nonnull)mesh x:(double)x y:(double)y z:(double)z turn:(BOOL * __nullable)turn {
    
    int i, n, k;
    double x1, y1, z1;
    double *nx, *ny, *nz;
    Element_t *element = NULL, *rightElement = NULL, *leftElement = NULL;
    Nodes_t *nodes = NULL;
    
    if (boundary->BoundaryInfo == NULL) return;
    
    nodes = mesh.getNodes;
    
    k = boundary->BoundaryInfo->Outbody;
    
    leftElement = boundary->BoundaryInfo->Left;
    
    if (leftElement != NULL) {
        rightElement = boundary->BoundaryInfo->Right;
        
        if (rightElement != NULL) {
            if (k > 0) {
                if (leftElement->BodyID == k) {
                    element = rightElement;
                } else {
                    element = leftElement;
                }
            } else {
                if (leftElement->BodyID > rightElement->BodyID) {
                    element = leftElement;
                } else {
                    element = rightElement;
                }
            }
        } else {
            element = leftElement;
        }
    } else {
        element = boundary->BoundaryInfo->Right;
    }
    
    if (element == NULL) return;
    
    n = element->Type.NumberOfNodes;
    
    nx = doublevec(0, n-1);
    ny = doublevec(0, n-1);
    nz = doublevec(0, n-1);
    
    for (i=0; i<n; i++) {
        nx[i] = nodes->x[element->NodeIndexes[i]];
        ny[i] = nodes->y[element->NodeIndexes[i]];
        nz[i] = nodes->z[element->NodeIndexes[i]];
    }
    
    x1 = 0.0; y1 = 0.0; z1 = 0.0;
    switch (element->Type.ElementCode / 100) {
        case 2:
        case 4:
        case 8:
            x1 = [self interpolateInElement:element nodalValues:nx evaluatedAt:0.0 andAt:0.0 andAt:0.0 withBasis:NULL];
            y1 = [self interpolateInElement:element nodalValues:ny evaluatedAt:0.0 andAt:0.0 andAt:0.0 withBasis:NULL];
            z1 = [self interpolateInElement:element nodalValues:nz evaluatedAt:0.0 andAt:0.0 andAt:0.0 withBasis:NULL];
            break;
        case 3:
            x1 = [self interpolateInElement:element nodalValues:nx evaluatedAt:1.0/3.0 andAt:1.0/3.0 andAt:0.0 withBasis:NULL];
            y1 = [self interpolateInElement:element nodalValues:ny evaluatedAt:1.0/3.0 andAt:1.0/3.0 andAt:0.0 withBasis:NULL];
            z1 = [self interpolateInElement:element nodalValues:nz evaluatedAt:1.0/3.0 andAt:1.0/3.0 andAt:0.0 withBasis:NULL];
            break;
        case 5:
            x1 = [self interpolateInElement:element nodalValues:nx evaluatedAt:1.0/4.0 andAt:1.0/4.0 andAt:1.0/4.0 withBasis:NULL];
            y1 = [self interpolateInElement:element nodalValues:ny evaluatedAt:1.0/4.0 andAt:1.0/4.0 andAt:1.0/4.0 withBasis:NULL];
            z1 = [self interpolateInElement:element nodalValues:nz evaluatedAt:1.0/4.0 andAt:1.0/4.0 andAt:1.0/4.0 withBasis:NULL];
            break;
        case 6:
            x1 = [self interpolateInElement:element nodalValues:nx evaluatedAt:0.0 andAt:0.0 andAt:1.0/3.0 withBasis:NULL];
            y1 = [self interpolateInElement:element nodalValues:ny evaluatedAt:0.0 andAt:0.0 andAt:1.0/3.0 withBasis:NULL];
            z1 = [self interpolateInElement:element nodalValues:nz evaluatedAt:0.0 andAt:0.0 andAt:1.0/3.0 withBasis:NULL];
            break;
        case 7:
            x1 = [self interpolateInElement:element nodalValues:nx evaluatedAt:1.0/3.0 andAt:1.0/3.0 andAt:0.0 withBasis:NULL];
            y1 = [self interpolateInElement:element nodalValues:ny evaluatedAt:1.0/3.0 andAt:1.0/3.0 andAt:0.0 withBasis:NULL];
            z1 = [self interpolateInElement:element nodalValues:nz evaluatedAt:1.0/3.0 andAt:1.0/3.0 andAt:0.0 withBasis:NULL];
            break;
        default:
            fatal("FEMElementDescription:checkNormalDirectionInBDElement", "Invalid element code for parent element.");
            break;
    }
    
    x1 = x1 - x;
    y1 = y1 - y;
    z1 = z1 - z;
    
    if (turn != NULL) *turn = NO;
    if (x1*normals[0] + y1*normals[1] + z1*normals[2] > 0) {
        if (element->BodyID != k) {
            vDSP_vnegD(normals, 1, normals, 1, 3);
            if (turn != NULL) *turn = YES;
        }
    } else if (element->BodyID == k) {
        vDSP_vnegD(normals, 1, normals, 1, 3);
        if (turn != NULL) *turn = YES;
    }
    
    free_dvector(nx, 0, n-1);
    free_dvector(ny, 0, n-1);
    free_dvector(nz, 0, n-1);
}

/***************************************************************************************
 
    Gives the normal vector of a boundary element.
    For non-curved elements, the normal vector does not depend on the local coordinate
    while otherwise it does. There are different uses of the function where some do 
    not have the luxury of knowing the local coordinates and hence the center point is
    used as default.
 
***************************************************************************************/
-(void)normalVectorForBDElement:(Element_t * __nonnull)boundary boundaryNodes:(Nodes_t * __nonnull)nodes mesh:(FEMMesh * __nonnull)mesh paraU:(double * __nullable)u0 paraV:(double * __nullable)v0 check:(BOOL * __nullable)check normals:(double * __nonnull)normals {
    
    double u, v, auu, auv, avv, detA, x=0.0, y=0.0, z=0.0;
    double dxdu, dxdv, dydu, dydv, dzdu, dzdv;
    BOOL doCheck;
    
    switch (boundary->Type.dimension) {
        case 0:
            normals[0] = 1.0;
            normals[1] = 0.0;
            normals[2] = 0.0;
            break;
        case 1:
            if (u0 != NULL) {
                u = *u0;
            } else u = 0.0;
            
            dxdu = [self firstDerivative1DInElement:boundary nodalValues:nodes->x evalutationPoint:u];
            dydu = [self firstDerivative1DInElement:boundary nodalValues:nodes->y evalutationPoint:u];
            
            detA = dxdu*dxdu + dydu*dydu;
            if (detA <= 0.0) {
                memset( normals, 0.0, 3*sizeof(double) );
                return;
            }
            detA = 1.0 / sqrt(detA);
            normals[0] = -dydu * detA;
            normals[1] = dxdu * detA;
            normals[2] = 0.0;
            break;
        case 2:
            if (u0 != NULL) {
                u = *u0;
                v = *v0;
            } else {
                if (boundary->Type.ElementCode / 100 == 3) {
                    u = 1.0 / 3.0;
                    v = 1.0 / 3.0;
                } else {
                    u = 0.0;
                    v = 0.0;
                }
            }
            
            dxdu = [self firstDerivativeU2DInElement:boundary nodalValues:nodes->x evaluatedAt:u andAt:v];
            dydu = [self firstDerivativeU2DInElement:boundary nodalValues:nodes->y evaluatedAt:u andAt:v];
            dzdu = [self firstDerivativeU2DInElement:boundary nodalValues:nodes->z evaluatedAt:u andAt:v];
            
            dxdv = [self firstDerivativeV2DInElement:boundary nodalValues:nodes->x evaluatedAt:u andAt:v];
            dydv = [self firstDerivativeV2DInElement:boundary nodalValues:nodes->y evaluatedAt:u andAt:v];
            dzdv = [self firstDerivativeV2DInElement:boundary nodalValues:nodes->z evaluatedAt:u andAt:v];
            
            auu = dxdu*dxdu + dydu*dydu + dzdu*dzdu;
            auv = dxdu*dxdv + dydu*dydv + dzdu*dzdv;
            avv = dxdv*dxdv + dydv*dydv + dzdv*dzdv;
            
            detA = 1.0 / sqrt(auu*avv - auv*auv);
            
            normals[0] = (dydu * dzdv - dydv * dzdu) * detA;
            normals[1] = (dxdv * dzdu - dxdu * dzdv) * detA;
            normals[2] = (dxdu * dydv - dxdv * dydu) * detA;
            break;
        default:
            fatal("FEMElementDescription:nomalVectorForBDElement", "Invalid dimension for determining normal.");
            break;
    }
    
    doCheck = NO;
    if (check != NULL) doCheck = *check;
    
    if (doCheck == YES) {
        switch (boundary->Type.ElementCode / 100) {
            case 1:
                x = nodes->x[0];
                y = nodes->x[0];
                z = nodes->z[0];
                break;
            case 2:
            case 4:
                x = [self interpolateInElement:boundary nodalValues:nodes->x evaluatedAt:0.0 andAt:0.0 andAt:0.0 withBasis:NULL];
                y = [self interpolateInElement:boundary nodalValues:nodes->y evaluatedAt:0.0 andAt:0.0 andAt:0.0 withBasis:NULL];
                z = [self interpolateInElement:boundary nodalValues:nodes->z evaluatedAt:0.0 andAt:0.0 andAt:0.0 withBasis:NULL];
                break;
            case 3:
                x = [self interpolateInElement:boundary nodalValues:nodes->x evaluatedAt:1.0/3.0 andAt:1.0/3.0 andAt:0.0 withBasis:NULL];
                y = [self interpolateInElement:boundary nodalValues:nodes->y evaluatedAt:1.0/3.0 andAt:1.0/3.0 andAt:0.0 withBasis:NULL];
                z = [self interpolateInElement:boundary nodalValues:nodes->z evaluatedAt:1.0/3.0 andAt:1.0/3.0 andAt:0.0 withBasis:NULL];
                break;
        }
        [self checkNormalDirectionInBDElement:boundary forNormals:normals mesh:mesh x:x y:y z:z turn:NULL];
    }
}

/***************************************************************************************************
 
    Convert global coordinates x, y, z, inside element to local coordinates u, v, w, of the element
    TODO: Needs support for p elements
 
    Method corresponds to Elmer from git on October 27 2015
 
***************************************************************************************************/
-(void)globalToLocalFromElement:(Element_t * __nonnull)element elementNodes:(Nodes_t * __nonnull)nodes localU:(double * __nonnull)u localV:(double * __nonnull)v localW:(double * __nonnull)w x:(double)x y:(double)y z:(double)z model:(FEMModel * __nonnull)model {
    
    int i, j, n;
    int const maxIter = 50;
    double r=0.0, s=0.0, t=0.0, buffer[3], delta[3], prevdelta[3], J[3][3], JJ[2][2], bf[3], det, acc, err=0.0, sum, yy[3];
    BOOL converged;
    FEMUtilities *utilities;

    *u = 0.0;
    *v = 0.0;
    *w = 0.0;
    if (element->Type.dimension == 0) return;
    
    utilities = [[FEMUtilities alloc] init];
    
    n = element->Type.NumberOfNodes;
    acc = DBL_EPSILON;
    converged = NO;
    
    double alpha = 0.5;
    for (i=1; i<=maxIter; i++) {
        
        r = [self interpolateInElement:element nodalValues:nodes->x evaluatedAt:*u andAt:*v andAt:*w withBasis:NULL] - x;
        s = [self interpolateInElement:element nodalValues:nodes->y evaluatedAt:*u andAt:*v andAt:*w withBasis:NULL] - y;
        t = [self interpolateInElement:element nodalValues:nodes->z evaluatedAt:*u andAt:*v andAt:*w withBasis:NULL] - z;
        
        err = pow(r, 2.0) + pow(s, 2.0) + pow(t, 2.0);
        
        if (err < acc) {
            converged = YES;
            break;
        }
        
        memcpy(prevdelta, delta, sizeof(delta));
        memset( delta, 0.0, sizeof(delta) );
        memset( *J, 0.0, (3*3)*sizeof(double) );
        
        switch (element->Type.dimension) {
            case 1:
                J[0][0] = [self firstDerivative1DInElement:element nodalValues:nodes->x evalutationPoint:*u];
                J[1][0] = [self firstDerivative1DInElement:element nodalValues:nodes->y evalutationPoint:*u];
                J[2][0] = [self firstDerivative1DInElement:element nodalValues:nodes->z evalutationPoint:*u];
                
                det = 0.0;
                for (j=0; j<3; j++) {
                    det = det + pow(J[j][0], 2.0);
                }
                delta[0] = (r * J[0][0] + s * J[1][0] + t * J[2][0]) / det;
                break;
                
            case 2:
                J[0][0] = [self firstDerivativeU2DInElement:element nodalValues:nodes->x evaluatedAt:*u andAt:*v];
                J[0][1] = [self firstDerivativeV2DInElement:element nodalValues:nodes->x evaluatedAt:*u andAt:*v];
                J[1][0] = [self firstDerivativeU2DInElement:element nodalValues:nodes->y evaluatedAt:*u andAt:*v];
                J[1][1] = [self firstDerivativeV2DInElement:element nodalValues:nodes->y evaluatedAt:*u andAt:*v];
                
                switch (model.dimension) {
                    case 3:
                        J[2][0] = [self firstDerivativeU2DInElement:element nodalValues:nodes->z evaluatedAt:*u andAt:*v];
                        J[2][1] = [self firstDerivativeV2DInElement:element nodalValues:nodes->z evaluatedAt:*u andAt:*v];
                        
                        delta[0] = r;
                        delta[1] = s;
                        delta[2] = t;
                        memset( yy, 0.0, sizeof(yy) );
                        cblas_dgemv(CblasRowMajor, CblasTrans, 3, 3, 1.0, (double *)J, 3, delta, 1, 0.0, yy, 1);
                        r = yy[0];
                        s = yy[1];
                        memset( *JJ, 0.0, (2*2)*sizeof(double) );
                        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 2, 2, 3, 1.0, (double *)J, 3, (double *)J, 3, 0.0, (double *)JJ, 2);
                        delta[2] = 0.0;
                        break;
                }
                bf[0] = r;
                bf[1] = s;
                bf[2] = 0.0;
                [utilities solveLinearSystem2x2:JJ afterSolve:delta rightHandSide:bf];
                break;
                
            case 3:
                J[0][0] = [self firstDerivativeU3DInElement:element nodalValues:nodes->x evaluatedAt:*u andAt:*v andAt:*w];
                J[0][1] = [self firstDerivativeV3DInElement:element nodalValues:nodes->x evaluatedAt:*u andAt:*v andAt:*w];
                J[0][2] = [self firstDerivativeW3DInElement:element nodalValues:nodes->x evaluatedAt:*u andAt:*v andAt:*w];
                
                J[1][0] = [self firstDerivativeU3DInElement:element nodalValues:nodes->y evaluatedAt:*u andAt:*v andAt:*w];
                J[1][1] = [self firstDerivativeV3DInElement:element nodalValues:nodes->y evaluatedAt:*u andAt:*v andAt:*w];
                J[1][2] = [self firstDerivativeW3DInElement:element nodalValues:nodes->y evaluatedAt:*u andAt:*v andAt:*w];
                
                J[2][0] = [self firstDerivativeU3DInElement:element nodalValues:nodes->z evaluatedAt:*u andAt:*v andAt:*w];
                J[2][1] = [self firstDerivativeV3DInElement:element nodalValues:nodes->z evaluatedAt:*u andAt:*v andAt:*w];
                J[2][2] = [self firstDerivativeW3DInElement:element nodalValues:nodes->z evaluatedAt:*u andAt:*v andAt:*w];
                
                bf[0] = r;
                bf[1] = s;
                bf[2] = t;
                [utilities solveLinearSystem3x3:J afterSolve:delta rightHandSide:bf];
                break;
        }
        
        if (i > 10) {
            // If some values is suggested over and over again, then exit.
            // This may be a sign that the node is off-plane and can not
            // be described within the element.
            memset(buffer, 0.0, sizeof(buffer));
            vDSP_vsubD(prevdelta, 1, delta, 1, buffer, 1, 3);
            vDSP_vabsD(buffer, 1, buffer, 1, 3);
            vDSP_sveD(buffer, 1, &sum, 3);
            if (sum < acc) break;
            
            // Use sloppier criteria when iteration still unsuccesfull
            if (i > 20) {
                if (sum < sqrt(acc)) break;
            }
            
            // If the iteration does not proceed try with some relaxation
            vDSP_vsmulD(delta, 1, &alpha, delta, 1, 3);
        }
        
        *u = *u - delta[0];
        *v = *v - delta[1];
        *w = *w - delta[2];
    }
    
    if (converged == NO) {
        if (err > sqrt(acc)) {
            if (i > maxIter) {
                fprintf(stdout, "FEMElementDescription:globalToLocalFromElement: did not converge.\n");
                fprintf(stdout, "rst, %d %f, %f, %f\n", i, r, s, t);
                fprintf(stdout, "err, %d, %f, %f, %f\n", i, err, acc, sqrt(acc));
                fprintf(stdout, "delta, %f %f %f\n", delta[0], delta[1], delta[2]);
                fprintf(stdout, "prevdelta, %f %f %f\n", prevdelta[0], prevdelta[1], prevdelta[2]);
                fprintf(stdout, "code, %d\n", element->Type.ElementCode);
                fprintf(stdout, "x: %f\n", x);
                for (i=0; i<element->Type.NumberOfNodes; i++) {
                    fprintf(stdout, "%f ", nodes->x[i]);
                }
                fprintf(stdout, "\n");
                fprintf(stdout, "y: %f\n", y);
                for (i=0; i<element->Type.NumberOfNodes; i++) {
                    fprintf(stdout, "%f ", nodes->y[i]);
                }
                fprintf(stdout, "\n");
                fprintf(stdout, "y: %f\n", z);
                for (i=0; i<element->Type.NumberOfNodes; i++) {
                    fprintf(stdout, "%f ", nodes->z[i]);
                }
                fprintf(stdout, "\n");
            } else {
                fprintf(stdout, "FEMElementDescription:globalToLocalFromElement: node may be out of element.");
                fprintf(stdout, "rst, %f, %f, %f, %f\n", r, s, t, DBL_EPSILON);
            }
        }
    }
}

/************************************************************************************************
 
    Given element structure, return value of a quantity x given at element nodes at local
    cooidinate points (u,v) inside the element. Element basis functions are used to compute the
    value. Used for 2d elements
 
    Element_t *element  ->  element structure
    double *x           ->  nodal values of the quantity whose partial derivative is required
    double u, v         ->  points at which to evaluate the partial derivative
 
    Return y = x(u,v)
 
************************************************************************************************/
-(double)interpolate2DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v {
    
    int i, n;
    int *p, *q;
    double y, s;
    double *coeff;
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (x[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                s = s + coeff[i] * pow(u, p[i]) * pow(v, q[i]);
            }
            y = y + s * x[n];
        }
    }
    
    return y;
}

/**********************************************************************************************
 
    Given element structure, return value of a quantity x given at element nodes at local
    cooidinate points (u,v,w) inside the element. Element basis functions are used to compute
    the value. Used for 3d elements
 
    Element_t *element  ->  element structure
    double *x           ->  nodal values of the quantity whose partial derivative is required
    double u, v, w      ->  points at which to evaluate the partial derivative
 
    Return y = x(u,v,w)
 
**********************************************************************************************/
-(double)interpolate3DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w {
    
    int i, n;
    int *p, *q, *r;
    double y, s;
    double *coeff;
    
    if (element->Type.ElementCode == 605) {
        s = 0.0;
        if (w == 1.0) w = 1.0 - 1.0e-12;
        s = 1.0 / (1-w);
        
        y = 0.0;
        y = y + x[0] * ( (1.0-u) * (1.0-v) - w + u*v*w * s) / 4.0;
        y = y + x[1] * ( (1.0+u) * (1.0-v) - w - u*v*w * s) / 4.0;
        y = y + x[2] * ( (1.0+u) * (1.0+v) - w + u*v*w * s) / 4.0;
        y = y + x[3] * ( (1.0-u) * (1.0+v) - w - u*v*w * s) / 4.0;
        y = y + x[4] * w;
        return y;
    } else if (element->Type.ElementCode == 613) {
        if (w == 1.0) w = 1.0 - 1.0e-12;
        s = 1.0 / (1-w);
        
        y = 0.0;
        y = y + x[0] * (-u-v-1.0) * ( (1.0-u) * (1.0-v) - w + u*v*w * s ) / 4.0;
        y = y + x[1] * (u-v-1.0) * ( (1.0+u) * (1.0-v) - w - u*v*w * s ) / 4.0;
        y = y + x[2] * (u+v-1.0) * ( (1.0+u) * (1.0+v) - w + u*v*w * s ) / 4.0;
        y = y + x[3] * (-u+v-1.0) * ( (1.0-u) * (1.0+v) - w - u*v*w * s ) / 4.0;
        y = y + x[4] * w*(2.0*w-1.0);
        y = y + x[5] * (1.0+u-w) * (1.0-u-w) * (1.0-v-w) * s /2.0;
        y = y + x[6] * (1.0+v-w) * (1.0-v-w) * (1.0+u-w) * s /2.0;
        y = y + x[7] * (1.0+u-w) * (1.0-u-w) * (1.0+v-w) * s /2.0;
        y = y + x[8] * (1.0+v-w) * (1.0-v-w) * (1.0-u-w) * s /2.0;
        y = y + x[9] * w * (1.0-u-w) * (1.0-v-w) * s;
        y = y + x[10] * w * (1.0+u-w) * (1.0-v-w) * s;
        y = y + x[11] * w * (1.0+u-w) * (1.0+v-w) * s;
        y = y + x[12] * w * (1.0-u-w) * (1.0+v-w) * s;
        return y;
    }
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (x[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            r = element->Type.BasisFunctions[n].r;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                s = s + coeff[i] * pow(u, p[i]) * pow(v, q[i]) * pow(w, r[i]);
            }
            y = y + s * x[n];
        }
    }
    
    return y;
}

/************************************************************
 
    Method corresponds to Elmer from git on October 27 2015

************************************************************/
-(void)invertMatrix3x3:(double[][3])G inverted:(double[][3])GI detG:(double)detG {
    
    double s = 1.0 / detG;
    
    GI[0][0] = s * ( G[1][1]*G[2][2] - G[2][1]*G[1][2] );
    GI[1][0] = -s * ( G[1][0]*G[2][2] - G[2][0]*G[1][2] );
    GI[2][0] = s * ( G[1][0]*G[2][1] - G[2][0]*G[1][1] );
    
    GI[0][1] = -s * ( G[0][1]*G[2][2] - G[2][1]*G[0][2] );
    GI[1][1] = s * ( G[0][0]*G[2][2] - G[2][0]*G[0][2] );
    GI[2][1] = -s * ( G[0][0]*G[2][1] - G[2][0]*G[0][1] );
    
    GI[0][2] = s * ( G[0][1]*G[1][2] - G[1][1]*G[0][2] );
    GI[1][2] = -s * ( G[0][0]*G[1][2] - G[1][0]*G[0][2] );
    GI[2][2] = s * ( G[0][0]*G[1][1] - G[1][0]*G[0][1] );
}

/****************************************************************************************************************
 
    This subroutine contains an older design for providing edge element basis functions of the lowest-degree. 
    Obtaining optimal accuracy with these elements may require that the element map is affine, while the edge 
    basis functions given by the newer design (the function EdgeElementInfo) should also work on general meshes.
 
    Method corresponds to Elmer from git on October 27 2015
 
****************************************************************************************************************/
-(void)getEdgeBasisElement:(Element_t * __nonnull)element wBasis:(double * __nonnull * __nonnull)wBasis rotWBasis:(double * __nonnull * __nonnull)rotWBasis basis:(double * __nonnull)basis dBasisdx:(double * __nonnull * __nonnull)dBasisdx {
 
    int j, k, n=0, nk, nj, **edgeMap = NULL, size=0;
    double base=0, curlBasis[8][3], dBase[3], dTriBase[3][3], detF=0, detG, du[3], dudx[3][3], dudx1[3], dudx2[3], dudx3[3], edgeBasis[8][3], f[3][3], g[3][3], g1[3], g2[3], g3[3], rBase[3], sum, tbase[3], triBase[3], u=0, v=0, w=0;
    
    //TODO: add support for parallel run
    
    if (element->Type.BasisFunctionDegree > 1) {
        fprintf(stderr, "FEMElementDescription:getEdgeBasisElement: can only handle linear elements, sorry.\n");
        fatal("FEMElementDescription:getEdgeBasisElement", "Saino will abort the simulation now...");
    }
    
    switch (element->Type.ElementCode/100) {
        case 4:
        case 7:
        case 8:
            n = element->Type.NumberOfNodes;
            u = cblas_ddot(n, basis, 1, element->Type.NodeU, 1);
            v = cblas_ddot(n, basis, 1, element->Type.NodeV, 1);
            w = cblas_ddot(n, basis, 1, element->Type.NodeW, 1);
            
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *dBasisdx, 3, element->Type.NodeU, 1, 0.0, dudx1, 1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *dBasisdx, 3, element->Type.NodeV, 1, 0.0, dudx2, 1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *dBasisdx, 3, element->Type.NodeW, 1, 0.0, dudx3, 1);
            
            triBase[0] = 1.0 - u - v;
            triBase[1] = u;
            triBase[2] = v;
            
            for (int i=0; i<3; i++) {
                dudx[0][i] = dudx1[i];
                dudx[1][i] = dudx2[i];
                dudx[2][i] = dudx3[i];
            }

            for (int i=0; i<3; i++) {
                dTriBase[0][i] = -dudx[0][i] - dudx[1][i];
                dTriBase[1][i] = dudx[0][i];
                dTriBase[2][i] = dudx[1][i];
            }
            break;
            
        case 6:
            n = element->Type.NumberOfNodes;
            u = cblas_ddot(n, basis, 1, element->Type.NodeU, 1);
            v = cblas_ddot(n, basis, 1, element->Type.NodeV, 1);
            w = cblas_ddot(n, basis, 1, element->Type.NodeW, 1);
            
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *dBasisdx, 3, element->Type.NodeU, 1, 0.0, g1, 1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *dBasisdx, 3, element->Type.NodeV, 1, 0.0, g2, 1);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *dBasisdx, 3, element->Type.NodeW, 1, 0.0, g3, 1);
            
            for (int i=0; i<3; i++) {
                g[0][i] = g1[i];
                g[1][i] = g2[i];
                g[2][i] = g3[i];
            }
            detG = g[0][0] * ( g[1][1]*g[2][2] - g[1][2]*g[2][1]) + g[0][1] * ( g[1][2]*g[2][0] - g[1][0]*g[2][2] ) +
                   g[0][2] * ( g[1][0]*g[2][1] - g[1][1]*g[2][0]);
            detF = 1.0 / detG;
            [self invertMatrix3x3:g inverted:f detG:detG];
            
            // The basis functions spanning the reference element space
            // and their curl with respect to the local coordinates
            edgeBasis[0][0] = (1.0 - v - w) / 4.0;
            edgeBasis[0][1] = 0.0;
            edgeBasis[0][2] = (u * (-1.0 + v + w)) / (4.0 * (-1.0 + w));
            curlBasis[0][0] = u / (4.0 * (-1.0 + w));
            curlBasis[0][1] = -(-2.0 + v + 2.0 * w) / (4.0 * (-1.0 + w));
            curlBasis[0][2] = 0.25;
            
            edgeBasis[1][0] = 0.0;
            edgeBasis[1][1] = (1.0 + u - w) /4.0;
            edgeBasis[1][2] = (v * (1.0 + u - w)) / (4.0 - 4.0 * w);
            curlBasis[1][0] = (2.0 + u - 2.0 * w) / (4.0 - 4.0 * w);
            curlBasis[1][1] = v / (4.0 * (-1.0 + w));
            curlBasis[1][2] = 0.25;
            
            edgeBasis[2][0] = (1.0 + v - w) / 4.0;
            edgeBasis[2][1] = 0.0;
            edgeBasis[2][2] = (u * (1.0 + v - w)) / (4.0 - 4.0 * w);
            curlBasis[2][0] = u / (4.0 - 4.0 * w);
            curlBasis[2][1] = (2.0 + v - 2.0 * w) / (4.0 * (-1.0 + w));
            curlBasis[2][2] = -0.25;
            
            edgeBasis[3][0] = 0.0;
            edgeBasis[3][1] = (1.0 - u - w) / 4.0;
            edgeBasis[3][2] = (v * (-1.0 + u + w)) / (4.0 * (-1.0 + w));
            curlBasis[3][0] = (-2.0 + u + 2.0 * w) / (4.0 * (-1.0 + w));
            curlBasis[3][1] = v / (4.0 - 4.0 * w);
            curlBasis[3][2] = -0.25;
            
            edgeBasis[4][0] = (w * (-1.0 + v + w)) / (4.0 * (-1.0 + w));
            edgeBasis[4][1] = (w * (-1.0 + u + w)) / (4.0 * (-1.0 + w));
            edgeBasis[4][2] = ( -((1.0 + v)*pow((-1.0 + w), 2.0)) + u * (v - pow((-1.0 + w), 2.0) - 2.0 * v * w) ) /
                              (4.0 * pow((-1.0 + w), 2.0));
            curlBasis[4][0] = -(-1.0 + u + w) / (2.0 * (-1.0 + w));
            curlBasis[4][1] = (-1.0 + v + w) / (2.0 * (-1.0 + w));
            curlBasis[4][2] = 0.0;
            
            edgeBasis[5][0] = -(w * (-1.0 + v + w)) / (4.0 * (-1.0 + w));
            edgeBasis[5][1] = (w * (-1.0 - u + w)) / (4.0 * (-1.0 + w));
            edgeBasis[5][2] = ( -((-1.0 + v)*pow((-1.0 + w), 2.0)) + u * (pow((-1.0 + w), 2.0) + v * (-1.0 + 2.0 * w)) ) /
                              (4.0 * pow((-1.0 + w), 2.0));
            curlBasis[5][0] = (1.0 + u - w) / (2.0 * (-1.0 + w));
            curlBasis[5][1] = -(-1.0 + v + w) / (2.0 * (-1.0 + w));
            curlBasis[5][2] = 0.0;
            
            edgeBasis[6][0] = ((1.0 + v - w) * w) / (4.0 * (-1.0 + w));
            edgeBasis[6][1] = ((1.0 + u - w) * w) / (4.0 * (-1.0 + w));
            edgeBasis[6][2] = ( (1.0 + v)*pow((-1.0 + w), 2.0) + u * (v + pow((-1.0 + w), 2.0) - 2.0 * v * w) ) /
                              (4.0 * pow((-1.0 + w), 2.0));
            curlBasis[6][0] = (1.0 + u - w) / (2.0 - 2.0 * w);
            curlBasis[6][1] = (1.0 + v - w) / (2.0 * (-1.0 + w));
            curlBasis[6][2] = 0.0;
            
            edgeBasis[7][0] = (w * (-1.0 - v + w)) / (4.0 * (-1.0 + w));
            edgeBasis[7][1] = -(w * (-1.0 + u + w)) / (4.0 * (-1.0 + w));
            edgeBasis[7][2] = ( (1.0 + v)*pow((-1.0 + w), 2.0) - u * (v + pow((-1.0 + w), 2.0) -2.0 * v * w) ) /
                              (4.0 * pow((-1.0 + w), 2.0));
            curlBasis[7][0] = (-1.0 + u + w) / (2.0 * (-1.0 + w));
            curlBasis[7][1] = (1.0 + v - w) / (2.0 - 2.0 * w);
            curlBasis[7][2] = 0.0;
            break;
    }
    
    edgeMap = [self getEdgeMap:element->Type.ElementCode/100];
    switch (element->Type.ElementCode/100) {
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
    for (int i=0; i<size; i++) {
        j = edgeMap[i][0]; k = edgeMap[i][1];
        
        nj = element->NodeIndexes[j];
        //TODO: add support for parallel run
        nk = element->NodeIndexes[k];
        //TODO: add support for parallel run
        
        switch (element->Type.ElementCode/100) {
            case 3:
            case 5:
                for (int l=0; l<3; l++) {
                    wBasis[i][l] = basis[j]*dBasisdx[k][l] - basis[k]*dBasisdx[j][l];
                }
                rotWBasis[i][0] = 2.0 * ( dBasisdx[j][1] * dBasisdx[k][2] - dBasisdx[j][2] * dBasisdx[k][1] );
                rotWBasis[i][1] = 2.0 * ( dBasisdx[j][2] * dBasisdx[k][0] - dBasisdx[j][0] * dBasisdx[k][2] );
                rotWBasis[i][2] = 2.0 * ( dBasisdx[j][0] * dBasisdx[k][1] - dBasisdx[j][1] * dBasisdx[k][0] );
                break;
                
            case 6:
                // Create the referential description of basis functions and their
                // spatial curl on the physical element via applying the Piola transform
                for (k=0; k<3; k++) {
                    sum = 0.0;
                    for (int l=0; l<3; l++) {
                        sum = sum + g[l][k] * edgeBasis[i][l];
                    }
                    wBasis[i][k] = sum;
                }
                for (k=0; k<3; k++) {
                    sum = 0.0;
                    for (int l=0; l<3; l++) {
                        sum = sum + f[k][l] * curlBasis[i][l];
                    }
                    rotWBasis[i][k] = 1.0 / detF * sum;
                }
                break;
                
            case 7:
                switch (i) {
                    case 0:
                        j = 0; k = 1;
                        base = (1.0 - w) / 2.0;
                        for (int l=0; l<3; l++) {
                            dBase[l] = -dudx[2][l] / 2.0;
                        }
                        break;
                    case 1:
                        j = 1; k = 2;
                        base = (1.0 - w) / 2.0;
                        for (int l=0; l<3; l++) {
                            dBase[l] = -dudx[2][l] / 2.0;
                        }
                        break;
                    case 2:
                        j = 2; k = 0;
                        base = (1.0 - w) / 2.0;
                        for (int l=0; l<3; l++) {
                            dBase[l] = -dudx[2][l] / 2.0;
                        }
                        break;
                    case 3:
                        j = 0; k = 1;
                        base = (1.0 + w) / 2.0;
                        for (int l=0; l<3; l++) {
                            dBase[l] = dudx[2][l] / 2.0;
                        }
                        break;
                    case 4:
                        j = 1; k = 2;
                        base = (1.0 + w) / 2.0;
                        for (int l=0; l<3; l++) {
                            dBase[l] = dudx[2][l] / 2.0;
                        }
                        break;
                    case 5:
                        j = 2; k = 0;
                        base = (1.0 + w) / 2.0;
                        for (int l=0; l<3; l++) {
                            dBase[l] = dudx[2][l] / 2.0;
                        }
                        break;
                    case 6:
                        base = triBase[0];
                        for (int l=0; l<3; l++) {
                            dBase[l] = dTriBase[0][l];
                            du[l] = dudx[2][l] / 2.0;
                        }
                        break;
                    case 7:
                        base = triBase[1];
                        for (int l=0; l<3; l++) {
                            dBase[l] = dTriBase[1][l];
                            du[l] = dudx[2][l] / 2.0;
                        }
                        break;
                    case 8:
                        base = triBase[2];
                        for (int l=0; l<3; l++) {
                            dBase[l] = dTriBase[2][l];
                            du[l] = dudx[2][l] / 2.0;
                        }
                        break;
                }
                
                if (i<=5) {
                    for (int l=0; l<3; l++) {
                        tbase[l] = triBase[j]*dTriBase[k][l] - triBase[k]*dTriBase[j][l];
                    }
                    rBase[0] = 2.0 * base * (dTriBase[j][1]*dTriBase[k][2] - dTriBase[k][1]*dTriBase[j][2]) +
                               dBase[1]*tbase[2] - dBase[2]*tbase[1];
                    rBase[1] = 2.0 * base * (dTriBase[j][2]*dTriBase[k][0] - dTriBase[k][2]*dTriBase[j][0]) +
                               dBase[2]*tbase[0] - dBase[0]*tbase[2];
                    rBase[2] = 2.0 * base * (dTriBase[j][0]*dTriBase[k][1] - dTriBase[k][0]*dTriBase[j][1]) +
                               dBase[0]*tbase[1] - dBase[1]*tbase[0];
                    
                    for (int l=0; l<3; l++) {
                        rotWBasis[i][l] = rBase[l];
                        wBasis[i][l] = tbase[l] * base;
                    }
                } else {
                    for (int l=0; l<3; l++) {
                        wBasis[i][l] = base * du[l];
                    }
                    rotWBasis[i][0] = (dBase[1]*du[2] - dBase[2]*du[1]);
                    rotWBasis[i][1] = (dBase[2]*du[0] - dBase[0]*du[2]);
                    rotWBasis[i][2] = (dBase[0]*du[1] - dBase[1]*du[0]);
                }
                break;
                
            case 4:
                switch (i) {
                    case 0:
                         for (int l=0; l<3; l++) {
                             du[l] = dudx[0][l];
                             dBase[l] = -dudx[1][l]*(1.0 - w) - (1.0 - v)*dudx[2][l];
                         }
                        base = (1.0 - v) * (1.0 - w);
                        break;
                    case 1:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[1][l];
                            dBase[l] = dudx[0][l]*(1.0 - w) - (1.0 + u)*dudx[2][l];
                        }
                        base = (1.0 + u) * (1.0 - w);
                        break;
                    case 2:
                        for (int l=0; l<3; l++) {
                            du[l] = -dudx[0][l];
                            dBase[l] = dudx[1][l]*(1.0 - w) - (1.0 + v)*dudx[2][l];
                        }
                        base = (1.0 + v) * (1.0 - w);
                        break;
                    case 3:
                        for (int l=0; l<3; l++) {
                            du[l] = -dudx[1][l];
                            dBase[l] = -dudx[0][l]*(1.0 - w) - (1.0 - u)*dudx[2][l];
                        }
                        base = (1.0 - u) * (1.0 - w);
                        break;
                }
                for (int l=0; l<3; l++) {
                    wBasis[i][l] = base * du[l] / n;
                }
                rotWBasis[i][0] = (dBase[1]*du[2] - dBase[2]*du[1]) / n;
                rotWBasis[i][1] = (dBase[2]*du[0] - dBase[0]*du[2]) / n;
                rotWBasis[i][2] = (dBase[0]*du[1] - dBase[1]*du[0]) / n;
                break;
                
            case 8:
                switch (i) {
                    case 0:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[0][l];
                            dBase[l] = -dudx[1][l]*(1.0 - w) - (1.0 - v)*dudx[2][l];
                        }
                        base = (1.0 - v) * (1.0 - w);
                        break;
                    case 1:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[1][l];
                            dBase[l] = dudx[0][l]*(1.0 - w) - (1.0 + u)*dudx[2][l];
                        }
                        base = (1.0 + u) * (1.0 - w);
                        break;
                    case 2:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[0][l];
                            dBase[l] = dudx[1][l]*(1.0 - w) - (1.0 + v)*dudx[2][l];
                        }
                        base = (1.0 + v) * (1.0 - w);
                        break;
                    case 3:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[1][l];
                            dBase[l] = -dudx[0][l]*(1.0 - w) - (1.0 - u)*dudx[2][l];
                        }
                        base = (1.0 - u) * (1.0 - w);
                        break;
                    case 4:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[0][l];
                            dBase[l] = -dudx[1][l]*(1.0 + w) + (1.0 - v)*dudx[2][l];
                        }
                        base = (1.0 - v) * (1.0 + w);
                        break;
                    case 5:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[1][l];
                            dBase[l] = dudx[0][l]*(1.0 + w) + (1.0 + u)*dudx[2][l];
                        }
                        base = (1.0 + u) * (1.0 + w);
                        break;
                    case 6:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[0][l];
                            dBase[l] = dudx[1][l]*(1.0 + w) + (1.0 + v)*dudx[2][l];
                        }
                        base = (1.0 + v) * (1.0 + w);
                        break;
                    case 7:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[1][l];
                            dBase[l] = -dudx[0][l]*(1.0 + w) + (1.0 - u)*dudx[2][l];
                        }
                        base = (1.0 - u) * (1.0 + w);
                        break;
                    case 8:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[2][l];
                            dBase[l] = -dudx[0][l]*(1.0 - v) - (1.0 - u)*dudx[1][l];
                        }
                        base = (1.0 - u) * (1.0 - v);
                        break;
                    case 9:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[2][l];
                            dBase[l] = dudx[0][l]*(1.0 - v) - (1.0 + u)*dudx[1][l];
                        }
                        base = (1.0 + u) * (1.0 - v);
                        break;
                    case 10:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[2][l];
                            dBase[l] = dudx[0][l]*(1.0 + v) + (1.0 + u)*dudx[1][l];
                        }
                        base = (1.0 + u) * (1.0 + v);
                        break;
                    case 11:
                        for (int l=0; l<3; l++) {
                            du[l] = dudx[2][l];
                            dBase[l] = -dudx[0][l]*(1.0 + v) + (1.0 - u)*dudx[1][l];
                        }
                        base = (1.0 - u) * (1.0 + v);
                        break;
                }
                
                for (int l=0; l<3; l++) {
                    wBasis[i][l] = base * du[l] / n;
                }
                rotWBasis[i][0] = (dBase[1]*du[2] - dBase[2]*du[1]) / n;
                rotWBasis[i][1] = (dBase[2]*du[0] - dBase[0]*du[2]) / n;
                rotWBasis[i][2] = (dBase[0]*du[1] - dBase[1]*du[0]) / n;
                break;
                
            default:
                fprintf(stderr, "FEMElementDescription:getEdgeBasisElement: not implemented for this element type: %d.\n", element->Type.ElementCode/100);
                fatal("FEMElementDescription:getEdgeBasisElement", "Saino will abort the simulation now...");
                break;
        }
        if (nk < nj) {
            for (int l=0; l<3; l++) {
                wBasis[i][l] = -wBasis[i][l];
                rotWBasis[i][l] = -rotWBasis[i][l];
            }
        }
    }
}

@end
