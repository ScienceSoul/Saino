//
//  ElementsDefinition.h
//  Saino
//
//  Created by Hakime Seddik on 29/05/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#ifndef Saino_ElementsDefinition_h
#define Saino_ElementsDefinition_h

typedef struct ElementDef {
    
    int dimension;
    char *topology;
    int code;
    int nodes;
    double nodeU[27];
    double nodeV[27];
    double nodeW[27];
    int basis[27];
    int gaussPoints[2];
    double stabilization;
        
} ElementDef;

//------------------------------
// Number of elements definition
//------------------------------
#define NUMBER_ELEMENT_DEFS 21

//---------------------------------------
// List of supported elements definition
//---------------------------------------
ElementDef _1nodePoint;
ElementDef _2nodePeriodic;
ElementDef _2nodeLine;
ElementDef _3nodeLine;
ElementDef _4nodeLine;
ElementDef _3nodeTriangle;
ElementDef _6nodeTriangle;
ElementDef _10nodeTriangle;
ElementDef _4nodeQuadrilateral;
ElementDef _8nodeQuadrilateral;
ElementDef _9nodeQuadrilateral;
ElementDef _12nodeQuadrilateral;
ElementDef _16nodeQuadrilateral;
ElementDef _4nodeTetrahedron;
ElementDef _10nodeTetrahedron;
ElementDef _5nodePyramid;
ElementDef _13nodePyramid;
ElementDef _6nodeWedge;
ElementDef _8nodeOctahedron;
ElementDef _20nodeOctahedron;
ElementDef _27nodeOctahedron;

//---------------------
// 1 node point element
//---------------------
_1nodePoint.dimension = 1;
_1nodePoint.topology = "point";
_1nodePoint.code = 101;
_1nodePoint.nodes = 1;
_1nodePoint.nodeU[0] = 0.0;
_1nodePoint.basis[0] = 1;

_1nodePoint.gaussPoints[0] = 1;
_1nodePoint.stabilization = 0.0;

//-------------------------
// 2 nodes periodic element
//-------------------------
_2nodePeriodic.dimension = 1;
_2nodePeriodic.topology = "point";
_2nodePeriodic.code = 102;
_2nodePeriodic.nodes = 2;
_2nodePeriodic.nodeU[0] = 0.0;
_2nodePeriodic.nodeU[1] = 1.0;
_2nodePeriodic.basis[0] = 1;
_2nodePeriodic.basis[1] = 2;

_2nodePeriodic.gaussPoints[0] = 0;
_2nodePeriodic.stabilization = 0.0;

//---------------------
// 2 nodes line element
//---------------------
_2nodeLine.dimension = 1;
_2nodeLine.topology = "line";
_2nodeLine.code = 202;
_2nodeLine.nodes = 2;
//          1    2 
_2nodeLine.nodeU[0] = -1.0;
_2nodeLine.nodeU[1] = 1.0;
//          1  2   3   4
//          1  u  u^2 u^3
_2nodeLine.basis[0] = 1;
_2nodeLine.basis[1] = 2;

_2nodeLine.gaussPoints[0] = 2;
_2nodeLine.gaussPoints[1] = 4;
_2nodeLine.stabilization = 0.33333333333333333333333;

//---------------------
// 3 nodes line element
//---------------------
_3nodeLine.dimension = 1;
_3nodeLine.topology = "line";
_3nodeLine.code = 203;
_3nodeLine.nodes = 3;
//          1     2     3
_3nodeLine.nodeU[0] = -1.0;
_3nodeLine.nodeU[1] = 1.0;
_3nodeLine.nodeU[2] = 0.0;
//          1  2   3   4
//          1  u  u^2 u^3
_3nodeLine.basis[0] = 1;
_3nodeLine.basis[1] = 2;
_3nodeLine.basis[2] = 3;

_3nodeLine.gaussPoints[0] = 3;
_3nodeLine.gaussPoints[1] = 5;
_3nodeLine.stabilization = 0.1666666666666666666666

//---------------------
// 4 nodes line element
//---------------------
_4nodeLine.dimension = 1;
_4nodeLine.topology = "line";
_4nodeLine.code = 204;
_4nodeLine.nodes = 4;
//          1     2     3     4
_4nodeLine.nodeU[0] = -1.0;
_4nodeLine.nodeU[1] = 1.0;
_4nodeLine.nodeU[2] = -0.333333333333333333;
_4nodeLine.nodeU[3] = 0.3333333333333333333;
//          1  2   3   4
//          1  u  u^2 u^3
_4nodeLine.basis[0] = 1;
_4nodeLine.basis[1] = 2;
_4nodeLine.basis[2] = 3;
_4nodeLine.basis[3] = 4;

_4nodeLine.gaussPoints[0] = 4;
_4nodeLine.gaussPoints[1] = 6;
_4nodeLine.stabilization = 0.03333333333333333333333;

//-------------------------
// 3 nodes triangle element
//-------------------------
_3nodeTriangle.dimension = 2;
_3nodeTriangle.topology = "triangle";
_3nodeTriangle.code = 303;
_3nodeTriangle.nodes = 3;
//          1     2     3
_3nodeTriangle.nodeU[0] = 0.0;
_3nodeTriangle.nodeU[1] = 1.0;
_3nodeTriangle.nodeU[2] = 0.0;
_3nodeTriangle.nodeV[0] = 0.0;
_3nodeTriangle.nodeV[1] = 0.0;
_3nodeTriangle.nodeV[2] = 1.0;
//       1  2   3   4   5   6   7     8    9     10 
//       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
//
//      11      12     13 14    15     16
//      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
_3nodeTriangle.basis[0] = 1;
_3nodeTriangle.basis[1] = 2;
_3nodeTriangle.basis[2] = 5;

_3nodeTriangle.gaussPoints[0] = 3;
_3nodeTriangle.gaussPoints[1] = 20;
_3nodeTriangle.stabilization = 0.333333333333333333333;

//-------------------------
// 6 nodes triangle element
//-------------------------
_6nodeTriangle.dimension = 2;
_6nodeTriangle.topology = "triangle";
_6nodeTriangle.code = 306;
_6nodeTriangle.nodes = 6;
//          1     2     3     4     5     6
_6nodeTriangle.nodeU[0] = 0.0;
_6nodeTriangle.nodeU[1] = 1.0;
_6nodeTriangle.nodeU[2] = 0.0;
_6nodeTriangle.nodeU[3] = 0.5;
_6nodeTriangle.nodeU[4] = 0.5;
_6nodeTriangle.nodeU[5] = 0.0;
_6nodeTriangle.nodeV[0] = 0.0;
_6nodeTriangle.nodeV[1] = 0.0;
_6nodeTriangle.nodeV[2] = 1.0;
_6nodeTriangle.nodeV[3] = 0.0;
_6nodeTriangle.nodeV[4] = 0.5;
_6nodeTriangle.nodeV[5] = 0.5;
//       1  2   3   4   5   6   7     8    9     10 
//       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
//
//      11      12     13 14    15     16
//      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
_6nodeTriangle.basis[0] = 1;
_6nodeTriangle.basis[1] = 2;
_6nodeTriangle.basis[2] = 3;
_6nodeTriangle.basis[3] = 5;
_6nodeTriangle.basis[4] = 6;
_6nodeTriangle.basis[5] = 9;

_6nodeTriangle.gaussPoints[0] = 6;
_6nodeTriangle.gaussPoints[1] = 17;
_6nodeTriangle.stabilization = 0.041666666666666666;

//--------------------------
// 10 nodes triangle element
//--------------------------
_10nodeTriangle.dimension = 2;
_10nodeTriangle.topology = "triangle";
_10nodeTriangle.code = 310;
_10nodeTriangle,nodes = 10;
//          1     2     3     ......
_10nodeTriangle.nodeU[0] = 0.0;
_10nodeTriangle.nodeU[1] = 1.0;
_10nodeTriangle.nodeU[2] = 0.0;
_10nodeTriangle.nodeU[3] = 0.333333333333333333;
_10nodeTriangle.nodeU[4] = 0.666666666666666667;
_10nodeTriangle.nodeU[5] = 0.666666666666666667;
_10nodeTriangle.nodeU[6] = 0.333333333333333333;
_10nodeTriangle.nodeU[7] = 0.000000000000000000;
_10nodeTriangle.nodeU[8] = 0.000000000000000000;
_10nodeTriangle.nodeU[9] = 0.333333333333333333;
_10nodeTriangle.nodeV[0] = 0.0;
_10nodeTriangle.nodeV[1] = 0.0;
_10nodeTriangle.nodeV[2] = 1.0;
_10nodeTriangle.nodeV[3] = 0.000000000000000000;
_10nodeTriangle.nodeV[4] = 0.000000000000000000;
_10nodeTriangle.nodeV[5] = 0.333333333333333333;
_10nodeTriangle.nodeV[6] = 0.666666666666666667;
_10nodeTriangle.nodeV[7] = 0.666666666666666667;
_10nodeTriangle.nodeV[8] = 0.333333333333333333;
_10nodeTriangle.nodeV[9] = 0.333333333333333333;
//       1  2   3   4   5   6   7     8    9     10 
//       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
//
//      11      12     13 14    15     16
//      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
_10nodeTriangle.basis[0] = 1;
_10nodeTriangle.basis[1] = 2;
_10nodeTriangle.basis[2] = 3;
_10nodeTriangle.basis[3] = 4;
_10nodeTriangle.basis[4] = 5;
_10nodeTriangle.basis[5] = 6;
_10nodeTriangle.basis[6] = 7;
_10nodeTriangle.basis[7] = 9;
_10nodeTriangle.basis[8] = 10;
_10nodeTriangle.basis[9] = 13;

_10nodeTriangle.gaussPoints[0] = 25;
_10nodeTriangle.gaussPoints[1] = 64;
_10nodeTriangle.stabilization = 0.01341555597798937329;

//------------------------------
// 4 nodes quadrilateral element
//------------------------------
_4nodeQuadrilateral.dimension = 2;
_4nodeQuadrilateral.topology = "quadrilateral";
_4nodeQuadrilateral.code = 404;
_4nodeQuadrilateral.nodes = 4;
//          1     2     3     4
_4nodeQuadrilateral.nodeU[0] = -1.0;
_4nodeQuadrilateral.nodeU[1] = 1.0;
_4nodeQuadrilateral.nodeU[2] = 1.0;
_4nodeQuadrilateral.nodeU[3] = -1.0;
_4nodeQuadrilateral.nodeV[0] = -1.0;
_4nodeQuadrilateral.nodeV[1] = -1.0;
_4nodeQuadrilateral.nodeV[2] = 1.0;
_4nodeQuadrilateral.nodeV[3] = 1.0;
//       1  2   3   4   5   6   7     8    9     10 
//       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
//
//      11      12     13 14    15     16
//      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
_4nodeQuadrilateral.basis[0] = 1;
_4nodeQuadrilateral.basis[1] = 2;
_4nodeQuadrilateral.basis[2] = 5;
_4nodeQuadrilateral.basis[3] = 6;

_4nodeQuadrilateral.gaussPoints[0] = 4;
_4nodeQuadrilateral.gaussPoints[1] = 16;
_4nodeQuadrilateral.stabilization = 0.333333333333333333333;

//------------------------------
// 8 nodes quadrilateral element
//------------------------------
_8nodeQuadrilateral.dimension = 2;
_8nodeQuadrilateral.topology = "quadrilateral";
_8nodeQuadrilateral.code = 408;
_8nodeQuadrilateral.nodes = 8;
//          1     2     3     4     5     6      7      8
_8nodeQuadrilateral.nodeU[0] = -1.0;
_8nodeQuadrilateral.nodeU[1] = 1.0;
_8nodeQuadrilateral.nodeU[2] = 1.0;
_8nodeQuadrilateral.nodeU[3] = -1.0;
_8nodeQuadrilateral.nodeU[4] = 0.0;
_8nodeQuadrilateral.nodeU[5] = 1.0;
_8nodeQuadrilateral.nodeU[6] = 0.0;
_8nodeQuadrilateral.nodeU[7] = -1.0;
_8nodeQuadrilateral.nodeV[0] = -1.0;
_8nodeQuadrilateral.nodeV[1] = -1.0;
_8nodeQuadrilateral.nodeV[2] = 1.0;
_8nodeQuadrilateral.nodeV[3] = 1.0;
_8nodeQuadrilateral.nodeV[4] = -1.0;
_8nodeQuadrilateral.nodeV[5] = 0.0;
_8nodeQuadrilateral.nodeV[6] = 1.0;
_8nodeQuadrilateral.nodeV[7] = 0.0;
//       1  2   3   4   5   6   7     8    9     10 
//       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
//
//      11      12     13 14    15     16
//      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
_8nodeQuadrilateral.basis[0] = 1;
_8nodeQuadrilateral.basis[1] = 2;
_8nodeQuadrilateral.basis[2] = 3;
_8nodeQuadrilateral.basis[3] = 5;
_8nodeQuadrilateral.basis[4] = 6;
_8nodeQuadrilateral.basis[5] = 7;
_8nodeQuadrilateral.basis[6] = 9;
_8nodeQuadrilateral.basis[7] = 10;

_8nodeQuadrilateral.gaussPoints[0] = 9;
_8nodeQuadrilateral.gaussPoints[1] = 25;
_8nodeQuadrilateral.stabilization = 0.0833333333333333333333;

//------------------------------
// 9 nodes quadrilateral element
//------------------------------
_9nodeQuadrilateral.dimension = 2;
_9nodeQuadrilateral.topology = "quadrilateral";
_9nodeQuadrilateral.code = 409;
_9nodeQuadrilateral.nodes = 9;
//          1     2     3     4     5     6      7      8      9
_9nodeQuadrilateral.nodeU[0] = -1.0;
_9nodeQuadrilateral.nodeU[1] = 1.0;
_9nodeQuadrilateral.nodeU[2] = 1.0;
_9nodeQuadrilateral.nodeU[3] = -1.0;
_9nodeQuadrilateral.nodeU[4] = 0.0;
_9nodeQuadrilateral.nodeU[5] = 1.0;
_9nodeQuadrilateral.nodeU[6] = 0.0;
_9nodeQuadrilateral.nodeU[7] = -1.0;
_9nodeQuadrilateral.nodeU[8] = 0.0;
_9nodeQuadrilateral.nodeV[0] = -1.0;
_9nodeQuadrilateral.nodeV[1] = -1.0;
_9nodeQuadrilateral.nodeV[2] = 1.0;
_9nodeQuadrilateral.nodeV[3] = 1.0;
_9nodeQuadrilateral.nodeV[4] = -1.0;
_9nodeQuadrilateral.nodeV[5] = 0.0;
_9nodeQuadrilateral.nodeV[6] = 1.0;
_9nodeQuadrilateral.nodeV[7] = 0.0;
_9nodeQuadrilateral.nodeV[8] = 0.0;
//       1  2   3   4   5   6   7     8    9     10 
//       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
//
//      11      12     13 14    15     16
//      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
_9nodeQuadrilateral.basis[0] = 1;
_9nodeQuadrilateral.basis[1] = 2;
_9nodeQuadrilateral.basis[2] = 3;
_9nodeQuadrilateral.basis[3] = 5;
_9nodeQuadrilateral.basis[4] = 6;
_9nodeQuadrilateral.basis[5] = 7;
_9nodeQuadrilateral.basis[6] = 9;
_9nodeQuadrilateral.basis[7] = 10;
_9nodeQuadrilateral.basis[8] = 11;

_9nodeQuadrilateral.gaussPoints[0] = 9;
_9nodeQuadrilateral.gaussPoints[1] = 25;
_9nodeQuadrilateral.stabilization = 0.0833333333333333333333;

//-------------------------------
// 12 nodes quadrilateral element
//-------------------------------
_12nodeQuadrilateral.dimension = 2;
_12nodeQuadrilateral.topology = "quadrilateral";
_12nodeQuadrilateral.code = 412;
_12nodeQuadrilateral.nodes = 12;
//          1     2     3     ......
_12nodeQuadrilateral.nodeU[0] = -1.0;
_12nodeQuadrilateral.nodeU[1] = 1.0;
_12nodeQuadrilateral.nodeU[2] = 1.0;
_12nodeQuadrilateral.nodeU[3] = -1.0
_12nodeQuadrilateral.nodeU[4] = -0.3333333333333333333;
_12nodeQuadrilateral.nodeU[5] = 0.3333333333333333333;
_12nodeQuadrilateral.nodeU[6] = 1.0000000000000000000;
_12nodeQuadrilateral.nodeU[7] = 1.0000000000000000000;
_12nodeQuadrilateral.nodeU[8] = 0.3333333333333333333;
_12nodeQuadrilateral.nodeU[9] = -0.3333333333333333333;
_12nodeQuadrilateral.nodeU[10] = -1.0000000000000000000;
_12nodeQuadrilateral.nodeU[11] = -1.0000000000000000000;
_12nodeQuadrilateral.nodeV[0] = -1.0;
_12nodeQuadrilateral.nodeV[1] = -1.0;
_12nodeQuadrilateral.nodeV[2] = 1.0;
_12nodeQuadrilateral.nodeV[3] = 1.0;
_12nodeQuadrilateral.nodeV[4] = -1.0000000000000000000;
_12nodeQuadrilateral.nodeV[5] = -1.0000000000000000000;
_12nodeQuadrilateral.nodeV[6] = -0.3333333333333333333;
_12nodeQuadrilateral.nodeV[7] = 0.3333333333333333333;
_12nodeQuadrilateral.nodeV[8] = 1.0000000000000000000;
_12nodeQuadrilateral.nodeV[9] = 1.0000000000000000000;
_12nodeQuadrilateral.nodeV[10] = 0.3333333333333333333;
_12nodeQuadrilateral.nodeV[11] = -0.3333333333333333333;
//       1  2   3   4   5   6   7     8    9     10 
//       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
//
//      11      12     13 14    15     16
//      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
_12nodeQuadrilateral.basis[0] = 1;
_12nodeQuadrilateral.basis[1] = 2;
_12nodeQuadrilateral.basis[2] = 3;
_12nodeQuadrilateral.basis[3] = 4;
_12nodeQuadrilateral.basis[4] = 5;
_12nodeQuadrilateral.basis[5] = 6;
_12nodeQuadrilateral.basis[6] = 7;
_12nodeQuadrilateral.basis[7] = 8;
_12nodeQuadrilateral.basis[8] = 9;
_12nodeQuadrilateral.basis[9] = 10;
_12nodeQuadrilateral.basis[10] = 13;
_12nodeQuadrilateral.basis[11] = 14;

_12nodeQuadrilateral.gaussPoints[0] = 16;
_12nodeQuadrilateral.gaussPoints[1] = 36;
_12nodeQuadrilateral.stabilization = 0.0;

//-------------------------------
// 16 nodes quadrilateral element
//-------------------------------
_16nodeQuadrilateral.dimension = 2;
_16nodeQuadrilateral.topology = "quadrilateral";
_16nodeQuadrilateral.code = 416;
_16nodeQuadrilateral.nodes = 16;
//          1     2     3     ......
_16nodeQuadrilateral.nodeU[0] = -1.0;
_16nodeQuadrilateral.nodeU[1] = 1.0;
_16nodeQuadrilateral.nodeU[2] = 1.0;
_16nodeQuadrilateral.nodeU[3] = -1.0;
_16nodeQuadrilateral.nodeU[4] = -0.3333333333333333333;
_16nodeQuadrilateral.nodeU[5] = 0.3333333333333333333;
_16nodeQuadrilateral.nodeU[6] = 1.0000000000000000000;
_16nodeQuadrilateral.nodeU[7] = 1.0000000000000000000;
_16nodeQuadrilateral.nodeU[8] = 0.3333333333333333333;
_16nodeQuadrilateral.nodeU[9] = -0.3333333333333333333;
_16nodeQuadrilateral.nodeU[10] = -1.0000000000000000000;
_16nodeQuadrilateral.nodeU[11] = -1.0000000000000000000;
_16nodeQuadrilateral.nodeU[12] = -0.3333333333333333333;
_16nodeQuadrilateral.nodeU[13] = 0.3333333333333333333;
_16nodeQuadrilateral.nodeU[14] = 0.3333333333333333333;
_16nodeQuadrilateral.nodeU[15] = -0.3333333333333333333;
_16nodeQuadrilateral.nodeV[0] = -1.0;
_16nodeQuadrilateral.nodeV[1] = -1.0;
_16nodeQuadrilateral.nodeV[2] = 1.0;
_16nodeQuadrilateral.nodeV[3] = 1.0;
_16nodeQuadrilateral.nodeV[4] = -1.0000000000000000000;
_16nodeQuadrilateral.nodeV[5] = -1.0000000000000000000;
_16nodeQuadrilateral.nodeV[6] = -0.3333333333333333333;
_16nodeQuadrilateral.nodeV[7] = 0.3333333333333333333;
_16nodeQuadrilateral.nodeV[8] = 1.0000000000000000000;
_16nodeQuadrilateral.nodeV[9] = 1.0000000000000000000; 
_16nodeQuadrilateral.nodeV[10] = 0.3333333333333333333;
_16nodeQuadrilateral.nodeV[11] = -0.3333333333333333333;
_16nodeQuadrilateral.nodeV[12] = -0.3333333333333333333;
_16nodeQuadrilateral.nodeV[13] = -0.3333333333333333333;
_16nodeQuadrilateral.nodeV[14] = 0.3333333333333333333;
_16nodeQuadrilateral.nodeV[15] = 0.3333333333333333333;
//       1  2   3   4   5   6   7     8    9     10 
//       1  u  u^2  u^3 v  uv  u^2v  u^3v v^2   uv^2 
//
//      11      12     13 14    15     16
//      u^2v^2 u^3v^2 v^3 uv^3 u^2v^3 u^3v^3
_16nodeQuadrilateral.basis[0] = 1;
_16nodeQuadrilateral.basis[1] = 2;
_16nodeQuadrilateral.basis[2] = 3;
_16nodeQuadrilateral.basis[3] = 4;
_16nodeQuadrilateral.basis[4] = 5;
_16nodeQuadrilateral.basis[5] = 6;
_16nodeQuadrilateral.basis[6] = 7;
_16nodeQuadrilateral.basis[7] = 8;
_16nodeQuadrilateral.basis[8] = 9;
_16nodeQuadrilateral.basis[9] = 10;
_16nodeQuadrilateral.basis[10] = 11;
_16nodeQuadrilateral.basis[11] = 12;
_16nodeQuadrilateral.basis[12] = 13;
_16nodeQuadrilateral.basis[13] = 14;
_16nodeQuadrilateral.basis[14] = 15;
_16nodeQuadrilateral.basis[15] = 16;

_16nodeQuadrilateral.gaussPoints[0] = 16;
_16nodeQuadrilateral.gaussPoints[1] = 36;
_16nodeQuadrilateral.stabilization = 0.01766875890919188522;

//----------------------------
// 4 nodes tetrahedron element
//----------------------------
_4nodeTetrahedron.dimension = 3;
_4nodeTetrahedron.topology = "tetrahedron";
_4nodeTetrahedron.code = 504;
_4nodeTetrahedron.nodes = 4;
//          1     2     3     4
_4nodeTetrahedron.nodeU[0] = 0.0;
_4nodeTetrahedron.nodeU[1] = 1.0;
_4nodeTetrahedron.nodeU[2] = 0.0;
_4nodeTetrahedron.nodeU[3] = 0.0;
_4nodeTetrahedron.nodeV[0] = 0.0;
_4nodeTetrahedron.nodeV[1] = 0.0;
_4nodeTetrahedron.nodeV[2] = 1.0;
_4nodeTetrahedron.nodeV[3] = 0.0;
_4nodeTetrahedron.nodeW[0] = 0.0;
_4nodeTetrahedron.nodeW[1] = 0.0;
_4nodeTetrahedron.nodeW[2] = 0.0;
_4nodeTetrahedron.nodeW[3] = 1.0;
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
_4nodeTetrahedron.basis[0] = 1;
_4nodeTetrahedron.basis[1] = 2;
_4nodeTetrahedron.basis[2] = 5;
_4nodeTetrahedron.basis[3] = 17;

_4nodeTetrahedron.gaussPoints[0] = 4;
_4nodeTetrahedron.gaussPoints[1] = 11;
_4nodeTetrahedron.stabilization = 0.333333333333333333;

//-----------------------------
// 10 nodes tetrahedron element
//-----------------------------
_10nodeTetrahedron.dimension = 3;
_10nodeTetrahedron.topology = "tetrahedron";
_10nodeTetrahedron.code = 510;
_10nodeTetrahedron.nodes = 10;
//          1     2     3     4     5     6      7      8      9     10
_10nodeTetrahedron.nodeU[0] = 0.0;
_10nodeTetrahedron.nodeU[1] = 1.0;
_10nodeTetrahedron.nodeU[2] = 0.0;
_10nodeTetrahedron.nodeU[3] = 0.0;
_10nodeTetrahedron.nodeU[4] = 0.5;
_10nodeTetrahedron.nodeU[5] = 0.5;
_10nodeTetrahedron.nodeU[6] = 0.0;
_10nodeTetrahedron.nodeU[7] = 0.0;
_10nodeTetrahedron.nodeU[8] = 0.5;
_10nodeTetrahedron.nodeU[9] = 0.0;
_10nodeTetrahedron.nodeV[0] = 0.0;
_10nodeTetrahedron.nodeV[1] = 0.0;
_10nodeTetrahedron.nodeV[2] = 1.0;
_10nodeTetrahedron.nodeV[3] = 0.0;
_10nodeTetrahedron.nodeV[4] = 0.0;
_10nodeTetrahedron.nodeV[5] = 0.5;
_10nodeTetrahedron.nodeV[6] = 0.5;
_10nodeTetrahedron.nodeV[7] = 0.0;
_10nodeTetrahedron.nodeV[8] = 0.0;
_10nodeTetrahedron.nodeV[9] = 0.5;
_10nodeTetrahedron.nodeW[0] = 0.0;
_10nodeTetrahedron.nodeW[1] = 0.0;
_10nodeTetrahedron.nodeW[2] = 0.0;
_10nodeTetrahedron.nodeW[3] = 1.0;
_10nodeTetrahedron.nodeW[4] = 0.0;
_10nodeTetrahedron.nodeW[5] = 0.0;
_10nodeTetrahedron.nodeW[6] = 0.0;
_10nodeTetrahedron.nodeW[7] = 0.5;
_10nodeTetrahedron.nodeW[8] = 0.5;
_10nodeTetrahedron.nodeW[9] = 0.5;
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
_10nodeTetrahedron.basis[0] = 1;
_10nodeTetrahedron.basis[1] = 2;
_10nodeTetrahedron.basis[2] = 3;
_10nodeTetrahedron.basis[3] = 5;
_10nodeTetrahedron.basis[4] = 6;
_10nodeTetrahedron.basis[5] = 9;
_10nodeTetrahedron.basis[6] = 17;
_10nodeTetrahedron.basis[7] = 18;
_10nodeTetrahedron.basis[8] = 21;
_10nodeTetrahedron.basis[9] = 33;

_10nodeTetrahedron.gaussPoints[0] = 27;
_10nodeTetrahedron.gaussPoints[1] = 64;
_10nodeTetrahedron.stabilization = 0.0416666666666666666;

//------------------------
// 5 nodes pyramid element
//------------------------
_5nodePyramid.dimension = 3;
_5nodePyramid.topology = "pyramid";
_5nodePyramid.code = 605;
_5nodePyramid.nodes = 5;
//          1     2     3     4     5
_5nodePyramid.nodeU[0] = -1.0;
_5nodePyramid.nodeU[1] = 1.0;
_5nodePyramid.nodeU[2] = 1.0;
_5nodePyramid.nodeU[3] = -1.0
_5nodePyramid.nodeU[4] = 0.0;
_5nodePyramid.nodeV[0] = -1.0;
_5nodePyramid.nodeV[1] = -1.0;
_5nodePyramid.nodeV[2] = 1.0;
_5nodePyramid.nodeV[3] = 1.0;
_5nodePyramid.nodeV[4] = 0.0;
_5nodePyramid.nodeW[0] = 0.0;
_5nodePyramid.nodeW[1] = 0.0;
_5nodePyramid.nodeW[2] = 0.0;
_5nodePyramid.nodeW[3] = 0.0;
_5nodePyramid.nodeW[4] = 1.0;
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
_5nodePyramid.basis[0] = 1; 
_5nodePyramid.basis[1] = 2;
_5nodePyramid.basis[2] = 5;
_5nodePyramid.basis[3] = 6;
_5nodePyramid.basis[4] = 17;

_5nodePyramid.gaussPoints[0] = 8;
_5nodePyramid.gaussPoints[1] = 64;
_5nodePyramid.stabilization = 0.333333333333333333;

//-------------------------
// 13 nodes pyramid element
//-------------------------
_13nodePyramid.dimension = 3;
_13nodePyramid.topology = "pyramid";
_13nodePyramid.code = 613;
_13nodePyramid.nodes = 13;
//          1     2     3     4     5     6      7      8      9     10    11    12    13
_13nodePyramid.nodeU[0] = -1.0;
_13nodePyramid.nodeU[1] = 1.0;
_13nodePyramid.nodeU[2] = 1.0;
_13nodePyramid.nodeU[3] = -1.0;
_13nodePyramid.nodeU[4] = 0.0;
_13nodePyramid.nodeU[5] = 0.0;
_13nodePyramid.nodeU[6] = 1.0
_13nodePyramid.nodeU[7] = 0.0;
_13nodePyramid.nodeU[8] = -1.0;
_13nodePyramid.nodeU[9] = -0.5;
_13nodePyramid.nodeU[10] = 0.5;
_13nodePyramid.nodeU[11] = 0.5;
_13nodePyramid.nodeU[12] = -0.5;
_13nodePyramid.nodeV[0] = -1.0;
_13nodePyramid.nodeV[1] = -1.0;
_13nodePyramid.nodeV[2] = 1.0;
_13nodePyramid.nodeV[3] = 1.0;
_13nodePyramid.nodeV[4] = 0.0;
_13nodePyramid.nodeV[5] = -1.0;
_13nodePyramid.nodeV[6] = 0.0;
_13nodePyramid.nodeV[7] = 1.0;
_13nodePyramid.nodeV[8] = 0.0;
_13nodePyramid.nodeV[9] = -0.5:
_13nodePyramid.nodeV[10] = -0.5;
_13nodePyramid.nodeV[11] = 0.5;
_13nodePyramid.nodeV[12] = 0.5;
_13nodePyramid.nodeW[0] = 0.0;
_13nodePyramid.nodeW[1] = 0.0;
_13nodePyramid.nodeW[2] = 0.0;
_13nodePyramid.nodeW[3] = 0.0;
_13nodePyramid.nodeW[4] = 1.0;
_13nodePyramid.nodeW[5] = 0.0;
_13nodePyramid.nodeW[6] = 0.0
_13nodePyramid.nodeW[7] = 0.0
_13nodePyramid.nodeW[8] = 0.0;
_13nodePyramid.nodeW[9] = 0.5;
_13nodePyramid.nodeW[10] = 0.5;
_13nodePyramid.nodeW[11] = 0.5;
_13nodePyramid.nodeW[12] = 0.5;
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
_13nodePyramid.basis[0] = 1;
_13nodePyramid.basis[1] = 2;
_13nodePyramid.basis[2] = 3;
_13nodePyramid.basis[3] = 5;
_13nodePyramid.basis[4] = 6;
_13nodePyramid.basis[5] = 7;
_13nodePyramid.basis[6] = 9;
_13nodePyramid.basis[7] = 10;
_13nodePyramid.basis[8] = 17;
_13nodePyramid.basis[9] = 18;
_13nodePyramid.basis[10] = 21;
_13nodePyramid.basis[11] = 22;
_13nodePyramid.basis[12] = 33;
_13nodePyramid.gaussPoints[0] = 27;
_13nodePyramid.gaussPoints[1] = 125;
_13nodePyramid.stabilization = 0.333333333333333333;

//----------------------
// 6 nodes wedge element
//----------------------
_6nodeWedge.dimension = 3;
_6nodeWedge.topology = "wedge";
_6nodeWedge.code = 706;
_6nodeWedge.nodes = 6;
//          1     2     3     4     5     6
_6nodeWedge.nodeU[0] = 0.0;
_6nodeWedge.nodeU[1] = 1.0;
_6nodeWedge.nodeU[2] = 0.0;
_6nodeWedge.nodeU[3] = 0.0;
_6nodeWedge.nodeU[4] = 1.0;
_6nodeWedge.nodeU[5] = 0.0;
_6nodeWedge.nodeV[0] = 0.0; 
_6nodeWedge.nodeV[1] = 0.0;
_6nodeWedge.nodeV[2] = 1.0;
_6nodeWedge.nodeV[3] = 0.0;
_6nodeWedge.nodeV[4] = 0.0;
_6nodeWedge.nodeV[5] = 1.0;
_6nodeWedge.nodeW[0] = 1.0;
_6nodeWedge.nodeW[1] = -1.0;
_6nodeWedge.nodeW[2] = -1.0;
_6nodeWedge.nodeW[3] = 1.0;
_6nodeWedge.nodeW[4] = 1.0;
_6nodeWedge.nodeW[5] = 1.0;
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
_6nodeWedge.basis[0] = 1;
_6nodeWedge.basis[1] = 2;
_6nodeWedge.basis[2] = 5;
_6nodeWedge.basis[3] = 17;
_6nodeWedge.basis[4] = 18;
_6nodeWedge.basis[5] = 21;

_6nodeWedge.gaussPoints[0] = 8;
_6nodeWedge.gaussPoints[1] = 64;
_6nodeWedge.stabilization = 0.333333333333333333;

//---------------------------
// 8 nodes octahedron element
//---------------------------
_8nodeOctahedron.dimension = 3;
_8nodeOctahedron.topology = "brick";
_8nodeOctahedron.code = 808;
_8nodeOctahedron.nodes = 8;
//          1     2     3     4     5     6     7     8
_8nodeOctahedron.nodeU[0] = -1.0; 
_8nodeOctahedron.nodeU[1] = 1.0;
_8nodeOctahedron.nodeU[2] = 1.0;
_8nodeOctahedron.nodeU[3] = -1.0;
_8nodeOctahedron.nodeU[4] = -1.0;
_8nodeOctahedron.nodeU[5] = 1.0;
_8nodeOctahedron.nodeU[6] = 1.0;
_8nodeOctahedron.nodeU[7] = -1.0;
_8nodeOctahedron.nodeV[0] = -1.0;
_8nodeOctahedron.nodeV[1] = -1.0;
_8nodeOctahedron.nodeV[2] = 1.0;
_8nodeOctahedron.nodeV[3] = 1.0;
_8nodeOctahedron.nodeV[4] = -1.0;
_8nodeOctahedron.nodeV[5] = -1.0;
_8nodeOctahedron.nodeV[6] = 1.0;
_8nodeOctahedron.nodeV[7] = 1.0;
_8nodeOctahedron.nodeW[0] = -1.0;
_8nodeOctahedron.nodeW[1] = -1.0;
_8nodeOctahedron.nodeW[2] = -1.0;
_8nodeOctahedron.nodeW[3] = -1.0;
_8nodeOctahedron.nodeW[4] = 1.0;
_8nodeOctahedron.nodeW[5] = 1.0;
_8nodeOctahedron.nodeW[6] = 1.0;
_8nodeOctahedron.nodeW[7] = 1.0;
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
_8nodeOctahedron.basis[0] = 1;
_8nodeOctahedron.basis[1] = 2;
_8nodeOctahedron.basis[2] = 5;
_8nodeOctahedron.basis[3] = 6;
_8nodeOctahedron.basis[4] = 17;
_8nodeOctahedron.basis[5] = 18;
_8nodeOctahedron.basis[6] = 21;
_8nodeOctahedron.basis[7] = 22;

_8nodeOctahedron.gaussPoints[0] = 8;
_8nodeOctahedron.gaussPoints[1] = 64;
_8nodeOctahedron.stabilization = 0.166666666666666666667;

//----------------------------
// 20 nodes octahedron element
//----------------------------
_20nodeOctahedron.dimension = 3;
_20nodeOctahedron.topology = "brick";
_20nodeOctahedron.code = 820;
_20nodeOctahedron.nodes = 20;
//          1     2     3     4     5     6      7      8      9     10    11    12    13    14    15   16   17    18    19    20
_20nodeOctahedron.nodeU[0] = -1.0;
_20nodeOctahedron.nodeU[1] = 1.0;
_20nodeOctahedron.nodeU[2] = 1.0;
_20nodeOctahedron.nodeU[3] = -1.0;
_20nodeOctahedron.nodeU[4] = -1.0;
_20nodeOctahedron.nodeU[5] = 1.0;
_20nodeOctahedron.nodeU[6] = 1.0;
_20nodeOctahedron.nodeU[7] = -1.0;
_20nodeOctahedron.nodeU[8] = 0.0;
_20nodeOctahedron.nodeU[9] = 1.0;
_20nodeOctahedron.nodeU[10] = 0.0;
_20nodeOctahedron.nodeU[11] = -1.0;
_20nodeOctahedron.nodeU[12] = -1.0;
_20nodeOctahedron.nodeU[13] = 1.0;
_20nodeOctahedron.nodeU[14] = 1.0;
_20nodeOctahedron.nodeU[15] = -1.0;
_20nodeOctahedron.nodeU[16] = 0.0;
_20nodeOctahedron.nodeU[17] = 1.0;
_20nodeOctahedron.nodeU[18] = 0.0;
_20nodeOctahedron.nodeU[19] = -1.0;
_20nodeOctahedron.nodeV[0] = -1.0;
_20nodeOctahedron.nodeV[1] = -1.0;
_20nodeOctahedron.nodeV[2] = 1.0;
_20nodeOctahedron.nodeV[3] = 1.0;
_20nodeOctahedron.nodeV[4] = -1.0;
_20nodeOctahedron.nodeV[5] = -1.0;
_20nodeOctahedron.nodeV[6] = 1.0;
_20nodeOctahedron.nodeV[7] = 1.0;
_20nodeOctahedron.nodeV[8] = -1.0;
_20nodeOctahedron.nodeV[9] = 0.0;
_20nodeOctahedron.nodeV[10] = 1.0;
_20nodeOctahedron.nodeV[11] = 0.0;
_20nodeOctahedron.nodeV[12] = -1.0;
_20nodeOctahedron.nodeV[13] = -1.0;
_20nodeOctahedron.nodeV[14] = 1.0;
_20nodeOctahedron.nodeV[15] = 1.0;
_20nodeOctahedron.nodeV[16] = -1.0;
_20nodeOctahedron.nodeV[17] = 0.0;
_20nodeOctahedron.nodeV[18] = 1.0;
_20nodeOctahedron.nodeV[19] = 0.0;
_20nodeOctahedron.nodeW[0] = -1.0;
_20nodeOctahedron.nodeW[1] = -1.0;
_20nodeOctahedron.nodeW[2] = -1.0;
_20nodeOctahedron.nodeW[3] = -1.0;
_20nodeOctahedron.nodeW[4] = 1.0;
_20nodeOctahedron.nodeW[5] = 1.0;
_20nodeOctahedron.nodeW[6] = 1.0;
_20nodeOctahedron.nodeW[7] = 1.0;
_20nodeOctahedron.nodeW[8] = -1.0;
_20nodeOctahedron.nodeW[9] = -1.0;
_20nodeOctahedron.nodeW[10] = -1.0;
_20nodeOctahedron.nodeW[11] = -1.0;
_20nodeOctahedron.nodeW[12] = 0.0;
_20nodeOctahedron.nodeW[13] = 0.0;
_20nodeOctahedron.nodeW[14] = 0.0;
_20nodeOctahedron.nodeW[15] = 0.0;
_20nodeOctahedron.nodeW[16] = 1.0;
_20nodeOctahedron.nodeW[17] = 1.0;
_20nodeOctahedron.nodeW[18] = 1.0;
_20nodeOctahedron.nodeW[19] = 1.0;
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
_20nodeOctahedron.basis[0] = 1;
_20nodeOctahedron.basis[1] = 2;
_20nodeOctahedron.basis[2] = 3;
_20nodeOctahedron.basis[3] = 5;
_20nodeOctahedron.basis[4] = 6;
_20nodeOctahedron.basis[5] = 7;
_20nodeOctahedron.basis[6] = 9;
_20nodeOctahedron.basis[7] = 10;
_20nodeOctahedron.basis[8] = 17;
_20nodeOctahedron.basis[9] = 18;
_20nodeOctahedron.basis[10] = 19;
_20nodeOctahedron.basis[11] = 21;
_20nodeOctahedron.basis[12] = 22;
_20nodeOctahedron.basis[13] = 23;
_20nodeOctahedron.basis[14] = 25;
_20nodeOctahedron.basis[15] = 26;
_20nodeOctahedron.basis[16] = 33;
_20nodeOctahedron.basis[17] = 34;
_20nodeOctahedron.basis[18] = 37;
_20nodeOctahedron.basis[19] = 38;

_20nodeOctahedron.gaussPoints[0] = 27;
_20nodeOctahedron.gaussPoints[1] = 125;
_20nodeOctahedron.stabilization = 0.08148148148148; 

//----------------------------
// 27 nodes octahedron element
//----------------------------
_27nodeOctahedron.dimension = 3;
_27nodeOctahedron.topology = "brick";
_27nodeOctahedron.code = 827;
_27nodeOctahedron.nodes = 27;
//          1     2     3     4     5     6      7      8      9     10    11    12    13    14    15   16   17    18    19    20    21    22    23    24    25
//          26    27
_27nodeOctahedron.nodeU[0] = -1.0;
_27nodeOctahedron.nodeU[1] = 1.0;
_27nodeOctahedron.nodeU[2] = 1.0;
_27nodeOctahedron.nodeU[3] = -1.0;
_27nodeOctahedron.nodeU[4] = -1.0;
_27nodeOctahedron.nodeU[5] = 1.0;
_27nodeOctahedron.nodeU[6] = 1.0;
_27nodeOctahedron.nodeU[7] = -1.0;
_27nodeOctahedron.nodeU[8] = 0.0;
_27nodeOctahedron.nodeU[9] = 1.0;
_27nodeOctahedron.nodeU[10] = 0.0;
_27nodeOctahedron.nodeU[11] = -1.0;
_27nodeOctahedron.nodeU[12] = -1.0;
_27nodeOctahedron.nodeU[13] = 1.0;
_27nodeOctahedron.nodeU[14] = 1.0;
_27nodeOctahedron.nodeU[15] = -1.0;
_27nodeOctahedron.nodeU[16] = 0.0;
_27nodeOctahedron.nodeU[17] = 1.0;
_27nodeOctahedron.nodeU[18] = 0.0;
_27nodeOctahedron.nodeU[19] = -1.0;
_27nodeOctahedron.nodeU[20] = 0.0;
_27nodeOctahedron.nodeU[21] = 1.0;
_27nodeOctahedron.nodeU[22] = 0.0;
_27nodeOctahedron.nodeU[23] = -1.0;
_27nodeOctahedron.nodeU[24] = 0.0;
_27nodeOctahedron.nodeU[25] = 0.0;
_27nodeOctahedron.nodeU[26] = 0.0;
_27nodeOctahedron.nodeV[0] = -1.0;
_27nodeOctahedron.nodeV[1] = -1.0;
_27nodeOctahedron.nodeV[2] = 1.0;
_27nodeOctahedron.nodeV[3] = 1.0;
_27nodeOctahedron.nodeV[4] = -1.0;
_27nodeOctahedron.nodeV[5] = -1.0;
_27nodeOctahedron.nodeV[6] = 1.0;
_27nodeOctahedron.nodeV[7] = 1.0;
_27nodeOctahedron.nodeV[8] = -1.0;
_27nodeOctahedron.nodeV[9] = 0.0;
_27nodeOctahedron.nodeV[10] = 1.0;
_27nodeOctahedron.nodeV[11] = 0.0;
_27nodeOctahedron.nodeV[12] = -1.0;
_27nodeOctahedron.nodeV[13] = -1.0;
_27nodeOctahedron.nodeV[14] = 1.0;
_27nodeOctahedron.nodeV[15] = 1.0;
_27nodeOctahedron.nodeV[16] = -1.0;
_27nodeOctahedron.nodeV[17] = 0.0;
_27nodeOctahedron.nodeV[18] = 1.0;
_27nodeOctahedron.nodeV[19] = 0.0;
_27nodeOctahedron.nodeV[20] = -1.0;
_27nodeOctahedron.nodeV[21] = 0.0;
_27nodeOctahedron.nodeV[22] = 1.0;
_27nodeOctahedron.nodeV[23] = 0.0;
_27nodeOctahedron.nodeV[24] = 0.0;
_27nodeOctahedron.nodeV[25] = 0.0;
_27nodeOctahedron.nodeV[26] = 0.0;
_27nodeOctahedron.nodeW[0] = -1.0;
_27nodeOctahedron.nodeW[1] = -1.0;
_27nodeOctahedron.nodeW[2] = -1.0;
_27nodeOctahedron.nodeW[3] = -1.0;
_27nodeOctahedron.nodeW[4] = 1.0;
_27nodeOctahedron.nodeW[5] = 1.0;
_27nodeOctahedron.nodeW[6] = 1.0;
_27nodeOctahedron.nodeW[7] = 1.0;
_27nodeOctahedron.nodeW[8] = -1.0;
_27nodeOctahedron.nodeW[9] = -1.0;
_27nodeOctahedron.nodeW[10] = -1.0;
_27nodeOctahedron.nodeW[11] = -1.0;
_27nodeOctahedron.nodeW[12] = 0.0;
_27nodeOctahedron.nodeW[13] = 0.0;
_27nodeOctahedron.nodeW[14] = 0.0;
_27nodeOctahedron.nodeW[15] = 0.0;
_27nodeOctahedron.nodeW[16] = 1.0;
_27nodeOctahedron.nodeW[17] = 1.0;
_27nodeOctahedron.nodeW[18] = 1.0;
_27nodeOctahedron.nodeW[19] = 1.0;
_27nodeOctahedron.nodeW[20] = 0.0;
_27nodeOctahedron.nodeW[21] = 0.0;
_27nodeOctahedron.nodeW[22] = 0.0;
_27nodeOctahedron.nodeW[23] = 0.0;
_27nodeOctahedron.nodeW[24] = -1.0;
_27nodeOctahedron.nodeW[25] = 1.0;
_27nodeOctahedron.nodeW[26] = 0.0;
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
_27nodeOctahedron.basis[0] = 1;
_27nodeOctahedron.basis[1] = 2;
_27nodeOctahedron.basis[2] = 3;
_27nodeOctahedron.basis[3] = 5;
_27nodeOctahedron.basis[4] = 6;
_27nodeOctahedron.basis[5] = 7;
_27nodeOctahedron.basis[6] = 9;
_27nodeOctahedron.basis[7] = 10;
_27nodeOctahedron.basis[8] = 11;
_27nodeOctahedron.basis[9] = 17;
_27nodeOctahedron.basis[10] = 18;
_27nodeOctahedron.basis[11] = 19;
_27nodeOctahedron.basis[12] = 21;
_27nodeOctahedron.basis[13] = 22;
_27nodeOctahedron.basis[14] = 23;
_27nodeOctahedron.basis[15] = 25;
_27nodeOctahedron.basis[16] = 26;
_27nodeOctahedron.basis[17] = 27;
_27nodeOctahedron.basis[18] = 33;
_27nodeOctahedron.basis[19] = 34;
_27nodeOctahedron.basis[20] = 35;
_27nodeOctahedron.basis[21] = 37;
_27nodeOctahedron.basis[22] = 38;
_27nodeOctahedron.basis[23] = 39;
_27nodeOctahedron.basis[24] = 41;
_27nodeOctahedron.basis[25] = 42;
_27nodeOctahedron.basis[26] = 43;

_27nodeOctahedron.gaussPoints[0] = 27;
_27nodeOctahedron.gaussPoints[1] = 125;
_27nodeOctahedron.stabilization = 0.08148148148148;

#endif
