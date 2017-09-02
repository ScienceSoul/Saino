//===----------------------------------------------------------------------===//
//  Constructors.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright (c) 2011 ScienceSoul. All rights reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

#import <float.h>
#import <complex.h>
#import <stdbool.h>

#ifndef CONSTRUCTORS_H
#define CONSTRUCTORS_H

#define MAX_NAME_LEN 128
#define MAX_ELEMENT_NODES 256

#define AEPS 10.0 * DBL_EPSILON

enum {
    MATRIX_CRS = 1,
    MATRIX_BAND,
    MATRIX_SBAND,
    MATRIX_LIST
};

enum {
    SOLUTION_SOLVE_NEVER = -1,
    SOLUTION_SOLVE_ALWAYS,
    SOLUTION_SOLVE_AHEAD_ALL,
    SOLUTION_SOLVE_AHEAD_TIME,
    SOLUTION_SOLVE_AFTER_ALL,
    SOLUTION_SOLVE_AFTER_TIME,
    SOLUTION_SOLVE_AHEAD_SAVE,
    SOLUTION_SOLVE_AFTER_SAVE
};


enum  {
    LIST_TYPE_CONSTANT_SCALAR = 1,
    LIST_TYPE_CONSTANT_TENSOR,
    LIST_TYPE_VARIABLE_SCALAR,
    LIST_TYPE_VARIABLE_TENSOR,
    LIST_TYPE_LOGICAL,
    LIST_TYPE_STRING,
    LIST_TYPE_INTEGER,
    List_TYPE_BLOCK,
    LIST_TYPE_CONSTANT_SCALAR_STR,
    LIST_TYPE_CONSTANT_TENSIOR_STR,
    LIST_TYPE_VARIABLE_SCALAR_STR,
    LIST_TYPE_VARIABLE_TENSOR_STR
};

enum {
    VARIABLE_ON_NODES = 0,
    VARIABLE_ON_EDGES,
    VARIABLE_ON_FACES,
    VARIABLE_ON_NODES_ON_ELEMENTS
};

enum {
    cartesian = 1,
    cylindric,
    cylindric_symmetric,
    axis_symmetric,
    polar
};

enum {
    incompressible = 0,
    user_defined1,
    user_defined2,
    perfect_gas1,
    perfect_gas2,
    perfect_gas3,
    thermal
};

enum {
    SOLUTION_MODE_DEFAULT = 0,  // Normal DPE
    SOLUTION_MODE_AUXILIARY,    // No FEM machinery (SaveData)
    SOLUTION_MODE_ASSEMBLY,     // Coupled solution with single block
    SOLUTION_MODE_COUPLED,      // Coupled solution with multiple blocks
    SOLUTION_MODE_BLOCK,        // Block solution
    SOLUTION_MODE_GLOBAL,       // Lumped variables (no mesh)
    SOLUTION_MODE_MATRIXFREE    // Normal field, no matrix
};

enum {
    PROJECTOR_TYPE_DEFAULT = 0, // Unspecified constraint matrix
    PROJECTOR_TYPE_NODAL,       // Nodal projector
    PROJECTOR_TYPE_GALERKIN     // Galerkin projector
};

typedef struct {
    
    int numberOfNodes;
    double * _Nullable x;                  // First coordinate
    double * _Nullable y;                  // Second coordinate
    double * _Nullable z;                  // Third coordinate
    
} Nodes_t;

typedef struct {
    
    int n;
    int * _Nullable p, * _Nullable q, * _Nullable r;
    double * _Nullable coeff;
    
} BasisFunctions_t;

typedef struct ElementType_t {
    
    struct ElementType_t * _Nullable NextElementType;
    
    int ElementCode;                       // Numeric code for element
    int BasisFunctionDegree,               // Linear or quadratic
        NumberOfNodes,
        NumberOfEdges,
        NumberOfFaces,
        dimension;                         // 1=Line; 2=Surface; 3=Volume
    
    int GaussPoints,                       // Number of Gauss points to use
        GaussPoints2,
        GaussPoints0;
    
    double StabilizationMK;
    double * _Nullable NodeU, * _Nullable NodeV, * _Nullable NodeW;         // They have size of NumberOfNodes
    
    BasisFunctions_t * _Nullable BasisFunctions;
    
} ElementType_t;

typedef struct {
    
    int NumberOfFactors;
    int NumberOfImplicitFactors;
    int * _Nullable Elements;
    double * _Nullable Factors;
    int sizeElements, sizeFactors;
    
} Factors_t;

typedef struct {
    
    Factors_t * _Nullable GebhardtFactors;
    int Constraint;                         // Initialize it to 0 somewhere!! (Index counted from 1 to n)
    int Outbody;                            // Initialize to -1 somewhere!!
    struct Element_t * _Nullable Left;    // Initialize to NULL somewhere!!
    struct Element_t * _Nullable Right;   // Initialize to NULL somewhere!!
    
} BoundaryInfo_t;

typedef struct {
    
    int p;
    int TetraType;          // Type of p tetrahedron={0,1,2}
    bool isEdge;            // Is element an edge or face.
    int GaussPoints;        // Number of gauss points to use when using p elements
    bool PyramidQuadEdge;   // Is element an edge pyramid quad face.
    int LocalNumber;        // Local number of an edge or face for element on boundary
    
} PElementDefs_t;

typedef struct {
    float red;              // The red component of the color object, specified as a value from 0.0 to 1.0.
    float green;            // The green component of the color object, specified as a value from 0.0 to 1.0.
    float blue;             // The blue component of the color object, specified as a value from 0.0 to 1.0.
    int colorIndex;
} RGBColors;

typedef struct Element_t {
    
    ElementType_t Type;
    BoundaryInfo_t * _Nullable BoundaryInfo;
    PElementDefs_t * _Nullable Pdefs;                         // Initialize to NULL somewhere!!
    
    RGBColors color;                                           // Color of the element
    
    bool colored;
    bool copy;
    int * _Nullable NodeIndexes, * _Nullable EdgeIndexes,
        * _Nullable FaceIndexes, * _Nullable BubbleIndexes,
        * _Nullable DGIndexes;                                // Initialize that to NULL somewhere!!
    
    int sizeNodeIndexes, sizeEdgeIndexes,
        sizeFaceIndexes, sizeBubbleIndexes,
        sizeDGIndexes;
    int BodyID;
    int Splitted;
    int ElementIndex,                                          // Index counted from 1 to n
        GElementIndex,                                         // Index counted from 1 to n
        PartIndex, NDOFs, BDOFs, DGDOFs;
    
    double StabilizationMK, hK;
    
} Element_t;

typedef struct QuadrantPointer_t {
    struct Quadrant_t * _Nullable quadrant;
} QuadrantPointer_t;

typedef struct Quadrant_t {
    
    int nElementsInQuadrant;
    double size, minElementSize, boundingBox[6];
    int * _Nullable elements;
    struct QuadrantPointer_t * _Nullable childQuadrants;
    int numberOfchildQuadrants;
} Quadrant_t;

typedef struct {
    //TODO
} Header;

typedef struct ListMatrixEntry_t {
    
    int Index;
    double Value;
    struct ListMatrixEntry_t * _Nullable Next;
    
} ListMatrixEntry_t;

typedef struct {
    
    int sizeOfContainer;     // If ListMatrix_t is contained within an array, the size of this array is given by
                             // this variable
    int Degree, Level;
    ListMatrixEntry_t * _Nullable Head;
    
} ListMatrix_t;

typedef struct {
    ListMatrixEntry_t * _Nullable p;
} Stack_t;

typedef struct matrixArraysContainer {
    
    ListMatrix_t * _Nullable ListMatrix;
    
    int * _Nullable Perm, * _Nullable InvPerm, * _Nullable RowOwner;
    int * _Nullable GRows, * _Nullable GOrder;
    int * _Nullable Rows, * _Nullable Cols, * _Nullable Diag;
    int sizePerm, sizeInvPerm, sizeRowOwner;
    int sizeGRows, sizeGOrder;
    int sizeRows, sizeCols, sizeDiag;
    
    double * _Nullable RHS, * _Nullable BulkRHS, * _Nullable RHS_im, * _Nullable * _Nullable Force;
    int sizeRHS, sizeBulkRHS, sizeRHS_im, size1force, size2Force;
    
    double * _Nullable Values, * _Nullable ILUValues;
    double * _Nullable MassValues, * _Nullable DampValues, * _Nullable BulkValues, * _Nullable DiagScaling;
    int sizeValues, sizeILUValues, sizeMassValues, sizeDampValues, sizeBulkValues, sizeDiagScaling;
    
    int * _Nullable ILURows, * _Nullable ILUCols, * _Nullable ILUDiag;
    int sizeILURows, sizeILUCols, sizeILUDiag;
    
    // For complex system
    double complex * _Nullable CRHS, * _Nullable CForce;
    double complex * _Nullable CValues, * _Nullable CILUValues;
    double complex * _Nullable CMassValues, * _Nullable CDampValues;
    int sizeCRHS, sizeCForce, sizeCValues, sizeCILUValues,
        sizeCMassValues, sizeCDampValues;
    
    // For flux corrected transport
    double * _Nullable FCT_D;
    double * _Nullable MassValuesLumped;
    int sizeFct, sizeMassValuesLumped;
    
} matrixArraysContainer;

typedef struct variableArraysContainer {
    
    int * _Nullable Perm;
    double * _Nullable Values;
    
    double * _Nullable * _Nullable ComponentValues;                            // Used if the variable is a component of a another variable with dofs > 1.
                                                                                 // This is a 1D array of pointers.
    double * _Nullable * _Nullable PrevValues;
    
    double * _Nullable * _Nullable * _Nullable ComponentPrevValues;           // Used if the variable is a component of a another variable with dofs > 1.
                                                                                 // This is a 2D array of pointers.
    
    double * _Nullable * _Nullable SecondaryToValues;                          // Some variables may be defined so that their values point to a given column of
                                                                                 // another variable PrevValues. This is used for that particular case
    double * _Nullable * _Nullable ComponentSecondaryToValues;
    double * _Nullable PValues;
    double * _Nullable NonLinValues;
    double * _Nullable SteadyValues;
    double complex * _Nullable EigenValues;
    double complex * _Nullable * _Nullable EigenVectors;
    double complex * _Nullable * _Nullable * _Nullable ComponentEigenVectors;  // This is a 2D array of pointers.
    double complex * _Nullable * _Nullable CValues;                             // This is a 1D array of pointers
    bool * _Nullable lowerLimitActive;
    bool * _Nullable upperLimitActive;
    int sizePerm;
    int sizeValues;
    int sizeComponentValues;
    int size1PrevValues;
    int size2PrevValues;
    int size1ComponentPrevValues;
    int size2ComponentPrevValues;
    int sizeSecondaryToValues;
    int sizeComponentSecondaryToValues;
    int sizePValues;
    int sizeNonLinValues;
    int sizeSteadyValues;
    int sizeEigenValues;
    int size1EigenVectors;
    int size2EigenVectors;
    int size1ComponentEigenVectors;
    int size2ComponentEigenVectors;
    int sizeCValues;
    int sizeLowerLimitActive;
    int sizeUpperLimitActive;
    
} variableArraysContainer;

typedef struct solutionArraysContainer {
    
    int * _Nullable activeElements;
    int * _Nullable * _Nullable defDofs;
    int sizeActiveElements;
    int size1DefDofs;
    int size2DefDofs;
} solutionArraysContainer;

typedef struct valueListArraysContainer {
    
    int * _Nullable iValues;
    double * _Nullable tValues;
    double * _Nullable * _Nullable * _Nullable fValues;
    double * _Nullable cubicCoeff;
    int sizeTValues;
    int sizeFValues1;
    int sizeFValues2;
    int sizeFValues3;
    int sizeIValues;
    int sizeCubicCoeff;
    
} valueListArraysContainer;

typedef struct modelArraysContainer {
    
    int * _Nullable freeSurfaceNodes;
    int * _Nullable rowNonZeros;
    double * _Nullable boundaryCurvatures;
    int sizeFreeSurfaceNodes;
    int sizeRowNonZeros;
    int sizeBoundaryCurvatures;
    
} modelArraysContainer;

typedef struct listBuffer {
    
    int * _Nullable ivector;
    double * _Nullable vector;
    double * _Nullable * _Nullable matrix;
    double * _Nullable * _Nullable * _Nullable tensor;
    int m, n, p;
    
} listBuffer;

typedef struct HashEntry_t {
    
    int node1, node2, face, edge;
    struct HashEntry_t * _Nullable next;
     
} HashEntry_t;

typedef struct HashTable_t {
    
    HashEntry_t * _Nullable head;

} HashTable_t;

typedef struct RungeKutta_t {

    double * _Nullable k1;
    double * _Nullable k2;
    double * _Nullable k3;
    double * _Nullable k4;
    
} RungeKutta_t;

typedef struct Dimensions_t {
    
    int mat1;
    int mat2;
    int vec;
    
} Dimensions_t;

#endif

void initNodes(Nodes_t * _Nonnull nodes);
void initElements(Element_t * _Nonnull elements, int n);
void initBoundaryInfo(BoundaryInfo_t * _Nonnull boundaryInfo);
variableArraysContainer * _Nonnull allocateVariableContainer(void);
RungeKutta_t * _Nonnull allocateRungeKutta(int n);

