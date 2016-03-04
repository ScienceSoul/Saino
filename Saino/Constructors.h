//
//  Constructors.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

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
    double * __nullable x;                  // First coordinate
    double * __nullable y;                  // Second coordinate
    double * __nullable z;                  // Third coordinate
    
} Nodes_t;

typedef struct {
    
    int n;
    int * __nullable p, * __nullable q, * __nullable r;
    double * __nullable coeff;
    
} BasisFunctions_t;

typedef struct ElementType_t {
    
    struct ElementType_t * __nullable NextElementType;
    
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
    double * __nullable NodeU, * __nullable NodeV, * __nullable NodeW;         // They have size of NumberOfNodes
    
    BasisFunctions_t * __nullable BasisFunctions;
    
} ElementType_t;

typedef struct {
    
    int NumberOfFactors;
    int NumberOfImplicitFactors;
    int * __nullable Elements;
    double * __nullable Factors;
    int sizeElements, sizeFactors;
    
} Factors_t;

typedef struct {
    
    Factors_t * __nullable GebhardtFactors;
    int Constraint;                         // Initialize it to 0 somewhere!! (Index counted from 1 to n)
    int Outbody;                            // Initialize to -1 somewhere!!
    struct Element_t * __nullable Left;    // Initialize to NULL somewhere!!
    struct Element_t * __nullable Right;   // Initialize to NULL somewhere!!
    
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
    BoundaryInfo_t * __nullable BoundaryInfo;
    PElementDefs_t * __nullable Pdefs;                         // Initialize to NULL somewhere!!
    
    RGBColors color;                                           // Color of the element
    
    bool colored;
    bool copy;
    int * __nullable NodeIndexes, * __nullable EdgeIndexes,
        * __nullable FaceIndexes, * __nullable BubbleIndexes,
        * __nullable DGIndexes;                                // Initialize that to NULL somewhere!!
    
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
    struct Quadrant_t * __nullable quadrant;
} QuadrantPointer_t;

typedef struct Quadrant_t {
    
    int nElementsInQuadrant;
    double size, minElementSize, boundingBox[6];
    int * __nullable elements;
    struct QuadrantPointer_t * __nullable childQuadrants;
    int numberOfchildQuadrants;
} Quadrant_t;

typedef struct {
    //TODO
} Header;

typedef struct ListMatrixEntry_t {
    
    int Index;
    double Value;
    struct ListMatrixEntry_t * __nullable Next;
    
} ListMatrixEntry_t;

typedef struct {
    
    int sizeOfContainer;     // If ListMatrix_t is contained within an array, the size of this array is given by
                             // this variable
    int Degree, Level;
    ListMatrixEntry_t * __nullable Head;
    
} ListMatrix_t;

typedef struct {
    ListMatrixEntry_t * __nullable p;
} Stack_t;

typedef struct matrixArraysContainer {
    
    ListMatrix_t * __nullable ListMatrix;
    
    int * __nullable Perm, * __nullable InvPerm, * __nullable RowOwner;
    int * __nullable GRows, * __nullable GOrder;
    int * __nullable Rows, * __nullable Cols, * __nullable Diag;
    int sizePerm, sizeInvPerm, sizeRowOwner;
    int sizeGRows, sizeGOrder;
    int sizeRows, sizeCols, sizeDiag;
    
    double * __nullable RHS, * __nullable BulkRHS, * __nullable RHS_im, * __nullable * __nullable Force;
    int sizeRHS, sizeBulkRHS, sizeRHS_im, size1force, size2Force;
    
    double * __nullable Values, * __nullable ILUValues;
    double * __nullable MassValues, * __nullable DampValues, * __nullable BulkValues, * __nullable DiagScaling;
    int sizeValues, sizeILUValues, sizeMassValues, sizeDampValues, sizeBulkValues, sizeDiagScaling;
    
    int * __nullable ILURows, * __nullable ILUCols, * __nullable ILUDiag;
    int sizeILURows, sizeILUCols, sizeILUDiag;
    
    // For complex system
    double complex * __nullable CRHS, * __nullable CForce;
    double complex * __nullable CValues, * __nullable CILUValues;
    double complex * __nullable CMassValues, * __nullable CDampValues;
    int sizeCRHS, sizeCForce, sizeCValues, sizeCILUValues,
        sizeCMassValues, sizeCDampValues;
    
    // For flux corrected transport
    double * __nullable FCT_D;
    double * __nullable MassValuesLumped;
    int sizeFct, sizeMassValuesLumped;
    
} matrixArraysContainer;

typedef struct variableArraysContainer {
    
    int * __nullable Perm;
    double * __nullable Values;
    
    double * __nullable * __nullable ComponentValues;                            // Used if the variable is a component of a another variable with dofs > 1.
                                                                                 // This is a 1D array of pointers.
    double * __nullable * __nullable PrevValues;
    
    double * __nullable * __nullable * __nullable ComponentPrevValues;           // Used if the variable is a component of a another variable with dofs > 1.
                                                                                 // This is a 2D array of pointers.
    
    double * __nullable * __nullable SecondaryToValues;                          // Some variables may be defined so that their values point to a given column of
                                                                                 // another variable PrevValues. This is used for that particular case
    double * __nullable * __nullable ComponentSecondaryToValues;
    double * __nullable PValues;
    double * __nullable NonLinValues;
    double * __nullable SteadyValues;
    double complex * __nullable EigenValues;
    double complex * __nullable * __nullable EigenVectors;
    double complex * __nullable * __nullable * __nullable ComponentEigenVectors;  // This is a 2D array of pointers.
    double complex * __nullable * __nullable CValues;                             // This is a 1D array of pointers
    bool * __nullable lowerLimitActive;
    bool * __nullable upperLimitActive;
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
    
    int * __nullable activeElements;
    int * __nullable * __nullable defDofs;
    int sizeActiveElements;
    int size1DefDofs;
    int size2DefDofs;
} solutionArraysContainer;

typedef struct valueListArraysContainer {
    
    int * __nullable iValues;
    double * __nullable tValues;
    double * __nullable * __nullable * __nullable fValues;
    double * __nullable cubicCoeff;
    int sizeTValues;
    int sizeFValues1;
    int sizeFValues2;
    int sizeFValues3;
    int sizeIValues;
    int sizeCubicCoeff;
    
} valueListArraysContainer;

typedef struct modelArraysContainer {
    
    int * __nullable freeSurfaceNodes;
    int * __nullable rowNonZeros;
    double * __nullable boundaryCurvatures;
    int sizeFreeSurfaceNodes;
    int sizeRowNonZeros;
    int sizeBoundaryCurvatures;
    
} modelArraysContainer;

typedef struct listBuffer {
    
    int * __nullable ivector;
    double * __nullable vector;
    double * __nullable * __nullable matrix;
    double * __nullable * __nullable * __nullable tensor;
    int m, n, p;
    
} listBuffer;

typedef struct HashEntry_t {
    
    int node1, node2, face, edge;
    struct HashEntry_t * __nullable next;
     
} HashEntry_t;

typedef struct HashTable_t {
    
    HashEntry_t * __nullable head;

} HashTable_t;

typedef struct RungeKutta_t {

    double * __nullable k1;
    double * __nullable k2;
    double * __nullable k3;
    double * __nullable k4;
    
} RungeKutta_t;

typedef struct Dimensions_t {
    
    int mat1;
    int mat2;
    int vec;
    
} Dimensions_t;

#endif

void initNodes(Nodes_t * __nonnull nodes);
void initElements(Element_t * __nonnull elements, int n);
void initBoundaryInfo(BoundaryInfo_t * __nonnull boundaryInfo);
variableArraysContainer * __nonnull allocateVariableContainer(void);
RungeKutta_t * __nonnull allocateRungeKutta(int n);

