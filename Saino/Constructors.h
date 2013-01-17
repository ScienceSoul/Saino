//
//  Constructors.h
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#include <complex.h>
#include <stdbool.h>

#ifndef CONSTRUCTORS_H
#define CONSTRUCTORS_H

#define MAX_NAME_LEN 128
#define MAX_ELEMENT_NODES 256

enum {
    MATRIX_CRS = 1,
    MATRIX_BAND,
    MATRIX_SBAND,
    MATRIX_LIST
};

enum  {
    LIST_TYPE_CONSTANT_SCALAR = 1,
    LIST_TYPE_CONSTANT_TENSOR,
    LIST_TYPE_VARIABLE_SCALAR,
    LIST_TYPE_VARIABLE_TENSOR,
    LIST_TYPE_LOGICAL,
    LIST_TYPE_STRING,
    LIST_TYPE_INTEGER,
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

typedef struct {
    
    double *x;                   // First coordinate
    double *y;                   // Second coordinate
    double *z;                   // Third coordinate
    
} Nodes_t;

typedef struct {
    
    int n;
    int *p, *q, *r;
    double *coeff;
    
    
} BasisFunctions_t;

typedef struct QuadrantPointer_t {
    struct Quadrant_t *quadrant;
} QuadrantPointer_t;

typedef struct Quadrant_t {
    
    int nElementsInQuadrant;
    double size, minElementSize, boundingBox[6];
    int *elements;
    struct QuadrantPointer_t *childQuadrants;
    int numberOfchildQuadrants;
} Quadrant_t;

typedef struct ElementType_t {
    
    struct ElementType_t *NextElementType; // List of types
    
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
    double *NodeU, *NodeV, *NodeW;         // They have size of NumberOfNodes
    
    BasisFunctions_t *BasisFunctions;
    
} ElementType_t;

typedef struct {
    
    int NumberOfFactors;
    int NumberOfImplicitFactors;
    int *Elements;
    double *Factors;
    
} Factors_t;

typedef struct {
    
    Factors_t *GebhardtFactors;
    int Constraint;            // Initialize it to 0 somewhere!!
    int Outbody;               // Initialize to -1 somewhere!!
    struct Element_t *Left;    // Initialize to NULL somewhere!!
    struct Element_t *Right;   // Initialize to NULL somewhere!!
    
} BoundaryInfo_t;

typedef struct {

    int p;  
    int TetraType;          // Type of p tetrahedron={0,1,2}
    bool isEdge;            // Is element an edge or face.
    int GaussPoints;        // Number of gauss points to use when using p elements
    bool PyramidQuadEdge;   // Is element an edge pyramid quad face.
    int LocalNumber;        // Local number of an edge or face for element on boundary
    
} PElementDefs_t;

typedef struct Element_t {
    
    ElementType_t Type;
    BoundaryInfo_t *BoundaryInfo;
    PElementDefs_t *Pdefs;                         // Initialize to NULL somewhere!!
    
    int *NodeIndexes, *EdgeIndexes, *FaceIndexes, 
        *BubbleIndexes, *DGIndexes;                // Initialize that to NULL somewhere!!
    
    int sizeNodeIndexes, sizeEdgeIndexes,
        sizeFaceIndexes, sizeBubbleIndexes, sizeDGIndexes;
    int BodyID;
    int ElementIndex, PartIndex, NDOFs, BDOFs, DGDOFs;
    
    double StabilizationMK, hK;
    
    
} Element_t;

typedef struct {               // Boundary element
    
    int ID;                    // Boundary element identification number 
    int boundary;              // Part of the boundary where element is located
    int p1;                    // First parent element
    int p2;                    // Second parent element
    int code;                  // Element code
    int *NodeIndexes;          // Element nodes
    
} BDElement_t;

typedef struct {
    //TODO
} Header;

typedef struct ListMatrixEntry_t {
    
    int Index;
    double Value;
    struct ListMatrixEntry_t *Next;
    
} ListMatrixEntry_t;

typedef struct {
    
    int sizeOfContainer;     // If ListMatrix_t is contained within an array, the size of this array is given by
                             // this variable
    int Degree, Level;
    ListMatrixEntry_t *Head;
    
} ListMatrix_t;

typedef struct {
    ListMatrixEntry_t *p;
} Stack_t;

typedef struct matrixArraysContainer {
    
    ListMatrix_t *ListMatrix;
    
    int *Perm, *InvPerm, *RowOwner;
    int *GRows, *GOrder;
    int *Rows, *Cols, *Diag;
    int sizePerm, sizeInvPerm, sizeRowOwner;
    int sizeGRows, sizeGOrder;
    int sizeRows, sizeCols, sizeDiag;
    
    double *RHS, *BulkRHS, *RHS_im, **Force;
    int sizeRHS, sizeBulkRHS, sizeRHS_im, size1force, size2Force;
    
    double *Values, *ILUValues;
    double *MassValues, *DampValues, *BulkValues;
    int sizeValues, sizeILUValues, sizeMassValues, sizeDampValues, sizeBulkValues;
    
    int *ILURows, *ILUCols, *ILUDiag;
    int sizeILURows, sizeILUCols, sizeILUDiag;
    
    // For complex system
    double complex *CRHS, *CForce;
    double complex *CValues, *CILUValues;
    double complex *CMassValues, *CDampValues;
    int sizeCRHS, sizeCForce, sizeCValues, sizeCILUValues,
        sizeCMassValues, sizeCDampValues;
    
} matrixArraysContainer;

typedef struct variableArraysContainer {
    
    int *Perm;
    double *Values;
    double **PrevValues;
    double *PValues;
    double *NonLinValues;
    double *SteadyValues;
    double complex *EigenValues;
    double complex **EigenVectors;
    int sizePerm;
    int sizeValues;
    int size1PrevValues;
    int size2PrevValues;
    int sizePValues;
    int sizeNonLinValues;
    int sizeSteadyValues;
    
} variableArraysContainer;

typedef struct solutionArraysContainer {
    
    bool **ntZeroingDone;
    int *activeElements;
    int **ntElement;
    int *boundaryReorder;
    int **defDofs;
    double **boundaryNormals;
    double **boundaryTangent1;
    double **boundaryTangent2;
    int sizeBoundaryReorder;
    int size1boundaryNormals;
    int size2boundaryNormals;
    int size1DefDofs;
    int size2DefDofs;
    
} solutionArraysContainer;

typedef struct valueListArraysContainer {
    
    int *iValues;
    double *tValues;
    double ***fValues;
    int sizeTValues;
    int sizeFValues1;
    int sizeFValues2;
    int sizeFValues3;
    int sizeIValues;
    
} valueListArraysContainer;

typedef struct modelArraysContainer {
    
    int *freeSurfaceNodes;
    int *rowNonZeros;
    double *boundaryCurvatures;
    int sizeFreeSurfaceNodes;
    int sizeRowNonZeros;
    int sizeBoundaryCurvatures;
    
} modelArraysContainer;

typedef struct listBuffer {
    
    int *ivector;
    double *vector;
    double **matrix;
    double ***tensor;
    int m, n, p;
    
} listBuffer;

typedef struct HashEntry_t {
    
    int node1, node2, face, edge;
    struct HashEntry_t *next;
     
} HashEntry_t;

typedef struct HashTable_t {
    
    HashEntry_t *head;

} HashTable_t;



#endif

Nodes_t *nodesvec(long nl, long nh);
Element_t *elementsvec(long nl, long nh);
BDElement_t *BDElementVec(long nl, long nh);