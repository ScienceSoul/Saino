//
//  FEMSolution.h
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMesh.h"

#import "Constructors.h"
#import "memory.h"
#import "Utils.h"

@interface FEMSolution : NSObject {
    
    NSMutableDictionary *solutionInfo;
    
    int normalTangentialNOFNodes;
    int coordinateSystemDimension;
    
    int *boundaryReorder;
    double **boundaryNormals, **boundaryTangent1, **boundaryTangent2;
    
    int sizeBoundaryReorder, size1boundaryNormals, size2boundaryNormals;
    
    int timeOrder, doneTime, order, nofEigenValues;
    int lineBeforeProc, linAfterProc;
    double alpha, beta, dt;
    
    int solutionComputeWhen;
    
    int numberOfActiveElements;
    int *acticeElements, *defDofs;
    int sizeDefDofs;
    
    ValueList_t values;
    Matrix_t matrix;
    Variable_t variable;
    
    Variable_t *exportedVariable;
    
    FEMMesh *mesh;
    
}

#pragma mark Setters and getters for FEMSolution variables

-(id)solutionInfoForKey:(NSString *)key;
-(int)normalTangentialNOFNodes;
-(int)coordinateSystemDimension;
-(int)boundaryReorder:(int)i;
-(double)boundaryNormals:(int)i: (int)j;
-(double)boundaryTangent1:(int)i: (int)j;
-(double)boundaryTangent2:(int)i: (int)j;
-(int)size1boundaryNormals;
-(int)size2boundaryNormals;
-(int)timeOrder;
-(int)doneTime;
-(int)order;
-(int)nofEigenValues;
-(double)alpha;
-(double)beta;
-(double)dt;
-(int)defDofs:(int)i;
-(int)sizeOfDefDofs;

-(void)setNormalTangentialNOFNodes:(int)i;
-(void)setCoordinateSystemDimension:(int)i;
-(void)setBoundaryReorder:(int)i: (int)n;
-(void)setBoundaryNormals:(int)i: (int)j: (double)n;
-(void)setBoundaryTangent1:(int)i: (int)j: (double)n;
-(void)setBoundaryTangent2:(int)i: (int)j: (double)n;
-(void)setSize1boundaryNormals:(int)n;
-(void)setSize2boundaryNormals:(int)n;
-(void)setTimeOrder:(int)n;
-(void)setDoneTime:(int)n;
-(void)setOrder:(int)n;
-(void)setAlpha:(double)n;
-(void)setBeta:(double)n;
-(void)setDt:(double)n;
-(void)setDefDofs:(int)i: (int)n;
-(void)setSizeOfDefDofs:(int)n;


#pragma mark Setters and getters for Matrix structure

-(int)matrixNumberOfRows;
-(int)matrixSubband;
-(int)matrixFormat;
-(int)matrixSolveCount;
-(int)matrixOrdered;
-(int)matrixLumped;
-(int)matrixComplex;
-(int)matrixRows:(int)i;
-(int)matrixCols:(int)i;
-(int)matrixDiag:(int)i;
-(double)matrixRHS:(int)i;
-(double)matrixBulkRHS:(int)i;
-(double)matrixForce:(int)i: (int)j;
-(double)matrixValues:(int)i;
-(double)matrixILUValues:(int)i;
-(double)matrixMassValues:(int)i;
-(double)matrixDampValues:(int)i;
-(double)matrixBulkValues:(int)i;
-(int)matrixILURows:(int)i;
-(int)matrixILUCols:(int)i;
-(int)matrixILUDiag:(int)i;
-(double complex)matrixCILUValues:(int)i;

-(void)setMatrixNumberOfRows:(int)i;
-(void)setMatrixSubband:(int)i;
-(void)setMatrixFormat:(int)i;
-(void)setMatrixSolveCount:(int)i;
-(void)setMatrixOrdered:(int)i;
-(void)setMatrixRows:(int)i: (int)n;
-(void)setMatrixCols:(int)i: (int)n;
-(void)setMatrixDiag:(int)i: (int)n;
-(void)setMatrixRHS:(int)i: (double)n;
-(void)setMatrixBulkRHS:(int)i: (double)n;
-(void)setMatrixForce:(int)i: (int)j: (double)n;
-(void)setMatrixValues:(int)i: (double)n;
-(void)setMatrixILUValues:(int)i: (double)n;
-(void)setMatrixMassValues:(int)i: (double)n;
-(void)setMatrixDampValues:(int)i: (double)n;
-(void)setMatrixBulkValues:(int)i: (double)n;
-(void)setMatrixILURows:(int)i: (int)n;
-(void)setMatrixILUCols:(int)i: (int)n;
-(void)setMatrixILUDiag:(int)i: (int)n;
-(void)setMatrixCILUValues:(int)i: (double complex)n;

#pragma mark Setters and getters for Variable structure

-(char *)variableName;
-(int)variableDofs;
-(int)variablePerm:(int)i;
-(double)variableNorm;
-(double)variablePrevNorm;
-(double)variableNonLinChange;
-(double)variableSteadyChange;
-(int)variableNonLinConverged;
-(int)variableSteadyConverged;
-(int)variableNonLinIter;
-(double)variableValues:(int)i;
-(double)variablePrevValues:(int)i: (int)j;
-(double)variableNonLinValues:(int)i;

-(void)setVariableName:(char *)string;
-(void)setVariableDofs:(int)i;
-(void)setVariablePerm:(int)i: (int)n;
-(void)setVariableNorm:(double)n;
-(void)setVariableNonLinChange:(double)n;
-(void)setVariableSteadyChange:(double)n;
-(void)setVariableNonLinConverged:(int)n;
-(void)setVariableSteadyConverged:(int)n;
-(void)setVariableNonLinIter:(int)n;
-(void)setVariableValues:(int)i: (double)n;
-(void)setVariablePrevValues:(int)i: (int)j: (double)n;
-(void)setVariableNonLinValues:(int)i: (double)n;


#pragma mark Sizes

-(int)sizeOfBoundaryReorder;
-(int)matrixSizeOfRows;
-(int)matrixSizeOfCols;
-(int)matrixSizeOfDiag;
-(int)matrixSizeOfRHS;
-(int)matrixSizeOfBulkRHS;
-(int)matrixSize1OfForce;
-(int)matrixSize2OfForce;
-(int)matrixSizeOfValues;
-(int)matrixSizeOfILUValues;
-(int)matrixSizeOfMassValues;
-(int)matrixSizeOfDampValues;
-(int)matrixSizeOfBulkValues;
-(int)matrixSizeOfILURows;
-(int)matrixSizeOFILUCols;
-(int)matrixSizeOfILUDiag;
-(int)matrixSizeOfCILUValues;
-(int)variableSizeOfValues;
-(int)variableSizeOfNonLinValues;

-(void)setSizeOfBoundaryReorder:(int)i;
-(void)setMatrixSizeOfRows:(int)i;
-(void)setMatrixSizeOfCols:(int)i;
-(void)setMatrixSizeOfDiag:(int)i;
-(void)setMatrixSizeOfRHS:(int)i;
-(void)setMatrixSizeOfBulkRHS:(int)i;
-(void)setMatrixSize1OfForce:(int)i;
-(void)setMatrixSize2OfForce:(int)i;
-(void)setMatrixSizeOfValues:(int)i;
-(void)setMatrixSizeOfILUValues:(int)i;
-(void)setMatrixSizeOfMassValues:(int)i;
-(void)setMatrixSizeOfDampValues:(int)i;
-(void)setMatrixSizeOfBulkValues:(int)i;
-(void)setMatrixSizeOfILURows:(int)i;
-(void)setMatrixSizeOfILUCols:(int)i;
-(void)setMatrixSizeOfILUDiag:(int)i;
-(void)setMatrixSizeOfCILUValues:(int)i;
-(void)setVariableSizeOfValues:(int)i;
-(void)setVariableSizeOfNonLinValues: (int)i;

#pragma Test associativity

-(BOOL)isAssociatedMatrixILUValues;
-(BOOL)isAssociatedMatrixMassValues;
-(BOOL)isAssociatedMatrixDampValues;
-(BOOL)isAssociatedMatrixILURows;
-(BOOL)isAssociatedMatrixILUCols;
-(BOOL)isAssociatedMatrixILUDiag;
-(BOOL)isAssociatedMatrixCILUValues;
-(BOOL)isAssociatedVariablePrevValues;
-(BOOL)isAssociatedVariableNonLinValues;
-(BOOL)isAssociatedVariableSteadyValues;


#pragma mark Methods returning pointers

-(int *)returnPointerToDefDofs;
-(FEMMesh *)returnPointerToMesh;

-(double *)matrixReturnPointerToRHS;
-(double *)matrixReturnPointerToBulkRHS;
-(double *)matrixReturnPointerToBulkValues;

-(double *)variableReturnPointerToValues;
-(double *)variableReturnPointerToNonLinValues;
-(double *)variableReturnPointerToSteadyValues;

-(Variable_t *)meshReturnPointerToVariables;

#pragma mark Initializations

-(void)initializeILU:(int)ilun;
-(void)initializeCILU:(int)ilun;
-(void)ilutWorkspaceCheck:(int)i: (int)n;
-(void)ilutComplexWorkspaceCheck:(int)i: (int)n;

#pragma mark Allocations

-(void)allocateMatrixBulkRHS:(int)n;
-(void)allocateMatrixBulkValues:(int)n;
-(void)allocateMatrixILUValues:(int)n;
-(void)allocateMatrixILURows:(int)n;
-(void)allocateMatrixILUCols:(int)n;
-(void)allocateMatrixILUDiag:(int)n;
-(void)allocateMatrixCILUValues:(int)n;
-(void)allocateVariableNonLinValues:(int)n;

-(void)freeMatrixBulkRHS:(int)n;
-(void)freeMatrixBulkValues:(int)n;
-(void)freeMatrixILUValues:(int)n;
-(void)freeMatrixILURows:(int)n;
-(void)freeMatrixILUCols:(int)n;
-(void)freeMatrixILUDiag:(int)n;
-(void)freeMatrixCILUValues:(int)n;
-(void)freeVariableNonLinValues:(int)n;



@end
