//
//  FEMSolution.m
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMSolution.h"

@interface FEMSolution ()

-(int)initializeILU1:(Matrix_t *)mat: (int)n;

@end

@implementation FEMSolution

#pragma mark Private methods...

-(int)initializeILU1:(Matrix_t *)mat: (int)n {
    
    int i, j, k, l, nonZeros, rowMin, rowMax;
    int *C;
    
    mat->ILURows = intvec(0, (n+1)-1);
    mat->ILUDiag = intvec(0, n-1);
    
    if (mat->ILURows == NULL || mat->ILUDiag == NULL) {
        errorfunct("initializeILU1", "Memory allocation error.");
    }
    
    // Count fills row by row   
    C = intvec(0, n-1);
    for (i=0; i<n; i++) {
        C[i] = 0;
    }
    nonZeros = mat->Rows[(n+1)-1]-1;
    
    for (i=0; i<n; i++) {
        for (k=mat->Rows[i]; k<=mat->Rows[i+1]-1; k++) {
            C[mat->Cols[k]] = 1;
        }
        
        for (k=mat->Cols[mat->Rows[i]]; k<=i-1; k++) {
            if (C[k] != 0) {
                for (l=mat->Diag[k]+1; l<=mat->Rows[k+1]; l++) {
                    j = mat->Cols[l];
                    if (C[j] == 0) nonZeros = nonZeros + 1;
                }
            }
        }
        
        for (k=mat->Rows[i]; k<=mat->Rows[i+1]-1; k++) {
            C[mat->Cols[k]] = 0;
        }
    }
    
    mat->ILUCols = intvec(0, nonZeros-1);
    if (mat->ILUCols == NULL) {
        errorfunct("initializeILU1", "Memory allocation error.");
    }
    
    // Update row nonzero structures
    for (i=0; i<n; i++) {
        C[i] = 0;
    }
    mat->ILURows[0] = 0;
    for (i=0; i<n; i++) {
        for (k=mat->Rows[i]; k<=mat->Rows[i+1]-1; k++) {
            C[mat->Cols[k]] = 1;
        }
        
        rowMin = mat->Cols[ mat->Rows[i] ];
        rowMax = mat->Cols[ mat->Rows[i+1]-1 ];
        
        for (k=rowMin ; k<=i-1; k++) {
            if (C[k] == 1) {
                for (l=mat->Diag[k]+1; l<=mat->Rows[k+1]-1; l++) {
                    j = mat->Cols[l];
                    if (C[j] == 0) {
                        C[j] = 2;
                        rowMax = max(rowMax, j);
                    }
                }
            }
        }
        
        j = mat->ILURows[i] - 1;
        for (k=rowMin; k<=rowMax; k++) {
            if (C[k] > 0) {
                j = j + 1;
                C[k] = 0;
                mat->ILUCols[j] = k;
                if (k == i) mat->ILUDiag[i] = j;
            }
        }
        
        mat->ILURows[i+1] = j + 1;
    }
    free_ivector(C, 0, n-1);
    
    return nonZeros;
    
}

#pragma mark Public methods...

#pragma mark Initializations

-(void)initializeILU:(int)ilun {
    
    int i;
    int n;
    
    Matrix_t A1;
    
    n = matrix.NumberOfRows;
    
    if (ilun == 0) {
        matrix.ILURows = matrix.Rows;
        matrix.ILUCols = matrix.Cols;
        matrix.ILUDiag = matrix.Diag;
        
        matrix.sizeILURows = matrix.sizeRows;
        matrix.sizeILUCols = matrix.sizeCols;
        matrix.sizeILUDiag = matrix.sizeDiag;
    } else {

        matrix.sizeILUCols = [self initializeILU1:&matrix :n];
        matrix.sizeILURows = n+1;
        matrix.sizeILUDiag = n;
        
        if (ilun > 1) {
            
            for (i=0; i<ilun-1; i++) {
                
                A1.Cols = matrix.ILUCols;
                A1.Rows = matrix.ILURows;
                A1.Diag = matrix.ILUDiag;
                
                [self initializeILU1:&A1 :n];
                
                matrix.ILUCols = A1.ILUCols;
                matrix.ILURows = A1.ILURows;
                matrix.ILUDiag = A1.ILUDiag;
                
                A1.Cols = NULL;
                A1.Rows = NULL;
                A1.Diag = NULL;
            }
        }
    }
    matrix.ILUValues = doublevec(0, matrix.ILURows[(n+1)-1]-1);
    if (matrix.ILUValues == NULL) {
        errorfunct("initializeILU", "Memory allocation error.");
    }
    matrix.sizeILUValues = matrix.ILURows[(n+1)-1]-1;
}

-(void)initializeCILU:(int)ilun {
    
    int i, j, k;
    int n;
    
    Matrix_t A1;
    
    n = matrix.NumberOfRows;
    
    A1.NumberOfRows = n/2;
    A1.Rows = intvec(0, (n/2+1)-1);
    A1.Diag = intvec(0, (n/2)-1);
    A1.Cols = intvec(0, (matrix.sizeCols / 4)-1);
    
    A1.Rows[0] = 0 ;
    k = 0;
    
    for (i=0; i<n; i+=2) {
        for (j=matrix.Rows[i]; j<=matrix.Rows[i+1]-1; j+=2) {
            k = k + 1;
            A1.Cols[k] = (matrix.Cols[j]+1) / 2;
            if (matrix.Cols[j] == i) A1.Diag[(i+1)/2] = k;
        }
        A1.Rows[(i+1)/2+1] = k + 1;
    }
    
    if (ilun == 0) {
        matrix.ILURows = matrix.Rows;
        matrix.ILUCols = matrix.Cols;
        matrix.ILUDiag = matrix.Diag;
        
        matrix.sizeILURows = matrix.sizeRows;
        matrix.sizeILUCols = matrix.sizeCols;
        matrix.sizeILUDiag = matrix.sizeDiag;
    } else {
        
        matrix.sizeILUCols = [self initializeILU1:&A1 :n/2];
        matrix.ILUCols = A1.ILUCols;
        matrix.ILURows = A1.ILURows;
        matrix.ILUDiag = A1.ILUDiag;
        
        matrix.sizeILURows = (n/2+1);
        matrix.sizeILUDiag = n/2;
        
        free_ivector(A1.Rows, 0, (n/2+1)-1);
        free_ivector(A1.Diag, 0, (n/2)-1);
        free_ivector(A1.Cols, 0, (matrix.sizeCols / 4)-1);
        A1.Rows = NULL;
        A1.Diag = NULL;
        A1.Cols = NULL;
        
        if (ilun > 1) {
            
            for (i=0; i<ilun-1; i++) {
                
                A1.Cols = matrix.ILUCols;
                A1.Rows = matrix.ILURows;
                A1.Diag = matrix.ILUDiag;
                
                [self initializeILU1:&A1 :n/2];
                
                matrix.ILUCols = A1.ILUCols;
                matrix.ILURows = A1.ILURows;
                matrix.ILUDiag = A1.ILUDiag;
                
                A1.Cols = NULL;
                A1.Rows = NULL;
                A1.Diag = NULL;
            }
        }
    }
    matrix.CILUValues = cdoublevec(0, matrix.ILURows[(n/2+1)-1]);
    if (matrix.CILUValues == NULL) {
        errorfunct("initializeCILU", "Memory allocation error.");
    }
    matrix.sizeCILUValues = matrix.ILURows[(n/2+1)-1];
}

-(void)ilutWorkspaceCheck:(int)i: (int)n {
    
    int j, k;
    int *iWork;
    double *cWork;
    
    k = matrix.ILURows[i+1] + min(0.75*matrix.ILURows[i+1], ((n-1)-i)*(1.0*n));
    
    iWork = intvec(0, k-1);
    if (iWork == NULL) {
        errorfunct("ilutWorkspaceCheck", "Memory allocation error.");
    }
    for (j=0; j<=matrix.ILURows[i+1]-1; j++) {
        iWork[j] = matrix.ILUCols[j];
    }
    [self freeMatrixILUCols:matrix.sizeILUCols];
    
    cWork = doublevec(0, k-1);
    if (cWork == NULL) {
        errorfunct("ilutWorkspaceCheck", "Memory allocation error.");
    }
    for (j=0; j<=matrix.ILURows[i+1]-1; j++) {
        cWork[j] = matrix.ILUValues[j];
    }
    [self freeMatrixILUValues:matrix.sizeILUValues];
    
    matrix.sizeILUCols = k;
    matrix.sizeILUValues = k;
    
    matrix.ILUCols = iWork;
    matrix.ILUValues = cWork;
    iWork = NULL;
    cWork = NULL;
    
}

-(void)ilutComplexWorkspaceCheck:(int)i: (int)n {
    
    int j, k;
    int *iWork;
    double complex *cWork;
    
    k = matrix.ILURows[i+1] + min(0.75*matrix.ILURows[i+1], ((n-1)-i)*(1.0*n));
    
    iWork = intvec(0, k-1);
    if (iWork == NULL) {
        errorfunct("ilutComplexWorkspaceCheck", "Memory allocation error.");
    }
    for (j=0; j<=matrix.ILURows[i+1]-1; j++) {
        iWork[j] = matrix.ILUCols[j];
    }
    [self freeMatrixILUCols:matrix.sizeILUCols];
    
    cWork = cdoublevec(0, k-1);
    if (cWork == NULL) {
        errorfunct("ilutComplexWorkspaceCheck", "Memory allocation error.");
    }
    for (j=0; j<=matrix.ILURows[i+1]-1; j++) {
        cWork[j] = matrix.CILUValues[j];
    }
    [self freeMatrixILUValues:matrix.sizeILUValues];
    
    matrix.sizeILUCols = k;
    matrix.sizeCILUValues = k;
    
    matrix.ILUCols = iWork;
    matrix.CILUValues = cWork;
    iWork = NULL;
    cWork = NULL;
    
}

#pragma mark Setters and getters for FEMSolution variables

-(id)solutionInfoForKey:(NSString *)key {
    
    return [solutionInfo objectForKey:key];
}

-(NSString *)normalTangentialName {
    
    return normalTangentialName;
}

-(int)simulationID {
    
    return simulationID;
}

-(int)normalTangentialNOFNodes {
    
    return normalTangentialNOFNodes;
}

-(int)coordinateSystemDimension {
    
    return coordinateSystemDimension;
}

-(int)ntElement:(int)i :(int)j {
    
    return ntElement[i][j];
}

-(BOOL)ntZeroingDone:(int)i :(int)j {
    
    return ntZeroingDone[i][j];
}

-(int)boundaryReorder:(int)i {
    
    return boundaryReorder[i];
}

-(double)boundaryNormals:(int)i :(int)j {
    
    return boundaryNormals[i][j];
}

-(double)boundaryTangent1:(int)i :(int)j {
    
    return boundaryTangent1[i][j];
}

-(double)boundaryTangent2:(int)i :(int)j {
    
    return boundaryTangent2[i][j];
}

-(int)size1boundaryNormals {
    
    return size1boundaryNormals;
}

-(int)size2boundaryNormals {
    
    return size2boundaryNormals;
}

-(int)timeOrder {
    
    return timeOrder;
}

-(int)doneTime {
    
    return doneTime;
}

-(int)order {
    
    return order;
}

-(int)nofEigenValues {
    
    return nofEigenValues;
}

-(double)alpha {
    
    return alpha;
}

-(double)beta {
    
    return beta;
}

-(double)dt {
    
    return dt;
}

-(int)defDofs:(int)i {
    
    return defDofs[i];
}

-(int)sizeOfDefDofs {
    
    return sizeDefDofs;
}

-(void)setNormalTangentialName:(NSString *)string {
    
    normalTangentialName = [NSString stringWithString:string];
}

-(void)setSimulationID:(int)i {
    
    simulationID = i;
}

-(void)setNormalTangentialNOFNodes:(int)i {
    
    normalTangentialNOFNodes = i;
}

-(void)setCoordinateSystemDimension:(int)i {
    
    coordinateSystemDimension = i;
}

-(void)setNtElement:(int)i :(int)j :(int)n {
    
    ntElement[i][j] = n;
}

-(void)setNtZeroingDone:(int)i :(int)j :(BOOL)n {
    
    ntZeroingDone[i][j] = n;
}

-(void)setBoundaryReorder:(int)i :(int)n {
    
    boundaryReorder[i] = n;
}

-(void)setBoundaryNormals:(int)i :(int)j :(double)n {
    
    boundaryNormals[i][j] = n;
}

-(void)setBoundaryTangent1:(int)i :(int)j :(double)n {
    
    boundaryTangent1[i][j] = n;
}

-(void)setBoundaryTangent2:(int)i :(int)j :(double)n {
    
    boundaryTangent2[i][j] = n;
}

-(void)setSize1boundaryNormals:(int)n {
    
    size1boundaryNormals = n;
}

-(void)setSize2boundaryNormals:(int)n {
    
    size2boundaryNormals = n;
}

-(void)setTimeOrder:(int)n {
    
    timeOrder = n;
}

-(void)setDoneTime:(int)n {
    
    doneTime = n;
}

-(void)setOrder:(int)n {
    
    order = n;
}

-(void)setAlpha:(double)n {
    
    alpha = n;
}

-(void)setBeta:(double)n {
    
    beta = n;
}

-(void)setDt:(double)n {
    
    dt = n;
}

-(void)setDefDofs:(int)i :(int)n {
    
    defDofs[i] = n;
}

-(void)setSizeOfDefDofs:(int)n {
    
    sizeDefDofs = n;
}

#pragma mark Setters and getters for Matrix structure

-(int)matrixNumberOfRows {
    
    return matrix.NumberOfRows;
}

-(int)matrixSubband {
    
    return matrix.Subband;
}

-(int)matrixFormat {
    
    return matrix.Format;
}

-(int)matrixSolveCount {
    
    return matrix.SolveCount;
}

-(int)matrixOrdered {
    
    return matrix.Ordered;
}

-(int)matrixLumped {
    
    return matrix.Lumped;
}

-(int)matrixSymmetric {
    
    return matrix.Symmetric;
}

-(int)matrixComplex {
    
    return matrix.Complex;
}

-(double)matrixRHS:(int)i {
    
    return matrix.RHS[i];
}

-(double)matrixBulkRHS:(int)i {
    
    return matrix.BulkRHS[i];
}

-(double)matrixForce:(int)i :(int)j {
    
    return matrix.Force[i][j];
}

-(int)matrixRows:(int)i {
    
    return matrix.Rows[i];
}

-(int)matrixCols:(int)i {
    
    return matrix.Cols[i];
}

-(int)matrixDiag:(int)i {
    
    return matrix.Diag[i];
}

-(double)matrixValues:(int)i {
    
    return matrix.Values[i];
}

-(double)matrixILUValues:(int)i {
    
    return matrix.ILUValues[i];
}

-(double)matrixMassValues:(int)i {
    
    return matrix.MassValues[i];
}

-(double)matrixDampValues:(int)i {
    
    return matrix.DampValues[i];
}

-(double)matrixBulkValues:(int)i {
    
    return matrix.BulkValues[i];
}

-(int)matrixILURows:(int)i {
    
    return matrix.ILURows[i];
}

-(int)matrixILUCols:(int)i {
    
    return matrix.ILUCols[i];
}

-(int)matrixILUDiag:(int)i {
    
    return matrix.Diag[i];
}

-(double complex)matrixCILUValues:(int)i {
    
    return matrix.CILUValues[i];
}

-(void)setMatrixNumberOfRows:(int)i {
    
    matrix.NumberOfRows = i;
}

-(void)setMatrixSubband:(int)i {
    
    matrix.Subband = i;
}

-(void)setMatrixFormat:(int)i {
    
    matrix.Format = i;
}

-(void)setMatrixSolveCount:(int)i {
    
    matrix.SolveCount = i;
}

-(void)setMatrixOrdered:(int)i {
    
    matrix.Ordered = i;
}

-(void)setMatrixLumped:(int)i {
    
    matrix.Lumped = i;
}

-(void)setMatrixSymmetric:(int)i {
    
    matrix.Symmetric = i;
}

-(void)setMatrixRows:(int)i :(int)n {
    
    matrix.Rows[i] = n;
}

-(void)setMatrixCols:(int)i :(int)n {
    
    matrix.Cols[i] = n;
}

-(void)setMatrixDiag:(int)i :(int)n {
    
    matrix.Diag[i] = n;
}

-(void)setMatrixRHS:(int)i :(double)n {
    
    matrix.RHS[i] = n; 
}

-(void)setMatrixBulkRHS:(int)i :(double)n {
    
    matrix.BulkRHS[i] = n;
}

-(void)setMatrixForce:(int)i :(int)j :(double)n {
    
    matrix.Force[i][j] = n;
}

-(void)setMatrixValues:(int)i :(double)n {
    
    matrix.Values[i] = n;
}

-(void)setMatrixILUValues:(int)i :(double)n {
    
    matrix.ILUValues[i] = n;
}

-(void)setMatrixMassValues:(int)i :(double)n {
    
    matrix.MassValues[i] = n;
}

-(void)setMatrixDampValues:(int)i :(double)n {
    
    matrix.DampValues[i] = n;
}

-(void)setMatrixBulkValues:(int)i :(double)n {
    
    matrix.BulkValues[i] = n;
}

-(void)setMatrixILURows:(int)i :(int)n {
    
    matrix.ILURows[i] = n;
}

-(void)setMatrixILUCols:(int)i :(int)n {
    
    matrix.ILUCols[i] = n;
}

-(void)setMatrixILUDiag:(int)i :(int)n {
    
    matrix.ILUDiag[i] = n;
}

-(void)setMatrixCILUValues:(int)i :(_Complex double)n {
    
    matrix.CILUValues[i] = n;
}

#pragma mark Setters and getters for Variable structure

-(char *)variableName {
    
    return variable.Name;
}

-(int)variableDofs {
    
    return variable.Dofs;
}

-(int)variablePerm:(int)i {
    
    return variable.Perm[i];
}

-(double)variableNorm {
    
    return variable.Norm;
}

-(double)variablePrevNorm {
    
    return variable.PrevNorm;
}

-(double)variableNonLinChange {
    
    return variable.NonLinChange;
}

-(double)variableSteadyChange {
    
    return variable.SteadyChange;
}

-(int)variableNonLinConverged {
    
    return variable.NonLinConverged;
}

-(int)variableSteadyConverged {
    
    return variable.SteadyConverged;
}

-(int)variableNonLinIter {
    
    return variable.NonLinIter;
}

-(double)variableValues:(int)i {
    
    return variable.Values[i];
}

-(double)variablePrevValues:(int)i :(int)j {
    
    return variable.PrevValues[i][j];
    
}

-(double)variableNonLinValues:(int)i {
    
    return variable.NonLinValues[i];
}

-(void)setVariableName:(char *)string {
    
    variable.Name = string;
}

-(void)setVariableDofs:(int)i {
    
    variable.Dofs = i;
}

-(void)setVariablePerm:(int)i :(int)n {
    
    variable.Perm[i] = n;
}

-(void)setVariableNorm:(double)n {
    
    variable.Norm = n;
}

-(void)setVariableNonLinChange:(double)n {
    
    variable.NonLinChange = n;
}

-(void)setVariableSteadyChange:(double)n {
    
    variable.SteadyChange = n;
}

-(void)setVariableNonLinConverged:(int)n {
    
    variable.NonLinConverged = n;
}

-(void)setVariableSteadyConverged:(int)n {
    
    variable.SteadyConverged = n;
}

-(void)setVariableNonLinIter:(int)n {
    
    variable.NonLinIter = n;
}

-(void)setVariableValues:(int)i :(double)n {
    
    variable.Values[i] = n;
}

-(void)setVariablePrevValues:(int)i :(int)j :(double)n {
    
    variable.PrevValues[i][j] = n;
}

-(void)setVariableNonLinValues:(int)i :(double)n {
    
    variable.NonLinValues[i] = n;
}

#pragma mark Sizes

-(int)sizeOfBoundaryReorder {
    
    return sizeBoundaryReorder;
}

-(int)matrixSizeOfRows {
    
    return matrix.sizeRows;
}

-(int)matrixSizeOfCols {
    
    return matrix.sizeCols;
}

-(int)matrixSizeOfDiag {
    
    return matrix.sizeDiag;
}

-(int)matrixSizeOfRHS {
    
    return matrix.sizeRHS;
}

-(int)matrixSizeOfBulkRHS {
    
    return matrix.sizeBulkRHS;
}

-(int)matrixSize1OfForce {
    
    return matrix.size1force;
}

-(int)matrixSize2OfForce {
    
    return matrix.size2Force;
}

-(int)matrixSizeOfValues {
    
    return matrix.sizeValues;
}

-(int)matrixSizeOfILUValues {
    
    return matrix.sizeILUValues;
}

-(int)matrixSizeOfMassValues {
    
    return matrix.sizeMassValues;
}

-(int)matrixSizeOfDampValues {
    
    return matrix.sizeDampValues;
}

-(int)matrixSizeOfBulkValues {
    
    return matrix.sizeBulkValues;
}

-(int)matrixSizeOfILURows {
    
    return matrix.sizeILURows;
}

-(int)matrixSizeOFILUCols {
    
    return matrix.sizeILUCols;
}

-(int)matrixSizeOfILUDiag {
    
    return matrix.sizeILUDiag;
}

-(int)matrixSizeOfCILUValues {
    
    return matrix.sizeCILUValues;
}

-(int)variableSizeOfPerm {
    
    return variable.sizePerm;
}

-(int)variableSizeOfValues {
    
    return variable.sizeValues;
}

-(int)variableSizeOfNonLinValues {
    
    return variable.sizeNonLinValues;
}

-(void)setSizeOfBoundaryReorder:(int)i {
    
    sizeBoundaryReorder = i;
}

-(void)setMatrixSizeOfRows:(int)i {
    
    matrix.sizeRows = i;
}

-(void)setMatrixSizeOfCols:(int)i {
    
    matrix.sizeCols = i;
}

-(void)setMatrixSizeOfDiag:(int)i {
    
    matrix.sizeDiag = i;
}

-(void)setMatrixSizeOfRHS:(int)i {
    
    matrix.sizeRHS = i;
}

-(void)setMatrixSizeOfBulkRHS:(int)i {
    
    matrix.sizeBulkRHS = i;
}

-(void)setMatrixSize1OfForce:(int)i {
    
    matrix.size1force = i;
}

-(void)setMatrixSize2OfForce:(int)i {
    
    matrix.size2Force = i;
}

-(void)setMatrixSizeOfValues:(int)i {
    
    matrix.sizeValues = i;
}

-(void)setMatrixSizeOfILUValues:(int)i {
    
    matrix.sizeILUValues = i;
}

-(void)setMatrixSizeOfMassValues:(int)i {
    
    matrix.sizeMassValues = i;
}

-(void)setMatrixSizeOfDampValues:(int)i {
    
    matrix.sizeDampValues = i;
}

-(void)setMatrixSizeOfBulkValues:(int)i {
    
    matrix.sizeBulkValues = i;
}

-(void)setMatrixSizeOfILURows:(int)i {
    
    matrix.sizeILURows = i;
}

-(void)setMatrixSizeOfILUCols:(int)i {
    
    matrix.sizeILUCols = i;
}

-(void)setMatrixSizeOfILUDiag:(int)i {
    
    matrix.sizeILUDiag = i;
}

-(void)setMatrixSizeOfCILUValues:(int)i {
    
    matrix.sizeCILUValues = i;
}

-(void)setVariableSizeOfPerm:(int)i {
    
    variable.sizePerm = i;
}

-(void)setVariableSizeOfValues:(int)i {
    
    variable.sizeValues = i;
}

-(void)setVariableSizeOfNonLinValues:(int)i {
    
    variable.sizeNonLinValues = i;
}

#pragma mark Test associativity

-(BOOL)isAssociatedMatrixDiag {

    if (matrix.Diag != NULL) {
        return YES;
    } else {
        return NO;
    }
    
}

-(BOOL)isAssociatedMatrixILUValues {
    
    if (matrix.ILUValues != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedMatrixMassValues {
    
    if (matrix.MassValues != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedMatrixDampValues {
    
    if (matrix.DampValues != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedMatrixILURows {
    
    if (matrix.ILURows != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedMatrixILUCols {
    
    if (matrix.ILUCols != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedMatrixILUDiag {
    
    if (matrix.ILUDiag != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedMatrixCILUValues {
    
    if (matrix.CILUValues != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedVariablePrevValues {
    
    if (variable.PrevValues != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedVariableNonLinValues {
    
    if (variable.NonLinValues != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedVariableSteadyValues {
    
    if (variable.SteadyValues != NULL) {
        return YES;
    } else {
        return NO;
    }
}


#pragma mark Methods returning pointers

-(int *)returnPointerToDefDofs {
    
    return defDofs;
}

-(FEMMesh *)returnPointerToMesh {
    
    return mesh;
}

-(Matrix_t *)matrixReturnPointerToConstraintMatrix {

    return matrix.ConstraintMatrix;
}

-(double *)matrixReturnPointerToRHS {
    
    return matrix.RHS;
}

-(double *)matrixReturnPointerToBulkRHS {
    
    return matrix.BulkRHS;
}

-(double *)matrixReturnPointerToBulkValues {
    
    return matrix.BulkValues;
}

-(double *)variableReturnPointerToValues {
    
    return variable.Values;
}

-(double *)variableReturnPointerToNonLinValues {
    
    return variable.NonLinValues;
}

-(double *)variableReturnPointerToSteadyValues {
    
    return variable.SteadyValues;
}

-(Variable_t *)meshReturnPointerToVariables {
    
    return [mesh returnPointerToVariables];
}

#pragma mark Methods assigning pointers

-(void)matrixAssignConstraintMatrix:(Matrix_t *)a {
    
    matrix.ConstraintMatrix = a;
}


#pragma mark Allocations and disallocations

-(void)allocateMatrixBulkRHS:(int)n {
    
    matrix.BulkRHS = doublevec(0, n-1);
}

-(void)allocateMatrixBulkValues:(int)n {
    
    matrix.BulkValues = doublevec(0, n-1);
}

-(void)allocateMatrixILUValues:(int)n {
    
    matrix.ILUValues = doublevec(0, n-1);
}

-(void)allocateMatrixILURows:(int)n {
    
    matrix.ILURows = intvec(0, n-1);
}

-(void)allocateMatrixILUCols:(int)n {
    
    matrix.ILUCols = intvec(0, n-1);
}

-(void)allocateMatrixILUDiag:(int)n {
    
    matrix.ILUDiag = intvec(0, n-1);
}

-(void)allocateMatrixCILUValues:(int)n {
    
    matrix.CILUValues = cdoublevec(0, n-1);
}

-(void)allocateVariableNonLinValues:(int)n {
    
    variable.NonLinValues = doublevec(0, n-1);
}

-(void)freeMatrixBulkRHS:(int)n {
    
    free_dvector(matrix.BulkRHS, 0, n-1);
    matrix.BulkRHS = NULL;
}

-(void)freeMatrixBulkValues:(int)n {
    
    free_dvector(matrix.BulkValues, 0, n-1);
    matrix.BulkValues = NULL;
}

-(void)freeMatrixILUValues:(int)n {
    
    free_dvector(matrix.ILUValues, 0, n-1);
    matrix.ILUValues = NULL;
}

-(void)freeMatrixILURows:(int)n {
    
    free_ivector(matrix.ILURows, 0, n-1);
    matrix.ILURows = NULL;
}

-(void)freeMatrixILUCols:(int)n {
    
    free_ivector(matrix.ILUCols, 0, n-1);
    matrix.ILUCols = NULL;
}

-(void)freeMatrixILUDiag:(int)n {
    
    free_ivector(matrix.ILUDiag, 0, n-1);
    matrix.ILUDiag = NULL;
}

-(void)freeMatrixCILUValues:(int)n {
    
    free_cdvector(matrix.CILUValues, 0, n-1);
    matrix.CILUValues = NULL;
}

-(void)freeVariableNonLinValues:(int)n {

    free_dvector(variable.NonLinValues, 0, n-1);
    variable.NonLinValues = NULL;
}


@end
