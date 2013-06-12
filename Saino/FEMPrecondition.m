//
//  FEMPreconditioners.m
//  Saino
//
//  Created by Hakime Seddik on 27/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <float.h>
#import "FEMPrecondition.h"
#import "Utils.h"

@interface FEMPrecondition ()

-(int)FEMPrecondition_initializeILU1:(matrixArraysContainer *)containers :(int)n;
-(void)FEMPrecondition_ilutWorkspaceCheckMatrix:(FEMMatrix *)matrix atIndex:(int)i numberOfRows:(int)n;
-(void)FEMPrecondition_ilutComplexWorkspaceCheckMatrix:(FEMMatrix *)matrix atIndex:(int)i numberOfRows:(int)n;
-(void)FEMPrecondition_computeIlutMatrix:(FEMMatrix *)matrix numberOfRows:(int)n tolerance:(int)tol;
-(void)FEMPrecondition_computeComplexIlutMatrix:(FEMMatrix *)matrix numberOfRows:(int)n tolerance:(int)tol;
-(void)FEMPrecondition_LUSolveSystemSize:(int)n matrix:(FEMMatrix *)matrix rightHandSide:(double *)b;
-(void)FEMPrecondition_ComplexLUSolveSystemSize:(int)n matrix:(FEMMatrix *)matrix rightHandSide:(double complex *)b;

@end

@implementation FEMPrecondition

#pragma mark Private methods

-(int)FEMPrecondition_initializeILU1:(matrixArraysContainer *)containers :(int)n {
    
    int i, j, k, l, nonZeros, rowMin, rowMax;
    int *C;
    
    containers->ILURows = intvec(0, (n+1)-1);
    containers->ILUDiag = intvec(0, n-1);
    
    if (containers->ILURows == NULL || containers->ILUDiag == NULL) {
        errorfunct("FEMPrecondition_initializeILU1", "Memory allocation error.");
    }
    
    // Count fills row by row
    C = intvec(0, n-1);
    for (i=0; i<n; i++) {
        C[i] = 0;
    }
    nonZeros = containers->Rows[(n+1)-1]-1;
    
    for (i=0; i<n; i++) {
        for (k=containers->Rows[i]; k<=containers->Rows[i+1]-1; k++) {
            C[containers->Cols[k]] = 1;
        }
        
        for (k=containers->Cols[containers->Rows[i]]; k<=i-1; k++) {
            if (C[k] != 0) {
                for (l=containers->Diag[k]+1; l<=containers->Rows[k+1]; l++) {
                    j = containers->Cols[l];
                    if (C[j] == 0) nonZeros = nonZeros + 1;
                }
            }
        }
        
        for (k=containers->Rows[i]; k<=containers->Rows[i+1]-1; k++) {
            C[containers->Cols[k]] = 0;
        }
    }
    
    containers->ILUCols = intvec(0, nonZeros-1);
    if (containers->ILUCols == NULL) {
        errorfunct("FEMPrecondition_initializeILU1", "Memory allocation error.");
    }
    
    // Update row nonzero structures
    for (i=0; i<n; i++) {
        C[i] = 0;
    }
    containers->ILURows[0] = 0;
    for (i=0; i<n; i++) {
        for (k=containers->Rows[i]; k<=containers->Rows[i+1]-1; k++) {
            C[containers->Cols[k]] = 1;
        }
        
        rowMin = containers->Cols[ containers->Rows[i] ];
        rowMax = containers->Cols[ containers->Rows[i+1]-1 ];
        
        for (k=rowMin ; k<=i-1; k++) {
            if (C[k] == 1) {
                for (l=containers->Diag[k]+1; l<=containers->Rows[k+1]-1; l++) {
                    j = containers->Cols[l];
                    if (C[j] == 0) {
                        C[j] = 2;
                        rowMax = max(rowMax, j);
                    }
                }
            }
        }
        
        j = containers->ILURows[i] - 1;
        for (k=rowMin; k<=rowMax; k++) {
            if (C[k] > 0) {
                j = j + 1;
                C[k] = 0;
                containers->ILUCols[j] = k;
                if (k == i) containers->ILUDiag[i] = j;
            }
        }
        
        containers->ILURows[i+1] = j + 1;
    }
    free_ivector(C, 0, n-1);
    
    return nonZeros;
}

-(void)FEMPrecondition_ilutWorkspaceCheckMatrix:(FEMMatrix *)matrix atIndex:(int)i numberOfRows:(int)n {
    
    int j, k;
    int *iWork = NULL;
    double *cWork = NULL;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    k = matContainers->ILURows[i+1] + min(0.75*matContainers->ILURows[i+1], ((n-1)-i)*(1.0*n));
    
    iWork = intvec(0, k-1);
    if (iWork == NULL) {
        errorfunct("FEMPrecondition_ilutWorkspaceCheckMatrix", "Memory allocation error.");
    }
    for (j=0; j<=matContainers->ILURows[i+1]-1; j++) {
        iWork[j] = matContainers->ILUCols[j];
    }
    free_ivector(matContainers->ILUCols, 0, matContainers->sizeILUCols-1);
    matContainers->ILUCols = NULL;
    
    
    cWork = doublevec(0, k-1);
    if (cWork == NULL) {
        errorfunct("FEMPrecondition_ilutWorkspaceCheckMatrix", "Memory allocation error.");
    }
    for (j=0; j<=matContainers->ILURows[i+1]-1; j++) {
        cWork[j] = matContainers->ILUValues[j];
    }
    free_dvector(matContainers->ILUValues, 0, matContainers->sizeILUValues-1);
    matContainers->ILUValues = NULL;
    
    matContainers->sizeILUCols = k;
    matContainers->sizeILUValues = k;
    
    matContainers->ILUCols = iWork;
    matContainers->ILUValues = cWork;
}

-(void)FEMPrecondition_ilutComplexWorkspaceCheckMatrix:(FEMMatrix *)matrix atIndex:(int)i numberOfRows:(int)n {
    
    int j, k;
    int *iWork = NULL;
    double complex *cWork = NULL;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    k = matContainers->ILURows[i+1] + min(0.75*matContainers->ILURows[i+1], ((n-1)-i)*(1.0*n));
    
    iWork = intvec(0, k-1);
    if (iWork == NULL) {
        errorfunct("FEMPrecondition_ilutComplexWorkspaceCheckMatrix", "Memory allocation error.");
    }
    for (j=0; j<=matContainers->ILURows[i+1]-1; j++) {
        iWork[j] = matContainers->ILUCols[j];
    }
    free_ivector(matContainers->ILUCols, 0, matContainers->sizeILUCols-1);
    matContainers->ILUCols = NULL;
    
    cWork = cdoublevec(0, k-1);
    if (cWork == NULL) {
        errorfunct("FEMPrecondition_ilutComplexWorkspaceCheckMatrix", "Memory allocation error.");
    }
    for (j=0; j<=matContainers->ILURows[i+1]-1; j++) {
        cWork[j] = matContainers->CILUValues[j];
    }
    free_cdvector(matContainers->CILUValues, 0, matContainers->sizeCILUValues-1);
    matContainers->CILUValues = NULL;
    
    matContainers->sizeILUCols = k;
    matContainers->sizeCILUValues = k;
    
    matContainers->ILUCols = iWork;
    matContainers->CILUValues = cWork;
}

-(void)FEMPrecondition_computeIlutMatrix:(FEMMatrix *)matrix numberOfRows:(int)n tolerance:(int)tol {
    
    int i, j, k, l, rowMin, rowMax;
    bool *C;
    const double WORKN = 128;
    double norma, cptime, ttime, t;
    double *S;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    ttime = cputime();
    cptime = 0.0;
    
    matContainers->ILURows = intvec(0, (n+1)-1);
    matContainers->ILUDiag = intvec(0, n-1);
    if (matContainers->ILURows == NULL || matContainers->ILUDiag == NULL) {
        errorfunct("FEMPrecondition_computeIlutMatrix", "Memory allocation error.");
    }
    matContainers->sizeILURows = n+1;
    matContainers->sizeILUDiag = n;
    
    matContainers->ILUCols = intvec(0, (WORKN*n)-1);
    matContainers->ILUValues = doublevec(0, (WORKN*n)-1);
    if (matContainers->ILUCols == NULL || matContainers->ILUValues == NO) {
        errorfunct("FEMPrecondition_computeIlutMatrix", "Memory allocation error.");
    }
    matContainers->sizeILUCols = WORKN*n;
    matContainers->sizeILUValues = WORKN*n;
    
    // The factorization row by row
    matContainers->ILURows[0] = 0;
    C = boolvec(0, n-1);
    S = doublevec(0, n-1);
    memset( C, false, n*sizeof(bool) );
    memset( S, 0.0, n*sizeof(double) );
    
    for (i=0; i<n; i++) {
        
        // Convert the current row to full from for speed,
        // only flagging the nonzero entries
        for (k=matContainers->Rows[i]; k<=matContainers->Rows[i+1]-1; k++) {
            C[matContainers->Cols[k]] = true;
            S[matContainers->Cols[k]] = matContainers->Values[k];
        }
        
        // Check for bandwidth for speed, bandwidth optimization,
        // it helps a lot.
        rowMin = matContainers->Cols[matContainers->Rows[i]];
        rowMax = matContainers->Cols[matContainers->Rows[i+1]-1];
        
        // Here is the factorization part for the current row;
        for (k=rowMin; k<=i-1; k++) {
            if (C[k] == true) {
                if (fabs(matContainers->ILUValues[matContainers->ILUDiag[k]]) > AEPS) S[k] = S[k] / matContainers->ILUValues[matContainers->ILUDiag[k]];
                
                for (l=matContainers->ILUDiag[k]+1; l<=matContainers->ILURows[k+1]-1; l++) {
                    j = matContainers->ILUCols[l];
                    if (C[j] == false) {
                        C[j] = true;
                        rowMax = max(rowMax, j);
                    }
                    S[j] = S[j] - S[k] * matContainers->ILUValues[l];
                }
            }
        }
        
        // This is the ILUT part, drop element ILU(i,j), if
        // ABS(ILU(i,j)) <= NORM(A(i.:))*tol
        norma = 0.0;
        for (k=matContainers->Rows[i]; k<=matContainers->Rows[i+1]-1; k++) {
            norma = norma + pow(fabs(matContainers->Values[k]), 2.0);
        }
        norma = sqrt(norma);
        
        j =matContainers->ILURows[i]-1;
        for (k=rowMin; k<=rowMax; k++) {
            if (C[k] == true) {
                if (fabs(S[k]) >= tol*norma || k == i) {
                    j = j + 1;
                    matContainers->ILUCols[j] = k;
                    matContainers->ILUValues[j] = S[k];
                    if (k == i) matContainers->ILUDiag[i] = j;;
                }
                S[k] = 0.0;
                C[k] = false;
            }
        }
        matContainers->ILURows[i+1] = j+1;
        
        // Preparations for the next row
        if (i < n-1) {
            
            // Check if still enough workspace
            if (matContainers->sizeILUCols < matContainers->ILURows[i+1]+n) {
                
                t = cputime();
                [self FEMPrecondition_ilutWorkspaceCheckMatrix:matrix atIndex:i numberOfRows:n];
                cptime = cptime + (cputime() - t);
                
            }
        }
    }
    
    // Prescale the diagonal for the LU solver
    for (i=0; i<n; i++) {
        
        if (fabs(matContainers->ILUValues[matContainers->ILUDiag[i]]) < AEPS) {
            matContainers->ILUValues[matContainers->ILUDiag[i]] = 1.0;
        } else {
            matContainers->ILUValues[matContainers->ILUDiag[i]] = 1.0 / matContainers->ILUValues[matContainers->ILUDiag[i]];
        }
    }
    free_bvector(C, 0, n-1);
    free_dvector(S, 0, n-1);
}

-(void)FEMPrecondition_computeComplexIlutMatrix:(FEMMatrix *)matrix numberOfRows:(int)n tolerance:(int)tol {
    
    int i, j, k, l, rowMin, rowMax;
    bool *C;
    const double WORKN = 128;
    double norma;
    double complex *S;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    matContainers->ILURows = intvec(0, (n+1)-1);
    matContainers->ILUDiag = intvec(0, n-1);
    if (matContainers->ILURows == NULL || matContainers->ILUDiag == NULL) {
        errorfunct("FEMPrecondition_computeComplexIlutMatrix", "Memory allocation error.");
    }
    matContainers->sizeILURows = n+1;
    matContainers->sizeILUDiag = n;
    
    matContainers->ILUCols = intvec(0, (WORKN*n)-1);
    matContainers->CILUValues = cdoublevec(0, (WORKN*n)-1);
    if (matContainers->ILUCols == NULL || matContainers->CILUValues == NULL) {
        errorfunct("FEMPrecondition_computeComplexIlutMatrix", "Memory allocation error.");
    }
    matContainers->sizeILUCols = WORKN*n;
    matContainers->sizeCILUValues = WORKN*n;
    
    // The factorization row by row
    matContainers->ILURows[0] = 0;
    C = boolvec(0, n-1);
    S = cdoublevec(0, n-1);
    memset( C, false, n*sizeof(bool) );
    memset( S, 0.0, n*sizeof(double complex) );
    
    for (i=0; i<n; i++) {
        
        // Convert the current row to full from for speed,
        // only flagging the nonzero entries
        for (k=matContainers->Rows[2*i-1]; k<=matContainers->Rows[2*i]-1; k+=2) {
            C[(matContainers->Cols[k]+1) / 2] = true;
            S[(matContainers->Cols[k]+1) / 2] = matContainers->Values[k] + (-matContainers->Values[k+1] * I);
        }
        
        // Check for bandwidth for speed, bandwidth optimization
        // helps a lot.
        rowMin = (matContainers->Cols[matContainers->Rows[2*i-1]] + 1) / 2;
        rowMax = (matContainers->Cols[matContainers->Rows[2*i]-1] + 1) / 2;
        
        // Here is the factorization part for the current row;
        for (k=rowMin; k<=i-1; k++) {
            if (C[k] == true) {
                if (fabs(matContainers->CILUValues[matContainers->ILUDiag[k]]) > AEPS) S[k] = S[k] / matContainers->CILUValues[matContainers->ILUDiag[k]];
                
                for (l=matContainers->ILUDiag[k]+1; l<=matContainers->ILURows[k+1]-1; l++) {
                    j = matContainers->ILUCols[l];
                    if (C[j] == false) {
                        C[j] = true;
                        rowMax = max(rowMax, j);
                    }
                    S[j] = S[j] - S[k] * matContainers->CILUValues[l];
                }
            }
        }
        
        // This is the ILUT part, drop element ILU(i,j), if
        // ABS(ILU(i,j)) <= NORM(A(i.:))*tol
        norma = 0.0;
        for (k=matContainers->Rows[2*i-1]; k<=matContainers->Rows[2*i]-1; k++) {
            norma = norma + pow(matContainers->Values[k], 2.0) + pow(matContainers->Values[k+1], 2.0);
        }
        norma = sqrt(norma);
        
        j = matContainers->ILURows[i]-1;
        for (k=rowMin; k<=rowMax; k++) {
            if (C[k] == true) {
                if (fabs(S[k]) >= tol*norma || k == i) {
                    j = j + 1;
                    matContainers->ILUCols[j] = k;
                    matContainers->CILUValues[j] = S[k];
                    if (k == i) matContainers->ILUDiag[i] = j;
                }
                S[k] = 0.0 + 0.0 * I;
                C[k] = false;
            }
        }
        matContainers->ILURows[i+1] = j+1;
        
        // Preparations for the next row
        if (i < n-1) {
            
            // Check if still enough workspace
            if (matContainers->sizeILUCols < matContainers->ILURows[i+1]+n) {
                [self FEMPrecondition_ilutComplexWorkspaceCheckMatrix:matrix atIndex:i numberOfRows:n];
            }
        }
    }
    
    // Prescale the diagonal for the LU solver
    for (i=0; i<n; i++) {
        
        if (fabs(matContainers->CILUValues[matContainers->ILUDiag[i]]) < AEPS) {
            matContainers->CILUValues[matContainers->ILUDiag[i]] = 1.0;
        } else {
            matContainers->CILUValues[matContainers->ILUDiag[i]] = 1.0 / matContainers->CILUValues[matContainers->ILUDiag[i]];
        }
    }
    free_bvector(C, 0, n-1);
    free_cdvector(S, 0, n-1);
}

/*******************************************************************************************
 
    Description:
    Solve a system  (Ax=b) after factorization A=LUD has been done. This method is
    meant as a part of a preconditioner for an iterative solver. Real version
 
        Arguments:
 
        int n                    -> Size of the system.
 
        FEMMatrix *matrix        -> Input matrix.
 
        double *b                -> On entry the RHS vector, on exit the solution vector.
 
*******************************************************************************************/
-(void)FEMPrecondition_LUSolveSystemSize:(int)n matrix:(FEMMatrix *)matrix rightHandSide:(double *)b {
    
    int i, j;
    double s;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    // If no ILU provided do diagonal solve
    if (matContainers->ILUValues == NULL) {
        for (i=0; i<n; i++) {
            b[i] = b[i] / matContainers->Values[matContainers->Diag[i]];
        }
        matContainers = NULL;
        return;
    }
    
    // Forward substitute (solve z from Lz = b)
    for (i=0; i<n; i++) {
        s = b[i];
        for (j=matContainers->ILURows[i]; j<=matContainers->ILUDiag[i]-1; j++) {
            s = s - matContainers->ILUValues[j] * b[matContainers->ILUCols[j]];
        }
        b[i] = s;
    }
    
    // Backward substitute (solve x from UDx = z)
    for (i=n-1; i>=0; i--) {
        s = b[i];
        for (j=matContainers->ILUDiag[i]+1; j<=matContainers->ILURows[i+1]-1; j++) {
            s = s - matContainers->ILUValues[j] * b[matContainers->ILUCols[j]];
        }
        b[i] = matContainers->ILUValues[matContainers->ILUDiag[i]] * s;
    }
}

/*******************************************************************************************
 
    Description:
    Solve a system  (Ax=b) after factorization A=LUD has been done. This method is
    meant as a part of a preconditioner for an iterative solver. Complex version
 
    Arguments:
 
        int n                  -> Size of the system.
 
        FEMMatrix *matrix      -> Input matrix.
 
        double complex *b      -> On entry the RHS vector, on exit the solution vector.
 
*******************************************************************************************/
-(void)FEMPrecondition_ComplexLUSolveSystemSize:(int)n matrix:(FEMMatrix *)matrix rightHandSide:(double complex *)b {
    
    int i, j;
    double complex s, x;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    // If no ILU provided do diagonal solve
    if (matContainers->CILUValues == NULL) {
        for (i=0; i<n/2; i++) {
            
            x = matContainers->Values[matContainers->Diag[2*i-1]] + (-matContainers->Values[matContainers->Diag[2*i-1]+1] * I);
            b[i] = b[i] / x;
        }
        matContainers = NULL;
        return;
    }
    
    // Forward substitute (solve z from Lz = b)
    for (i=0; i<n; i++) {
        s = b[i];
        for (j=matContainers->ILURows[i]; j<=matContainers->ILUDiag[i]-1; j++) {
            s = s - matContainers->CILUValues[j] * b[matContainers->ILUCols[j]];
        }
        b[i] = s;
    }
    
    // Backward substitute (solve x from UDx = z)
    for (i=n-1; i>=0; i--) {
        s = b[i];
        for (j=matContainers->ILUDiag[i]+1; j<=matContainers->ILURows[i+1]-1; j++) {
            s = s - matContainers->CILUValues[j] * b[matContainers->ILUCols[j]];
        }
        b[i] = matContainers->CILUValues[matContainers->ILUDiag[i]] * s;
    }
}

#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
    }
    
    return self;
}


/*******************************************************************************************
 
    Description: Diagonal preconditioning of a CRS format matrix. Matrix is accessed from
    the FEMSolution class.
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix.
 
        double u               ->  Resulting approximate solution after preconditioning
 
        double v               -> Given right-hand-side
 
        int ipar               -> Input stucture holding info from the HUTIter iterative solver
 
*******************************************************************************************/
-(void)CRSDiagPreconditionMatrix:(FEMMatrix *)matrix afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar {

    int i, j, n;
    int *range1;
    double *range2;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    n = matrix.numberOfRows;  
    
    if (matrix.isOrdered == NO) {
        for (i=0; i<n; i++) {
            range1 = intvec(0, (matContainers->RHS[i+1]-matContainers->RHS[i])-1);
            range2 = doublevec(0, (matContainers->RHS[i+1]-matContainers->RHS[i])-1);
            for (j=matContainers->RHS[i]; j<=matContainers->RHS[i+1]-1; j++) {
                range1[j] = matContainers->Cols[j];
                range2[j] = matContainers->Values[j];
            }
            sort(matContainers->RHS[i+1]-matContainers->RHS[i], range1-1, range2-1);
            for (j=matContainers->RHS[i]; j<=matContainers->RHS[i+1]-1; j++) {
                matContainers->Cols[j] = range1[j];
                matContainers->Values[j] = range2[j];
            }
            free_ivector(range1, 0, (matContainers->RHS[i+1]-matContainers->RHS[i])-1);
            free_dvector(range2, 0, (matContainers->RHS[i+1]-matContainers->RHS[i])-1);
        }
        for (i=0; i<n; i++) {
            for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                if (matContainers->Cols[j] == i) {
                    matContainers->Diag[i] = j;
                    break;
                }
            }
        }
        matrix.ordered = YES;
    }
    
    for (i=0; i<n; i++) {
        if (fabs(matContainers->Values[matContainers->Diag[i]]) > AEPS) {
            u[i] = v[i] / matContainers->Values[matContainers->Diag[i]];
        } else {
            u[i] = v[i];
        }
    }
}


/*******************************************************************************************
 
    Description: Diagonal preconditioning of a CRS format matrix. Matrix is accessed from
    the FEMSolution class.
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix.
 
        double u               ->  Resulting approximate solution after preconditioning
 
        double v               -> Given right-hand-side
 
        int ipar               -> Input stucture holding info from the HUTIter iterative solver
 
*******************************************************************************************/
-(void)CRSComplexDiagPreconditionMatrix:(FEMMatrix *)matrix afterPrecondition:(double complex *)u rightHandSide:(double complex *)v info:(int *)ipar {
    
    int i, j, n;
    double complex A;
    int *range1;
    double *range2;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    n = matrix.numberOfRows; 
    
    if (matrix.isOrdered == NO) {
        for (i=0; i<n; i++) {
            range1 = intvec(0, (matContainers->RHS[i+1]-matContainers->RHS[i])-1);
            range2 = doublevec(0, (matContainers->RHS[i+1]-matContainers->RHS[i])-1);
            for (j=matContainers->RHS[i]; j<=matContainers->RHS[i+1]-1; j++) {
                range1[j] = matContainers->Cols[j];
                range2[j] = matContainers->Values[j];
            }
            sort(matContainers->RHS[i+1]-matContainers->RHS[i], range1-1, range2-1);
            for (j=matContainers->RHS[i]; j<=matContainers->RHS[i+1]-1; j++) {
                matContainers->Cols[j] = range1[j];
                matContainers->Values[j] = range2[j];
            }
            free_ivector(range1, 0, (matContainers->RHS[i+1]-matContainers->RHS[i])-1);
            free_dvector(range2, 0, (matContainers->RHS[i+1]-matContainers->RHS[i])-1);
        }
        for (i=0; i<n; i++) {
            for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                if (matContainers->Cols[j] == i) {
                    matContainers->Diag[i] = j;
                    break;
                }
            }
        }
        matrix.ordered = YES;
    }
    
    for (i=0; n/2; i++) {
        A = matContainers->Values[matContainers->Diag[2*i-1]] + (-matContainers->Values[matContainers->Diag[2*i-1]+1] * I);
        u[i] = v[i] / A;
    }
}

/*******************************************************************************************
    Description:
    Pics the block diagonal entries form matrix solution.matrix to build matrix B.
 
        FEMMatrix *matrix      -> Input matrix.
        FEMMatrix *B           ->  The block diagonal matrix
        int blocks             ->  Number of blocks used in the decomposition
 
*******************************************************************************************/

-(void)CRSBlockDiagonalMatrix:(FEMMatrix *)matrix blockDiagMatrix:(FEMMatrix *)B numberOfBlocks:(int)blocks {
    
    int n;
    int i, k, l, kb;
    matrixArraysContainer *matContainers = NULL, *bContainers = NULL;
    
    if (blocks <= 1) return;
    
    matContainers = matrix.getContainers;
    bContainers = B.getContainers;
    
    n = matrix.numberOfRows;
    B.numberOfRows = n;
    
    kb = 0;
    for (i=0; i<n; i++) {
        for (k=matContainers->Rows[i]; k<=matContainers->Rows[i+1]-1; k++) {
            l = matContainers->Cols[k];
            if ((i % blocks) == (l % blocks)) kb++;
        }
    }
    
    bContainers->Rows = intvec(0, (n+1)-1);
    bContainers->sizeRows = n+1;
    bContainers->Cols = intvec(0, kb-1);
    bContainers->sizeCols = kb;
    bContainers->Values = doublevec(0, kb-1);
    bContainers->sizeValues = kb;
    bContainers->Diag = intvec(0, n-1);
    bContainers->sizeDiag = n;
    
    kb = 0;
    for (i=0; i<n; i++) {
        bContainers->Rows[i] = kb;
        for (k=matContainers->Rows[i]; k<=matContainers->Rows[i+1]-1; k++) {
            l = matContainers->Cols[k];
            if ((i % blocks) == (l % blocks)) {
                bContainers->Values[kb] = matContainers->Values[k];
                bContainers->Cols[kb] = matContainers->Cols[k];
                if (bContainers->Cols[kb] == i) bContainers->Diag[i] = kb;
                kb++;
            }
        }
    }
    bContainers->Rows[(n+1)-1] = kb;
}

-(void)initializeILUMatrix:(FEMMatrix *)matrix numberOfRows:(int)ilun {
    
    int i, n, m;
    
    FEMMatrix *a1;
    matrixArraysContainer *matContainers = NULL;
    matrixArraysContainer *a1Containers = NULL;
    
    a1 = [[FEMMatrix alloc] init];
    matContainers = matrix.getContainers;
    a1Containers = a1.getContainers;
    
    n = matrix.numberOfRows;
    
    if (ilun == 0) {
        matContainers->ILURows = matContainers->Rows;
        matContainers->ILUCols = matContainers->Cols;
        matContainers->ILUDiag = matContainers->Diag;
        
        matContainers->sizeILURows = matContainers->sizeRows;
        matContainers->sizeILUCols = matContainers->sizeCols;
        matContainers->sizeILUDiag = matContainers->sizeDiag;
    } else {
        
        matContainers->sizeILUCols = [self FEMPrecondition_initializeILU1:matContainers :n];
        matContainers->sizeILURows = n+1;
        matContainers->sizeILUDiag = n;
        
        if (ilun > 1) {
            
            for (i=0; i<ilun-1; i++) {
                
                a1Containers->Cols = matContainers->ILUCols;
                a1Containers->Rows = matContainers->ILURows;
                a1Containers->Diag = matContainers->ILUDiag;
                
                m = [self FEMPrecondition_initializeILU1:a1Containers :n];
                
                matContainers->ILUCols = a1Containers->ILUCols;
                matContainers->ILURows = a1Containers->ILURows;
                matContainers->ILUDiag = a1Containers->ILUDiag;
                
                a1Containers->Cols = NULL;
                a1Containers->Rows = NULL;
                a1Containers->Diag = NULL;
            }
        }
    }
    matContainers->ILUValues = doublevec(0, matContainers->ILURows[(n+1)-1]-1);
    if (matContainers->ILUValues == NULL) {
        errorfunct("initializeILUMatrix", "Memory allocation error.");
    }
    matContainers->sizeILUValues = matContainers->ILURows[(n+1)-1]-1;
}

-(void)initializeCILUMatrix:(FEMMatrix *)matrix :(int)ilun {
    
    int i, j, k;
    int n;
    
    FEMMatrix *a1;
    matrixArraysContainer *matContainers = NULL;
    matrixArraysContainer *a1Containers = NULL;
    
    a1 = [[FEMMatrix alloc] init];
    matContainers = matrix.getContainers;
    a1Containers = a1.getContainers;
    
    n = matrix.numberOfRows;
    
    a1.numberOfRows = n/2;
    a1Containers->Rows = intvec(0, (n/2+1)-1);
    a1Containers->Diag = intvec(0, (n/2)-1);
    a1Containers->Cols = intvec(0, (matContainers->sizeCols / 4)-1);
    
    a1Containers->Rows[0] = 0 ;
    k = 0;
    
    for (i=0; i<n; i+=2) {
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j+=2) {
            k = k + 1;
            a1Containers->Cols[k] = (matContainers->Cols[j]+1) / 2;
            if (matContainers->Cols[j] == i) a1Containers->Diag[(i+1)/2] = k;
        }
        a1Containers->Rows[(i+1)/2+1] = k + 1;
    }
    
    if (ilun == 0) {
        matContainers->ILURows = matContainers->Rows;
        matContainers->ILUCols = matContainers->Cols;
        matContainers->ILUDiag = matContainers->Diag;
        
        matContainers->sizeILURows = matContainers->sizeRows;
        matContainers->sizeILUCols = matContainers->sizeCols;
        matContainers->sizeILUDiag = matContainers->sizeDiag;
    } else {
        
        matContainers->sizeILUCols = [self FEMPrecondition_initializeILU1:a1Containers :n/2];
        matContainers->ILUCols = a1Containers->ILUCols;
        matContainers->ILURows = a1Containers->ILURows;
        matContainers->ILUDiag = a1Containers->ILUDiag;
        
        matContainers->sizeILURows = (n/2+1);
        matContainers->sizeILUDiag = n/2;
        
        free_ivector(a1Containers->Rows, 0, (n/2+1)-1);
        free_ivector(a1Containers->Diag, 0, (n/2)-1);
        free_ivector(a1Containers->Cols, 0, (matContainers->sizeCols / 4)-1);
        a1Containers->Rows = NULL;
        a1Containers->Diag = NULL;
        a1Containers->Cols = NULL;
        
        if (ilun > 1) {
            
            for (i=0; i<ilun-1; i++) {
                
                a1Containers->Cols = matContainers->ILUCols;
                a1Containers->Rows = matContainers->ILURows;
                a1Containers->Diag = matContainers->ILUDiag;
                
                k = [self FEMPrecondition_initializeILU1:a1Containers :n/2];
                
                matContainers->ILUCols = a1Containers->ILUCols;
                matContainers->ILURows = a1Containers->ILURows;
                matContainers->ILUDiag = a1Containers->ILUDiag;
                
                a1Containers->Cols = NULL;
                a1Containers->Rows = NULL;
                a1Containers->Diag = NULL;
            }
        }
    }
    matContainers->CILUValues = cdoublevec(0, matContainers->ILURows[(n/2+1)-1]);
    if (matContainers->CILUValues == NULL) {
        errorfunct("initializeCILUMatrix", "Memory allocation error.");
    }
    matContainers->sizeCILUValues = matContainers->ILURows[(n/2+1)-1];
}


/*******************************************************************************************
 
    Description:
    Builds an incomplete ILU(n) factorization for an iterative solver precondioner.
    Real matrix version.
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix
 
        int ilun               -> Order of fills allowed 0-9
 
        Return Value           -> A BOOL whether or not the factorization succeeded.
 
*******************************************************************************************/
-(BOOL)CRSIncompleteLUMatrix:(FEMMatrix *)matrix fillsOrder:(int)ilun {
    
    int i, j, k, l, n;
    double t, *S;
    bool *C;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    errorfunct("CRSIncompleteLUMatrix", "ILU (Real). Starting factorization with order:", ilun);
    t = cputime();
    
    n = matrix.numberOfRows;
    
    if (matContainers->ILUValues == NULL) {
        
        [self initializeILUMatrix:matrix numberOfRows:ilun];
    }
    
    // Allocate space for storing one full row
    C = boolvec(0, n-1);
    S = doublevec(0, n-1);
    memset( C, false, n*sizeof(bool) );
    memset( S, 0.0, n*sizeof(double) );
    
    // The factorization row by row
    for (i=0; i<n; i++) {
        
        // Convert current row to full form for speed, 
        // only flagging the nonzeros entries
        for (k=matContainers->Rows[i]; k<=matContainers->Rows[i+1]-1; k++) {
            S[matContainers->Cols[k]] = matContainers->Values[k];
        }
        for (k=matContainers->ILURows[i]; k<=matContainers->ILURows[i+1]-1; k++) {
            C[matContainers->ILUCols[k]] = true;
        }
        
        // This is the factorization part for the current row
        for (k=matContainers->ILUCols[matContainers->ILURows[i]]; k<=i-1; k++) {
            if (C[k] == true) {
                if (fabs(matContainers->ILUValues[matContainers->ILUDiag[k]]) > AEPS) S[k] = S[k] / matContainers->ILUValues[matContainers->ILUDiag[k]];
                
                for (l=matContainers->ILUDiag[k]; l<=matContainers->ILURows[k+1]-1; l++) {
                    j = matContainers->ILUCols[l];
                    if (C[j] == true) S[j] = S[j] - S[k] * matContainers->ILUValues[l];
                }
            }
        }
        
        // Convert the row back to CRS format
        for (k=matContainers->ILURows[i]; k<=matContainers->ILURows[i+1]-1; k++) {
            if (C[matContainers->ILUCols[k]] == true) {
                matContainers->ILUValues[k] = S[matContainers->ILUCols[k]];
                S[matContainers->ILUCols[k]] = 0.0;
                C[matContainers->ILUCols[k]] = false;
            }
        }
    }
    
    free_bvector(C, 0, n-1);
    free_dvector(S, 0, n-1);
    
    // Prescale the diagonal for the LU solver
    for (i=0; i<n; i++) {
        if (fabs(matContainers->ILUValues[matContainers->ILUDiag[i]]) < AEPS) {
            matContainers->ILUValues[matContainers->ILUDiag[i]] = 1.0;
        } else {
            matContainers->ILUValues[matContainers->ILUDiag[i]] = 1.0 / matContainers->ILUValues[matContainers->ILUDiag[i]];
        }
    }
    
    warnfunct("CRSIncompleteLUMatrix", "ILU (Real), NOF nonzeros", matContainers->ILURows[(n+1)-1]);
    warnfunct("CRSIncompleteLUMatrix", "ILU (Real), Filling (%):", floor(matContainers->ILURows[(n+1)-1]) * (100.0 / matContainers->Rows[(n+1)-1]));
    warnfunct("CRSIncompleteLUMatrix", "ILU (Real), Factorization ready at (s):", cputime() - t);
    
    return YES;
}

/*******************************************************************************************
 
    Description:
    Builds an incomplete ILU(n) factorization for an iterative solver precondioner.
    Complex matrix version.
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix
 
        int ilun               -> Order of fills allowed 0-9
 
        Return Value           -> A BOOL whether or not the factorization succeeded.
 
*******************************************************************************************/
-(BOOL)CRSComplexIncompleteLUMatrix:(FEMMatrix *)matrix fillsOrder: (int)ilun {
    
    int i, j, k, l, n;
    double t;
    double complex *S;
    bool *C;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    errorfunct("CRSComplexIncompleteLUMatrix", "ILU (Complex), Starting factorization with order:", ilun);
    t = cputime();
    
    n = matrix.numberOfRows;
    
    if (matContainers->CILUValues == NULL) {
        
        [self initializeCILUMatrix:matrix :ilun];
    }
    
    // Allocate space for storing one full row
    C = boolvec(0, n/2-1);
    S = cdoublevec(0, n/2-1);
    memset( C, false, n*sizeof(bool) );
    memset( S, 0.0, n*sizeof(double complex) );
    
    // The factorization row by row
    for (i=0; i<n/2; i++) {
        
        // Convert current row to full form for speed, 
        // only flagging the nonzeros entries
        for (k=matContainers->ILURows[i]; k<=matContainers->ILURows[i+1]-1; k++) {
            C[matContainers->ILUCols[k]] = true;
        }
        for (k=matContainers->Rows[2*i-1]; k<=matContainers->Rows[2*i]-1; k+=2) {
            S[(matContainers->Cols[k]+1)/2] = matContainers->Values[k] + (-matContainers->Values[k+1] * I);
        }
        
        // This is the factorization part for the current row
        for (k=matContainers->ILUCols[matContainers->ILURows[i]]; k<=i-1; k++) {
            if (C[k] == true) {
                if (fabs(matContainers->CILUValues[matContainers->ILUDiag[k]]) > AEPS) S[k] = S[k] / matContainers->CILUValues[matContainers->ILUDiag[k]];
                
                for (l=matContainers->ILUDiag[k]+1; l<=matContainers->ILURows[k+1]-1; l++) {
                    j = matContainers->ILUCols[l];
                    if (C[j] == 1) S[j] = S[j] - S[k] * matContainers->CILUValues[l];
                }
            }
        }
        
        // Convert the row back to CRS format
        for (k=matContainers->ILURows[i]; k<=matContainers->ILURows[i+1]-1; k++) {
            if (C[matContainers->ILUCols[k]] == true) {
                matContainers->CILUValues[k] = S[matContainers->ILUCols[k]];
                S[matContainers->ILUCols[k]] = 0.0;
                C[matContainers->ILUCols[k]] = false;
            }
        }
    }
    
    free_bvector(C, 0,  n/2-1);
    free_cdvector(S, 0, n/2-1);
    
    // Prescale the diagonal for the LU solver
    for (i=0; i<n/2; i++) {
        if (fabs(matContainers->CILUValues[matContainers->ILUDiag[i]]) < AEPS) {
            matContainers->CILUValues[matContainers->ILUDiag[i]] = 1.0;
        } else {
            matContainers->CILUValues[matContainers->ILUDiag[i]] = 1.0 / matContainers->CILUValues[matContainers->ILUDiag[i]];
        }
    }
    
    warnfunct("CRSComplexIncompleteLUMatrix", "ILU (Complex), NOF nonzeros", matContainers->ILURows[(n/2+1)-1]);
    warnfunct("CRSComplexIncompleteLUMatrix", "ILU (Complex), Filling (%):", floor(matContainers->ILURows[(n/2+1)-1]) * (400.0 / matContainers->Rows[(n+1)-1]));
    warnfunct("CRSComplexIncompleteLUMatrix", "ILU (Complex), Factorization ready at (s):", cputime() - t);
    
    return YES;
}

/*******************************************************************************************
 
    Description:
    Builds an incomplete (ILUT) factorization for an iterative solver preconditioner.
    Real matrix version
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix
 
        int tol                -> Drop tolerance: if ILUT(i,j) <= NORM(A(i.;))*tol
                                  the value is dropped.
 
        Return Value           -> A BOOL whether or not the factorization succeeded.
 
*******************************************************************************************/
-(BOOL)CRSIlutMatrix:(FEMMatrix *)matrix dropTolerance:(int)tol {
    
    int n;
    double t;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    warnfunct("CRSIlutMatrix", "Starting factorization:");
    t = cputime();
    
    n = matrix.numberOfRows;
    
    if (matContainers->ILUValues != NULL) {
        free_ivector(matContainers->ILURows, 0, matContainers->sizeILURows-1);
        free_ivector(matContainers->ILUDiag, 0, matContainers->sizeILUDiag-1);
        free_ivector(matContainers->ILUCols, 0, matContainers->sizeILUCols-1);
        free_dvector(matContainers->ILUValues, 0, matContainers->sizeILUValues-1);
    }
    
    [self FEMPrecondition_computeIlutMatrix:matrix numberOfRows:n tolerance:tol];
    
    warnfunct("CRSIlutMatrix", "ILU(T) (Real), NOF nonzeros", matContainers->ILURows[(n+1)-1]);
    warnfunct("CRSIlutMatrix", "ILU(T) (Real), Filling (%)", floor(matContainers->ILURows[(n+1)-1]) * (100.0 / matContainers->Rows[(n+1)-1]));
    warnfunct("CRSIlutMatrix", "ILU(T) (Real), Factorization ready at (s):", cputime() - t);
    
    return YES;
}

/*******************************************************************************************
 
    Description:
    Builds an incomplete (ILUT) factorization for an iterative solver preconditioner.
    Complex matrix version
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix
 
        int tol                -> Drop tolerance: if ILUT(i,j) <= NORM(A(i.;))*tol
                                  the value is dropped.
 
        Return Value           -> A BOOL whether or not the factorization succeeded.
 
*******************************************************************************************/
-(BOOL)CRSComplexIlutMatrix:(FEMMatrix *)matrix dropTolerance:(int)tol {
    
    int n;
    double t;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    warnfunct("CRSComplexIlutMatrix", "Starting factorization:");
    t = cputime();
    
    n = matrix.numberOfRows / 2;
    
    if (matContainers->CILUValues != NULL) {
        free_ivector(matContainers->ILURows, 0, matContainers->sizeILURows-1);
        free_ivector(matContainers->ILUDiag, 0, matContainers->sizeILUDiag-1);
        free_ivector(matContainers->ILUCols, 0, matContainers->sizeILUCols-1);
        free_cdvector(matContainers->CILUValues, 0, matContainers->sizeCILUValues-1);
    }
    
    [self FEMPrecondition_computeComplexIlutMatrix:matrix numberOfRows:n tolerance:tol];
    
    warnfunct("CRSComplexIlutMatrix", "ILU(T) (Complex), NOF nonzeros", matContainers->ILURows[(n+1)-1]);
    warnfunct("CRSComplexIlutMatrix", "ILU(T) (Complex), Filling (%)", floor(matContainers->ILURows[(n+1)-1]) * (400.0 / matContainers->Rows[(2*n+1)-1]));
    warnfunct("CRSComplexIlutMatrix", "ILU(T) (Complex), Factorization ready at (s):", cputime() - t);
    
    return YES;
}

/*******************************************************************************************
 
    Description:
    Incomplete factorization preconditioner solver for a CRS format matrix.
    Matrix is accessed from the solution class. Real matrix version.
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix
 
        double *u              -> Solution vector
 
        double *v              -> Right-hand-side vector
 
        int ipar               -> Structure holding info from the HUTIter iterative solver.
 
*******************************************************************************************/
-(BOOL)CRSLuPreconditionMatrix:(FEMMatrix *)matrix afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar {
    
    memcpy(u, v, ipar[2]*sizeof(double));
    
    [self FEMPrecondition_LUSolveSystemSize:ipar[2] matrix:matrix rightHandSide:u];
    
    return YES;
}

/*****************************************************************************************************************
 
    Description:
    Incomplete factorization preconditioner solver for a CRS format matrix.
    Matrix is accessed from the solution class. Real matrix version.
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix
 
        double complex *u      -> Solution vector
 
        double complex *v      -> Right-hand-side vector
 
        int ipar               -> Structure holding info from the HUTIter iterative solver.
 
****************************************************************************************************************/
-(BOOL)CRSComplexLuPreconditionMatrix:(FEMMatrix *)matrix afterPrecondition:(double complex *)u rightHandSide:(double complex *)v info:(int *)ipar {
    
    memcpy(u, v, ipar[2]*sizeof(double));
    
    [self FEMPrecondition_ComplexLUSolveSystemSize:ipar[2] matrix:matrix rightHandSide:u];

    return YES;
}

/*******************************************************************************************
 
    Description:
    Matrix vector product (v = Au) for a matrix given in CRS format. The matrix
    is accessed from the solution class. Real version.
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix
 
        double *u              -> Vector to multiply
 
        double *v              -> Result vector
 
        int ipar               -> Structure holding info from the HUTIter iterative solver.
 
*******************************************************************************************/
-(void)CRSMatrixVectorProduct:(FEMMatrix *)matrix vector:(double *)u result:(double *)v info:(int *)ipar {
    
    int i, j, n;
    double s;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    n = matrix.numberOfRows;
    
    if (ipar[5] == 0) {
        
        for (i=0; i<n; i++) {
            s = 0.0;
            for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                s = s + matContainers->Values[j] * u[matContainers->Cols[j]];
            }
            v[i] = s;
        }
    } else {
        memset( v, 0.0, n*sizeof(double) );
        for (i=0; i<n; i++) {
            s = u[i];
            for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                v[matContainers->Cols[j]] = v[matContainers->Cols[j]] + s * matContainers->Values[j];
            }
        }
    }
}

/*******************************************************************************************
 
    Description:
    Matrix vector product (v = Au) for a matrix given in CRS format. The matrix
    is accessed from the solution class. Complex version.
 
    Arguments:
 
        FEMMatrix *matrix      -> Input matrix
 
        double complex *u      -> Vector to multiply
 
        double complex *v      -> Result vector
 
        int ipar               -> Structure holding info from the HUTIter iterative solver.
 
*******************************************************************************************/
-(void)CRSComplexMatrixVectorProduct:(FEMMatrix *)matrix vector:(double complex *)u result:(double complex *)v info:(int *)ipar {
    
    int i, j, n;
    double complex s, rsum;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    n = ipar[2];
    
    if (ipar[5] == 0) {
        
        for (i=0; i<n; i++) {
            rsum = 0.0 + 0.0 * I;
            for (j=matContainers->Rows[2*i-1]; j<=matContainers->Rows[2*i]-1; i+=2) {
                s = matContainers->Values[j] + (-matContainers->Values[j+1] * I);
                rsum = rsum + s * u[(matContainers->Cols[j]+1)/2];
            }
            v[i] = rsum;
        }
    } else {
        for (i=0; i<n; i++) {
            v[i] = 0.0 + 0.0 * I; 
        }
        for (i=0; i<n; i++) {
            rsum = u[i];
            for (j=matContainers->Rows[2*i-1]; j<=matContainers->Rows[2*i]-1; j+=2) {
                s = matContainers->Values[j] + (-matContainers->Values[j+1] * I);
                v[(matContainers->Cols[j]+1)/2] = v[(matContainers->Cols[j]+1)/2] + s * rsum;
            }
        }
    }
}

-(void)CRSPCondDummyMatrix:(FEMMatrix *)matrix afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar {
    
    memcpy(u, v, ipar[2]*sizeof(double));
}

@end
