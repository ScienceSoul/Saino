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

static double AEPS = 10.0 * DBL_EPSILON;

@interface FEMPrecondition ()

-(void)FEMPrecondition_computeIlutInSolution:(FEMSolution *)solution numberOfRows:(int)n tolerance:(int)tol;
-(void)FEMPrecondition_computeComplexIlutInSolution:(FEMSolution *)solution numberOfRows:(int)n tolerance:(int)tol;
-(void)FEMPrecondition_LUSolveSystemSize:(int)n solution:(FEMSolution *)solution rightHandSide:(double *)b;
-(void)FEMPrecondition_ComplexLUSolveSystemSize:(int)n solution:(FEMSolution *)solution rightHandSide:(double complex *)b;

@end

@implementation FEMPrecondition

#pragma mark Private methods

-(void)FEMPrecondition_computeIlutInSolution:(FEMSolution *)solution numberOfRows:(int)n tolerance:(int)tol {
    
    int i, j, k, l, rowMin, rowMax;
    bool *C;
    const double WORKN = 128;
    double norma, cptime, ttime, t;
    double *S;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    ttime = cputime();
    cptime = 0.0;
    
    matContainers->ILURows = intvec(0, (n+1)-1);
    matContainers->ILUDiag = intvec(0, n-1);
    if (matContainers->ILURows == NULL || matContainers->ILUDiag == NULL) {
        errorfunct("computeILUT", "Memory allocation error.");
    }
    matContainers->sizeILURows = n+1;
    matContainers->sizeILUDiag = n;
    
    matContainers->ILUCols = intvec(0, (WORKN*n)-1);
    matContainers->ILUValues = doublevec(0, (WORKN*n)-1);
    if (matContainers->ILUCols == NULL || matContainers->ILUValues == NO) {
        errorfunct("computeILUT", "Memory allocation error.");
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
                [solution ilutWorkspaceCheckAtIndex:i numberOfRows:n];
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

-(void)FEMPrecondition_computeComplexIlutInSolution:(FEMSolution *)solution numberOfRows:(int)n tolerance:(int)tol {
    
    int i, j, k, l, rowMin, rowMax;
    bool *C;
    const double WORKN = 128;
    double norma;
    double complex *S;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    matContainers->ILURows = intvec(0, (n+1)-1);
    matContainers->ILUDiag = intvec(0, n-1);
    if (matContainers->ILURows == NULL || matContainers->ILUDiag == NULL) {
        errorfunct("computeILUT", "Memory allocation error.");
    }
    matContainers->sizeILURows = n+1;
    matContainers->sizeILUDiag = n;
    
    matContainers->ILUCols = intvec(0, (WORKN*n)-1);
    matContainers->CILUValues = cdoublevec(0, (WORKN*n)-1);
    if (matContainers->ILUCols == NULL || matContainers->CILUValues == NULL) {
        errorfunct("computeILUT", "Memory allocation error.");
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
                [solution ilutComplexWorkspaceCheckAtIndex:i numberOfRows:n];
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

-(void)FEMPrecondition_LUSolveSystemSize:(int)n solution:(FEMSolution *)solution rightHandSide:(double *)b {
/*******************************************************************************************
 
    Description:
    Solve a system  (Ax=b) after factorization A=LUD has been done. This method is
    meant as a part of a preconditioner for an iterative solver. Real version
 
    Arguments:
 
    int n                    -> Size of the system.
 
    FEMSolution *solution    -> Class holding input matrix.
 
    double *b                -> On entry the RHS vector, on exit the solution vector.
 
*******************************************************************************************/
    
    int i, j;
    double s;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
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

-(void)FEMPrecondition_ComplexLUSolveSystemSize:(int)n solution:(FEMSolution *)solution rightHandSide:(double complex *)b {
/*******************************************************************************************
 
    Description:
    Solve a system  (Ax=b) after factorization A=LUD has been done. This method is
    meant as a part of a preconditioner for an iterative solver. Complex version
 
    Arguments:
 
    int n                  -> Size of the system.
 
    FEMSolution *solution  -> Class holding input matrix.
 
    double complex *b      -> On entry the RHS vector, on exit the solution vector.
 
*******************************************************************************************/
    
    int i, j;
    double complex s, x;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
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


-(void)CRSDiagPreconditionInSolution:(FEMSolution *)solution afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar {
/*******************************************************************************************
                                                                             
    Description: Diagonal preconditioning of a CRS format matrix. Matrix is accessed from   
    the FEMSolution class.
 
    Arguments:
       
        double u  ->  Resulting approximate solution after preconditioning

        double v  -> Given right-hand-side

        int ipar  -> Input stucture holding info from the HUTIter iterative solver
 
*******************************************************************************************/
    int i, j, n;
    int *range1;
    double *range2;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    n = solution.matrix.numberOfRows;  
    
    if (solution.matrix.isOrdered == NO) {
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
        solution.matrix.ordered = YES;
    }
    
    for (i=0; i<n; i++) {
        if (fabs(matContainers->Values[matContainers->Diag[i]]) > AEPS) {
            u[i] = v[i] / matContainers->Values[matContainers->Diag[i]];
        } else {
            u[i] = v[i];
        }
    }
}

-(void)CRSComplexDiagPreconditionInSolution:(FEMSolution *)solution afterPrecondition:(double complex *)u rightHandSide:(double complex *)v info:(int *)ipar {
/*******************************************************************************************
 
    Description: Diagonal preconditioning of a CRS format matrix. Matrix is accessed from   
    the FEMSolution class.
 
    Arguments:
 
        double u  ->  Resulting approximate solution after preconditioning
 
        double v  -> Given right-hand-side
 
        int ipar  -> Input stucture holding info from the HUTIter iterative solver

*******************************************************************************************/
    
    int i, j, n;
    double complex A;
    int *range1;
    double *range2;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    n = solution.matrix.numberOfRows; 
    
    if (solution.matrix.isOrdered == NO) {
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
        solution.matrix.ordered = YES;
    }
    
    for (i=0; n/2; i++) {
        A = matContainers->Values[matContainers->Diag[2*i-1]] + (-matContainers->Values[matContainers->Diag[2*i-1]+1] * I);
        u[i] = v[i] / A;
    }
}

-(void)CRSBlockDiagonalInSolution:(FEMSolution *)solution blockDiagMatrix:(FEMMatrix *) B numberOfBlocks:(int)blocks {
/*******************************************************************************************
    Description: 
        Pics the block diagonal entries form matrix solution.matrix to build matrix B.
 
        FEMMatrix *B  ->  The block diagonal matrix
        int blocks    ->  Number of blocks used in the decomposition

*******************************************************************************************/
    
    int n;
    int i, k, l, kb;
    matrixArraysContainer *matContainers = NULL, *bContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    bContainers = B.getContainers;
    
    n = solution.matrix.numberOfRows;
    B.numberOfRows = n;
    
    kb = 0;
    for (i=0; i<n; i++) {
        for (k=matContainers->Rows[i]; k<=matContainers->Rows[i+1]-1; k++) {
            l = matContainers->Cols[k];
            if ((i % blocks) == (l % blocks)) kb = kb + 1;
        }
    }
    
    // From here we need that B.Rows, B.Cols, B.Values and B.Diag have allocation in memory.
    // This should be taken care by the caller of this method.
    
    kb = 0;
    for (i=0; i<n; i++) {
        bContainers->Rows[i] = kb;
        for (k=matContainers->Rows[i]; k<=matContainers->Rows[i+1]-1; k++) {
            l = matContainers->Cols[k];
            if ((i % blocks) == (l % blocks)) {
                bContainers->Values[kb] = matContainers->Values[k];
                bContainers->Cols[kb] = matContainers->Cols[k];
                if (bContainers->Cols[kb] == i) bContainers->Diag[i] = kb;
                kb = kb + 1;
            }
        }
    }
    bContainers->Rows[(n-1)+1] = kb;
}

-(BOOL)CRSIncompleteLUInSolution:(FEMSolution *)solution fillsOrder:(int)ilun {
/*******************************************************************************************
 
    Description: 
        Builds an incomplete ILU(n) factorization for an iterative solver precondioner.
        Real matrix version.
 
    Arguments:
 
    FEMSolution *solution  -> Class holding input matrix, will also hold
                              the factorization on exit.
 
    int ilun               -> Order of fills allowed 0-9
 
    Return Value           -> A BOOL whether or not the factorization succeeded.
 
 *******************************************************************************************/
    
    int i, j, k, l, n;
    double t, *S;
    bool *C;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    errorfunct("CRSIncompleteLUInSolution", "ILU (Real). Starting factorization with order:", ilun);
    t = cputime();
    
    n = solution.matrix.numberOfRows;
    
    if (matContainers->ILUValues == NULL) {
        
        [solution initializeILU:ilun];
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
    
    warnfunct("CRSIncompleteLUInSolution", "ILU (Real), NOF nonzeros", matContainers->ILURows[(n+1)-1]);
    warnfunct("CRSIncompleteLUInSolution", "ILU (Real), Filling (%):", floor(matContainers->ILURows[(n+1)-1]) * (100.0 / matContainers->Rows[(n+1)-1]));
    warnfunct("CRSIncompleteLUInSolution", "ILU (Real), Factorization ready at (s):", cputime() - t);
    
    return YES;
}

-(BOOL)CRSComplexIncompleteLUInSolution:(FEMSolution *)solution fillsOrder: (int)ilun {
/*******************************************************************************************
 
    Description: 
        Builds an incomplete ILU(n) factorization for an iterative solver precondioner.
        Complex matrix version.
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix, will also hold
                                  the factorization on exit.
 
        int ilun               -> Order of fills allowed 0-9
 
        Return Value           -> A BOOL whether or not the factorization succeeded.
 
*******************************************************************************************/
    
    int i, j, k, l, n;
    double t;
    double complex *S;
    bool *C;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    errorfunct("CRSComplexIncompleteLUInSolution", "ILU (Complex), Starting factorization with order:", ilun);
    t = cputime();
    
    n = solution.matrix.numberOfRows;
    
    if (matContainers->CILUValues == NULL) {
        
        [solution initializeCILU:ilun];
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
    
    warnfunct("CRSComplexIncompleteLUInSolution", "ILU (Complex), NOF nonzeros", matContainers->ILURows[(n/2+1)-1]);
    warnfunct("CRSComplexIncompleteLUInSolution", "ILU (Complex), Filling (%):", floor(matContainers->ILURows[(n/2+1)-1]) * (400.0 / matContainers->Rows[(n+1)-1]));
    warnfunct("CRSComplexIncompleteLUInSolution", "ILU (Complex), Factorization ready at (s):", cputime() - t);
    
    return YES;
}

-(BOOL)CRSIlutInSolution:(FEMSolution *)solution dropTolerance:(int)tol {
/*******************************************************************************************
 
    Description: 
        Builds an incomplete (ILUT) factorization for an iterative solver preconditioner.
        Real matrix version
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix, will also hold
                                  the factorization on exit.
 
        int tol                -> Drop tolerance: if ILUT(i,j) <= NORM(A(i.;))*tol
                                  the value is dropped.
 
        Return Value           -> A BOOL whether or not the factorization succeeded.
 
*******************************************************************************************/
    
    int n; 
    double t;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    warnfunct("CRS_ICRSIlutInSolutionLUT", "Starting factorization:");
    t = cputime();
    
    n = solution.matrix.numberOfRows;
    
    if (matContainers->ILUValues != NULL) {
        free_ivector(matContainers->ILURows, 0, matContainers->sizeILURows-1);
        free_ivector(matContainers->ILUDiag, 0, matContainers->sizeILUDiag-1);
        free_ivector(matContainers->ILUCols, 0, matContainers->sizeILUCols-1);
        free_dvector(matContainers->ILUValues, 0, matContainers->sizeILUValues-1);
    }
    
    [self FEMPrecondition_computeIlutInSolution:solution numberOfRows:n tolerance:tol];
    
    warnfunct("CRSIlutInSolution", "ILU(T) (Real), NOF nonzeros", matContainers->ILURows[(n+1)-1]);
    warnfunct("CRSIlutInSolution", "ILU(T) (Real), Filling (%)", floor(matContainers->ILURows[(n+1)-1]) * (100.0 / matContainers->Rows[(n+1)-1]));
    warnfunct("CRSIlutInSolution", "ILU(T) (Real), Factorization ready at (s):", cputime() - t);
    
    return YES;
}

-(BOOL)CRSComplexIlutInSolution:(FEMSolution *)solution dropTolerance:(int)tol {
/*******************************************************************************************
 
    Description: 
        Builds an incomplete (ILUT) factorization for an iterative solver preconditioner.
        Complex matrix version
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix, will also hold
        the factorization on exit.
 
        int tol                -> Drop tolerance: if ILUT(i,j) <= NORM(A(i.;))*tol
                                  the value is dropped.
 
        Return Value           -> A BOOL whether or not the factorization succeeded.
 
*******************************************************************************************/
    
    int n; 
    double t;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    warnfunct("CRSComplexIlutInSolution", "Starting factorization:");
    t = cputime();
    
    n = solution.matrix.numberOfRows / 2;
    
    if (matContainers->CILUValues != NULL) {
        free_ivector(matContainers->ILURows, 0, matContainers->sizeILURows-1);
        free_ivector(matContainers->ILUDiag, 0, matContainers->sizeILUDiag-1);
        free_ivector(matContainers->ILUCols, 0, matContainers->sizeILUCols-1);
        free_cdvector(matContainers->CILUValues, 0, matContainers->sizeCILUValues-1);
    }
    
    [self FEMPrecondition_computeComplexIlutInSolution:solution numberOfRows:n tolerance:tol];
    
    warnfunct("CRSComplexIlutInSolution", "ILU(T) (Complex), NOF nonzeros", matContainers->ILURows[(n+1)-1]);
    warnfunct("CRSComplexIlutInSolution", "ILU(T) (Complex), Filling (%)", floor(matContainers->ILURows[(n+1)-1]) * (400.0 / matContainers->Rows[(2*n+1)-1]));
    warnfunct("CRSComplexIlutInSolution", "ILU(T) (Complex), Factorization ready at (s):", cputime() - t);
    
    return YES;
}

-(BOOL)CRSLuPreconditionInSolution:(FEMSolution *)solution afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar {
/*******************************************************************************************
 
    Description: 
        Incomplete factorization preconditioner solver for a CRS format matrix.
        Matrix is accessed from the solution class. Real matrix version.
 
 Arguments:

    FEMSolution *solution  -> Class holding input matrix.
 
    double *u              -> Solution vector

    double *v              -> Right-hand-side vector

    int ipar               -> Structure holding info from the HUTIter iterative solver.

*******************************************************************************************/

    int i;
    
    for (i=0; i<ipar[2]; i++) {
        u[i] = v[i];
    }
    
    [self FEMPrecondition_LUSolveSystemSize:ipar[2] solution:solution rightHandSide:u];
    
    return YES;
    
}

-(BOOL)CRSComplexLuPreconditionInSolution:(FEMSolution *)solution afterPrecondition:(double complex *)u rightHandSide:(double complex *)v info:(int *)ipar {
/*****************************************************************************************************************
 
    Description: 
        Incomplete factorization preconditioner solver for a CRS format matrix.
        Matrix is accessed from the solution class. Real matrix version.
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix.
 
        double complex *u      -> Solution vector
 
        double complex *v      -> Right-hand-side vector
 
        int ipar               -> Structure holding info from the HUTIter iterative solver.
 
 ****************************************************************************************************************/
    
    int i;
    
    for (i=0; i<ipar[2]; i++) {
        u[i] = v[i];
    }
    
    [self FEMPrecondition_ComplexLUSolveSystemSize:ipar[2] solution:solution rightHandSide:u];

    return YES;
    
}

-(void)CRSMatrixVectorProdInSolution:(FEMSolution *)solution multiplyVector:(double *)u resultVector:(double *)v info:(int *)ipar {
/*******************************************************************************************
 
    Description: 
        Matrix vector product (v = Au) for a matrix given in CRS format. The matrix 
        is accessed from the solution class. Real version.
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix.
 
        double *u              -> Vector to multiply

        double *v              -> Result vector
 
        int ipar               -> Structure holding info from the HUTIter iterative solver.
 
*******************************************************************************************/
    
    int i, j, n;
    double s;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    n = solution.matrix.numberOfRows;
    
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

-(void)CRSComplexMatrixVectorProdInSolution:(FEMSolution *)solution multiplyVector:(double complex *)u resultVector:(double complex *)v info:(int *)ipar {
/*******************************************************************************************
 
    Description: 
        Matrix vector product (v = Au) for a matrix given in CRS format. The matrix 
        is accessed from the solution class. Complex version.
 
    Arguments:
 
    FEMSolution *solution  -> Class holding input matrix.
 
    double complex *u      -> Vector to multiply
 
    double complex *v      -> Result vector
 
    int ipar               -> Structure holding info from the HUTIter iterative solver.
 
*******************************************************************************************/
    
    int i, j, n;
    double complex s, rsum;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
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

-(void)CRSPCondDummyInSolution:(FEMSolution *)solution afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar {
    
    int i;
    
    for (i=0; i<ipar[2]; i++) {
        u[i] = v[i];
    }
}

@end
