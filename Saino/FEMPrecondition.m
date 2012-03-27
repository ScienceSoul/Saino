//
//  FEMPreconditioners.m
//  Saino
//
//  Created by Hakime Seddik on 27/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMPrecondition.h"

#import <float.h>

static double AEPS = 10.0 * DBL_EPSILON;

@interface FEMPrecondition ()

// Diagonal preconditioning
-(void)CRS_DiagPrecondition:(FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar;
-(void)CRS_ComplexDiagPrecondition:(FEMSolution *)solution: (double complex *)u: (double complex *)v: (int *)ipar;

-(void)CRS_BlockDiagonal:(FEMSolution *)solution: (Matrix_t) B: (int)blocks;

// ILU(n) preconditioning
-(BOOL)CRS_IncompleteLU:(FEMSolution *)solution: (int)ilun;
-(BOOL)CRS_ComplexIncompleteLU:(FEMSolution *)solution: (int)ilun;

// ILU(T) preconditioning
-(void)computeILUT:(FEMSolution *)solution: (int)n: (int)tol;
-(BOOL)CRS_ILUT:(FEMSolution *)solution: (int)tol;
-(void)computeComplexILUT:(FEMSolution *)solution: (int)n: (int)tol;
-(BOOL)CRS_ComplexILUT:(FEMSolution *)solution: (int)tol;

// LU Solve
-(void)CRS_LUSolve:(int)n: (FEMSolution *)solution: (double *)b;
-(void)CRS_ComplexLUSolve:(int)n: (FEMSolution *)solution: (double complex *)b;
-(BOOL)CRS_LUPrecondition: (FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar;
-(BOOL)CRS_ComplexLUPrecondition: (FEMSolution *)solution: (double complex *)u: (double complex *)v: (int *)ipar;

// Matrix-vector product
-(void)CRS_MatrixVectorProd:(FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar;
-(void)CRS_ComplexMatrixVectorProd:(FEMSolution *)solution: (double complex *)u: (double complex *)v: (int *)ipar;

// Dummy method when preconditioning is not needed
-(void)CRS_pcond_dummy:(FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar;



@end

@implementation FEMPrecondition

#pragma mark Private methods

-(void)CRS_DiagPrecondition:(FEMSolution *)solution :(double *)u :(double *)v :(int *)ipar {
/*******************************************************************************************
                                                                             
    Description: Diagonal preconditioning of a CRS format matrix. Matrix is accessed from   
    the FEMSolution class.
 
    Arguments:
       
        double u, v;
 
        int ipar    -> Input stucture holding info from the HUTIter iterative solver
 
*******************************************************************************************/
    int i, j, n;
    
    int *range1;
    double *range2;
    
    n = [solution matrixNumberOfRows];  
    
    if ([solution matrixOrdered] == 0) {
        for (i=0; i<n; i++) {
            range1 = intvec(0, ([solution matrixRHS:i+1]-[solution matrixRHS:i])-1);
            range2 = doublevec(0, ([solution matrixRHS:i+1]-[solution matrixRHS:i])-1);
            for (j=[solution matrixRHS:i]; j<=[solution matrixRHS:i+1]-1; j++) {
                range1[j] = [solution matrixCols:j];
                range2[j] = [solution matrixValues:j];
            }
            sort([solution matrixRHS:i+1]-[solution matrixRHS:i], range1-1, range2-1);
            for (j=[solution matrixRHS:i]; j<=[solution matrixRHS:i+1]-1; j++) {
                [solution setMatrixCols:j :range1[j]];
                [solution setMatrixValues:j :range2[j]];
            }
            free_ivector(range1, 0, ([solution matrixRHS:i+1]-[solution matrixRHS:i])-1);
            free_dvector(range2, 0, ([solution matrixRHS:i+1]-[solution matrixRHS:i])-1);
        }
        for (i=0; i<n; i++) {
            for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                if ([solution matrixCols:j] == i) {
                    [solution setMatrixDiag:i :j];
                    break;
                }
            }
        }
        [solution setMatrixOrdered:1];
    }
    
    for (i=0; i<n; i++) {
        if (fabs([solution matrixValues:[solution matrixDiag:i]]) > AEPS) {
            u[i] = v[i] / [solution matrixValues:[solution matrixDiag:i]];
        } else {
            u[i] = v[i];
        }
    }
    
}

-(void)CRS_ComplexDiagPrecondition:(FEMSolution *)solution: (double complex *)u: (double complex *)v: (int *)ipar {
/*******************************************************************************************
 
    Description: Diagonal preconditioning of a CRS format matrix. Matrix is accessed from   
    the FEMSolution class.
 
    Arguments:
 
    double complex u, v;
 
    int ipar    -> Input stucture holding info from the HUTIter iterative solver
 
 *******************************************************************************************/
    
    int i, j, n;
    
    double complex A;
    
    int *range1;
    double *range2;
    
    n = [solution matrixNumberOfRows]; 
    
    if ([solution matrixOrdered] == 0) {
        for (i=0; i<n; i++) {
            range1 = intvec(0, ([solution matrixRHS:i+1]-[solution matrixRHS:i])-1);
            range2 = doublevec(0, ([solution matrixRHS:i+1]-[solution matrixRHS:i])-1);
            for (j=[solution matrixRHS:i]; j<=[solution matrixRHS:i+1]-1; j++) {
                range1[j] = [solution matrixCols:j];
                range2[j] = [solution matrixValues:j];
            }
            sort([solution matrixRHS:i+1]-[solution matrixRHS:i], range1-1, range2-1);
            for (j=[solution matrixRHS:i]; j<=[solution matrixRHS:i+1]-1; j++) {
                [solution setMatrixCols:j :range1[j]];
                [solution setMatrixValues:j :range2[j]];
            }
            free_ivector(range1, 0, ([solution matrixRHS:i+1]-[solution matrixRHS:i])-1);
            free_dvector(range2, 0, ([solution matrixRHS:i+1]-[solution matrixRHS:i])-1);
        }
        for (i=0; i<n; i++) {
            for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                if ([solution matrixCols:j] == i) {
                    [solution setMatrixDiag:i :j];
                    break;
                }
            }
        }
        [solution setMatrixOrdered:1];
    }
    
    for (i=0; n/2; i++) {
        A = [solution matrixValues:[solution matrixDiag:2*i-1]] + (-[solution matrixValues:[solution matrixDiag:2*i-1]+1] * I);
        u[i] = v[i] / A;
    }
    
}

-(void)CRS_BlockDiagonal:(FEMSolution *)solution: (Matrix_t) B: (int)blocks {
/*******************************************************************************************
    Description: 
        Pics the block diagonal entries form matrix Solution.matrix to build matrix B.
 
 *******************************************************************************************/
    
    int n;
    int i, k, l, kb;
    
    n = [solution matrixNumberOfRows];
    B.NumberOfRows = n;
    
    kb = 0;
    for (i=0; i<n; i++) {
        for (k=[solution matrixRows:i]; k<=[solution matrixRows:i+1]-1; k++) {
            l = [solution matrixCols:k];
            if ((i % blocks) == (l % blocks)) kb = kb + 1;
        }
    }
    
    // From here we need that B.Rows, B.Cols, B.Values and B.Diag have allocation in memory.
    // This should be taken care by the caller of this method.
    
    kb = 0;
    for (i=0; i<n; i++) {
        B.Rows[i] = kb;
        for (k=[solution matrixRows:i]; k<=[solution matrixRows:i+1]-1; k++) {
            l = [solution matrixCols:k];
            if ((i % blocks) == (l % blocks)) {
                B.Values[kb] = [solution matrixValues:k];
                B.Cols[kb] = [solution matrixCols:k];
                if (B.Cols[kb] == i) B.Diag[i] = kb;
                kb = kb + 1;
            }
        }
    }
    B.Rows[(n-1)+1] = kb;
    
}

-(BOOL)CRS_IncompleteLU:(FEMSolution *)solution: (int)ilun {
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
    int *C;                      //Boolean
    
    errorfunct("CRS_IncompleteLU", "ILU (Real). Starting factorization with order:", ilun);
    t = cputime();
    
    n = [solution matrixNumberOfRows];
    
    if ([solution isAssociatedMatrixILUValues] == NO) {
        
        [solution initializeILU:ilun];
    }
    
    // Allocate space for storing one full row
    C = intvec(0, n-1);
    S = doublevec(0, n-1);
    for (i=0; i<n; i++) {
        C[i] = 0;
        S[i] = 0.0;
    }
    
    // The factorization row by row
    for (i=0; i<n; i++) {
        
        // Convert current row to full form for speed, 
        // only flagging the nonzeros entries
        for (k=[solution matrixRows:i]; k<=[solution matrixRows:i+1]-1; k++) {
            S[[solution matrixCols:k]] = [solution matrixValues:k];
        }
        for (k=[solution matrixILURows:i]; k<=[solution matrixILURows:i+1]-1; k++) {
            C[[solution matrixILUCols:k]] = 1;
        }
        
        // This is the factorization part for the current row
        for (k=[solution matrixILUCols:[solution matrixILURows:i]]; k<=i-1; k++) {
            if (C[k] == 1) {
                if (fabs([solution matrixILUValues:[solution matrixILUDiag:k]]) > AEPS) S[k] = S[k] / [solution matrixILUValues:[solution matrixILUDiag:k]];
                
                for (l=[solution matrixILUDiag:k]; l<=[solution matrixILURows:k+1]-1; l++) {
                    j = [solution matrixILUCols:l];
                    if (C[j] == 1) S[j] = S[j] - S[k] * [solution matrixILUValues:l];
                }
            }
        }
        
        // Convert the row back to CRS format
        for (k=[solution matrixILURows:i]; k<=[solution matrixILURows:i+1]-1; k++) {
            if (C[[solution matrixILUCols:k]] == 1) {
                [solution setMatrixILUValues:k :S[[solution matrixILUCols:k]]];
                S[[solution matrixILUCols:k]] = 0.0;
                C[[solution matrixILUCols:k]] = 0;
            }
        }
    }
    
    free_ivector(C, 0, n-1);
    free_dvector(S, 0, n-1);
    
    // Prescale the diagonal for the LU solver
    for (i=0; i<n; i++) {
        if (fabs([solution matrixILUValues:[solution matrixILUDiag:i]]) < AEPS) {
            [solution setMatrixILUValues:[solution matrixILUDiag:i] :1.0];
        } else {
            [solution setMatrixILUValues:[solution matrixILUDiag:i] :(1.0 / [solution matrixILUValues:[solution matrixILUDiag:i]])];
        }
    }
    
    warnfunct("CRS_IncompleteLU", "ILU (Real), NOF nonzeros", [solution matrixILURows:(n+1)-1]);
    warnfunct("CRS_IncompleteLU", "ILU (Real), Filling (%):", floor([solution matrixILURows:(n+1)-1]) * (100.0 / [solution matrixRows:(n+1)-1]));
    warnfunct("CRS_IncompleteLU", "ILU (Real), Factorization ready at (s):", cputime() - t);
    
    return YES;
    
}

-(BOOL)CRS_ComplexIncompleteLU:(FEMSolution *)solution: (int)ilun {
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
    int *C;                      //Boolean
    
    errorfunct("CRS_IncompleteLU", "ILU (Complex), Starting factorization with order:", ilun);
    t = cputime();
    
    n = [solution matrixNumberOfRows];
    
    if ([solution isAssociatedMatrixCILUValues] == NO) {
        
        [solution initializeCILU:ilun];
    }
    
    // Allocate space for storing one full row
    C = intvec(0, n/2-1);
    S = cdoublevec(0, n/2-1);
    for (i=0; i<n; i++) {
        C[i] = 0;
        S[i] = 0.0;
    }
    
    // The factorization row by row
    for (i=0; i<n/2; i++) {
        
        // Convert current row to full form for speed, 
        // only flagging the nonzeros entries
        for (k=[solution matrixILURows:i]; k<=[solution matrixILURows:i+1]-1; k++) {
            C[[solution matrixILUCols:k]] = 1;
        }
        for (k=[solution matrixRows:2*i-1]; k<=[solution matrixRows:2*i]-1; k+=2) {
            S[([solution matrixCols:k]+1)/2] = [solution matrixValues:k] + (-[solution matrixValues:k+1] * I);
        }
        
        // This is the factorization part for the current row
        for (k=[solution matrixILUCols:[solution matrixILURows:i]]; k<=i-1; k++) {
            if (C[k] == 1) {
                if (fabs([solution matrixCILUValues:[solution matrixILUDiag:k]]) > AEPS) S[k] = S[k] / [solution matrixCILUValues:[solution matrixILUDiag:k]];
                
                for (l=[solution matrixILUDiag:k]+1; l<=[solution matrixILURows:k+1]-1; l++) {
                    j = [solution matrixILUCols:l];
                    if (C[j] == 1) S[j] = S[j] - S[k] * [solution matrixCILUValues:l];
                }
            }
        }
        
        // Convert the row back to CRS format
        for (k=[solution matrixILURows:i]; k<=[solution matrixILURows:i+1]-1; k++) {
            if (C[[solution matrixILUCols:k]] == 1) {
                [solution setMatrixCILUValues:k :S[[solution matrixILUCols:k]]];
                S[[solution matrixILUCols:k]] = 0.0;
                C[[solution matrixILUCols:k]] = 0;
            }
        }
    }
    
    free_ivector(C, 0, n/2-1);
    free_cdvector(S, 0, n/2-1);
    
    // Prescale the diagonal for the LU solver
    for (i=0; i<n/2; i++) {
        if (fabs([solution matrixCILUValues:[solution matrixILUDiag:i]]) < AEPS) {
            [solution setMatrixCILUValues:[solution matrixILUDiag:i] :1.0];
        } else {
            [solution setMatrixCILUValues:[solution matrixILUDiag:i] :1.0 / [solution matrixCILUValues:[solution matrixILUDiag:i]]];
        }
    }
    
    warnfunct("CRS_IncompleteLU", "ILU (Complex), NOF nonzeros", [solution matrixILURows:(n/2+1)-1]);
    warnfunct("CRS_IncompleteLU", "ILU (Complex), Filling (%):", floor([solution matrixILURows:(n/2+1)-1]) * (400.0 / [solution matrixRows:(n+1)-1]));
    warnfunct("CRS_IncompleteLU", "ILU (Complex), Factorization ready at (s):", cputime() - t);
    
    return YES;
    
}

-(void)computeILUT:(FEMSolution *)solution :(int)n :(int)tol {
    
    int i, j, k, l, rowMin, rowMax;
    int *C;                         //Boolean
    const double WORKN = 128;
    double norma, cptime, ttime, t;
    double *S;
    
    ttime = cputime();
    cptime = 0.0;
    
    [solution allocateMatrixILURows:n+1];
    [solution allocateMatrixILUDiag:n];
    if ([solution isAssociatedMatrixILURows] == NO || [solution isAssociatedMatrixILUDiag] == NO) {
        errorfunct("computeILUT", "Memory allocation error.");
    }
    [solution setMatrixSizeOfILURows:n+1];
    [solution setMatrixSizeOfILUDiag:n];
    
    [solution allocateMatrixILUCols:WORKN*n];
    [solution allocateMatrixILUValues:WORKN*n];
    if ([solution isAssociatedMatrixILUCols] == NO || [solution isAssociatedMatrixILUValues] == NO) {
        errorfunct("computeILUT", "Memory allocation error.");
    }
    [solution setMatrixSizeOfILUCols:WORKN*n];
    [solution setMatrixSizeOfILUValues:WORKN*n];
    
    // The factorization row by row
    [solution setMatrixILURows:0 :0];
    C = intvec(0, n);
    S = doublevec(0, n);
    for (i=0; i<n; i++) {
        S[i] = 0.0;
        C[i] = 0;
    }
    
    for (i=0; i<n; i++) {
        
        // Convert the current row to full from for speed, 
        // only flagging the nonzero entries
        for (k=[solution matrixRows:i]; k<=[solution matrixRows:i+1]-1; k++) {
            C[[solution matrixCols:k]] = 1;
            S[[solution matrixCols:k]] = [solution matrixValues:k];
        }
        
        // Check for bandwidth for speed, bandwidth optimization
        // helps a lot.
        rowMin = [solution matrixCols:[solution matrixRows:i]];
        rowMax = [solution matrixCols:[solution matrixRows:i+1]-1];
        
        // Here is the factorization part for the current row;
        for (k=rowMin; k<=i-1; k++) {
            if (C[k] == 1) {
                if (fabs([solution matrixILUValues:[solution matrixILUDiag:k]]) > AEPS) S[k] = S[k] / [solution matrixILUValues:[solution matrixILUDiag:k]];
                
                for (l=[solution matrixILUDiag:k]+1; l<=[solution matrixILURows:k+1]-1; l++) {
                    j = [solution matrixILUCols:l];
                    if (C[j] == 0) {
                        C[j] = 1;
                        rowMax = max(rowMax, j);
                    }
                    S[j] = S[j] - S[k] * [solution matrixILUValues:l];
                }
            }
        }
        
        // This is the ILUT part, drop element ILU(i,j), if
        // ABS(ILU(i,j)) <= NORM(A(i.:))*tol
        norma = 0.0;
        for (k=[solution matrixRows:i]; k<=[solution matrixRows:i+1]-1; k++) {
            norma = norma + pow(fabs([solution matrixValues:k]), 2.0);
        }
        norma = sqrt(norma);
        
        j = [solution matrixILURows:i]-1;
        for (k=rowMin; k<=rowMax; k++) {
            if (C[k] == 1) {
                if (fabs(S[k]) >= tol*norma || k == i) {
                    j = j + 1;
                    [solution setMatrixILUCols:j :k];
                    [solution setMatrixILUValues:j :S[k]];
                    if (k == i) [solution setMatrixILUDiag:i :j];
                }
                S[k] = 0.0;
                C[k] = 0;
            }
        }
        [solution setMatrixILURows:i+1 :j+1];
        
        // Preparations for the next row
        if (i < n-1) {
            
            // Check if still enough workspace
            if ([solution matrixSizeOFILUCols] < [solution matrixILURows:i+1]+n) {
                
                t = cputime();
                [solution ilutWorkspaceCheck:i :n];
                cptime = cptime + (cputime() - t);
                
            }
        }
    }
    
    // Prescale the diagonal for the LU solver
    for (i=0; i<n; i++) {
        
        if (fabs([solution matrixILUValues:[solution matrixILUDiag:i]]) < AEPS) {
            [solution setMatrixILUValues:[solution matrixILUDiag:i] :1.0];
        } else {
            
            [solution setMatrixILUValues:[solution matrixILUDiag:i] :(1.0 / [solution matrixILUValues:[solution matrixILUDiag:i]])];
        }
    }
    
}

-(BOOL)CRS_ILUT:(FEMSolution *)solution: (int)tol {
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
    
    warnfunct("CRS_ILUT", "Starting factorization:");
    t = cputime();
    
    n = [solution matrixNumberOfRows];
    
    if ([solution isAssociatedMatrixILUValues] == YES) {
        [solution freeMatrixILURows:[solution matrixSizeOfILURows]];
        [solution freeMatrixILUDiag:[solution matrixSizeOfILUDiag]];
        [solution freeMatrixILUCols:[solution matrixSizeOFILUCols]];
        [solution freeMatrixILUValues:[solution matrixSizeOfILUValues]];
    }
    
    [self computeILUT:solution :n :tol];
    
    warnfunct("CRS_ILUT", "ILU(T) (Real), NOF nonzeros", [solution matrixILURows:(n+1)-1]);
    warnfunct("CRS_ILUT", "ILU(T) (Real), Filling (%)", floor([solution matrixILURows:(n+1)-1]) * (100.0 / [solution matrixRows:(n+1)-1]));
    warnfunct("CRS_ILUT", "ILU(T) (Real), Factorization ready at (s):", cputime() - t);
    
    return YES;
    
}

-(void)computeComplexILUT:(FEMSolution *)solution: (int)n: (int)tol {
    
    int i, j, k, l, rowMin, rowMax;
    int *C;                         //Boolean
    const double WORKN = 128;
    double norma;
    double complex *S;
    
    [solution allocateMatrixILURows:n+1];
    [solution allocateMatrixILUDiag:n];
    if ([solution isAssociatedMatrixILURows] == NO || [solution isAssociatedMatrixILUDiag] == NO) {
        errorfunct("computeILUT", "Memory allocation error.");
    }
    [solution setMatrixSizeOfILURows:n+1];
    [solution setMatrixSizeOfILUDiag:n];
    
    [solution allocateMatrixILUCols:WORKN*n];
    [solution allocateMatrixCILUValues:WORKN*n];
    if ([solution isAssociatedMatrixILUCols] == NO || [solution isAssociatedMatrixCILUValues] == NO) {
        errorfunct("computeILUT", "Memory allocation error.");
    }
    [solution setMatrixSizeOfILUCols:WORKN*n];
    [solution setMatrixSizeOfCILUValues:WORKN*n];
    
    // The factorization row by row
    [solution setMatrixILURows:0 :0];
    C = intvec(0, n);
    S = cdoublevec(0, n);
    for (i=0; i<n; i++) {
        S[i] = 0.0 + 0.0 * I;
        C[i] = 0;
    }
    
    for (i=0; i<n; i++) {
        
        // Convert the current row to full from for speed, 
        // only flagging the nonzero entries
        for (k=[solution matrixRows:2*i-1]; k<=[solution matrixRows:2*i]-1; k+=2) {
            C[([solution matrixCols:k]+1) / 2] = 1;
            S[([solution matrixCols:k]+1) / 2] = [solution matrixValues:k] + (-[solution matrixValues:k+1] * I);
        }
        
        // Check for bandwidth for speed, bandwidth optimization
        // helps a lot.
        rowMin = ([solution matrixCols:[solution matrixRows:2*i-1]] + 1) / 2;
        rowMax = ([solution matrixCols:[solution matrixRows:2*i]-1] + 1) / 2;
        
        // Here is the factorization part for the current row;
        for (k=rowMin; k<=i-1; k++) {
            if (C[k] == 1) {
                if (fabs([solution matrixCILUValues:[solution matrixILUDiag:k]]) > AEPS) S[k] = S[k] / [solution matrixCILUValues:[solution matrixILUDiag:k]];
                
                for (l=[solution matrixILUDiag:k]+1; l<=[solution matrixILURows:k+1]-1; l++) {
                    j = [solution matrixILUCols:l];
                    if (C[j] == 0) {
                        C[j] = 1;
                        rowMax = max(rowMax, j);
                    }
                    S[j] = S[j] - S[k] * [solution matrixCILUValues:l];
                }
            }
        }
        
        // This is the ILUT part, drop element ILU(i,j), if
        // ABS(ILU(i,j)) <= NORM(A(i.:))*tol
        norma = 0.0;
        for (k=[solution matrixRows:2*i-1]; k<=[solution matrixRows:2*i]-1; k++) {
            norma = norma + pow([solution matrixValues:k], 2.0) + pow([solution matrixValues:k+1], 2.0);
        }
        norma = sqrt(norma);
        
        j = [solution matrixILURows:i]-1;
        for (k=rowMin; k<=rowMax; k++) {
            if (C[k] == 1) {
                if (fabs(S[k]) >= tol*norma || k == i) {
                    j = j + 1;
                    [solution setMatrixILUCols:j :k];
                    [solution setMatrixCILUValues:j :S[k]];
                    if (k == i) [solution setMatrixILUDiag:i :j];
                }
                S[k] = 0.0 + 0.0 * I;
                C[k] = 0;
            }
        }
        [solution setMatrixILURows:i+1 :j+1];
        
        // Preparations for the next row
        if (i < n-1) {
            
            // Check if still enough workspace
            if ([solution matrixSizeOFILUCols] < [solution matrixILURows:i+1]+n) {
                
                [solution ilutComplexWorkspaceCheck:i :n];
            }
        }
    }
    
    // Prescale the diagonal for the LU solver
    for (i=0; i<n; i++) {
        
        if (fabs([solution matrixCILUValues:[solution matrixILUDiag:i]]) < AEPS) {
            [solution setMatrixCILUValues:[solution matrixILUDiag:i] :1.0];
        } else {
            
            [solution setMatrixCILUValues:[solution matrixILUDiag:i] :(1.0 / [solution matrixCILUValues:[solution matrixILUDiag:i]])];
        }
    }
    
}

-(BOOL)CRS_ComplexILUT:(FEMSolution *)solution: (int)tol {
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
    
    warnfunct("CRS_ComplexILUT", "Starting factorization:");
    t = cputime();
    
    n = [solution matrixNumberOfRows] / 2;
    
    if ([solution isAssociatedMatrixCILUValues] == YES) {
        [solution freeMatrixILURows:[solution matrixSizeOfILURows]];
        [solution freeMatrixILUDiag:[solution matrixSizeOfILUDiag]];
        [solution freeMatrixILUCols:[solution matrixSizeOFILUCols]];
        [solution freeMatrixCILUValues:[solution matrixSizeOfCILUValues]];
    }
    
    [self computeComplexILUT:solution :n :tol];
    
    warnfunct("CRS_ComplexILUT", "ILU(T) (Complex), NOF nonzeros", [solution matrixILURows:(n+1)-1]);
    warnfunct("CRS_ComplexILUT", "ILU(T) (Complex), Filling (%)", floor([solution matrixILURows:(n+1)-1]) * (400.0 / [solution matrixRows:(2*n+1)-1]));
    warnfunct("CRS_ComplexILUT", "ILU(T) (Complex), Factorization ready at (s):", cputime() - t);
    
    return YES;
    
}

-(void)CRS_LUSolve:(int)n: (FEMSolution *)solution: (double *)b {
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
    
    // If no ILU provided do diagonal solve
    if ([solution isAssociatedMatrixILUValues] == NO) {
        for (i=0; i<n; i++) {
            b[i] = b[i] / [solution matrixValues:[solution matrixDiag:i]];
        }
        return;
    }
    
    // Forward substitute (solve z from Lz = b)
    for (i=0; i<n; i++) {
        s = b[i];
        for (j=[solution matrixILURows:i]; j<=[solution matrixILUDiag:i]-1; j++) {
            s = s - [solution matrixILUValues:j] * b[[solution matrixILUCols:j]];
        }
        b[i] = s;
    }
    
    // Backward substitute (solve x from UDx = z)
    for (i=n-1; i>=0; i--) {
        s = b[i];
        for (j=[solution matrixILUDiag:i]+1; j<=[solution matrixILURows:i+1]-1; j++) {
            s = s - [solution matrixILUValues:j] * b[[solution matrixILUCols:j]];
        }
        b[i] = [solution matrixILUValues:[solution matrixILUDiag:i]] * s;
    }
    
}

-(void)CRS_ComplexLUSolve:(int)n: (FEMSolution *)solution: (double complex *)b {
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
    
    // If no ILU provided do diagonal solve
    if ([solution isAssociatedMatrixCILUValues] == NO) {
        for (i=0; i<n/2; i++) {
            
            x = [solution matrixValues:[solution matrixDiag:2*i-1]] + (-[solution matrixValues:[solution matrixDiag:2*i-1]+1] * I);
            b[i] = b[i] / x;
        }
        return;
    }
    
    // Forward substitute (solve z from Lz = b)
    for (i=0; i<n; i++) {
        s = b[i];
        for (j=[solution matrixILURows:i]; j<=[solution matrixILUDiag:i]-1; j++) {
            s = s - [solution matrixCILUValues:j] * b[[solution matrixILUCols:j]];
        }
        b[i] = s;
    }
    
    // Backward substitute (solve x from UDx = z)
    for (i=n-1; i>=0; i--) {
        s = b[i];
        for (j=[solution matrixILUDiag:i]+1; j<=[solution matrixILURows:i+1]-1; j++) {
            s = s - [solution matrixCILUValues:j] * b[[solution matrixILUCols:j]];
        }
        b[i] = [solution matrixCILUValues:[solution matrixILUDiag:i]] * s;
    }    
    
}

-(BOOL)CRS_LUPrecondition: (FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar {
/*******************************************************************************************
 
    Description: 
        Incomplete factorization preconditioner solver for a CRS format matrix.
        Matrix is accessed from the solution class. Real matrix version.
 
 Arguments:

    FEMSolution *solution  -> Class holding input matrix.
 
    double *u 
    
    double *v
 
    int ipar               -> Structure holding info from the HUTIter iterative solver.

*******************************************************************************************/

    int i;
    
    for (i=0; i<ipar[2]; i++) {
        u[i] = v[i];
    }
    
    [self CRS_LUSolve:ipar[2] :solution :u];
    
    return YES;
    
}

-(BOOL)CRS_ComplexLUPrecondition: (FEMSolution *)solution: (double complex *)u: (double complex *)v: (int *)ipar {
/*****************************************************************************************************************
 
    Description: 
        Incomplete factorization preconditioner solver for a CRS format matrix.
        Matrix is accessed from the solution class. Real matrix version.
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix.
 
        double complex *u 
 
        double complex *v
 
        int ipar               -> Structure holding info from the HUTIter iterative solver.
 
 ****************************************************************************************************************/
    
    int i;
    
    for (i=0; i<ipar[2]; i++) {
        u[i] = v[i];
    }
    
    [self CRS_ComplexLUSolve:ipar[2] :solution :u];

    return YES;
    
}

-(void)CRS_MatrixVectorProd:(FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar {
/*******************************************************************************************
 
    Description: 
        Matrix vector product (v = Au) for a matrix given in CRS format. The matrix 
        is accessed from the solution class. Real version.
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix.
 
        double *u 
 
        double *v
 
        int ipar               -> Structure holding info from the HUTIter iterative solver.
 
*******************************************************************************************/
    
    int i, j, n;
    double s;
    
    n = [solution matrixNumberOfRows];
    
    if (ipar[5] == 0) {
        
        for (i=0; i<n; i++) {
            s = 0.0;
            for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                s = s + [solution matrixValues:j] * u[[solution matrixCols:j]];
            }
            v[i] = s;
        }
    } else {
        for (i=0; i<n; i++) {
            v[i] = 0.0;
        }
        for (i=0; i<n; i++) {
            s = u[i];
            for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                v[[solution matrixCols:j]] = v[[solution matrixCols:j]] + s * [solution matrixValues:j];
            }
        }
    }
}

-(void)CRS_ComplexMatrixVectorProd:(FEMSolution *)solution: (double complex *)u: (double complex *)v: (int *)ipar {
/*******************************************************************************************
 
    Description: 
        Matrix vector product (v = Au) for a matrix given in CRS format. The matrix 
        is accessed from the solution class. Complex version.
 
    Arguments:
 
    FEMSolution *solution  -> Class holding input matrix.
 
    double complex *u 
 
    double complex *v
 
    int ipar               -> Structure holding info from the HUTIter iterative solver.
 
*******************************************************************************************/
    
    int i, j, n;
    double complex s, rsum;
    
    n = ipar[2];
    
    if (ipar[5] == 0) {
        
        for (i=0; i<n; i++) {
            rsum = 0.0 + 0.0 * I;
            for (j=[solution matrixRows:2*i-1]; j<=[solution matrixRows:2*i]-1; i+=2) {
                s = [solution matrixValues:j] + (-[solution matrixValues:j+1] * I);
                rsum = rsum + s * u[([solution matrixCols:j]+1)/2];
            }
            v[i] = rsum;
        }
    } else {
        for (i=0; i<n; i++) {
            v[i] = 0.0 + 0.0 * I; 
        }
        for (i=0; i<n; i++) {
            rsum = u[i];
            for (j=[solution matrixRows:2*i-1]; j<=[solution matrixRows:2*i]-1; j+=2) {
                s = [solution matrixValues:j] + (-[solution matrixValues:j+1] * I);
                v[([solution matrixCols:j]+1)/2] = v[([solution matrixCols:j]+1)/2] + s * rsum;
            }
        }
    }
}

-(void)CRS_pcond_dummy:(FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar {
    
    int i;
    
    for (i=0; i<ipar[2]; i++) {
        u[i] = v[i];
    }
}

#pragma mark Public methods

-(void)CRS_MatrixVectorMultiply:(FEMSolution *)solution: (double *)u: (double *)v {
    /*******************************************************************************************
     
     Description: 
     Matrix vector product (v = Au) for a matrix given in CRS format. The matrix 
     is accessed from the solution class. Real version.
     
     Arguments:
     
     FEMSolution *solution  -> Class holding input matrix.
     
     double *u 
     
     double *v
     
     *******************************************************************************************/
    
    int i, j, n;
    double rsum;
    
    n = [solution matrixNumberOfRows];
    
    for (i=0; i<n; i++) {
        rsum = 0.0;
        for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
            rsum = rsum + u[[solution matrixCols:j]] * [solution matrixValues:j];
        }
        v[i] = rsum;
    }
    
}

-(void)CRS_ComplexMatrixVectorMultiply:(FEMSolution *)solution: (double complex *)u: (double complex *)v {
    /*******************************************************************************************
     
     Description: 
     Matrix vector product (v = Au) for a matrix given in CRS format. The matrix 
     is accessed from the solution class. Complex version.
     
     Arguments:
     
     FEMSolution *solution  -> Class holding input matrix.
     
     double *u 
     
     double *v
     
     *******************************************************************************************/
    
    int i, j, n;
    double complex s, rsum;
    
    n = [solution matrixNumberOfRows] / 2;
    
    for (i=0; i<n; i++) {
        rsum = 0.0 + 0.0 * I;
        for (j=[solution matrixRows:2*i-1]; j<=[solution matrixRows:2*i]-1; j+=2) {
            s = [solution matrixValues:j] + (-[solution matrixValues:j+1] * I);
            rsum = rsum + s * u[([solution matrixCols:j]+1)/2];
        }
        v[i] = rsum;
    }
    
}




@end
