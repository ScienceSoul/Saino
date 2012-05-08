//
//  FEMMatrixBand.m
//  Saino
//
//  Created by Hakime Seddik on 23/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMatrixBand.h"

@interface FEMMatrixBand ()

@end


@implementation FEMMatrixBand

-(void)zeroRowInGlobal:(FEMSolution *)solution: (int)n {
    
    int j;
    
    if ([solution matrixFormat] == MATRIX_BAND) {
        for (j=max(0, n-[solution matrixSubband]); j<min([solution matrixNumberOfRows], n+[solution matrixSubband]); j++) {
            [solution setMatrixValues:(j*(3*[solution matrixSubband]+1) + (n)-(j)+2*[solution matrixSubband]) :0.0];
        }
    } else {
        for (j=max(0, n-[solution matrixSubband]); j<=n; j++) {
            [solution setMatrixValues:(j*([solution matrixSubband]+1) + (n)-(j)) :0.0];
        }
    }
    
}

-(void)setMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value {
/************************************************************************************************
    Set a given value to an element of a Band format Matrix
 
    Arguments:
        FEMSolution *solution   ->  solution class containing the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
 
************************************************************************************************/
    
    if ([solution matrixFormat] == MATRIX_BAND) {
        [solution setMatrixValues:(j*(3*[solution matrixSubband]+1) + (i)-(j)+2*[solution matrixSubband]) :value];
    } else {
        if (j <= i) [solution setMatrixValues:(j*([solution matrixSubband]+1) + (i)-(j)) :value];
    }
    
}

-(void)addToMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value {
/************************************************************************************************
    Add a given value to an element of a Band format matrix
 
    Arguments:
        FEMSolution *solution   ->  solution class containing the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
 
 ************************************************************************************************/
    int k;
    
    if ([solution matrixFormat] == MATRIX_BAND) {
        k = (j*(3*[solution matrixSubband]+1) + (i)-(j)+2*[solution matrixSubband]);
        [solution setMatrixValues:k :[solution matrixValues:k]+value];
    } else {
        k = (j*([solution matrixSubband]+1) + (i)-(j));
        if (j <= i) [solution setMatrixValues:k :[solution matrixValues:k]+value];
    }
    
}

-(void)glueLocalMatrixInGlobal:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes {
/*************************************************************************************************************
    
    Add a set of values (i.e., element stiffness matrix) to a Band format matrix
 
    Arguments:
 
        FEMSolution *solution   -> Solution class holding the global matrix
        double **matrix         -> (n x dofs) x (n x dofs) matrix holding the values to be
                                   added to the Band format matrix
        int n                   -> number of nodes in element
        int dofs                -> number of degrees of freemdom for one node
        int **indexes           -> Maps element node number to global (or partition) node number
                                   (to matrix rows and cols if dofs = 1)
 
************************************************************************************************************/
    
    int i, j, k, l, ind, row, col;
    
    if ([solution matrixFormat] == MATRIX_BAND) {
        for (i=0; i<n; i++) {
            for (k=1; k<=dofs; k++) {
                row = dofs * (indexes[i]+1) - k;
                for (j=0; j<n; j++) {
                    for (l=1; l<=dofs; l++) {
                        col = dofs * (indexes[j]+1) - l;
                        ind = ( (col-1)*(3*[solution matrixSubband]+1) + row - col + 2*[solution matrixSubband]+1 );
                        [solution setMatrixValues:ind :matrix[dofs*(i+1)-k][dofs*(j+1)-l]];
                    }
                }
            
            }
        }
    } else {
        for (i=0; i<n; i++) {
            for (k=1; k<=dofs; k++) {
                row = dofs * (indexes[i]+1) - k;
                for (j=0; j<n; j++) {
                    for (l=1; l<=dofs; l++) {
                        col = dofs * (indexes[j]+1) - l;
                        if (col <= row) {
                            ind = ( (col-1)*([solution matrixSubband]+1) + row - col + 1 );
                            [solution setMatrixValues:ind :matrix[dofs*(i+1)-k][dofs*(j+1)-l]];
                        }
                    }
                }
            }
        }
    }
    
}

-(void)sBand_setDirichlet:(FEMSolution *)solution: (int)n: (double)value {
/*************************************************************************************************************
 
    Set value of unknown x(n) to given value for symmetric band matrix. This is done by replacing the
    equation of the unknown by x(n) = value (i.e., zeroing the row of the unknown in the matrix and
    setting diagonal to identity). Also the repective column is set to zero (except for the diagonal)
    to preserve symmetry, while also substituting the rhs by rhs(i) = rhs(i) - A(i,n) * value.
 
    Arguments:
 
    FEMSolution *solution   -> Solution class holding the global matrix
    int n                   -> ordered number of the unknown (i.e., matrix row and column number)
    double value            -> value for the unknown
 
 ************************************************************************************************************/
    
    int j;
    
    for (j=max(0, n-[solution matrixSubband]); j<n-1; j++) {
        [solution setMatrixRHS:j :[solution matrixRHS:j]-value*[solution matrixValues:j*([solution matrixSubband]+1)+(n)-(j)]];
        [solution setMatrixValues:j*([solution matrixSubband]+1)+(n)-(j) :0.0];
    }
    
    for (j=n+1; j<min(n+[solution matrixSubband], [solution matrixNumberOfRows]); j++) {
        [solution setMatrixRHS:j :[solution matrixRHS:j]-value*[solution matrixValues:n*([solution matrixSubband]+1)+(j)-(n)]];
        [solution setMatrixValues:n*([solution matrixSubband]+1)+(j)-(n) :0.0];
    }
    
    [solution setMatrixRHS:n :value];
    [solution setMatrixValues:n*([solution matrixSubband]+1)+(n)-(n) :1.0];
    
}

-(void)zeroRowInMatrix:(Matrix_t *)a: (int)n {
    
    int j;
    
    if (a->Format == MATRIX_BAND) {
        for (j=max(0, n-a->Subband); j<min(a->NumberOfRows, n+a->Subband); j++) {
            a->Values[(j*(3*a->Subband+1) + (n)-(j)+2*a->Subband)] = 0.0;
        }
    } else {
        for (j=max(0, n-a->Subband); j<=n; j++) {
            a->Values[(j*(a->Subband+1) + (n)-(j))] = 0.0;
        }
    }
    
}

-(void)setMatrixElementInMatrix:(Matrix_t *)a: (int)i: (int)j: (double)value {
/************************************************************************************************
    Set a given value to an element of a Band format Matrix
 
    Arguments:
    Matrix_t *a             ->  the matrix
    int i, j                ->  row and column numbers respectively of the matrix element
    double value            ->  value to be set
 
 ************************************************************************************************/
    
    if (a->Format == MATRIX_BAND) {
        a->Values[(j*(3*a->Subband+1) + (i)-(j)+2*a->Subband)] = value;
    } else {
        if (j <= i) a->Values[(j*(a->Subband+1) + (i)-(j))] = value;
    }
    
}

@end
