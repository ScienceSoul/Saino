//
//  FEMMatrixBand.m
//  Saino
//
//  Created by Hakime Seddik on 23/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMatrixBand.h"
#import "Utils.h"

@interface FEMMatrixBand ()

@end


@implementation FEMMatrixBand

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
    }
    
    return self;
}

/********************************************************************************************************************
 
    Create the structures required for a Band format matrix
 
    int rows            ->  Number of rows in the matrix
    int subBand         ->  max(abs(Col-Diag(Row))) of the matrix
    BOOL symmetric      ->  Symmetric or not  
    BOOL allocateValues ->  Should the values arrays be allocated?
 
 *********************************************************************************************************************/
-(FEMMatrix * __nonnull)createMatrixWithNumberOfRows:(int)rows subBand:(int)subBand symmetric:(BOOL)symmetric allocateValues:(BOOL)allocateValues {
    
    FEMMatrix *matrix;
    matrixArraysContainer *matContainers = NULL;
    
    matrix = [[FEMMatrix alloc] init];
    
    matrix.subband = subBand;
    matrix.numberOfRows = rows;
    
    matContainers = matrix.getContainers;
    
    if (allocateValues == YES) {
        if (symmetric == YES) {
            matContainers->Values = doublevec(0, ((matrix.subband+1)*rows)-1);
        } else {
            matContainers->Values = doublevec(0, ((3*matrix.subband+1)*rows)-1);
        }
    }
    
    return matrix;
}

-(void)zeroRowInGlobal:(FEMSolution * __nonnull)solution numberOfRows:(int)n {
    
    int j;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    if (solution.matrix.format == MATRIX_BAND) {
        for (j=max(0, n-solution.matrix.subband); j<=min(solution.matrix.numberOfRows-1, n+solution.matrix.subband); j++) {
            matContainers->Values[(j*(3*solution.matrix.subband+1) + (n)-(j)+2*solution.matrix.subband)] = 0.0;
        }
    } else {
        for (j=max(0, n-solution.matrix.subband); j<=n; j++) {
            matContainers->Values[(j*(solution.matrix.subband+1) + (n)-(j))] = 0.0;
        }
    }
}

/**********************************************************************************************************
    Set a given value to an element of a Band format Matrix
 
    Arguments:
       
        FEMSolution *solution   ->  solution class containing the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
**********************************************************************************************************/
-(void)setElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value {
    
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    if (solution.matrix.format == MATRIX_BAND) {
        matContainers->Values[(j*(3*solution.matrix.subband+1) + (i)-(j)+2*solution.matrix.subband)] = value;
    } else {
        if (j <= i) matContainers->Values[(j*(solution.matrix.subband+1) + (i)-(j))] = value;
    }
}


/***********************************************************************************************************
    Add a given value to an element of a Band format matrix
 
    Arguments:
 
        FEMSolution *solution   ->  solution class containing the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
**********************************************************************************************************/
-(void)addToElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value {
    
    int k;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    if (solution.matrix.format == MATRIX_BAND) {
        k = (j*(3*solution.matrix.subband+1) + (i)-(j)+2*solution.matrix.subband);
        matContainers->Values[k] = matContainers->Values[k] + value;
    } else {
        if (j <= i) {
            k = (j*(solution.matrix.subband+1) + (i)-(j));
            matContainers->Values[k] = matContainers->Values[k] + value;
        }
    }
}

/*****************************************************************************************************************************************
 
    Add a set of values (i.e., element stiffness matrix) to a Band format matrix
 
    Arguments:
 
        FEMSolution *solution   -> Solution class holding the global matrix
        double **matrix         -> (n x dofs) x (n x dofs) matrix holding the values to be
                                   added to the Band format matrix
        int n                   -> number of nodes in element
        int dofs                -> number of degrees of freemdom for one node
        int **indexes           -> Maps element node number to global (or partition) node number
                                   (to matrix rows and cols if dofs = 1)
*****************************************************************************************************************************************/
-(void)glueLocalMatrix:(double * __nonnull * __nonnull)localMatrix inGlobal:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dofs:(int)dofs indexes:(int * __nonnull)indexes {
    
    int i, j, k, l, ind, row, col;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    if (solution.matrix.format == MATRIX_BAND) {
        for (i=0; i<n; i++) {
            for (k=1; k<=dofs; k++) {
                row = dofs * (indexes[i]+1) - k;
                for (j=0; j<n; j++) {
                    for (l=1; l<=dofs; l++) {
                        col = dofs * (indexes[j]+1) - l;
                        ind = ( (col-1)*(3*solution.matrix.subband+1) + row - col + 2*solution.matrix.subband );
                        matContainers->Values[ind] = localMatrix[dofs*(i+1)-k][dofs*(j+1)-l];
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
                            ind = ( (col-1)*(solution.matrix.subband+1) + row - col + 1 );
                            matContainers->Values[ind] = localMatrix[dofs*(i+1)-k][dofs*(j+1)-l];
                        }
                    }
                }
            }
        }
    }
}

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
-(void)sBand_setDirichlet:(FEMSolution * __nonnull)solution orderedNumber:(int)n value:(double)value {
    
    int j;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    for (j=max(0, n-solution.matrix.subband); j<=n-1; j++) {
        matContainers->RHS[j] = matContainers->RHS[j]-value*matContainers->Values[j*(solution.matrix.subband+1)+(n)-(j)];
        matContainers->Values[j*(solution.matrix.subband+1)+(n)-(j)] = 0.0;
    }
    
    for (j=n+1; j<=min(n+solution.matrix.subband, solution.matrix.numberOfRows-1); j++) {
        matContainers->RHS[j] = matContainers->RHS[j]-value*matContainers->Values[n*(solution.matrix.subband+1)+(j)-(n)];
        matContainers->Values[n*(solution.matrix.subband+1)+(j)-(n)] = 0.0;
    }
    
    matContainers->RHS[n] = value;
    matContainers->Values[n*(solution.matrix.subband+1)+(n)-(n)] = 1.0;
}

/*******************************************************************************************
 
    Description:
        Matrix vector product (v = Au) for a matrix given in band format. The matrix
        is accessed from the solution class. Real version.
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix.
        double *u              -> Vector to multiply
        double *v              -> Result vector
*******************************************************************************************/
-(void)matrixVectorMultiplyInGlobal:(FEMSolution * __nonnull)solution vector:(double * __nonnull)u result:(double * __nonnull)v {
    
    int i, j, n;
    double s;
    matrixArraysContainer *matContainers = NULL;

    n = solution.matrix.numberOfRows;
    matContainers = solution.matrix.getContainers;
    
    if (solution.matrix.format == MATRIX_BAND) {
        for (i=0; i<n; i++) {
            s = 0.0;
            for (j=max(0, i-solution.matrix.subband); j<=min(n-1, i+solution.matrix.subband); j++) {
                s = s + u[j] * matContainers->Values[((j)*(3*solution.matrix.subband+1) + (i)-(j)+2*solution.matrix.subband)];
            }
            v[i] = s;
        }
    } else {
        for (i=0; i<n; i++) {
            s = 0.0;
            for (j=max(0, i-solution.matrix.subband); j<=i; j++) {
                s = s + u[j] * matContainers->Values[(j)*(solution.matrix.subband+1) + (i)-(j)];
            }
            
            for (j=i+1; j<=min(i+solution.matrix.subband, solution.matrix.numberOfRows-1); j++) {
                s = s + u[j] + matContainers->Values[(i)*(solution.matrix.subband+1) + (j)-(i)];
            }
            v[i] = s;
        }
    }
}

-(void)zeroRowInMatrix:(FEMMatrix * __nonnull)matrix numberOfRows:(int)n {
    
    int j;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    if (matrix.format == MATRIX_BAND) {
        for (j=max(0, n-matrix.subband); j<=min(matrix.numberOfRows-1, n+matrix.subband); j++) {
            matContainers->Values[(j*(3*matrix.subband+1) + (n)-(j)+2*matrix.subband)] = 0.0;
        }
    } else {
        for (j=max(0, n-matrix.subband); j<=n; j++) {
            matContainers->Values[(j*(matrix.subband+1) + (n)-(j))] = 0.0;
        }
    }
}

/************************************************************************************************
    Set a given value to an element of a Band format Matrix
 
    Arguments:
 
        Matrix_t *a             ->  the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
************************************************************************************************/
-(void)setElementInMatrix:(FEMMatrix * __nonnull)matrix row:(int)i col:(int)j value:(double)value {
    
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    if (matrix.format == MATRIX_BAND) {
        matContainers->Values[(j*(3*matrix.subband+1) + (i)-(j)+2*matrix.subband)] = value;
    } else {
        if (j <= i) matContainers->Values[(j*(matrix.subband+1) + (i)-(j))] = value;
    }
}

/*******************************************************************************************
 
    Matrix vector product (v = Au) for a matrix given in band format.
 
    Arguments:
 
        FEMMatrix *a  -> input matrix.
        double *u     -> Vector to multiply
        double *v     -> Result vector
*******************************************************************************************/
-(void)matrixVectorMultiply:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)u result:(double * __nonnull)v {
    
    int i, j, n;
    double s;
    matrixArraysContainer *matContainers = NULL;
    
    n = matrix.numberOfRows;
    matContainers = matrix.getContainers;
    
    if (matrix.format == MATRIX_BAND) {
        for (i=0; i<n; i++) {
            s = 0.0;
            for (j=max(0, i-matrix.subband); j<=min(n-1, i+matrix.subband); j++) {
                s = s + u[j] * matContainers->Values[((j)*(3*matrix.subband+1) + (i)-(j)+2*matrix.subband)];
            }
            v[i] = s;
        }
    } else {
        for (i=0; i<n; i++) {
            s = 0.0;
            for (j=max(0, i-matrix.subband); j<=i; j++) {
                s = s + u[j] * matContainers->Values[(j)*(matrix.subband+1) + (i)-(j)];
            }
            
            for (j=i+1; j<=min(i+matrix.subband, matrix.numberOfRows-1); j++) {
                s = s + u[j] + matContainers->Values[(i)*(matrix.subband+1) + (j)-(i)];
            }
            v[i] = s;
        }
    }
}

@end
