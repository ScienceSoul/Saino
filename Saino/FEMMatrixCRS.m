//
//  FEMMatrixCRS.m
//  Saino
//
//  Created by Hakime Seddik on 23/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMatrixCRS.h"

@interface FEMMatrixCRS ()

-(int)CRS_Search:(int)n: (int *)array: (int)value;

@end


@implementation FEMMatrixCRS

#pragma mark Private methods

-(int)CRS_Search:(int)n: (int *)array: (int)value {
    
    int lower, upper, lou, index;
    
    index = 0;
    upper = n;
    lower = 0;
    
    // Handle the special case returns -1
    if (upper < 0) return -1;
    
    while (1) {
        if (array[lower] == value) {
            index = lower;
            break;
        } else if (array[upper] == value) {
            index = upper;
            break;
        }
        
        if ( (upper-lower) > 0 ) {
            lou = (upper+lower) >> 1;
            if (array[lou] < value) {
                lower = lou;
            } else {
                upper = lou;
            }
        } else {
            break;
        }
    }
    
    return index;
}

#pragma mark Public methods

-(void)zeroRowInGlobal:(FEMSolution *)solution: (int)n {
    
    int i;
    matrixArraysContainer *matContainers;
    
    matContainers = solution.matrix.getContainers;
    
    
    for (i=matContainers->Rows[n]; i<=matContainers->Rows[n+1]-1; i++) {
        matContainers->Values[i] = 0.0;
    }
    
    if (matContainers->MassValues != NULL) {
        if (matContainers->sizeMassValues == matContainers->sizeValues) {
            for (i=matContainers->Rows[n]; i<=matContainers->Rows[n+1]-1; i++) {
                matContainers->MassValues[i] = 0.0;
            }
        }
    }
    
    if (matContainers->DampValues != NULL) {
        if (matContainers->sizeDampValues == matContainers->sizeValues) {
            for (i=matContainers->Rows[n]; i<=matContainers->Rows[n+1]-1; i++) {
                matContainers->DampValues[i] = 0.0;
            }
        }
    }
    
    matContainers = NULL;
}

-(void)sortInGlobal:(FEMSolution *)solution: (BOOL *)alsoValues {
/************************************************************************************************
    Sort columns to ascending order for rows of a CRS format matrix
 
    Arguments:
        FEMSolution *solution   ->  solution class containing the matrix
        BOOL *alsoValues        ->  whether values are sorted
 
************************************************************************************************/
    
    int i, j, k, n;
    int *buffer1;
    double *buffer2;
    BOOL sortValues;
    matrixArraysContainer *matContainers;
    
    matContainers = solution.matrix.getContainers;
    
    sortValues = NO;
    if (alsoValues != NULL) sortValues = *alsoValues;
    
    n = solution.matrix.numberOfRows;
    
    if (solution.matrix.isOrdered == NO) {
        if (sortValues == YES) {
            buffer1 = intvec(0, matContainers->sizeValues-1);
            buffer2 = doublevec(0, matContainers->sizeValues-1);
            memset( buffer1, 0, (matContainers->sizeValues*sizeof(buffer1)) );
            memset( buffer2, 0.0, (matContainers->sizeValues*sizeof(buffer2)) );
            for (i=0; i<n; i++) {
                k = 0;
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    buffer1[k] = matContainers->Cols[j];
                    buffer2[k] = matContainers->Values[j];
                    k++;
                }
                sort(matContainers->Rows[i+1]-matContainers->Rows[i], buffer1-1, buffer2-1);
                k = 0;
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    matContainers->Cols[j] = buffer1[k];
                    matContainers->Values[j] = buffer2[k];
                    k++;
                }
            }
            free_ivector(buffer1, 0, matContainers->sizeValues-1);
            free_dvector(buffer2, 0, matContainers->sizeValues-1);
            
        } else {
            buffer1 = intvec(0, matContainers->sizeValues-1);
            memset( buffer1, 0, (matContainers->sizeValues*sizeof(buffer1)) );
            for (i=0; i<n; i++) {
                k = 0;
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    buffer1[k] = matContainers->Cols[j];
                    k++;
                }
                sort(matContainers->Rows[i+1]-matContainers->Rows[i], buffer1-1);
                k = 0;
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    matContainers->Cols[j] = buffer1[k];
                    k++;
                }
            }
            free_ivector(buffer1, 0, matContainers->sizeValues-1);
        }
        
        if (matContainers->Diag != NULL) {
            for (i=0; i<n; i++) {
                for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    if (matContainers->Cols[j] == i) {
                        matContainers->Diag[i] = j;
                        break;
                    }
                }
            }
        }
        
        solution.matrix.ordered = YES;
    }
    
    matContainers = NULL;
    
}

-(void)setMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value {
/************************************************************************************************
    Set a given value to an element of a CRS format Matrix
 
    Arguments:
        FEMSolution *solution   ->  solution class containing the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
 
************************************************************************************************/
    
    int ii, jj, k;
    int *buffer;
    matrixArraysContainer *matContainers;
    
    matContainers = solution.matrix.getContainers;
    
    if (matContainers->Diag == NULL || i != j || solution.matrix.isOrdered == NO) {
        jj = matContainers->Rows[i+1]-matContainers->Rows[i];
        buffer = intvec( 0, jj );
        memset( buffer, 0.0, ((jj+1)*sizeof(buffer)) );
        for (ii=matContainers->Rows[i]; ii<=matContainers->Rows[i+1]-1; ii++) {
            buffer[ii] = matContainers->Cols[ii];
        }
        k = [self CRS_Search:matContainers->Rows[i+1]-matContainers->Rows[i] :buffer :j];
        if (k < 0) {
            warnfunct("CRS:setMatrixElement", "Trying to set value to non existent element:");
            printf("%d %d %f\n", i, j, value);
            return;
        }
        k = k + matContainers->Rows[i];
        free_ivector(buffer, 0, jj);
    } else {
        k = matContainers->Diag[i];
    }
    matContainers->Values[k] = value;

    matContainers = NULL;
}

-(void)addToMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value {
/************************************************************************************************
    Add a given value to an element of a CRS format Matrix
 
    Arguments:
        FEMSolution *solution   ->  solution class containing the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be added
 
************************************************************************************************/
    
    int k, ii, jj;
    int *buffer;
    matrixArraysContainer *matContainers;
    
    matContainers = solution.matrix.getContainers;
    
    if (matContainers->Diag == NULL || i != j || solution.matrix.isOrdered == NO) {
        jj = matContainers->Rows[i+1]-matContainers->Rows[i];
        buffer = intvec( 0, jj );
        memset( buffer, 0.0, ((jj+1)*sizeof(buffer)) );
        for (ii=matContainers->Rows[i]; ii<=matContainers->Rows[i+1]-1; ii++) {
            buffer[ii] = matContainers->Cols[ii];
        }
        k = [self CRS_Search:matContainers->Rows[i+1]-matContainers->Rows[i] :buffer :j];
        if (k < 0 && value != 0) warnfunct("addToMatrixElement", "Trying to add value to non existent element:");
        printf("%d %d %f\n", i, j, value);
        if (k < 0) return;
        k = k + matContainers->Rows[i];
        free_ivector(buffer, 0, jj);
    } else {
        k = matContainers->Diag[i];
    }
    matContainers->Values[k] = matContainers->Values[k]+value;
    
    matContainers = NULL;
}


-(void)glueLocalMatrixInGlobal:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes {
/*************************************************************************************************************
   
    Add a set of values (i.e., element stiffness matrix) to a CRS format matrix
 
    Arguments:
    
        FEMSolution *solution   -> Solution class holding the global matrix
        double **matrix         -> (n x dofs) x (n x dofs) matrix holding the values to be
                                   added to the CRS format matrix
        int n                   -> number of nodes in element
        int dofs                -> number of degrees of freemdom for one node
        int **indexes           -> Maps element node number to global (or partition) node number
                                   (to matrix rows and cols if dofs = 1)
 
************************************************************************************************************/
    
    int i, j, k, l, c, row, col;
    matrixArraysContainer *matContainers;
    
    matContainers = solution.matrix.getContainers;
    
    if (dofs == 1) {
        for (i=0; i<n; i++) {
            row = indexes[i];
            if (row < 0) continue;
            for (j=0; j<n; j++) {
                col = indexes[j];
                if (col < 0) continue;
                if (col >= row) {
                    for (c=matContainers->Diag[row]; c<=matContainers->Rows[row+1]-1; c++) {
                        if (matContainers->Cols[c] == col) {
                            matContainers->Values[c] = matContainers->Values[c] + matrix[i][j];
                            break;
                        }
                    }
                } else {
                    for (c=matContainers->Rows[row]; c<=matContainers->Diag[row]-1; c++) {
                        if (matContainers->Cols[c] == col) {
                            matContainers->Values[c] = matContainers->Values[c] + matrix[i][j];
                            break;
                        }
                    }
                }
            }
        }
    } else {
        for (i=0; i<n; i++) {
            for (k=1; k<=dofs; k++) {
                if (indexes[i] < 0) continue;
                row = dofs * (indexes[i]+1) - k;
                for (j=0; j<n; j++) {
                    for (l=1; l<=dofs; l++) {
                        if (indexes[j] < 0) continue;
                        col = dofs * (indexes[j]+1) - l;
                        if (col >= row) {
                            for (c=matContainers->Diag[row]; c<=matContainers->Rows[row+1]-1; c++) {
                                if (matContainers->Cols[c] == col) {
                                    matContainers->Values[c] = matContainers->Values[c] + matrix[dofs*(i+1)-k][dofs*(j+1)-l];
                                }
                            }
                        } else {
                            for (c=matContainers->Rows[row]; c<=matContainers->Diag[row]-1; c++) {
                                if (matContainers->Cols[c] == col) {
                                    matContainers->Values[c] = matContainers->Values[c] + matrix[dofs*(i+1)-k][dofs*(j+1)-l];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    matContainers = NULL;
}

-(void)setSymmetricDirichletInGlobal:(FEMSolution *)solution: (int)n: (double)val {
    
    int i, j, k, l, m, k1, k2;
    int *buffer;
    BOOL isMass, isDamp;
    matrixArraysContainer *matContainers;
    
    matContainers = solution.matrix.getContainers;
    
    isMass = (matContainers->MassValues != NULL) ? YES : NO;
    isMass = (isMass == YES && matContainers->sizeMassValues == matContainers->sizeValues) ? YES : NO;
    isDamp = (matContainers->DampValues != NULL) ? YES : NO;
    isDamp = (isDamp == YES && matContainers->sizeDampValues == matContainers->sizeValues) ? YES : NO;
    
    for (l=matContainers->Rows[n]; l<=matContainers->Rows[n+1]-1; l++) {
        i = matContainers->Cols[l];
        if (i == n) continue;
        
        if (n > i) {
            k1 = matContainers->Diag[i]+1;
            k2 = matContainers->Rows[i+1]-1;
        } else {
            k1 = matContainers->Rows[i];
            k2 = matContainers->Diag[i]-1;
        }
        
        k = k2 - k1;
        if (k <= 30) {
            for (j=k1; j<=k2; j++) {
                if (matContainers->Cols[j] == n) {
                    matContainers->Rows[i] = matContainers->RHS[i] - matContainers->Values[j] * val;
                    matContainers->Values[j] = 0.0;
                    if (isMass == YES) matContainers->MassValues[j] = 0.0;
                    if (isDamp == YES) matContainers->DampValues[j] = 0.0;
                    break;
                } else if (matContainers->Cols[j] > n) {
                    break;
                }
            }
        } else {
            buffer = intvec(k1, k2);
            m = 0;
            for (i=k1; i<=k2; i++) {
                buffer[m] = matContainers->Cols[i];
                m++;
            }
            j = [self CRS_Search:k1 :buffer :n];
            if (j >= 0) {
                j = j + k1;
                matContainers->RHS[i] = matContainers->RHS[i] - matContainers->Values[j] * val;
                matContainers->Values[j] = 0.0;
                if (isMass == YES) matContainers->MassValues[j] = 0.0;
                if (isDamp == YES) matContainers->DampValues[j] = 0.0;
            }
            free_ivector(buffer, k1, k2);
        }
    }
    
    [self zeroRowInGlobal:solution :n];
    matContainers->RHS[n] = val;
    matContainers->Values[matContainers->Diag[n]] = 1.0;
    
    matContainers = NULL;
}

-(void)zeroRowInMatrix:(FEMMatrix *)a: (int)n {
    
    int i;
    matrixArraysContainer *aContainers;
    
    aContainers = a.getContainers;
    
    for (i=aContainers->Rows[n]; i<=aContainers->Rows[n+1]-1; i++) {
        aContainers->Values[i] = 0.0;
    }
    
    if (aContainers->MassValues != NULL) {
        if (aContainers->sizeMassValues == aContainers->sizeValues) {
            for (i=aContainers->Rows[n]; i<=aContainers->Rows[n+1]-1; i++) {
                aContainers->MassValues[i] = 0.0;
            }
        }
    }
    
    if (aContainers->DampValues != NULL) {
        if (aContainers->sizeDampValues == aContainers->sizeValues) {
            for (i=aContainers->Rows[n]; i<=aContainers->Rows[n+1]-1; i++) {
                aContainers->DampValues[i] = 0.0;
            }
        }
    }
    
    aContainers = NULL;
}

-(void)sortInMatrix:(FEMMatrix *)a: (BOOL *)alsoValues {
/************************************************************************************************
    Sort columns to ascending order for rows of a CRS format matrix
 
    Arguments:
    Matrix_t *a             ->  solution class containing the matrix
    BOOL *alsoValues        ->  whether values are sorted
 
************************************************************************************************/
    
    int i, j, k, n;
    int *buffer1;
    double *buffer2;
    BOOL sortValues;
    matrixArraysContainer *aContainers;
    
    aContainers = a.getContainers;
    
    sortValues = NO;
    if (alsoValues != NULL) sortValues = *alsoValues;
    
    n = a.numberOfRows;
    
    if (a.isOrdered == NO) {
        if (sortValues == YES) {
            buffer1 = intvec(0, aContainers->sizeValues-1);
            buffer2 = doublevec(0, aContainers->sizeValues-1);
            memset( buffer1, 0, (aContainers->sizeValues*sizeof(buffer1)) );
            memset( buffer2, 0.0, (aContainers->sizeValues*sizeof(buffer2)) );
            for (i=0; i<n; i++) {
                k = 0;
                for (j=aContainers->Rows[i]; j<=aContainers->Rows[i+1]-1; j++) {
                    buffer1[k] = aContainers->Cols[j];
                    buffer2[k] = aContainers->Values[j];
                    k++;
                }
                sort(aContainers->Rows[i+1]-aContainers->Rows[i], buffer1-1, buffer2-1);
                k = 0;
                for (j=aContainers->Rows[i]; j<=aContainers->Rows[i+1]-1; j++) {
                    aContainers->Cols[j] = buffer1[k];
                    aContainers->Values[j] = buffer2[k];
                    k++;
                }
            }
            free_ivector(buffer1, 0, aContainers->sizeValues-1);
            free_dvector(buffer2, 0, aContainers->sizeValues-1);
            
        } else {
            buffer1 = intvec(0, aContainers->sizeValues-1);
            memset( buffer1, 0, (aContainers->sizeValues*sizeof(buffer1)) );
            for (i=0; i<n; i++) {
                k = 0;
                for (j=aContainers->Rows[i]; j<=aContainers->Rows[i+1]-1; j++) {
                    buffer1[k] = aContainers->Cols[j];
                    k++;
                }
                sort(aContainers->Rows[i+1]-aContainers->Rows[i], buffer1-1);
                k = 0;
                for (j=aContainers->Rows[i]; j<=aContainers->Rows[i+1]-1; j++) {
                    aContainers->Cols[j] = buffer1[k];
                    k++;
                }
            }
            free_ivector(buffer1, 0, aContainers->sizeValues-1);
        }
        
        if (aContainers->Diag != NULL) {
            for (i=0; i<n; i++) {
                for (j=aContainers->Rows[i]; j<=aContainers->Rows[i+1]-1; j++) {
                    if (aContainers->Cols[j] == i) {
                        aContainers->Diag[i] = j;
                        break;
                    }
                }
            }
        }
        
        a.ordered = YES;
    }
    
    aContainers = NULL;
}

-(void)setMatrixElementInMatrix:(FEMMatrix *)a: (int)i: (int)j: (double)value {
/************************************************************************************************
    Set a given value to an element of a CRS format Matrix
 
    Arguments:
        Matrix_t *a             ->  the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
 
************************************************************************************************/
    
    int ii, jj, k;
    int *buffer;
    matrixArraysContainer *aContainers;
    
    aContainers = a.getContainers;
    
    if (aContainers->Diag != NULL || i != j || a.isOrdered == NO) {
        jj = aContainers->Rows[i+1]-aContainers->Rows[i];
        buffer = intvec( 0, jj );
        memset( buffer, 0.0, ((jj+1)*sizeof(buffer)) );
        for (ii=aContainers->Rows[i]; ii<=aContainers->Rows[i+1]-1; ii++) {
            buffer[ii] = aContainers->Cols[ii];
        }
        k = [self CRS_Search:aContainers->Rows[i+1]-aContainers->Rows[i] :buffer :j];
        if (k < 0) {
            warnfunct("CRS:setMatrixElement", "Trying to set value to non existent element:");
            printf("%d %d %f\n", i, j, value);
            return;
        }
        k = k + aContainers->Rows[i];
        free_ivector(buffer, 0, jj);
    } else {
        k = aContainers->Diag[i];
    }
    
    aContainers->Values[k] = value;
    
    aContainers = NULL;
}

@end
