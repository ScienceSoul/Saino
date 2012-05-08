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
    
    for (i=[solution matrixRows:n]; i<=[solution matrixRows:n+1]-1; i++) {
        [solution setMatrixValues:i :0.0];
    }
    
    if ([solution isAssociatedMatrixMassValues] == YES) {
        if ([solution matrixSizeOfMassValues] == [solution matrixSizeOfValues]) {
            for (i=[solution matrixRows:n]; i<=[solution matrixRows:n+1]-1; i++) {
                [solution setMatrixMassValues:i :0.0];
            }
        }
    }
    
    if ([solution isAssociatedMatrixDampValues] == YES) {
        if ([solution matrixSizeOfDampValues] == [solution matrixSizeOfValues]) {
            for (i=[solution matrixRows:n]; i<=[solution matrixRows:n+1]-1; i++) {
                [solution setMatrixDampValues:i :0.0];
            }
        }
    }
    
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
    
    sortValues = NO;
    if (alsoValues != NULL) sortValues = *alsoValues;
    
    n = [solution matrixNumberOfRows];
    
    if ([solution matrixOrdered] == 0) {
        if (sortValues == YES) {
            buffer1 = intvec(0, [solution matrixSizeOfValues]-1);
            buffer2 = doublevec(0, [solution matrixSizeOfValues]-1);
            memset( buffer1, 0, ([solution matrixSizeOfValues]*sizeof(buffer1)) );
            memset( buffer2, 0.0, ([solution matrixSizeOfValues]*sizeof(buffer2)) );
            for (i=0; i<n; i++) {
                k = 0;
                for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                    buffer1[k] = [solution matrixCols:j];
                    buffer2[k] = [solution matrixValues:j];
                    k++;
                }
                sort([solution matrixRows:i+1]-[solution matrixRows:i], buffer1-1, buffer2-1);
                k = 0;
                for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                    [solution setMatrixCols:j :buffer1[k]];
                    [solution setMatrixValues:j :buffer2[k]];
                    k++;
                }
            }
            free_ivector(buffer1, 0, [solution matrixSizeOfValues]-1);
            free_dvector(buffer2, 0, [solution matrixSizeOfValues]-1);
            
        } else {
            buffer1 = intvec(0, [solution matrixSizeOfValues]-1);
            memset( buffer1, 0, ([solution matrixSizeOfValues]*sizeof(buffer1)) );
            for (i=0; i<n; i++) {
                k = 0;
                for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                    buffer1[k] = [solution matrixCols:j];
                    k++;
                }
                sort([solution matrixRows:i+1]-[solution matrixRows:i], buffer1-1);
                k = 0;
                for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                    [solution setMatrixCols:j :buffer1[k]];
                    k++;
                }
            }
            free_ivector(buffer1, 0, [solution matrixSizeOfValues]-1);
        }
        
        if ([solution isAssociatedMatrixDiag] == YES) {
            for (i=0; i<n; i++) {
                for (j=[solution matrixRows:i]; j<=[solution matrixRows:i+1]-1; j++) {
                    if ([solution matrixCols:j] == i) {
                        [solution setMatrixDiag:i :j];
                        break;
                    }
                }
            }
        }
        
        [solution setMatrixOrdered:1];
    }
    
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
    
    if ([solution isAssociatedMatrixDiag] == NO || i != j || [solution matrixOrdered] == 0) {
        jj = [solution matrixRows:i]-[solution matrixRows:i+1];
        buffer = intvec( 0, jj );
        memset( buffer, 0.0, ((jj+1)*sizeof(buffer)) );
        for (ii=[solution matrixRows:i]; ii<=[solution matrixRows:i+1]-1; ii++) {
            buffer[ii] = [solution matrixCols:ii];
        }
        k = [self CRS_Search:[solution matrixRows:i+1]-[solution matrixRows:i] :buffer :j];
        if (k < 0) {
            warnfunct("CRS:setMatrixElement", "Trying to set value to non existent element:");
            printf("%d %d %f\n", i, j, value);
            return;
        }
        k = k + [solution matrixRows:i];
        free_ivector(buffer, 0, jj);
    } else {
        k = [solution matrixDiag:i];
    }
    
    [solution setMatrixValues:k :value];
    
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
    
    if ([solution isAssociatedMatrixDiag] == NO || i != j || [solution matrixOrdered] == 0) {
        jj = [solution matrixRows:i]-[solution matrixRows:i+1];
        buffer = intvec( 0, jj );
        memset( buffer, 0.0, ((jj+1)*sizeof(buffer)) );
        for (ii=[solution matrixRows:i]; ii<=[solution matrixRows:i+1]-1; ii++) {
            buffer[ii] = [solution matrixCols:ii];
        }
        k = [self CRS_Search:[solution matrixRows:i+1]-[solution matrixRows:i] :buffer :j];
        if (k < 0 && value != 0) warnfunct("addToMatrixElement", "Trying to add value to non existent element:");
        printf("%d %d %f\n", i, j, value);
        if (k < 0) return;
        k = k + [solution matrixRows:i];
        free_ivector(buffer, 0, jj);
    } else {
        k = [solution matrixDiag:i];
    }
    [solution setMatrixValues:k :[solution matrixValues:k]+value];
    
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
    
    if (dofs == 1) {
        for (i=0; i<n; i++) {
            row = indexes[i];
            if (row < 0) continue;
            for (j=0; j<n; j++) {
                col = indexes[j];
                if (col < 0) continue;
                if (col >= row) {
                    for (c=[solution matrixDiag:row]; c<=[solution matrixRows:row+1]-1; c++) {
                        if ([solution matrixCols:c] == col) {
                            [solution setMatrixValues:c :[solution matrixValues:c] + matrix[i][j]];
                            break;
                        }
                    }
                } else {
                    for (c=[solution matrixRows:row]; c<=[solution matrixDiag:row]-1; c++) {
                        if ([solution matrixCols:c] == col) {
                            [solution setMatrixValues:c :[solution matrixValues:c] + matrix[i][j]];
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
                            for (c=[solution matrixDiag:row]; c<=[solution matrixRows:row+1]-1; c++) {
                                if ([solution matrixCols:c] == col) {
                                    [solution setMatrixValues:c :[solution matrixValues:c] + matrix[dofs*(i+1)-k][dofs*(j+1)-l]];
                                }
                            }
                        } else {
                            for (c=[solution matrixRows:row]; c<=[solution matrixDiag:row]-1; c++) {
                                if ([solution matrixCols:c] == col) {
                                    [solution setMatrixValues:c :[solution matrixValues:c] + matrix[dofs*(i+1)-k][dofs*(j+1)-l]];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
}

-(void)setSymmetricDirichletInGlobal:(FEMSolution *)solution: (int)n: (double)val {
    
    int i, j, k, l, m, k1, k2;
    int *buffer;
    BOOL isMass, isDamp;
    
    isMass = ([solution isAssociatedMatrixMassValues] == YES) ? YES : NO;
    isMass = (isMass == YES && [solution matrixSizeOfMassValues] == [solution matrixSizeOfValues]) ? YES : NO;
    isDamp = ([solution isAssociatedMatrixDampValues] == YES) ? YES : NO;
    isDamp = (isDamp == YES && [solution matrixSizeOfDampValues] == [solution matrixSizeOfValues]) ? YES : NO;
    
    for (l=[solution matrixRows:n]; l<=[solution matrixRows:n+1]-1; l++) {
        i = [solution matrixCols:l];
        if (i == n) continue;
        
        if (n > i) {
            k1 = [solution matrixDiag:i]+1;
            k2 = [solution matrixRows:i+1]-1;
        } else {
            k1 = [solution matrixRows:i];
            k2 = [solution matrixDiag:i]-1;
        }
        
        k = k2 - k1;
        if (k <= 30) {
            for (j=k1; j<=k2; j++) {
                if ([solution matrixCols:j] == n) {
                    [solution setMatrixRHS:i :[solution matrixRHS:i] - [solution matrixValues:j] * val];
                    [solution setMatrixValues:j :0.0];
                    if (isMass == YES) [solution setMatrixMassValues:j :0.0];
                    if (isDamp == YES) [solution setMatrixDampValues:j :0.0];
                    break;
                } else if ([solution matrixCols:j] > n) {
                    break;
                }
            }
        } else {
            buffer = intvec(k1, k2);
            m = 0;
            for (i=k1; i<=k2; i++) {
                buffer[m] = [solution matrixCols:i];
                m++;
            }
            j = [self CRS_Search:k1 :buffer :n];
            if (j >= 0) {
                j = j + k1;
                [solution setMatrixRHS:i :[solution matrixRHS:i] - [solution matrixValues:j] * val];
                [solution setMatrixValues:j :0.0];
                if (isMass == YES) [solution setMatrixMassValues:j :0.0];
                if (isDamp == YES) [solution setMatrixDampValues:j :0.0];
            }
            free_ivector(buffer, k1, k2);
        }
    }
    
    [self zeroRowInGlobal:solution :n];
    [solution setMatrixRHS:n :val];
    [solution setMatrixValues:[solution matrixDiag:n] : 1.0];
    
}

-(void)zeroRowInMatrix:(Matrix_t *)a: (int)n {
    
    int i;
    
    for (i=a->Rows[n]; i<=a->Rows[n+1]-1; i++) {
        a->Values[i] = 0.0;
    }
    
    if (a->MassValues != NULL) {
        if (a->sizeMassValues == a->sizeValues) {
            for (i=a->Rows[n]; i<=a->Rows[n+1]-1; i++) {
                a->MassValues[i] = 0.0;
            }
        }
    }
    
    if (a->DampValues != NULL) {
        if (a->sizeDampValues == a->sizeValues) {
            for (i=a->Rows[n]; i<=a->Rows[n+1]-1; i++) {
                a->DampValues[i] = 0.0;
            }
        }
    }

}

-(void)sortInMatrix:(Matrix_t *)a: (BOOL *)alsoValues {
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
    
    sortValues = NO;
    if (alsoValues != NULL) sortValues = *alsoValues;
    
    n = a->NumberOfRows;
    
    if (a->Ordered == 0) {
        if (sortValues == YES) {
            buffer1 = intvec(0, a->sizeValues-1);
            buffer2 = doublevec(0, a->sizeValues-1);
            memset( buffer1, 0, (a->sizeValues*sizeof(buffer1)) );
            memset( buffer2, 0.0, (a->sizeValues*sizeof(buffer2)) );
            for (i=0; i<n; i++) {
                k = 0;
                for (j=a->Rows[i]; j<=a->Rows[i+1]-1; j++) {
                    buffer1[k] = a->Cols[j];
                    buffer2[k] = a->Values[j];
                    k++;
                }
                sort(a->Rows[i+1]-a->Rows[i], buffer1-1, buffer2-1);
                k = 0;
                for (j=a->Rows[i]; j<=a->Rows[i+1]-1; j++) {
                    a->Cols[j] = buffer1[k];
                    a->Values[j] = buffer2[k];
                    k++;
                }
            }
            free_ivector(buffer1, 0, a->sizeValues-1);
            free_dvector(buffer2, 0, a->sizeValues-1);
            
        } else {
            buffer1 = intvec(0, a->sizeValues-1);
            memset( buffer1, 0, (a->sizeValues*sizeof(buffer1)) );
            for (i=0; i<n; i++) {
                k = 0;
                for (j=a->Rows[i]; j<=a->Rows[i+1]-1; j++) {
                    buffer1[k] = a->Cols[j];
                    k++;
                }
                sort(a->Rows[i+1]-a->Rows[i], buffer1-1);
                k = 0;
                for (j=a->Rows[i]; j<=a->Rows[i+1]-1; j++) {
                    a->Cols[j] = buffer1[k];
                    k++;
                }
            }
            free_ivector(buffer1, 0, a->sizeValues-1);
        }
        
        if (a->Diag != NULL) {
            for (i=0; i<n; i++) {
                for (j=a->Rows[i]; j<=a->Rows[i+1]-1; j++) {
                    if (a->Cols[j] == i) {
                        a->Diag[i] = j;
                        break;
                    }
                }
            }
        }
        
        a->Ordered = 1;
    }
    
}

-(void)setMatrixElementInMatrix:(Matrix_t *)a: (int)i: (int)j: (double)value {
/************************************************************************************************
    Set a given value to an element of a CRS format Matrix
 
    Arguments:
        Matrix_t *a             ->  the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
 
************************************************************************************************/
    
    int ii, jj, k;
    int *buffer;
    
    if (a->Diag != NULL || i != j || a->Ordered == 0) {
        jj = a->Rows[i]-a->Rows[i+1];
        buffer = intvec( 0, jj );
        memset( buffer, 0.0, ((jj+1)*sizeof(buffer)) );
        for (ii=a->Rows[i]; ii<=a->Rows[i+1]-1; ii++) {
            buffer[ii] = a->Cols[ii];
        }
        k = [self CRS_Search:a->Rows[i+1]-a->Rows[i] :buffer :j];
        if (k < 0) {
            warnfunct("CRS:setMatrixElement", "Trying to set value to non existent element:");
            printf("%d %d %f\n", i, j, value);
            return;
        }
        k = k + a->Rows[i];
        free_ivector(buffer, 0, jj);
    } else {
        k = a->Diag[i];
    }
    
    a->Values[k] = value;

}

@end
