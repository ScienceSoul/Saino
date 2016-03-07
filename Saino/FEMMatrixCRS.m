//
//  FEMMatrixCRS.m
//  Saino
//
//  Created by Hakime Seddik on 23/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Accelerate/Accelerate.h>
#import "FEMMatrixCRS.h"
#import "Utils.h"

@interface FEMMatrixCRS ()

-(int)FEMMatrixCRS_SearchWithinLength:(int)n inArray:(int * __nonnull)array theValue:(int)value;

@end


@implementation FEMMatrixCRS

#pragma mark Private methods

-(int)FEMMatrixCRS_SearchWithinLength:(int)n inArray:(int * __nonnull)array theValue:(int)value {
    
    int lower, upper, lou, index;
    
    index = -1;
    upper = n;
    lower = 1;
    
    // Handle the special case returns -1
    if (upper == 0) return index;
    
    while (1) {
        if (array[lower-1] == value) {
            index = lower-1;
            break;
        } else if (array[upper-1] == value) {
            index = upper-1;
            break;
        }
        
        if ( (upper-lower) > 1 ) {
            lou = (upper+lower) >> 1;
            if (array[lou-1] < value) {
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

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
    }
    
    return self;
}

/********************************************************************************************************************
 
    Create the structures required for a CRS format matrix
    
    int rows            ->  Number of rows in the matrix
    int totalNonZeros   ->  Total number of nonzero entries in the matrix
    int *rowNonZeros    ->  Number of nonzero entries in the rows of the matrix
    int degreesFreedom  ->  N degrees of freedom
    int *reorder        ->  Permutation index for bandwidth reduction
    int sizeOfReorder   ->  Size of permutation array
    BOOL allocateValues ->  Should the values arrays be allocated?
 
*********************************************************************************************************************/
-(FEMMatrix * __nonnull)createMatrixWithNumberOfRows:(int)rows totalNonZeros:(int)totalNonZeros rowNonZeros:(int * __nonnull)rowNonZeros degreesFreedom:(int)degreesFreedom reorder:(int * __nonnull)reorder sizeOfReorder:(int)sizeOfReorder allocateValues:(BOOL)allocateValues {
    
    int i, j, k;
    FEMMatrix *matrix;
    matrixArraysContainer *matContainers = NULL;
    
    matrix = [[FEMMatrix alloc] init];
    
    k = degreesFreedom*degreesFreedom*totalNonZeros;
    
    matContainers = matrix.getContainers;
    matContainers->Rows = intvec(0, (rows+1)-1);
    matContainers->Diag = intvec(0, rows-1);
    matContainers->Cols = intvec(0, k-1);
    matContainers->sizeRows = rows+1;
    matContainers->sizeDiag = rows;
    matContainers->sizeCols = k;
    
    if (allocateValues == YES) {
        matContainers->Values = doublevec(0, k-1);
        matContainers->sizeValues = k;
    }

    j = 0;
    for (i=0; i<sizeOfReorder; i++) {
        if (reorder[i] >= 0) {
            matContainers->Diag[reorder[i]] = j;
            j++;
        }
    }
    
    matrix.numberOfRows = rows;
    matContainers->Rows[0] = 0;
    for (i=1; i<rows; i++) {
        j = matContainers->Diag[(i-1)/degreesFreedom];
        matContainers->Rows[i] = matContainers->Rows[i-1] + degreesFreedom*rowNonZeros[j];
    }
    j = matContainers->Diag[(rows-1)/degreesFreedom];
    matContainers->Rows[rows] = matContainers->Rows[rows-1] + degreesFreedom*rowNonZeros[j];
    
    // Because those arrays contains indexes ranging from 0 to n, we initialize them at -1
    memset( matContainers->Cols, -1, matContainers->sizeCols*sizeof(int) );
    memset( matContainers->Diag, -1, matContainers->sizeDiag*sizeof(int) );
    
    matrix.ordered = NO;
    
    return matrix;
}

-(void)zeroRowInGlobal:(FEMSolution * __nonnull)solution numberOfRows:(int)n {
    
    int i;
    matrixArraysContainer *matContainers = NULL;
    
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
}

/************************************************************************************************
    Sort columns to ascending order for rows of a CRS format matrix
 
    Arguments:
        FEMSolution *solution   ->  solution class containing the matrix
        BOOL *alsoValues        ->  whether values are sorted
************************************************************************************************/
-(void)sortGlobal:(FEMSolution * __nonnull)solution alsoValues:(BOOL * __nullable)alsoValues {
    
    int i, j, n;
    BOOL sortValues;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    sortValues = NO;
    if (alsoValues != NULL) sortValues = *alsoValues;
    
    n = solution.matrix.numberOfRows;
    
    if (solution.matrix.isOrdered == NO) {
        if (sortValues == YES) {
            for (i=0; i<n; i++) {
                sort(matContainers->Rows[i+1]-matContainers->Rows[i],  matContainers->Cols+(matContainers->Rows[i]-1), matContainers->Values+(matContainers->Rows[i]-1));
                }
        } else {
            for (i=0; i<n; i++) {
                vDSP_vsort((float *)matContainers->Cols+matContainers->Rows[i], matContainers->Rows[i+1]-matContainers->Rows[i], 1);
            }
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
}

/*********************************************************************************************************
    Set a given value to an element of a CRS format Matrix
 
    Arguments:
        FEMSolution *solution   ->  solution class containing the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
*********************************************************************************************************/
-(void)setElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value {
    
    int k;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    if (matContainers->Diag == NULL || i != j || solution.matrix.isOrdered == NO) {
        k = [self FEMMatrixCRS_SearchWithinLength:matContainers->Rows[i+1]-matContainers->Rows[i] inArray:matContainers->Cols+matContainers->Rows[i] theValue:j];
        if (k < 0) {
            NSLog(@"FEMMatrixCRS:setMatrixElementInGlobal: trying to set value to non existent element: %d %d %f.\n", i, j, value);
            return;
        }
        k = k + matContainers->Rows[i];
    } else {
        k = matContainers->Diag[i];
    }
    matContainers->Values[k] = value;
}

/************************************************************************************************************
    Add a given value to an element of a CRS format Matrix
 
    Arguments:
        FEMSolution *solution   ->  solution class containing the matrix
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be added
************************************************************************************************************/
-(void)addToElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value {
    
    int k;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    if (matContainers->Diag == NULL || i != j || solution.matrix.isOrdered == NO) {
        k = [self FEMMatrixCRS_SearchWithinLength:matContainers->Rows[i+1]-matContainers->Rows[i] inArray:matContainers->Cols+matContainers->Rows[i] theValue:j];
        if (k < 0 && value != 0) NSLog(@"FEMMatrixCRS:addToMatrixElementInGlobal: trying to add value to non existent element: %d %d %f.\n", i, j, value);
        if (k < 0) return;
        k = k + matContainers->Rows[i];
    } else {
        k = matContainers->Diag[i];
    }
    matContainers->Values[k] = matContainers->Values[k]+value;
}


/******************************************************************************************************************************************
 
    Add a set of values (i.e., element stiffness matrix) to a CRS format matrix
 
    Arguments:
 
        FEMSolution *solution  -> Solution class holding the global matrix
        double **matrix        -> (n x dofs) x (n x dofs) matrix holding the values to be
                                  added to the CRS format matrix
        int n                  -> number of nodes in element
        int dofs               -> number of degrees of freemdom for one node
        int *indexes           -> Maps element node number to global (or partition) node number
                                  (to matrix rows and cols if dofs = 1)
*******************************************************************************************************************************************/
-(void)glueLocalMatrix:(double * __nonnull * __nonnull)localMatrix inGlobal:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dofs:(int)dofs indexes:(int * __nonnull)indexes {
    
    int i, j, k, l, c, row, col;
    matrixArraysContainer *matContainers = NULL;
    
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
                            matContainers->Values[c] = matContainers->Values[c] + localMatrix[i][j];
                            break;
                        }
                    }
                } else {
                    for (c=matContainers->Rows[row]; c<=matContainers->Diag[row]-1; c++) {
                        if (matContainers->Cols[c] == col) {
                            matContainers->Values[c] = matContainers->Values[c] + localMatrix[i][j];
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
                                    matContainers->Values[c] = matContainers->Values[c] + localMatrix[dofs*(i+1)-k][dofs*(j+1)-l];
                                    break;
                                }
                            }
                        } else {
                            for (c=matContainers->Rows[row]; c<=matContainers->Diag[row]-1; c++) {
                                if (matContainers->Cols[c] == col) {
                                    matContainers->Values[c] = matContainers->Values[c] + localMatrix[dofs*(i+1)-k][dofs*(j+1)-l];
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/*************************************************************************************************************
 
    When the Dirichlet conditions are set by zeroing the row except for setting the diagonal entry to one,
    the matrix symmetry is broken. This routine maintains the symmetric structure of the matrix equation.
 
    Arguments:
 
        FEMMatrix *matrix    -> solution class holding the global matrix
        int n                -> index of the dofs to be fixed
        double value         -> Dirichlet value to be set
************************************************************************************************************/
-(void)setSymmetricDirichletInGlobal:(FEMSolution * __nonnull)solution atIndex:(int)n value:(double)value {
    
    int i, j, k, l, k1, k2;
    BOOL isMass, isDamp;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    isMass = (matContainers->MassValues != NULL) ? YES : NO;
    if (isMass == YES) isMass = (isMass == YES && matContainers->sizeMassValues == matContainers->sizeValues) ? YES : NO;
    isDamp = (matContainers->DampValues != NULL) ? YES : NO;
    if (isDamp == YES) isDamp = (isDamp == YES && matContainers->sizeDampValues == matContainers->sizeValues) ? YES : NO;
    
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
        
        k = k2 - k1 + 1;
        if ((k+1) <= 30) {
            for (j=k1; j<=k2; j++) {
                if (matContainers->Cols[j] == n) {
                    matContainers->RHS[i] = matContainers->RHS[i] - matContainers->Values[j] * value;
                    matContainers->Values[j] = 0.0;
                    if (isMass == YES) matContainers->MassValues[j] = 0.0;
                    if (isDamp == YES) matContainers->DampValues[j] = 0.0;
                    break;
                } else if (matContainers->Cols[j] > n) {
                    break;
                }
            }
        } else {
            j = [self FEMMatrixCRS_SearchWithinLength:k inArray:matContainers->Cols+k1 theValue:n];
            if (j >= 0) {
                j = j + k1;
                matContainers->RHS[i] = matContainers->RHS[i] - matContainers->Values[j] * value;
                matContainers->Values[j] = 0.0;
                if (isMass == YES) matContainers->MassValues[j] = 0.0;
                if (isDamp == YES) matContainers->DampValues[j] = 0.0;
            }
        }
    }
    
    [self zeroRowInGlobal:solution numberOfRows:n];
    matContainers->RHS[n] = value;
    matContainers->Values[matContainers->Diag[n]] = 1.0;
}

/*******************************************************************************************
 
    Description:
        Matrix vector product (v = Au) for a matrix given in CRS format. The matrix
        is accessed from the solution class. Real version.
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix.
        double *u              -> Vector to multiply
        double *v              -> Result vector
*******************************************************************************************/
-(void)matrixVectorMultiplyInGlobal:(FEMSolution * __nonnull)solution vector:(double * __nonnull)u result:(double * __nonnull)v {
    
    int i, j, n;
    double rsum;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    n = solution.matrix.numberOfRows;
    
    //TODO: Can we optimize the loop below?
    for (i=0; i<n; i++) {
        rsum = 0.0;
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            rsum = rsum + u[matContainers->Cols[j]] * matContainers->Values[j];
        }
        v[i] = rsum;
    }
}

/*******************************************************************************************
 
    Description:
        Matrix vector product (v = Au) for a matrix given in CRS format. The matrix
        is accessed from the solution class. Complex version.
 
    Arguments:
 
        FEMSolution *solution  -> Class holding input matrix.
        double *u              -> Vector to multiply
        double *v              -> Result vector
*******************************************************************************************/
-(void)complexMatrixVectorMultiplyInGlobal:(FEMSolution * __nonnull)solution vector:(double complex * __nonnull)u result:(double complex * __nonnull)v {
    
    int i, j, n;
    double complex s, rsum;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = solution.matrix.getContainers;
    
    n = solution.matrix.numberOfRows / 2;
    
    //TODO: Can we optimize the loop below
    for (i=0; i<n; i++) {
        rsum = 0.0 + 0.0 * I;
        for (j=matContainers->Rows[2*i-1]; j<=matContainers->Rows[2*i]-1; j+=2) {
            s = matContainers->Values[j] + (-matContainers->Values[j+1] * I);
            rsum = rsum + s * u[(matContainers->Cols[j]+1)/2];
        }
        v[i] = rsum;
    }
}

/*********************************************************************************************
    Makes an algebraic lower order scheme assuming steady state advection-diffusion equation.
    This can be applied together with flux corrected transport (FCT) scheme.
 
    Dimitri Kuzmin (2008): "Explicit and implicit FEM-FCT algorithms with flux linearization"
    The current implementation is really kust a starting point
 
    This method takes as argument a solution from which the matrix is used  
    or an individual matrix
*********************************************************************************************/
-(void)fctlLowOrderInSolution:(FEMSolution * __nullable)solution orMatrix:(FEMMatrix * __nullable)matrix {
    
    int i, j, k, k2, n;
    double aij, aji, aii, dij;
    BOOL found;
    matrixArraysContainer *matContainers = NULL;
    
    if (solution == nil && matrix == nil) {
        NSLog(@"FEMMatrixCRS:fctlLowOrderInSolution: no matrix available. At least one the method argumens should be non-nil.\n");
        return;
    }
    
    NSLog(@"FEMMatrixCRS:fctlLowOrderInSolution: making low order FCT correction to matrix.\n");
    
    if (solution != nil) {
        matContainers = solution.matrix.getContainers;
        n = solution.matrix.numberOfRows;
    } else {
        matContainers = matrix.getContainers;
        n = matrix.numberOfRows;
    }
    
    if (matContainers->FCT_D == NULL) {
        matContainers->FCT_D = doublevec(0, matContainers->sizeValues-1);
        matContainers->sizeFct = matContainers->sizeValues;
    }
    memset( matContainers->FCT_D, 0.0, matContainers->sizeFct*sizeof(double) );
    
    if (matContainers->BulkValues == NULL) {
        matContainers->BulkValues = doublevec(0, matContainers->sizeValues-1);
        matContainers->sizeBulkValues = matContainers->sizeValues;
    }
    memcpy(matContainers->BulkValues, matContainers->Values, matContainers->sizeValues*sizeof(double));
    
    for (i=0; i<n; i++) {
        for (k=matContainers->Rows[i]; k<=matContainers->Rows[i+1]-1; k++) {
            j = matContainers->Cols[k];
            
            // Go through the lower symmetric part (i,j) and find the corresponding entry (j,i)
            if (i >= j) continue;
            
            // First find entry (j,i)
            found = NO;
            for (k2=matContainers->Rows[j]; k2<=matContainers->Rows[j+1]-1; k2++) {
                if (matContainers->Cols[k2] == i) {
                    found = YES;
                    break;
                }
            }
            
            if (found == NO) {
                NSLog(@"FEMMatrixCRS:fctlLowOrderInSolution: entry not found, matrix might not be symmetric.\n");
                continue;
            }
            
            // Equation (30) in Kuzmin's paper
            // Modified so that it also work similarly to A and -A
            aij = matContainers->Values[k];
            aji = matContainers->Values[k2];
            aii = matContainers->Values[matContainers->Diag[i]];
            if (aii < 0.0) {
                dij = max(-aij, -aji, 0.0);
            } else {
                dij = min(-aij, -aji, 0.0);
            }
            
            if (/* DISABLES CODE */ (NO)) {
                NSLog(@"FEMMatrixCRS:fctlLowOrderInSolution: ij: %d %d %d %d.\n", i, j, matContainers->Cols[k2], matContainers->Cols[k]);
                NSLog(@"FEMMatrixCRS:fctlLowOrderInSolution: Diag: %d %d.\n", matContainers->Cols[matContainers->Diag[i]], matContainers->Cols[matContainers->Diag[j]]);
                NSLog(@"FEMMatrixCRS:fctlLowOrderInSolution: A: %f %f %f %f.\n", aij, aji, aii, dij);
            }
            
            // Equation (32) in Kuzmin's paper
            if (fabs(dij) > 0.0) {
                matContainers->FCT_D[k] = matContainers->FCT_D[k] + dij;
                matContainers->FCT_D[k2] = matContainers->FCT_D[k2] + dij;
                matContainers->FCT_D[matContainers->Diag[i]] = matContainers->FCT_D[matContainers->Diag[i]] - dij;
                matContainers->FCT_D[matContainers->Diag[j]] = matContainers->FCT_D[matContainers->Diag[j]] -  dij;
            }
        }
    }
    
    for (i=0; i<matContainers->sizeValues; i++) {
        matContainers->Values[i] = matContainers->Values[i] + matContainers->FCT_D[i];
    }
}

/*************************************************************************************************************
 
    Add a set of values (i.e., element stiffness matrix) to a CRS format matrix
 
        Arguments:
 
        FEMMatrix *matrix       -> The matrix class
        double **localMatrix    -> (n x dofs) x (n x dofs) matrix holding the values to be
                                   added to the CRS format matrix
        int n                   -> number of nodes in element
        int dofs                -> number of degrees of freemdom for one node
        int *indexes            -> Maps element node number to global (or partition) node number
                                   (to matrix rows and cols if dofs = 1)
************************************************************************************************************/
-(void)glueLocalMatrix:(double * __nonnull * __nonnull)localMatrix inMatrix:(FEMMatrix * __nonnull)matrix numberOfNodes:(int)numberOfNodes dofs:(int)dofs indexes:(int * __nonnull)indexes {
    
    int i, j, k, l, c, row, col;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    if (dofs == 1) {
        for (i=0; i<numberOfNodes; i++) {
            row = indexes[i];
            if (row < 0) continue;
            for (j=0; j<numberOfNodes; j++) {
                col = indexes[j];
                if (col < 0) continue;
                if (col >= row) {
                    for (c=matContainers->Diag[row]; c<=matContainers->Rows[row+1]-1; c++) {
                        if (matContainers->Cols[c] == col) {
                            matContainers->Values[c] = matContainers->Values[c] + localMatrix[i][j];
                            break;
                        }
                    }
                } else {
                    for (c=matContainers->Rows[row]; c<=matContainers->Diag[row]-1; c++) {
                        if (matContainers->Cols[c] == col) {
                            matContainers->Values[c] = matContainers->Values[c] + localMatrix[i][j];
                            break;
                        }
                    }
                }
            }
        }
    } else {
        for (i=0; i<numberOfNodes; i++) {
            for (k=1; k<=dofs; k++) {
                if (indexes[i] < 0) continue;
                row = dofs * (indexes[i]+1) - k;
                for (j=0; j<numberOfNodes; j++) {
                    for (l=1; l<=dofs; l++) {
                        if (indexes[j] < 0) continue;
                        col = dofs * (indexes[j]+1) - l;
                        if (col >= row) {
                            for (c=matContainers->Diag[row]; c<=matContainers->Rows[row+1]-1; c++) {
                                if (matContainers->Cols[c] == col) {
                                    matContainers->Values[c] = matContainers->Values[c] + localMatrix[dofs*(i+1)-k][dofs*(j+1)-l];
                                }
                            }
                        } else {
                            for (c=matContainers->Rows[row]; c<=matContainers->Diag[row]-1; c++) {
                                if (matContainers->Cols[c] == col) {
                                    matContainers->Values[c] = matContainers->Values[c] + localMatrix[dofs*(i+1)-k][dofs*(j+1)-l];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/*****************************************************************************************
 
    Fill-in the column number to a CRS format matrix (values are not affected in any way)
 
    Arguments:
 
        FEMMatrix *a -> the matrix class
        int i        -> row number of the matrix element
        int j        -> column number of the matrix element
*****************************************************************************************/
-(void)makeMatrixIndex:(FEMMatrix * __nonnull)matrix row:(int)i col:(int)j {
    
    int k, n;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    n = matContainers->Rows[i];
    for (k=matContainers->Rows[i]; k<=matContainers->Rows[i+1]-1; k++) {
        if (matContainers->Cols[k] == j) {
            return;
        } else if (matContainers->Cols[k] < 0) {
            n = k;
            break;
        }
    }
    
    if (matContainers->Cols[n] >= 0) {
        NSLog(@"FEMMatrixCRS:makeMatrixIndex: trying to access non-existent column: %d, %d.\n", n, matContainers->Cols[n]);
        fatal("FEMMatrixCRS:makeMatrixIndex", "Programm terminating now...\n");
    }
    
    matContainers->Cols[n] = j;
}

-(void)zeroRowInMatrix:(FEMMatrix * __nonnull)matrix numberOfRows:(int)n {
    
    int i;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
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
}

/************************************************************************************************
    Sort columns to ascending order for rows of a CRS format matrix
 
    Arguments:
        Matrix_t *a             ->  the matrix class
        BOOL *alsoValues        ->  whether values are sorted
************************************************************************************************/
-(void)sortMatrix:(FEMMatrix * __nonnull)matrix alsoValues:(BOOL * __nullable)alsoValues {
        
    matrixArraysContainer *matContainers = matrix.getContainers;
    
    BOOL sortValues = NO;
    if (alsoValues != NULL) sortValues = *alsoValues;
    
    int n = matrix.numberOfRows;
    
    if (matrix.isOrdered == NO) {
        if (sortValues == YES) {
            for (int i=0; i<n; i++) {
                sort(matContainers->Rows[i+1]-matContainers->Rows[i],  matContainers->Cols+(matContainers->Rows[i]-1), matContainers->Values+(matContainers->Rows[i]-1));
            }
        } else {
            for (int i=0; i<n; i++) {
                vDSP_vsort((float *)matContainers->Cols+matContainers->Rows[i], matContainers->Rows[i+1]-matContainers->Rows[i], 1);
            }
        }
        
        if (matContainers->Diag != NULL) {
            for (int i=0; i<n; i++) {
                for (int j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
                    if (matContainers->Cols[j] == i) {
                        matContainers->Diag[i] = j;
                        break;
                    }
                }
            }
        }
        
        matrix.ordered = YES;
    }
}

/************************************************************************************************
    Set a given value to an element of a CRS format Matrix
 
    Arguments:
        Matrix_t *a             ->  the matrix class
        int i, j                ->  row and column numbers respectively of the matrix element
        double value            ->  value to be set
 
************************************************************************************************/
-(void)setElementInMatrix:(FEMMatrix * __nonnull)matrix row:(int)i col:(int)j value:(double)value {
    
    int k;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    if (matContainers->Diag != NULL || i != j || matrix.isOrdered == NO) {
        k = [self FEMMatrixCRS_SearchWithinLength:matContainers->Rows[i+1]-matContainers->Rows[i] inArray:matContainers->Cols+matContainers->Rows[i] theValue:j];
        if (k < 0) {
            NSLog(@"FEMMatrixCRS:setMatrixElementInMatrix: trying to set value to non existent element: %d %d %f.\n", i, j, value);
            return;
        }
        k = k + matContainers->Rows[i];
    } else {
        k = matContainers->Diag[i];
    }
    
    matContainers->Values[k] = value;
}

-(void)applyProjector:(FEMMatrix * __nonnull)pMatrix values:(double * __nonnull)u permutation:(int * __nullable)uperm values:(double * __nonnull)v permutation:(int * __nullable)vperm transpose:(BOOL * __nullable)trans {
    
    int i, j, k, l, n;
    matrixArraysContainer *containers = NULL;
    BOOL ltrans, any;
    
    ltrans = NO;
    if (trans != NULL) {
        ltrans = *trans;
    }
    
    n = pMatrix.numberOfRows;
    containers = pMatrix.getContainers;
    
    if (uperm != NULL && vperm != NULL) {
        if (ltrans == YES) {
            for (i=0; i<n; i++) {
                k = uperm[i];
                if (k >= 0) {
                    for (j=containers->Rows[i]; j<=containers->Rows[i+1]-1; j++) {
                        l = vperm[containers->Cols[j]];
                        if (l >= 0) v[l] = v[l] + u[k] * containers->Values[j];
                    }
                }
            }
        } else {
            for (i=0; i<n; i++) {
                l = vperm[i];
                if (l >= 0) {
                    any = NO;
                    for (j=containers->Rows[i]; j<=containers->Rows[i+1]-1; j++) {
                        if (containers->Values[j] != 0.0) {
                            any = YES;
                            break;
                        }
                    }
                    if (any == YES) v[l] = 0.0;
                }
            }
            
            for (i=0; i<n; i++) {
                l = vperm[i];
                if (l >= 0) {
                    for (j=containers->Rows[i]; j<=containers->Rows[i+1]-1; j++) {
                        k = uperm[containers->Cols[j]];
                        if (k >= 0) v[l] = v[l] + u[k] * containers->Values[j];
                    }
                }
            }
        }
    } else {
        if (ltrans == YES) {
            for (i=0; i<n; i++) {
                for (j=containers->Rows[i]; j<=containers->Rows[i+1]-1; j++) {
                    v[containers->Cols[j]] = v[containers->Cols[j]] + u[i] * containers->Values[j];
                }
            }
        } else {
            for (i=0; i<n; i++) {
                for (j=containers->Rows[i]; j<=containers->Rows[i+1]-1; j++) {
                    v[i] = v[i] + u[containers->Cols[j]] * containers->Values[j];
                }
            }
        }
    }
}

/*******************************************************************************************
 
    Description:
    Matrix vector product (v = Au) for a matrix given in CRS format.
 
    Arguments:
 
        FEMMatrix *matrix  -> the input matrix class.
        double *u          -> Vector to multiply
        double *v          -> Result vector
 
*******************************************************************************************/
-(void)matrixVectorMultiply:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)u result:(double * __nonnull)v {
    
    int i, j, n;
    double rsum;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    n = matrix.numberOfRows;
    
    //TODO: Can we optimize the loop below?
    for (i=0; i<n; i++) {
        rsum = 0.0;
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            rsum = rsum + u[matContainers->Cols[j]] * matContainers->Values[j];
        }
        v[i] = rsum;
    }
}

/*******************************************************************************************
 
    Description:
    Matrix vector product (v = Au) for a matrix given in CRS format.
 
    Arguments:
 
        FEMMatrix *matrix  -> the input matrix class.
        double *u          -> Vector to multiply
        double *v          -> Result vector
 
*******************************************************************************************/
-(void)complexMatrixVectorMultiply:(FEMMatrix * __nonnull)matrix vector:(double complex * __nonnull)u result:(double complex * __nonnull)v {
    
    int i, j, n;
    double complex s, rsum;
    matrixArraysContainer *matContainers = NULL;
    
    matContainers = matrix.getContainers;
    
    n = matrix.numberOfRows / 2;
    
    //TODO: Can we optimize the loop below
    for (i=0; i<n; i++) {
        rsum = 0.0 + 0.0 * I;
        for (j=matContainers->Rows[2*i]; j<=matContainers->Rows[2*i+1]-1; j+=2) {
            s = matContainers->Values[j] + (-matContainers->Values[j+1] * I);
            rsum = rsum + s * u[(matContainers->Cols[j]+1)/2];
        }
        v[i] = rsum;
    }
}

-(void)zeroMatrix:(FEMMatrix * __nonnull)matrix {
    
    matrixArraysContainer *matContainers = NULL;
    matContainers = matrix.getContainers;
    if (matContainers == NULL) return;
    
    if (matContainers->Values != NULL) memset(matContainers->Values, 0.0, matContainers->sizeValues*sizeof(double) );
}

@end
