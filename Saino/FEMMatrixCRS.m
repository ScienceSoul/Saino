//
//  FEMMatrixCRS.m
//  Saino
//
//  Created by Hakime Seddik on 23/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMatrixCRS.h"

@interface FEMMatrixCRS ()

@end

@implementation FEMMatrixCRS

-(void)CRS_glueLocalMatrix:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes {
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

@end
