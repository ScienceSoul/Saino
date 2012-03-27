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

-(void)BAND_glueLocalMatrix:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes {
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

@end
