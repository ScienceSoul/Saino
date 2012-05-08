//
//  FEMUtilities.m
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMUtilities.h"

@implementation FEMUtilities

-(Matrix_t *)allocateMatrix {
    
    Matrix_t *matrix;
    
    matrix = (Matrix_t*) malloc( sizeof(Matrix_t) * 1 );
    
    matrix->Format = MATRIX_CRS;
    
    matrix->Child = NULL;
    matrix->Parent = NULL;
    matrix->EMatrix = NULL;
    matrix->ConstraintMatrix = NULL;
    
    matrix->Perm = NULL;
    matrix->InvPerm = NULL;
    
    matrix->Cols = NULL;
    matrix->Rows = NULL;
    matrix->Diag = NULL;
    matrix->GRows = NULL;
    
    matrix->RHS = NULL;
    matrix->Force = NULL;
    matrix->RHS_im = NULL;
    
    matrix->Values = NULL;
    matrix->ILUValues = NULL;
    matrix->MassValues = NULL;
    matrix->DampValues = NULL;
    
    matrix->BulkRHS = NULL;
    matrix->BulkValues = NULL;
    
    matrix->ILUCols = NULL;
    matrix->ILURows = NULL;
    matrix->ILUDiag = NULL;
    
    matrix->CRHS = NULL;
    matrix->CForce = NULL;
    
    matrix->RowOwner = NULL;
    
    matrix->CValues = NULL;
    matrix->CILUValues = NULL;
    matrix->CMassValues = NULL;
    matrix->CDampValues = NULL;
    
    matrix->GRows = NULL;
    matrix->GOrder = NULL;
    
    matrix->UMFPack_Numeric = 0;
    
    matrix->Lumped = 0;
    matrix->Ordered = 0;
    matrix->Complex = 0;
    matrix->Symmetric = 0;
    matrix->SolveCount = 0;
    matrix->NumberOfRows = 0;
    
    return matrix;
    
}

-(void)zeroTheNumberOfRows:(int)n inMatrix:(Matrix_t *)a {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    if (a->Format == MATRIX_CRS) {
        crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix zeroRowInMatrix:a :n];
        
    } else if (a->Format == MATRIX_LIST) {
        
        // To do's: implement the zeroRow method for list matrix.
        
    } else if (a->Format == MATRIX_BAND || a->Format == MATRIX_SBAND) {
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix zeroRowInMatrix:a :n];
        
    }
    
}

-(void)setMatrixElement:(Matrix_t *)a: (int)i: (int)j: (double)value {
    
    FEMMatrixCRS *crsMatrix;
    FEMMatrixBand *bandMatrix;
    
    if (a->Format == MATRIX_CRS) {
        crsMatrix = [[FEMMatrixCRS alloc] init];
        [crsMatrix setMatrixElementInMatrix:a :i :j :value];
        
    } else if (a->Format == MATRIX_LIST) {
        // To do's: implement the setMatrixElement method for list matrix.
        
    } else if (a->Format == MATRIX_BAND || a->Format == MATRIX_SBAND) {
        bandMatrix = [[FEMMatrixBand alloc] init];
        [bandMatrix setMatrixElementInMatrix:a :i :j :value];
    }
    
}

-(double)interpolateCurve:(double *)tValues: (double *)fValues: (double)t: (int)n {
    
    int i;
    double f;
    
    for (i=0; i<n; i++) {
        if (tValues[i] >= t) break;
    }
    if (i > n-1) i = n-1;
    if (i < 1) i = 1;
    
    f = (t - tValues[i-1]) / (tValues[i] - tValues[i-1]);
    f = (1.0 - f)*fValues[i-1] + f*fValues[i];
    
    return f;
}

-(void)solveLinearSystem2x2:(double **)a :(double *)x :(double *)b {
    
    double detA;
    
    detA = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    
    if (detA == 0.0) {
        errorfunct("solveLinearSystem2x2", "Singular matrix, bad!!!");
        return;
    }
    
    detA = 1.0 / detA;
    x[0] = detA * ( a[1][1] * b[0] - a[0][1] * b[1] );
    x[1] = detA * ( a[0][0] * b[1] - a[1][0] * b[0] );
    
}

-(void)solveLinearSystem3x3:(double **)a: (double *)x: (double *)b {
    
    double **c, *y, *g, s, t, q;
    
    c = doublematrix(0, 1, 0, 1);
    y = doublevec(0, 1);
    g = doublevec(0, 1);
    
    if ( (fabs(a[0][0]) > fabs(a[0][1])) && (fabs(a[0][0]) > fabs(a[0][2])) ) {
        q = 1.0 / a[0][0];
        s = q * a[1][0];
        t = q * a[2][0];
        c[0][0] = a[1][1] - s * a[0][1];
        c[0][1] = a[1][2] - s * a[0][2];
        c[1][0] = a[2][1] - t * a[0][1];
        c[1][1] = a[2][2] - t * a[0][2];
        
        g[0] = b[1] - s * b[0];
        g[1] = b[2] - t * b[0];
        [self solveLinearSystem2x2:c :y :g];
        
        x[1] = y[0];
        x[2] = y[1];
        x[0] = q * ( b[0] - a[0][1] * x[1] - a[0][2] * x[2] );
    } else if (fabs(a[0][1]) > fabs(a[0][2])) {
        q = 1.0 / a[0][1];
        s = q * a[1][1];
        t = q * a[2][1];
        c[0][0] = a[1][0] - s * a[0][0];
        c[0][1] = a[1][2] - s * a[0][2];
        c[1][0] = a[2][0] - t * a[0][0];
        c[1][1] = a[2][2] - t * a[0][2];
        
        g[0] = b[1] - s * b[0];
        g[1] = b[2] - t * b[0];
        [self solveLinearSystem2x2:c :y :g];
        
        x[0] = y[0];
        x[2] = y[1];
        x[1] = q * ( b[0] - a[0][0] * x[0] - a[0][2] * x[2] );
    } else {
        q = 1.0 / a[0][2];
        s = q * a[1][2];
        t = q * a[2][2];
        c[0][0] = a[1][0] - s * a[0][0];
        c[0][1] = a[1][1] - s * a[0][1];
        c[1][0] = a[2][0] - t * a[0][0];
        c[1][1] = a[2][1] - t * a[0][1];
        
        g[0] = b[1] - s * b[0];
        g[1] = b[2] - t * b[0];
        [self solveLinearSystem2x2:c :y :g];
        
        x[0] = y[0];
        x[1] = y[1];
        x[2] = q * ( b[0] - a[0][0] * x[0] - a[0][1] * x[1] );
    }
    
    free_dmatrix(c, 0, 1, 0, 1);
    free_dvector(y, 0, 1);
    free_dvector(g, 0, 1);
    
}

@end
