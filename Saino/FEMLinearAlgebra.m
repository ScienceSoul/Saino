//
//  FEMLinearAlgebra.m
//  Saino
//
//  Created by Hakime Seddik on 01/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMLinearAlgebra.h"

#import "memory.h"
#import "Utils.h"

@implementation FEMLinearAlgebra

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
    }
    
    return self;
}

-(void)LUDecompositionOfMatrix:(double **)a ofSize:(int)n resultPivot:(int *)pivot {
    
    int i, j, k, l;
    double swap;
    
    for (i=0; i<n; i++) {
        j = i;
        for (k=i+1; k<n; k++) {
            if (fabs(a[i][k]) > fabs(a[i][j])) j = k;
        }
        
        if (fabs(a[i][j]) == 0.0) {
            NSLog(@"FEMLinearAlgebra:LUDecompositionOfMatrix: matrix is singular.\n");
            return;
        }
        
        pivot[i] = j;
        
        if (j != i) {
            for (k=0; k<=i; k++) {
                swap = a[k][j];
                a[k][j] = a[k][i];
                a[k][i] = swap;
            }
        }
        
        for (k=i+1; k<n; k++) {
            a[i][k] = a[i][k] / a[i][i];
        }
        
        for (k=i+1; k<n; k++) {
            if (j != i) {
                swap = a[k][i];
                a[k][i] = a[k][j];
                a[k][j] = swap;
            }
            
            for (l=i+1; l<n; l++) {
                a[k][l] = a[k][l] - a[k][i] * a[i][l];
            }
        }
    }
    pivot[n-1] = n-1;
    if (fabs(a[n-1][n-1]) == 0.0) {
        NSLog(@"FEMLinearAlgebra:LUDecompositionOfMatrix: matrix is (at least almost) singular.\n");
    }
    
}

-(void)invertMatrix:(double **)a ofSize:(int)n {
    
    int i, j, k;
    double s;
    int *pivot;
    
    pivot = intvec(0, n-1);
    
    [self LUDecompositionOfMatrix:a ofSize:n resultPivot:pivot];
    
    for (i=0; i<n; i++) {
        if (fabs(a[i][i]) == 0.0) {
            NSLog(@"FEMLinearAlgebra:invertMatrix: matrix is singular.\n");
            return;
        }
        a[i][i] = 1.0 / a[i][i];
    }
    
    for (i=n-2; i>=0; i--) {
        for (j=n-1; j>=i+1; j--) {
            s = -a[i][j];
            for (k=i+1; k<=j-1; k++) {
                s = s - a[i][k]*a[k][j];
            }
            a[i][j] = s;
        }
    }
    
    for (i=n-2; i>=0; i--) {
        for (j=n-1; j>=i+1; j--) {
            s = 0.0;
            for (k=i+1; k<=j; k++) {
                s = s - a[j][k]*a[k][i];
            }
            a[j][i] = a[i][i]*s;
        }
    }
    
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            s = 0.0;
            for (k=max(i,j); k<n; k++) {
                if (k != i) {
                    s = s + a[i][k]*a[k][j];
                } else {
                    s = s + a[k][j];
                }
            }
            a[i][j] = s;
        }
    }
    
    for (i=n-1; i>=0; i--) {
        if (pivot[i] != i) {
            for (j=0; j<n; j++) {
                s = a[i][j];
                a[i][j] = a[pivot[i]][j];
                a[pivot[i]][j] = s;
            }
        }
    }
    
    free_ivector(pivot, 0, n-1);
}

@end
