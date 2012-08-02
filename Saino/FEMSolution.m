//
//  FEMSolution.m
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMSolution.h"

@interface FEMSolution ()

-(int)initializeILU1:(matrixArraysContainer *)containers: (int)n;

@end

@implementation FEMSolution

#pragma mark Private methods...

-(int)initializeILU1:(matrixArraysContainer *)containers: (int)n {
    
    int i, j, k, l, nonZeros, rowMin, rowMax;
    int *C;
    
    containers->ILURows = intvec(0, (n+1)-1);
    containers->ILUDiag = intvec(0, n-1);
    
    if (containers->ILURows == NULL || containers->ILUDiag == NULL) {
        errorfunct("initializeILU1", "Memory allocation error.");
    }
    
    // Count fills row by row   
    C = intvec(0, n-1);
    for (i=0; i<n; i++) {
        C[i] = 0;
    }
    nonZeros = containers->Rows[(n+1)-1]-1;
    
    for (i=0; i<n; i++) {
        for (k=containers->Rows[i]; k<=containers->Rows[i+1]-1; k++) {
            C[containers->Cols[k]] = 1;
        }
        
        for (k=containers->Cols[containers->Rows[i]]; k<=i-1; k++) {
            if (C[k] != 0) {
                for (l=containers->Diag[k]+1; l<=containers->Rows[k+1]; l++) {
                    j = containers->Cols[l];
                    if (C[j] == 0) nonZeros = nonZeros + 1;
                }
            }
        }
        
        for (k=containers->Rows[i]; k<=containers->Rows[i+1]-1; k++) {
            C[containers->Cols[k]] = 0;
        }
    }
    
    containers->ILUCols = intvec(0, nonZeros-1);
    if (containers->ILUCols == NULL) {
        errorfunct("initializeILU1", "Memory allocation error.");
    }
    
    // Update row nonzero structures
    for (i=0; i<n; i++) {
        C[i] = 0;
    }
    containers->ILURows[0] = 0;
    for (i=0; i<n; i++) {
        for (k=containers->Rows[i]; k<=containers->Rows[i+1]-1; k++) {
            C[containers->Cols[k]] = 1;
        }
        
        rowMin = containers->Cols[ containers->Rows[i] ];
        rowMax = containers->Cols[ containers->Rows[i+1]-1 ];
        
        for (k=rowMin ; k<=i-1; k++) {
            if (C[k] == 1) {
                for (l=containers->Diag[k]+1; l<=containers->Rows[k+1]-1; l++) {
                    j = containers->Cols[l];
                    if (C[j] == 0) {
                        C[j] = 2;
                        rowMax = max(rowMax, j);
                    }
                }
            }
        }
        
        j = containers->ILURows[i] - 1;
        for (k=rowMin; k<=rowMax; k++) {
            if (C[k] > 0) {
                j = j + 1;
                C[k] = 0;
                containers->ILUCols[j] = k;
                if (k == i) containers->ILUDiag[i] = j;
            }
        }
        
        containers->ILURows[i+1] = j + 1;
    }
    free_ivector(C, 0, n-1);
    
    return nonZeros;
    
}

#pragma mark Public methods...

@synthesize matrix = _matrix;
@synthesize variable = _variable;
@synthesize exportedVariable = _exportedVariable;
@synthesize mesh = _mesh;

#pragma mark Initializations

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here

    }
    
    return self;
}


-(void)initializeILU:(int)ilun {
    
    int i, n, m;
    
    FEMMatrix *a1;
    matrixArraysContainer *matContainers;
    matrixArraysContainer *a1Containers;
    
    a1 = [[FEMMatrix alloc] init];
    matContainers = self.matrix.getContainers;
    a1Containers = a1.getContainers;
    
    n = self.matrix.numberOfRows;
    
    if (ilun == 0) {
        matContainers->ILURows = matContainers->Rows;
        matContainers->ILUCols = matContainers->Cols;
        matContainers->ILUDiag = matContainers->Diag;
        
        matContainers->sizeILURows = matContainers->sizeRows;
        matContainers->sizeILUCols = matContainers->sizeCols;
        matContainers->sizeILUDiag = matContainers->sizeDiag;
    } else {

        matContainers->sizeILUCols = [self initializeILU1:matContainers :n];
        matContainers->sizeILURows = n+1;
        matContainers->sizeILUDiag = n;
        
        if (ilun > 1) {
            
            for (i=0; i<ilun-1; i++) {
                
                a1Containers->Cols = matContainers->ILUCols;
                a1Containers->Rows = matContainers->ILURows;
                a1Containers->Diag = matContainers->ILUDiag;
                
                m = [self initializeILU1:a1Containers :n];
                
                matContainers->ILUCols = a1Containers->ILUCols;
                matContainers->ILURows = a1Containers->ILURows;
                matContainers->ILUDiag = a1Containers->ILUDiag;
                
                a1Containers->Cols = NULL;
                a1Containers->Rows = NULL;
                a1Containers->Diag = NULL;
            }
        }
    }
    matContainers->ILUValues = doublevec(0, matContainers->ILURows[(n+1)-1]-1);
    if (matContainers->ILUValues == NULL) {
        errorfunct("initializeILU", "Memory allocation error.");
    }
    matContainers->sizeILUValues = matContainers->ILURows[(n+1)-1]-1;
    
    matContainers = NULL;
    a1Containers = NULL;
}

-(void)initializeCILU:(int)ilun {
    
    int i, j, k;
    int n;
    
    FEMMatrix *a1;
    matrixArraysContainer *matContainers;
    matrixArraysContainer *a1Containers;
    
    a1 = [[FEMMatrix alloc] init];
    matContainers = self.matrix.getContainers;
    a1Containers = a1.getContainers;
    
    n = self.matrix.numberOfRows;
    
    a1.numberOfRows = n/2;
    a1Containers->Rows = intvec(0, (n/2+1)-1);
    a1Containers->Diag = intvec(0, (n/2)-1);
    a1Containers->Cols = intvec(0, (matContainers->sizeCols / 4)-1);
    
    a1Containers->Rows[0] = 0 ;
    k = 0;
    
    for (i=0; i<n; i+=2) {
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j+=2) {
            k = k + 1;
            a1Containers->Cols[k] = (matContainers->Cols[j]+1) / 2;
            if (matContainers->Cols[j] == i) a1Containers->Diag[(i+1)/2] = k;
        }
        a1Containers->Rows[(i+1)/2+1] = k + 1;
    }
    
    if (ilun == 0) {
        matContainers->ILURows = matContainers->Rows;
        matContainers->ILUCols = matContainers->Cols;
        matContainers->ILUDiag = matContainers->Diag;
        
        matContainers->sizeILURows = matContainers->sizeRows;
        matContainers->sizeILUCols = matContainers->sizeCols;
        matContainers->sizeILUDiag = matContainers->sizeDiag;
    } else {
        
        matContainers->sizeILUCols = [self initializeILU1:a1Containers :n/2];
        matContainers->ILUCols = a1Containers->ILUCols;
        matContainers->ILURows = a1Containers->ILURows;
        matContainers->ILUDiag = a1Containers->ILUDiag;
        
        matContainers->sizeILURows = (n/2+1);
        matContainers->sizeILUDiag = n/2;
        
        free_ivector(a1Containers->Rows, 0, (n/2+1)-1);
        free_ivector(a1Containers->Diag, 0, (n/2)-1);
        free_ivector(a1Containers->Cols, 0, (matContainers->sizeCols / 4)-1);
        a1Containers->Rows = NULL;
        a1Containers->Diag = NULL;
        a1Containers->Cols = NULL;
        
        if (ilun > 1) {
            
            for (i=0; i<ilun-1; i++) {
                
                a1Containers->Cols = matContainers->ILUCols;
                a1Containers->Rows = matContainers->ILURows;
                a1Containers->Diag = matContainers->ILUDiag;
                
                k = [self initializeILU1:a1Containers :n/2];
                
                matContainers->ILUCols = a1Containers->ILUCols;
                matContainers->ILURows = a1Containers->ILURows;
                matContainers->ILUDiag = a1Containers->ILUDiag;
                
                a1Containers->Cols = NULL;
                a1Containers->Rows = NULL;
                a1Containers->Diag = NULL;
            }
        }
    }
    matContainers->CILUValues = cdoublevec(0, matContainers->ILURows[(n/2+1)-1]);
    if (matContainers->CILUValues == NULL) {
        errorfunct("initializeCILU", "Memory allocation error.");
    }
    matContainers->sizeCILUValues = matContainers->ILURows[(n/2+1)-1];
    
    matContainers = NULL;
    a1Containers = NULL;
}

-(void)ilutWorkspaceCheck:(int)i: (int)n {
    
    int j, k;
    int *iWork;
    double *cWork;
    matrixArraysContainer *matContainers;
    
    matContainers = self.matrix.getContainers;
    
    k = matContainers->ILURows[i+1] + min(0.75*matContainers->ILURows[i+1], ((n-1)-i)*(1.0*n));
    
    iWork = intvec(0, k-1);
    if (iWork == NULL) {
        errorfunct("ilutWorkspaceCheck", "Memory allocation error.");
    }
    for (j=0; j<=matContainers->ILURows[i+1]-1; j++) {
        iWork[j] = matContainers->ILUCols[j];
    }
    free_ivector(matContainers->ILUCols, 0, matContainers->sizeILUCols-1);
    matContainers->ILUCols = NULL;

    
    cWork = doublevec(0, k-1);
    if (cWork == NULL) {
        errorfunct("ilutWorkspaceCheck", "Memory allocation error.");
    }
    for (j=0; j<=matContainers->ILURows[i+1]-1; j++) {
        cWork[j] = matContainers->ILUValues[j];
    }
    free_dvector(matContainers->ILUValues, 0, matContainers->sizeILUValues-1);
    matContainers->ILUValues = NULL;
    
    matContainers->sizeILUCols = k;
    matContainers->sizeILUValues = k;
    
    matContainers->ILUCols = iWork;
    matContainers->ILUValues = cWork;
    
    iWork = NULL;
    cWork = NULL;
    matContainers = NULL;
    
}

-(void)ilutComplexWorkspaceCheck:(int)i: (int)n {
    
    int j, k;
    int *iWork;
    double complex *cWork;
    matrixArraysContainer *matContainers;
    
    matContainers = self.matrix.getContainers;
    
    k = matContainers->ILURows[i+1] + min(0.75*matContainers->ILURows[i+1], ((n-1)-i)*(1.0*n));
    
    iWork = intvec(0, k-1);
    if (iWork == NULL) {
        errorfunct("ilutComplexWorkspaceCheck", "Memory allocation error.");
    }
    for (j=0; j<=matContainers->ILURows[i+1]-1; j++) {
        iWork[j] = matContainers->ILUCols[j];
    }
    free_ivector(matContainers->ILUCols, 0, matContainers->sizeILUCols-1);
    matContainers->ILUCols = NULL;
    
    cWork = cdoublevec(0, k-1);
    if (cWork == NULL) {
        errorfunct("ilutComplexWorkspaceCheck", "Memory allocation error.");
    }
    for (j=0; j<=matContainers->ILURows[i+1]-1; j++) {
        cWork[j] = matContainers->CILUValues[j];
    }
    free_cdvector(matContainers->CILUValues, 0, matContainers->sizeCILUValues-1);
    matContainers->CILUValues = NULL;
    
    matContainers->sizeILUCols = k;
    matContainers->sizeCILUValues = k;
    
    matContainers->ILUCols = iWork;
    matContainers->CILUValues = cWork;
    
    iWork = NULL;
    cWork = NULL;
    matContainers = NULL;
    
}

#pragma mark Setters and getters for FEMSolution variables

-(id)solutionInfoForKey:(NSString *)key {
    
    return [solutionInfo objectForKey:key];
}

-(NSString *)normalTangentialName {
    
    return normalTangentialName;
}

-(int)simulationID {
    
    return simulationID;
}

-(int)normalTangentialNOFNodes {
    
    return normalTangentialNOFNodes;
}

-(int)coordinateSystemDimension {
    
    return coordinateSystemDimension;
}

-(int)ntElement:(int)i :(int)j {
    
    return ntElement[i][j];
}

-(BOOL)ntZeroingDone:(int)i :(int)j {
    
    return ntZeroingDone[i][j];
}

-(int)boundaryReorder:(int)i {
    
    return boundaryReorder[i];
}

-(double)boundaryNormals:(int)i :(int)j {
    
    return boundaryNormals[i][j];
}

-(double)boundaryTangent1:(int)i :(int)j {
    
    return boundaryTangent1[i][j];
}

-(double)boundaryTangent2:(int)i :(int)j {
    
    return boundaryTangent2[i][j];
}

-(int)size1boundaryNormals {
    
    return size1boundaryNormals;
}

-(int)size2boundaryNormals {
    
    return size2boundaryNormals;
}

-(int)timeOrder {
    
    return timeOrder;
}

-(int)doneTime {
    
    return doneTime;
}

-(int)order {
    
    return order;
}

-(int)nofEigenValues {
    
    return nofEigenValues;
}

-(double)alpha {
    
    return alpha;
}

-(double)beta {
    
    return beta;
}

-(double)dt {
    
    return dt;
}

-(int)defDofs:(int)i {
    
    return defDofs[i];
}

-(int)sizeOfDefDofs {
    
    return sizeDefDofs;
}

-(void)setNormalTangentialName:(NSString *)string {
    
    normalTangentialName = [NSString stringWithString:string];
}

-(void)setSimulationID:(int)i {
    
    simulationID = i;
}

-(void)setNormalTangentialNOFNodes:(int)i {
    
    normalTangentialNOFNodes = i;
}

-(void)setCoordinateSystemDimension:(int)i {
    
    coordinateSystemDimension = i;
}

-(void)setNtElement:(int)i :(int)j :(int)n {
    
    ntElement[i][j] = n;
}

-(void)setNtZeroingDone:(int)i :(int)j :(BOOL)n {
    
    ntZeroingDone[i][j] = n;
}

-(void)setBoundaryReorder:(int)i :(int)n {
    
    boundaryReorder[i] = n;
}

-(void)setBoundaryNormals:(int)i :(int)j :(double)n {
    
    boundaryNormals[i][j] = n;
}

-(void)setBoundaryTangent1:(int)i :(int)j :(double)n {
    
    boundaryTangent1[i][j] = n;
}

-(void)setBoundaryTangent2:(int)i :(int)j :(double)n {
    
    boundaryTangent2[i][j] = n;
}

-(void)setSize1boundaryNormals:(int)n {
    
    size1boundaryNormals = n;
}

-(void)setSize2boundaryNormals:(int)n {
    
    size2boundaryNormals = n;
}

-(void)setTimeOrder:(int)n {
    
    timeOrder = n;
}

-(void)setDoneTime:(int)n {
    
    doneTime = n;
}

-(void)setOrder:(int)n {
    
    order = n;
}

-(void)setAlpha:(double)n {
    
    alpha = n;
}

-(void)setBeta:(double)n {
    
    beta = n;
}

-(void)setDt:(double)n {
    
    dt = n;
}

-(void)setDefDofs:(int)i :(int)n {
    
    defDofs[i] = n;
}

-(void)setSizeOfDefDofs:(int)n {
    
    sizeDefDofs = n;
}

#pragma mark Sizes

-(int)sizeOfBoundaryReorder {
    
    return sizeBoundaryReorder;
}

-(void)setSizeOfBoundaryReorder:(int)i {
    
    sizeBoundaryReorder = i;
}

#pragma mark Methods returning pointers

-(int *)returnPointerToDefDofs {
    
    return defDofs;
}

@end
