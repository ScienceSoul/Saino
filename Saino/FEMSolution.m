//
//  FEMSolution.m
//  Saino
//
//  Created by Hakime Seddik on 18/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMSolution.h"
#import "Utils.h"

@interface FEMSolution ()

-(int)FEMSolution_initializeILU1:(matrixArraysContainer *)containers :(int)n;

@end

@implementation FEMSolution

#pragma mark Private methods...

-(int)FEMSolution_initializeILU1:(matrixArraysContainer *)containers :(int)n {
    
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
@synthesize simulationID = _simulationID;
@synthesize coordinateSystemDimension  = _coordinateSystemDimension;
@synthesize timeOrder = _timeOrder;
@synthesize doneTime = _doneTime;
@synthesize order = _order;
@synthesize nOfEigenValues = _nOfEigenValues;
@synthesize solutionSolveWhen = _solutionSolveWhen;
@synthesize solutionMode = _solutionMode;
@synthesize numberOfActiveElements = _numberOfActiveElements;
@synthesize multigridLevel = _multigridLevel;
@synthesize multiGridTotal = _multiGridTotal;
@synthesize multiGridSweep = _multiGridSweep;
@synthesize alpha = _alpha;
@synthesize beta = _beta;
@synthesize dt = _dt;
@synthesize multigridSolution = _multigridSolution;
@synthesize multigridEqualPlit = _multigridEqualPlit;
@synthesize plugInPrincipalClassInstance = _plugInPrincipalClassInstance;
@synthesize matrix = _matrix;
@synthesize variable = _variable;
@synthesize mesh = _mesh;
@synthesize solutionInfo = _solutionInfo;
@synthesize normalTangentialName = _normalTangentialName;
@synthesize exportedVariables = _exportedVariables;
@synthesize valuesList = _valuesList;

#pragma mark Initializations

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _plugInPrincipalClassInstance = nil;
        
        
        _containers = (solutionArraysContainer*)malloc(sizeof(solutionArraysContainer));
        _containers->ntZeroingDone = NULL;
        _containers->activeElements = NULL;
        _containers->ntElement = NULL;
        _containers->defDofs = NULL;
        _containers->size1DefDofs = 0;
        _containers->size2DefDofs = 0;
    }
    
    return self;
}

-(void)deallocation {
    free(_containers);
}

-(solutionArraysContainer*)getContainers {
    
    return _containers;
}

-(void)initializeILU:(int)ilun {
    
    int i, n, m;
    
    FEMMatrix *a1;
    matrixArraysContainer *matContainers = NULL;
    matrixArraysContainer *a1Containers = NULL;
    
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

        matContainers->sizeILUCols = [self FEMSolution_initializeILU1:matContainers :n];
        matContainers->sizeILURows = n+1;
        matContainers->sizeILUDiag = n;
        
        if (ilun > 1) {
            
            for (i=0; i<ilun-1; i++) {
                
                a1Containers->Cols = matContainers->ILUCols;
                a1Containers->Rows = matContainers->ILURows;
                a1Containers->Diag = matContainers->ILUDiag;
                
                m = [self FEMSolution_initializeILU1:a1Containers :n];
                
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
}

-(void)initializeCILU:(int)ilun {
    
    int i, j, k;
    int n;
    
    FEMMatrix *a1;
    matrixArraysContainer *matContainers = NULL;
    matrixArraysContainer *a1Containers = NULL;
    
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
        
        matContainers->sizeILUCols = [self FEMSolution_initializeILU1:a1Containers :n/2];
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
                
                k = [self FEMSolution_initializeILU1:a1Containers :n/2];
                
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
}

-(void)ilutWorkspaceCheckAtIndex:(int)i numberOfRows:(int)n {
    
    int j, k;
    int *iWork = NULL;
    double *cWork = NULL;
    matrixArraysContainer *matContainers = NULL;
    
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
}

-(void)ilutComplexWorkspaceCheckAtIndex:(int)i numberOfRows:(int)n {
    
    int j, k;
    int *iWork = NULL;
    double complex *cWork = NULL;
    matrixArraysContainer *matContainers = NULL;
    
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
}

@end
