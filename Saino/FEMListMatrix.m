//===----------------------------------------------------------------------===//
//  FEMListMatrix.m
//  Saino
//
//  Created by Seddik hakime on 07/01/13.
//  Copyright (c) 2013 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

#import "FEMListMatrix.h"
#import "FEMMatrixCRS.h"
#import "Utils.h"

@implementation FEMListMatrix

@synthesize listMatrixGrowth = _listMatrixGrowth;

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _listMatrixGrowth = 1000;
    }
    
    return self;
}

-(ListMatrix_t * __nonnull)allocateMatrix:(int)n {
    
    int i;
    ListMatrix_t *matrix = NULL;
    
    matrix = (ListMatrix_t*)malloc(sizeof(ListMatrix_t) * n );
    for (i=0; i<n; i++) {
        matrix[i].Head = NULL;
        matrix[i].Level = 0;
        matrix[i].Degree = 0;
        matrix[i].sizeOfContainer = n;
    }
    
    return matrix;
}

-(void)freeMatrix:(ListMatrix_t * __nonnull)list size:(int)n {
    
    int i;
    ListMatrixEntry_t *p = NULL, *p1 = NULL;
    
    if (list == NULL) return;
    
    for (i=0; i<n; i++) {
        p = list[i].Head;
        while (p != NULL) {
            p1 = p->Next;
            free(p);
            p = p1;
        }
    }
    free(list);
}

-(ListMatrix_t * __nonnull)enlargeMatrix:(ListMatrix_t * __nonnull)matrix toSize:(int)n {
    
    int i;
    ListMatrix_t *newMatrix = NULL;
    
    newMatrix = [self allocateMatrix:n];
    if (matrix != NULL) {
        for (i=0; i<matrix->sizeOfContainer; i++) {
            newMatrix[i] = matrix[i];
        }
    }
    free(matrix);
    
    return newMatrix;
}

-(ListMatrixEntry_t * __nonnull)getMatrixIndexInListMatrix:(ListMatrix_t * __nullable)list atIndex:(int)k1 andIndex:(int)k2 {
    
    ListMatrixEntry_t *cList = NULL, *entry = NULL, *prev = NULL;
    
    if (list == NULL) list = [self allocateMatrix:k1+1];
    
    if ((k1+1) > list->sizeOfContainer) {
        list = [self enlargeMatrix:list toSize:max((k1+1), (list->sizeOfContainer+self.listMatrixGrowth))];
    }
    
    cList = list[k1].Head;
    
    if (cList == NULL) {
        entry = (ListMatrixEntry_t*)malloc(sizeof(ListMatrixEntry_t));
        entry->Value = 0.0;
        entry->Index = k2;
        entry->Next = NULL;
        list[k1].Degree = 1;
        list[k1].Head = entry;
        return entry;
    }
    
    prev = NULL;
    while (cList != NULL) {
        if (cList->Index >= k2) break;
        prev = cList;
        cList = cList->Next;
    }
    
    if (cList != NULL) {
        if (cList->Index == k2) {
            entry = cList;
            return entry;
        }
    }
    
    entry = (ListMatrixEntry_t*)malloc(sizeof(ListMatrixEntry_t));
    entry->Value = 0.0;
    entry->Index = k2;
    entry->Next = cList;
    if (prev != NULL) {
        prev->Next = entry;
    } else {
        list[k1].Head = entry;
    }
    
    list[k1].Degree = list[k1].Degree + 1;
    
    return entry;
}

/*************************************************************
 
    Method corresponds to Elmer from git on October 27 2015

*************************************************************/
-(void)addToMatrixElement:(ListMatrix_t * __nonnull)list atIndex:(int)k1 andIndex:(int)k2 value:(double)value setValue:(BOOL * __nullable)setValue {
    
    ListMatrixEntry_t *entry = NULL;
    BOOL set=NO;
    
    if (setValue != NULL) set = *setValue;
    
    entry = [self getMatrixIndexInListMatrix:list atIndex:k1 andIndex:k2];
    if (set == YES) {
        entry->Value = value;
    } else {
        entry->Value = entry->Value + value;
    }
}

/*************************************************************************************
 
    Transfer the flexible list matrix to the more efficient CRS matrix that is
    used in most places of the code. The matrix structure can accomodate both forms.
 
    Method corresponds to Elmer from git on October 27 2015
 
*************************************************************************************/
-(void)convertToCRSMatrix:(FEMMatrix * __nonnull)matrix{
    
    int n, *rows = NULL, *cols = NULL, *diag = NULL;
    double *values = NULL;
    ListMatrixEntry_t *p = NULL;
    
    if (matrix.format != MATRIX_LIST) {
        fprintf(stdout, "FEMListMatrix:convertToCRSMatrix: the initial matrix type is not a list matrix.\n");
        return;
    }
    matrixArraysContainer *containers = matrix.getContainers;
    ListMatrix_t *l = containers->ListMatrix;
    
    if (l == NULL) {
        matrix.format = MATRIX_CRS;
        matrix.numberOfRows = 0;
        return;
    }
    
    for (n=l->sizeOfContainer-1; n>=0; n--) {
        if (l[n].Degree > 0) break;
    }
    n++;
    rows = intvec(0, (n+1)-1);
    containers->sizeRows = n+1;
    diag = intvec(0, n-1);
    containers->sizeDiag = n;
    memset(diag, -1, n*sizeof(int));
    rows[0] = 0;
    for (int i=0; i<n; i++) {
        rows[i+1] = rows[i] + l[i].Degree;
    }
    
    fprintf(stdout, "FEMListMatrix:convertToCRSMatrix: number of entries in CRS matrix: %d.\n", rows[n]);
    
    cols = intvec(0, rows[n]-1);
    containers->sizeCols = rows[n];
    values = doublevec(0, rows[n]-1);
    containers->sizeValues = rows[n];
    
    int j = 0;
    for (int i=0; i<n; i++) {
        p = l[i].Head;
        while (p != NULL) {
            cols[j] = p->Index;
            values[j] = p->Value;
            j++;
            p = p->Next;
        }
    }
    
    matrix.numberOfRows = n;
    containers->Rows = rows;
    containers->Diag = diag;
    containers->Cols = cols;
    containers->Values = values;
    
    matrix.ordered = NO;
    FEMMatrixCRS *crsMatrix = [[FEMMatrixCRS alloc] init];
    [crsMatrix sortMatrix:matrix alsoValues:NULL];
    
    [self freeMatrix:l size:l->sizeOfContainer];
    containers->ListMatrix = NULL;
    matrix.format = MATRIX_CRS;
    fprintf(stdout, "FEMListMatrix:convertToCRSMatrix: matrix format changed from List to CRS.\n");
}

@end
