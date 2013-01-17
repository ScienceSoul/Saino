//
//  FEMListMatrix.m
//  Saino
//
//  Created by Seddik hakime on 07/01/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMListMatrix.h"

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

-(ListMatrix_t *)allocateMatrix:(int)n {
    
    int i;
    ListMatrix_t *matrix;
    
    matrix = (ListMatrix_t*)malloc(sizeof(ListMatrix_t) * n );
    for (i=0; i<n; i++) {
        matrix[i].Head = NULL;
        matrix[i].Level = 0;
        matrix[i].Degree = 0;
        matrix[i].sizeOfContainer = n;
    }
    
    return matrix;
}

-(void)freeMatrix:(ListMatrix_t *)list size:(int)n {
    
    int i;
    ListMatrixEntry_t *p, *p1;
    
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

-(ListMatrix_t *)enlargeMatrix:(ListMatrix_t *)matrix toSize:(int)n {
    
    int i;
    ListMatrix_t *newMatrix;
    
    newMatrix = [self allocateMatrix:n];
    if (matrix != NULL) {
        for (i=0; i<matrix->sizeOfContainer; i++) {
            newMatrix[i] = matrix[i];
        }
    }
    free(matrix);
    
    return newMatrix;
}

-(ListMatrixEntry_t *)getMatrixIndexInListMatrix:(ListMatrix_t *)list atIndex:(int)k1 andIndex:(int)k2 {
    
    ListMatrixEntry_t *cList, *entry, *prev;
    
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

@end
