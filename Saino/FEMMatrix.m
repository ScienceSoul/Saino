//
//  FEMMatrix.m
//  Saino
//
//  Created by Hakime Seddik on 23/07/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMatrix.h"

@implementation FEMMatrix

@synthesize child = _child;
@synthesize parent = _parent;
@synthesize constraint = _constraint;
@synthesize ematrix = _ematrix;
@synthesize numberOfRows = _numberOfRows;
@synthesize subband = _subband;
@synthesize format = _format;
@synthesize solveCount = _solveCount;
@synthesize comm = _comm;
@synthesize umfPackNumeric = _umfPackNumeric;
@synthesize ordered = _ordered;
@synthesize lumped = _lumped;
@synthesize symmetric = _symmetric;
@synthesize complexMatrix = _complexMatrix;
@synthesize dgMatrix = _dgMatrix;

- (id)init
{
    self = [super init];
    if (self) {
        // NOTE: The classes below need to be assigned to allocated FEMMatrix classes when used
        _child = nil;
        _parent = nil;
        _constraint = nil; 
        _ematrix = nil;
        
        _format = MATRIX_CRS;
        _solveCount = 0;
        _numberOfRows = 0;
        _umfPackNumeric = 0;
        _lumped = NO;
        _ordered = NO;
        _complexMatrix = NO;
        _symmetric = NO;
        _dgMatrix = NO;
        
        _containers = (matrixArraysContainer*)malloc(sizeof(matrixArraysContainer) * 1 );
        _containers->Perm = NULL;
        _containers->InvPerm = NULL;
        
        _containers->Cols = NULL;
        _containers->Rows = NULL;
        _containers->Diag = NULL;
        _containers->GRows = NULL;
        
        _containers->RHS = NULL;
        _containers->Force = NULL;
        _containers->RHS_im = NULL;
        
        _containers->Values = NULL;
        _containers->ILUValues = NULL;
        _containers->MassValues = NULL;
        _containers->DampValues = NULL;
        
        _containers->BulkRHS = NULL;
        _containers->BulkValues = NULL;
        
        _containers->ILUCols = NULL;
        _containers->ILURows = NULL;
        _containers->ILUDiag = NULL;
        
        _containers->CRHS = NULL;
        _containers->CForce = NULL;
        
        _containers->RowOwner = NULL;
        
        _containers->CValues = NULL;
        _containers->CILUValues = NULL;
        _containers->CMassValues = NULL;
        _containers->CDampValues = NULL;
        
        _containers->GRows = NULL;
        _containers->GOrder = NULL;
    }
    
    return self;
}

-(void)deallocation {
    free(_containers);
}

-(matrixArraysContainer *)getContainers {
    
    return _containers;
}

@end
