//
//  FEMMatrix.m
//  Saino
//
//  Created by Hakime Seddik on 23/07/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMatrix.h"

#import "memory.h"

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
        
        _containers = (matrixArraysContainer*)malloc(sizeof(matrixArraysContainer));
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
    
    // TODO: Specific deallocations relative to direct solvers
    
    if (_containers->Perm != NULL) {
        free_ivector(_containers->Perm, 0, _containers->sizePerm-1);
        _containers->Perm = NULL;
    }
    
    if (_containers->InvPerm != NULL) {
        free_ivector(_containers->InvPerm, 0, _containers->sizeInvPerm-1);
        _containers->InvPerm = NULL;
    }
    
    if (_containers->Cols != NULL) {
        if (_containers->ILUCols == _containers->Cols) _containers->ILUCols = NULL;
        free_ivector(_containers->Cols, 0, _containers->sizeCols-1);
        _containers->Cols = NULL;
    }
    
    if (_containers->Rows != NULL) {
        if (_containers->ILURows == _containers->Rows) _containers->ILURows = NULL;
        free_ivector(_containers->Rows, 0, _containers->sizeRows-1);
        _containers->Rows = NULL;
    }
    
    if (_containers->Diag != NULL) {
        if (_containers->ILUDiag == _containers->Diag) _containers->ILUDiag = NULL;
        free_ivector(_containers->Diag, 0, _containers->sizeDiag-1);
        _containers->Diag = NULL;
    }
    
    if (_containers->RHS != NULL) {
        free_dvector(_containers->RHS, 0, _containers->sizeRHS-1);
        _containers->RHS = NULL;
    }
    
    if (_containers->Force != NULL) {
        free_dmatrix(_containers->Force, 0, _containers->size1force-1, 0, _containers->size2Force-1);
        _containers->Force = NULL;
    }
    
    if (_containers->RHS_im != NULL) {
        free_dvector(_containers->RHS_im, 0, _containers->sizeRHS_im-1);
        _containers->RHS_im = NULL;
    }
    
    if (_containers->Values != NULL) {
        free_dvector(_containers->Values, 0, _containers->sizeValues-1);
        _containers->Values = NULL;
    }
    
    if (_containers->MassValues != NULL) {
        free_dvector(_containers->MassValues, 0, _containers->sizeMassValues-1);
        _containers->MassValues = NULL;
    }
    
    if (_containers->DampValues != NULL) {
        free_dvector(_containers->DampValues, 0, _containers->sizeDampValues-1);
        _containers->DampValues = NULL;
    }
    
    if (_containers->ILUValues != NULL) {
        free_dvector(_containers->ILUValues, 0, _containers->sizeILUValues-1);
        _containers->ILUValues = NULL;
    }
    
    if (_containers->ILUCols != NULL) {
        free_ivector(_containers->ILUCols, 0, _containers->sizeILUCols-1);
        _containers->ILUCols = NULL;
    }
    
    if (_containers->ILURows != NULL) {
        free_ivector(_containers->ILURows, 0, _containers->sizeILURows-1);
        _containers->ILURows = NULL;
    }
    
    if (_containers->ILUDiag != NULL) {
        free_ivector(_containers->ILUDiag, 0, _containers->sizeILUDiag-1);
        _containers->ILUDiag = NULL;
    }
    
    if (_containers->CRHS != NULL) {
        free_cdvector(_containers->CRHS, 0, _containers->sizeCRHS-1);
        _containers->CRHS = NULL;
    }
    
    if (_containers->CForce != NULL) {
        free_cdvector(_containers->CForce, 0, _containers->sizeCForce-1);
        _containers->CForce = NULL;
    }
    
    if (_containers->CValues != NULL) {
        free_cdvector(_containers->CValues, 0, _containers->sizeCValues-1);
        _containers->CValues = NULL;
    }
    
    if (_containers->CILUValues != NULL) {
        free_cdvector(_containers->CILUValues, 0, _containers->sizeCILUValues-1);
        _containers->CILUValues = NULL;
    }
    
    if (_containers->CMassValues != NULL) {
        free_cdvector(_containers->CMassValues, 0, _containers->sizeCMassValues-1);
        _containers->CMassValues = NULL;
    }
    
    if (_containers->CDampValues != NULL) {
        free_cdvector(_containers->CDampValues, 0, _containers->sizeCDampValues-1);
        _containers->CDampValues = NULL;
    }
    
    if (_containers->GRows != NULL) {
        free_ivector(_containers->GRows, 0, _containers->sizeGRows-1);
        _containers->GRows = NULL;
    }
    
    if (_containers->RowOwner != NULL) {
        free_ivector(_containers->RowOwner, 0, _containers->sizeRowOwner-1);
        _containers->RowOwner = NULL;
    }
    
    if (_containers->GOrder != NULL) {
        free_ivector(_containers->GOrder, 0, _containers->sizeGOrder-1);
        _containers->GOrder = NULL;
    }
    
    free(_containers);
}

-(matrixArraysContainer *)getContainers {
    
    return _containers;
}

@end
