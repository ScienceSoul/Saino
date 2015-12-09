//
//  FEMMatrix.m
//  Saino
//
//  Created by Hakime Seddik on 23/07/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMatrix.h"

#import "memory.h"

@interface FEMMatrix ()
-(void)FEMMatrix_deallocateContainers:(matrixArraysContainer *)containers;
@end

@implementation FEMMatrix

-(void)FEMMatrix_deallocateContainers:(matrixArraysContainer *)containers {
    if (containers != NULL) {
        
        if (containers->Perm != NULL) {
            free_ivector(containers->Perm, 0, containers->sizePerm-1);
            containers->Perm = NULL;
        }
        
        if (containers->InvPerm != NULL) {
            free_ivector(containers->InvPerm, 0, containers->sizeInvPerm-1);
            containers->InvPerm = NULL;
        }
        
        if (containers->Cols != NULL) {
            if (containers->ILUCols == containers->Cols) containers->ILUCols = NULL;
            free_ivector(containers->Cols, 0, containers->sizeCols-1);
            containers->Cols = NULL;
        }
        
        if (containers->Rows != NULL) {
            if (containers->ILURows == containers->Rows) containers->ILURows = NULL;
            free_ivector(containers->Rows, 0, containers->sizeRows-1);
            containers->Rows = NULL;
        }
        
        if (containers->Diag != NULL) {
            if (containers->ILUDiag == containers->Diag) containers->ILUDiag = NULL;
            free_ivector(containers->Diag, 0, containers->sizeDiag-1);
            containers->Diag = NULL;
        }
        
        if (containers->RHS != NULL) {
            free_dvector(containers->RHS, 0, containers->sizeRHS-1);
            containers->RHS = NULL;
        }
        
        if (containers->Force != NULL) {
            free_dmatrix(containers->Force, 0, containers->size1force-1, 0, containers->size2Force-1);
            containers->Force = NULL;
        }
        
        if (containers->RHS_im != NULL) {
            free_dvector(containers->RHS_im, 0, containers->sizeRHS_im-1);
            containers->RHS_im = NULL;
        }
        
        if (containers->Values != NULL) {
            free_dvector(containers->Values, 0, containers->sizeValues-1);
            containers->Values = NULL;
        }
        
        if (containers->MassValues != NULL) {
            free_dvector(containers->MassValues, 0, containers->sizeMassValues-1);
            containers->MassValues = NULL;
        }
        
        if (containers->DampValues != NULL) {
            free_dvector(containers->DampValues, 0, containers->sizeDampValues-1);
            containers->DampValues = NULL;
        }
        
        if (containers->ILUValues != NULL) {
            free_dvector(containers->ILUValues, 0, containers->sizeILUValues-1);
            containers->ILUValues = NULL;
        }
        
        if (containers->ILUCols != NULL) {
            free_ivector(containers->ILUCols, 0, containers->sizeILUCols-1);
            containers->ILUCols = NULL;
        }
        
        if (containers->ILURows != NULL) {
            free_ivector(containers->ILURows, 0, containers->sizeILURows-1);
            containers->ILURows = NULL;
        }
        
        if (containers->ILUDiag != NULL) {
            free_ivector(containers->ILUDiag, 0, containers->sizeILUDiag-1);
            containers->ILUDiag = NULL;
        }
        
        if (containers->CRHS != NULL) {
            free_cdvector(containers->CRHS, 0, containers->sizeCRHS-1);
            containers->CRHS = NULL;
        }
        
        if (containers->CForce != NULL) {
            free_cdvector(containers->CForce, 0, containers->sizeCForce-1);
            containers->CForce = NULL;
        }
        
        if (containers->CValues != NULL) {
            free_cdvector(containers->CValues, 0, containers->sizeCValues-1);
            containers->CValues = NULL;
        }
        
        if (containers->CILUValues != NULL) {
            free_cdvector(containers->CILUValues, 0, containers->sizeCILUValues-1);
            containers->CILUValues = NULL;
        }
        
        if (containers->CMassValues != NULL) {
            free_cdvector(containers->CMassValues, 0, containers->sizeCMassValues-1);
            containers->CMassValues = NULL;
        }
        
        if (containers->CDampValues != NULL) {
            free_cdvector(containers->CDampValues, 0, containers->sizeCDampValues-1);
            containers->CDampValues = NULL;
        }
        
        if (containers->GRows != NULL) {
            free_ivector(containers->GRows, 0, containers->sizeGRows-1);
            containers->GRows = NULL;
        }
        
        if (containers->RowOwner != NULL) {
            free_ivector(containers->RowOwner, 0, containers->sizeRowOwner-1);
            containers->RowOwner = NULL;
        }
        
        if (containers->GOrder != NULL) {
            free_ivector(containers->GOrder, 0, containers->sizeGOrder-1);
            containers->GOrder = NULL;
        }
        
        free(containers);
        containers = NULL;
    }
}

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
@synthesize projectorBC = _projectorBC;
@synthesize projectorType = _projectorType;
@synthesize rhsScaling = _rhsScaling;
@synthesize cholesky = _cholesky;
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
        _projectorBC = 0;
        _projectorType = PROJECTOR_TYPE_DEFAULT;
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
        _containers->DiagScaling = NULL;
        
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
        
        _containers->FCT_D = NULL;
        _containers->MassValuesLumped = NULL;
    }
    
    return self;
}

-(void)deallocation {
    
    // TODO: Specific deallocations relative to direct solvers
    
    [self FEMMatrix_deallocateContainers:_containers];
    
    if (_ematrix != nil) {
        matrixArraysContainer *containers = _ematrix.getContainers;
        [self FEMMatrix_deallocateContainers:containers];
        _ematrix = nil;
    }
    if (_constraint != nil) {
         matrixArraysContainer *containers = _constraint.getContainers;
        [self FEMMatrix_deallocateContainers:containers];
        _constraint = nil;
    }
}

-(matrixArraysContainer *)getContainers {
    
    return _containers;
}

@end
