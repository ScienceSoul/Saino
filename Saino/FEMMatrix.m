//===----------------------------------------------------------------------===//
//  FEMMatrix.m
//  Saino
//
//  Created by Hakime Seddik on 23/07/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
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
        _containers->sizePerm = 0;
        _containers->sizeInvPerm = 0;
        _containers->sizeRowOwner = 0;
        _containers->sizeGRows = 0;
        _containers->sizeGOrder = 0;
        _containers->sizeRows = 0;
        _containers->sizeCols = 0;
        _containers->sizeDiag = 0;
        _containers->sizeRHS = 0;
        _containers->sizeBulkRHS = 0;
        _containers->sizeRHS_im = 0;
        _containers->size1force = 0;
        _containers->size2Force = 0;
        _containers->sizeValues = 0;
        _containers->sizeILUValues = 0;
        _containers->sizeMassValues = 0;
        _containers->sizeDampValues = 0;
        _containers->sizeBulkValues = 0;
        _containers->sizeDiagScaling = 0;
        _containers->sizeILURows = 0;
        _containers->sizeILUCols = 0;
        _containers->sizeILUDiag = 0;
        _containers->sizeCRHS = 0;
        _containers->sizeCForce = 0;
        _containers->sizeCValues = 0;
        _containers->sizeCILUValues = 0;
        _containers->sizeCMassValues = 0;
        _containers->sizeCDampValues = 0;
        _containers->sizeFct = 0;
        _containers->sizeCDampValues = 0;
        
        _containers->ListMatrix = NULL;
        
        _containers->Perm = NULL;
        _containers->InvPerm = NULL;
        _containers->RowOwner = NULL;
        _containers->GRows = NULL;
        _containers->GOrder = NULL;
        
        _containers->Rows = NULL;
        _containers->Cols = NULL;
        _containers->Diag = NULL;
        
        _containers->RHS = NULL;
        _containers->BulkRHS = NULL;
        _containers->RHS_im = NULL;
        _containers->Force = NULL;
        
        _containers->Values = NULL;
        _containers->ILUValues = NULL;
        _containers->MassValues = NULL;
        _containers->DampValues = NULL;
        _containers->BulkValues = NULL;
        _containers->DiagScaling = NULL;
        
        _containers->ILURows = NULL;
        _containers->ILUCols = NULL;
        _containers->ILUDiag = NULL;
        
        _containers->CRHS = NULL;
        _containers->CForce = NULL;
        _containers->CValues = NULL;
        _containers->CILUValues = NULL;
        _containers->CMassValues = NULL;
        _containers->CDampValues = NULL;
        
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

-(matrixArraysContainer * _Nonnull)getContainers {
    
    return _containers;
}

// Retruns an ummutable copy of a FEMMatrix class

-(FEMMatrix * _Nonnull)copyWithZone:(NSZone *)zone {
    
    FEMMatrix *matrix = [[FEMMatrix allocWithZone:zone] init];
    
    // We don't try to copy these objects from the original matrix.
    // Instead just assign them to nil since anyway they were only
    // referenced in the original matrix.
    matrix.child = nil;
    matrix.parent = nil;
    matrix.constraint = nil;
    matrix.ematrix = nil;
    
    matrix.numberOfRows = self.numberOfRows;
    matrix.subband = self.subband;
    matrix.format = self.format;
    matrix.solveCount = self.solveCount;
    matrix.comm = self.comm;
    matrix.umfPackNumeric = self.umfPackNumeric;
    matrix.projectorBC = self.projectorBC;
    matrix.projectorType = self.projectorType;
    matrix.rhsScaling = self.rhsScaling;
    matrix.cholesky = self.cholesky;
    matrix.ordered = self.ordered;
    matrix.lumped = self.lumped;
    matrix.symmetric = self.symmetric;
    matrix.complexMatrix = self.complexMatrix;
    matrix.dgMatrix = self.dgMatrix;
    
    matrixArraysContainer *matrixContainers = matrix.getContainers;
    
    if (_containers->ListMatrix != NULL) {
        matrixContainers->ListMatrix = (ListMatrix_t*)malloc(sizeof(ListMatrix_t));
        memcpy(matrixContainers->ListMatrix, _containers->ListMatrix, sizeof(ListMatrix_t));
    }
    
    if (_containers->Perm != NULL) {
        matrixContainers->Perm = intvec(0, _containers->sizePerm-1);
        matrixContainers->sizePerm = _containers->sizePerm;
        memcpy(matrixContainers->Perm, _containers->Perm, _containers->sizePerm*sizeof(int));
    }
    
    if (_containers->InvPerm != NULL) {
        matrixContainers->InvPerm = intvec(0, _containers->sizeInvPerm-1);
        matrixContainers->sizeInvPerm = _containers->sizeInvPerm;
        memcpy(matrixContainers->InvPerm, _containers->InvPerm, _containers->sizeInvPerm*sizeof(int));
    }
    
    if (_containers->RowOwner != NULL) {
        matrixContainers->RowOwner = intvec(0, _containers->sizeRowOwner-1);
        matrixContainers->sizeRowOwner = _containers->sizeRowOwner;
        memcpy(matrixContainers->RowOwner, _containers->RowOwner, _containers->sizeRowOwner*sizeof(int));
    }
    
    if (_containers->GRows != NULL) {
        matrixContainers->GRows = intvec(0, _containers->sizeGRows-1);
        matrixContainers->sizeGRows = _containers->sizeGRows;
        memcpy(matrixContainers->GRows, _containers->GRows, _containers->sizeGRows*sizeof(int));
    }
    
    if (_containers->GOrder != NULL) {
        matrixContainers->GOrder = intvec(0, _containers->sizeGOrder-1);
        matrixContainers->sizeGOrder = _containers->sizeGOrder;
        memcpy(matrixContainers->GOrder, _containers->GOrder, _containers->sizeGOrder*sizeof(int));
    }
    
    if (_containers->Rows != NULL) {
        matrixContainers->Rows = intvec(0, _containers->sizeRows-1);
        matrixContainers->sizeRows = _containers->sizeRows;
        memcpy(matrixContainers->Rows, _containers->Rows, _containers->sizeRows*sizeof(int));
    }
    
    if (_containers->Cols != NULL) {
        matrixContainers->Cols = intvec(0, _containers->sizeCols-1);
        matrixContainers->sizeCols = _containers->sizeCols;
        memcpy(matrixContainers->Cols, _containers->Cols, _containers->sizeCols*sizeof(int));
    }
    
    if (_containers->Diag != NULL) {
        matrixContainers->Diag = intvec(0, _containers->sizeDiag-1);
        matrixContainers->sizeDiag = _containers->sizeDiag;
        memcpy(matrixContainers->Diag, _containers->Diag, _containers->sizeDiag*sizeof(int));
    }
    
    if (_containers->RHS != NULL) {
        matrixContainers->RHS = doublevec(0, _containers->sizeRHS-1);
        matrixContainers->sizeRHS = _containers->sizeRHS;
        memcpy(matrixContainers->RHS, _containers->RHS, _containers->sizeRHS*sizeof(double));
    }
    
    if (_containers->BulkRHS != NULL) {
        matrixContainers->BulkRHS = doublevec(0, _containers->sizeBulkRHS-1);
        matrixContainers->sizeBulkRHS = _containers->sizeBulkRHS;
        memcpy(matrixContainers->BulkRHS, _containers->BulkRHS, _containers->sizeBulkRHS*sizeof(double));
    }
    
    if (_containers->RHS_im != NULL) {
        matrixContainers->RHS_im = doublevec(0, _containers->sizeRHS_im-1);
        matrixContainers->sizeRHS_im = _containers->sizeRHS_im;
        memcpy(matrixContainers->RHS_im, _containers->RHS_im, _containers->sizeRHS_im*sizeof(double));
    }
    
    if (_containers->Force != NULL) {
        matrixContainers->Force = doublematrix(0, _containers->size1force-1, 0, _containers->size2Force-1);
        matrixContainers->size1force = _containers->size1force;
        matrixContainers->size2Force = _containers->size2Force;
        memcpy(*matrixContainers->Force, *_containers->Force, (_containers->size1force*_containers->size2Force)*sizeof(double));
    }
    
    if (_containers->Values != NULL) {
        matrixContainers->Values = doublevec(0, _containers->sizeValues-1);
        matrixContainers->sizeValues = _containers->sizeValues;
        memcpy(matrixContainers->Values, _containers->Values, _containers->sizeValues*sizeof(double));
    }
    
    if (_containers->ILUValues != NULL) {
        matrixContainers->ILUValues = doublevec(0, _containers->sizeILUValues-1);
        matrixContainers->sizeILUValues = _containers->sizeILUValues;
        memcpy(matrixContainers->ILUValues, _containers->ILUValues, _containers->sizeILUValues*sizeof(double));
    }
    
    if (_containers->MassValues != NULL) {
        matrixContainers->MassValues = doublevec(0, _containers->sizeMassValues-1);
        matrixContainers->sizeMassValues = _containers->sizeMassValues;
        memcpy(matrixContainers->MassValues, _containers->MassValues, _containers->sizeMassValues*sizeof(double));
    }
    
    if (_containers->DampValues != NULL) {
        matrixContainers->DampValues = doublevec(0, _containers->sizeDampValues-1);
        matrixContainers->sizeDampValues = _containers->sizeDampValues;
        memcpy(matrixContainers->DampValues, _containers->DampValues, _containers->sizeDampValues*sizeof(double));
    }
    
    if (_containers->BulkValues != NULL) {
        matrixContainers->BulkValues = doublevec(0, _containers->sizeBulkValues-1);
        matrixContainers->sizeBulkValues = _containers->sizeBulkValues;
        memcpy(matrixContainers->BulkValues, _containers->BulkValues, _containers->sizeBulkValues*sizeof(double));
    }
    
    if (_containers->DiagScaling != NULL) {
        matrixContainers->DiagScaling = doublevec(0, _containers->sizeDiagScaling-1);
        matrixContainers->sizeDiagScaling = _containers->sizeDiagScaling;
        memcpy(matrixContainers->DiagScaling, _containers->DiagScaling, _containers->sizeDiagScaling*sizeof(double));
    }
    
    if (_containers->ILURows != NULL) {
        matrixContainers->ILURows = intvec(0, _containers->sizeILURows-1);
        matrixContainers->sizeILURows = _containers->sizeILURows;
        memcpy(matrixContainers->ILURows, _containers->ILURows, _containers->sizeILURows*sizeof(int));
    }
    
    if (_containers->ILUCols != NULL) {
        matrixContainers->ILUCols = intvec(0, _containers->sizeILUCols-1);
        matrixContainers->sizeILUCols = _containers->sizeILUCols;
        memcpy(matrixContainers->ILUCols, _containers->ILUCols, _containers->sizeILUCols*sizeof(int));
    }
    
    if (_containers->ILUDiag != NULL) {
        matrixContainers->ILUDiag = intvec(0, _containers->sizeILUDiag-1);
        matrixContainers->sizeILUDiag = _containers->sizeILUDiag;
        memcpy(matrixContainers->ILUDiag, _containers->ILUDiag, _containers->sizeILUDiag*sizeof(int));
    }
    
    if (_containers->CRHS != NULL) {
        matrixContainers->CRHS = cdoublevec(0, _containers->sizeCRHS-1);
        matrixContainers->sizeCRHS = _containers->sizeCRHS;
        memcpy(matrixContainers->CRHS, _containers->CRHS, _containers->sizeCRHS*sizeof(double complex));
    }
    
    if (_containers->CForce != NULL) {
        matrixContainers->CForce = cdoublevec(0, _containers->sizeCForce-1);
        matrixContainers->sizeCForce = _containers->sizeCForce;
        memcpy(matrixContainers->CForce, _containers->CForce, _containers->sizeCForce*sizeof(double complex));
    }
    
    if (_containers->CValues != NULL) {
        matrixContainers->CValues = cdoublevec(0, _containers->sizeCValues-1);
        matrixContainers->sizeCValues = _containers->sizeCValues;
        memcpy(matrixContainers->CValues, _containers->CValues, _containers->sizeCValues*sizeof(double complex));
    }
    
    if (_containers->CILUValues != NULL) {
        matrixContainers->CILUValues = cdoublevec(0, _containers->sizeCILUValues-1);
        matrixContainers->sizeCILUValues = _containers->sizeCILUValues;
        memcpy(matrixContainers->CILUValues, _containers->CILUValues, _containers->sizeCILUValues*sizeof(double complex));
    }
    
    if (_containers->CMassValues != NULL) {
        matrixContainers->CMassValues = cdoublevec(0, _containers->sizeCMassValues-1);
        matrixContainers->sizeCMassValues = _containers->sizeCMassValues;
        memcpy(matrixContainers->CMassValues, _containers->CMassValues, _containers->sizeCMassValues*sizeof(double complex));
    }
    
    if (_containers->CDampValues != NULL) {
        matrixContainers->CDampValues = cdoublevec(0, _containers->sizeCDampValues-1);
        matrixContainers->sizeCDampValues = _containers->sizeCDampValues;
        memcpy(matrixContainers->CDampValues, _containers->CDampValues, _containers->sizeCDampValues*sizeof(double complex));
    }
    
    if (_containers->FCT_D != NULL) {
        matrixContainers->FCT_D = doublevec(0, _containers->sizeFct-1);
        matrixContainers->sizeFct = _containers->sizeFct;
        memcpy(matrixContainers->FCT_D, _containers->FCT_D, _containers->sizeFct*sizeof(double));
    }
    
    if (_containers->MassValuesLumped != NULL) {
        matrixContainers->MassValuesLumped = doublevec(0, _containers->sizeMassValuesLumped-1);
        matrixContainers->sizeMassValuesLumped = _containers->sizeMassValuesLumped;
        memcpy(matrixContainers->MassValuesLumped, _containers->MassValuesLumped, _containers->sizeMassValuesLumped*sizeof(double));
    }
    
    return matrix;
}

@end
