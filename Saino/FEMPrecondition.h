//
//  FEMPreconditioners.h
//  Saino
//
//  Created by Hakime Seddik on 27/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMSolution.h"

#import "memory.h"
#import "Constructors.h"
#import "Utils.h"
#import "TimeProfile.h"


@interface FEMPrecondition : NSObject

// Diagonal preconditioning
-(void)CRS_DiagPrecondition:(FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar;
-(void)CRS_ComplexDiagPrecondition:(FEMSolution *)solution: (double complex *)u: (double complex *)v: (int *)ipar;

-(void)CRS_BlockDiagonal:(FEMSolution *)solution: (FEMMatrix *) B: (int)blocks;

// ILU(n) preconditioning
-(BOOL)CRS_IncompleteLU:(FEMSolution *)solution: (int)ilun;
-(BOOL)CRS_ComplexIncompleteLU:(FEMSolution *)solution: (int)ilun;

// ILU(T) preconditioning
-(BOOL)CRS_ILUT:(FEMSolution *)solution: (int)tol;
-(BOOL)CRS_ComplexILUT:(FEMSolution *)solution: (int)tol;

// LU Solve
-(BOOL)CRS_LUPrecondition: (FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar;
-(BOOL)CRS_ComplexLUPrecondition: (FEMSolution *)solution: (double complex *)u: (double complex *)v: (int *)ipar;

// Matrix-vector product
-(void)CRS_MatrixVectorProd:(FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar;
-(void)CRS_ComplexMatrixVectorProd:(FEMSolution *)solution: (double complex *)u: (double complex *)v: (int *)ipar;

// Dummy method when preconditioning is not needed
-(void)CRS_pcond_dummy:(FEMSolution *)solution: (double *)u: (double *)v: (int *)ipar;

// CRS Matrix-vector multiply
-(void)CRS_MatrixVectorMultiply:(FEMSolution *)solution: (double *)u: (double *)v;
-(void)CRS_ComplexMatrixVectorMultiply:(FEMSolution *)solution: (double complex *)u: (double complex *)v;

@end

