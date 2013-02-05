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
-(void)CRSDiagPreconditionInSolution:(FEMSolution *)solution afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar;
-(void)CRSComplexDiagPreconditionInSolution:(FEMSolution *)solution afterPrecondition:(double complex *)u rightHandSide:(double complex *)v info:(int *)ipar;

-(void)CRSBlockDiagonalInSolution:(FEMSolution *)solution blockDiagMatrix:(FEMMatrix *) B numberOfBlocks:(int)blocks;

// ILU(n) preconditioning
-(BOOL)CRSIncompleteLUInSolution:(FEMSolution *)solution fillsOrder:(int)ilun;
-(BOOL)CRSComplexIncompleteLUInSolution:(FEMSolution *)solution fillsOrder: (int)ilun;

// ILU(T) preconditioning
-(BOOL)CRSIlutInSolution:(FEMSolution *)solution dropTolerance:(int)tol;
-(BOOL)CRSComplexIlutInSolution:(FEMSolution *)solution dropTolerance:(int)tol;

// LU Solve
-(BOOL)CRSLuPreconditionInSolution:(FEMSolution *)solution afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar;
-(BOOL)CRSComplexLuPreconditionInSolution:(FEMSolution *)solution afterPrecondition:(double complex *)u rightHandSide:(double complex *)v info:(int *)ipar;

// Matrix-vector product
-(void)CRSMatrixVectorProdInSolution:(FEMSolution *)solution multiplyVector:(double *)u resultVector:(double *)v info:(int *)ipar;
-(void)CRSComplexMatrixVectorProdInSolution:(FEMSolution *)solution multiplyVector:(double complex *)u resultVector:(double complex *)v info:(int *)ipar;

// Dummy method when preconditioning is not needed
-(void)CRSPCondDummyInSolution:(FEMSolution *)solution afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar;

// CRS Matrix-vector multiply
-(void)CRSMatrixVectorMultiplyInSolution:(FEMSolution *)solution multiplyVector:(double *)u resultVector:(double *)v;
-(void)CRSComplexMatrixVectorMultiplyInSolution:(FEMSolution *)solution multiplyVector:(double complex *)u resultVector:(double complex *)v;

@end

