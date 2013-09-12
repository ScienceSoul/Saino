//
//  FEMPreconditioners.h
//  Saino
//
//  Created by Hakime Seddik on 27/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMatrix.h"
#import "TimeProfile.h"

@interface FEMPrecondition : NSObject

// Diagonal preconditioning
-(void)CRSDiagPreconditionMatrix:(FEMMatrix *)matrix afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar;
-(void)CRSComplexDiagPreconditionMatrix:(FEMMatrix *)matrix afterPrecondition:(double complex *)u rightHandSide:(double complex *)v info:(int *)ipar;

-(void)CRSBlockDiagonalMatrix:(FEMMatrix *)matrix blockDiagMatrix:(FEMMatrix *)B numberOfBlocks:(int)blocks;

// ILU(n) preconditioning
-(void)initializeILUMatrix:(FEMMatrix *)matrix numberOfRows:(int)ilun;
-(void)initializeCILUMatrix:(FEMMatrix *)matrix :(int)ilun;
-(BOOL)CRSIncompleteLUMatrix:(FEMMatrix *)matrix fillsOrder:(int)ilun;
-(BOOL)CRSComplexIncompleteLUMatrix:(FEMMatrix *)matrix fillsOrder: (int)ilun;

// ILU(T) preconditioning
-(BOOL)CRSIlutMatrix:(FEMMatrix *)matrix dropTolerance:(int)tol;
-(BOOL)CRSComplexIlutMatrix:(FEMMatrix *)matrix dropTolerance:(int)tol;

// LU Solve
-(BOOL)CRSLuPreconditionMatrix:(FEMMatrix *)matrix afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar;
-(BOOL)CRSComplexLuPreconditionMatrix:(FEMMatrix *)matrix afterPrecondition:(double complex *)u rightHandSide:(double complex *)v info:(int *)ipar;

// Matrix-vector product
-(void)CRSMatrixVectorProduct:(FEMMatrix *)matrix vector:(double *)u result:(double *)v info:(int *)ipar;
-(void)CRSComplexMatrixVectorProduct:(FEMMatrix *)matrix vector:(double complex *)u result:(double complex *)v info:(int *)ipar;

// Dummy method when preconditioning is not needed
-(void)CRSPCondDummyMatrix:(FEMMatrix *)matrix afterPrecondition:(double *)u rightHandSide:(double *)v info:(int *)ipar;
-(void)CRSPCondDummyComplexMatrix:(FEMMatrix *)matrix afterPrecondition:(double complex *)u rightHandSide:(double complex *)v info:(int *)ipar;

@end

