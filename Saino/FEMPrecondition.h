//===----------------------------------------------------------------------===//
//  FEMPreconditioners.h
//  Saino
//
//  Created by Hakime Seddik on 27/01/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
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

#import <Foundation/Foundation.h>
#import "FEMMatrix.h"
#import "TimeProfile.h"

@interface FEMPrecondition : NSObject

// Diagonal preconditioning
-(void)CRSDiagPreconditionMatrix:(FEMMatrix * _Nonnull)matrix afterPrecondition:(double * _Nonnull)u rightHandSide:(double * _Nonnull)v info:(int * _Nullable)ipar;
-(void)CRSComplexDiagPreconditionMatrix:(FEMMatrix * _Nonnull)matrix afterPrecondition:(double complex * _Nonnull)u rightHandSide:(double complex * _Nonnull)v info:(int * _Nullable)ipar;

-(void)CRSBlockDiagonalMatrix:(FEMMatrix * _Nonnull)matrix blockDiagMatrix:(FEMMatrix * _Nonnull)B numberOfBlocks:(int)blocks;

// ILU(n) preconditioning
-(void)initializeILUMatrix:(FEMMatrix * _Nonnull)matrix numberOfRows:(int)ilun;
-(void)initializeCILUMatrix:(FEMMatrix * _Nonnull)matrix :(int)ilun;
-(BOOL)CRSIncompleteLUMatrix:(FEMMatrix * _Nonnull)matrix fillsOrder:(int)ilun;
-(BOOL)CRSComplexIncompleteLUMatrix:(FEMMatrix * _Nonnull)matrix fillsOrder: (int)ilun;

// ILU(T) preconditioning
-(BOOL)CRSIlutMatrix:(FEMMatrix * _Nonnull)matrix dropTolerance:(int)tol;
-(BOOL)CRSComplexIlutMatrix:(FEMMatrix * _Nonnull)matrix dropTolerance:(int)tol;

// LU Solve
-(BOOL)CRSLuPreconditionMatrix:(FEMMatrix * _Nonnull)matrix afterPrecondition:(double * _Nonnull)u rightHandSide:(double * _Nonnull)v info:(int * _Nonnull)ipar;
-(BOOL)CRSComplexLuPreconditionMatrix:(FEMMatrix * _Nonnull)matrix afterPrecondition:(double complex * _Nonnull)u rightHandSide:(double complex * _Nonnull)v info:(int * _Nonnull)ipar;

// Matrix-vector product
-(void)CRSMatrixVectorProduct:(FEMMatrix * _Nonnull)matrix vector:(double * _Nonnull)u result:(double * _Nonnull)v info:(int * _Nonnull)ipar;
-(void)CRSComplexMatrixVectorProduct:(FEMMatrix * _Nonnull)matrix vector:(double complex * _Nonnull)u result:(double complex * _Nonnull)v info:(int * _Nonnull)ipar;

// Dummy method when preconditioning is not needed
-(void)CRSPCondDummyMatrix:(FEMMatrix * _Nonnull)matrix afterPrecondition:(double * _Nonnull)u rightHandSide:(double * _Nonnull)v info:(int * _Nonnull)ipar;
-(void)CRSPCondDummyComplexMatrix:(FEMMatrix * _Nonnull)matrix afterPrecondition:(double complex * _Nonnull)u rightHandSide:(double complex * _Nonnull)v info:(int * _Nonnull)ipar;

@end

