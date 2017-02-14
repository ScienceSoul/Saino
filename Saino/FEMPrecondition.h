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
-(void)CRSDiagPreconditionMatrix:(FEMMatrix * __nonnull)matrix afterPrecondition:(double * __nonnull)u rightHandSide:(double * __nonnull)v info:(int * __nullable)ipar;
-(void)CRSComplexDiagPreconditionMatrix:(FEMMatrix * __nonnull)matrix afterPrecondition:(double complex * __nonnull)u rightHandSide:(double complex * __nonnull)v info:(int * __nullable)ipar;

-(void)CRSBlockDiagonalMatrix:(FEMMatrix * __nonnull)matrix blockDiagMatrix:(FEMMatrix * __nonnull)B numberOfBlocks:(int)blocks;

// ILU(n) preconditioning
-(void)initializeILUMatrix:(FEMMatrix * __nonnull)matrix numberOfRows:(int)ilun;
-(void)initializeCILUMatrix:(FEMMatrix * __nonnull)matrix :(int)ilun;
-(BOOL)CRSIncompleteLUMatrix:(FEMMatrix * __nonnull)matrix fillsOrder:(int)ilun;
-(BOOL)CRSComplexIncompleteLUMatrix:(FEMMatrix * __nonnull)matrix fillsOrder: (int)ilun;

// ILU(T) preconditioning
-(BOOL)CRSIlutMatrix:(FEMMatrix * __nonnull)matrix dropTolerance:(int)tol;
-(BOOL)CRSComplexIlutMatrix:(FEMMatrix * __nonnull)matrix dropTolerance:(int)tol;

// LU Solve
-(BOOL)CRSLuPreconditionMatrix:(FEMMatrix * __nonnull)matrix afterPrecondition:(double * __nonnull)u rightHandSide:(double * __nonnull)v info:(int * __nonnull)ipar;
-(BOOL)CRSComplexLuPreconditionMatrix:(FEMMatrix * __nonnull)matrix afterPrecondition:(double complex * __nonnull)u rightHandSide:(double complex * __nonnull)v info:(int * __nonnull)ipar;

// Matrix-vector product
-(void)CRSMatrixVectorProduct:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)u result:(double * __nonnull)v info:(int * __nonnull)ipar;
-(void)CRSComplexMatrixVectorProduct:(FEMMatrix * __nonnull)matrix vector:(double complex * __nonnull)u result:(double complex * __nonnull)v info:(int * __nonnull)ipar;

// Dummy method when preconditioning is not needed
-(void)CRSPCondDummyMatrix:(FEMMatrix * __nonnull)matrix afterPrecondition:(double * __nonnull)u rightHandSide:(double * __nonnull)v info:(int * __nonnull)ipar;
-(void)CRSPCondDummyComplexMatrix:(FEMMatrix * __nonnull)matrix afterPrecondition:(double complex * __nonnull)u rightHandSide:(double complex * __nonnull)v info:(int * __nonnull)ipar;

@end

