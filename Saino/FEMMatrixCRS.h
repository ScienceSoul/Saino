//===----------------------------------------------------------------------===//
//  FEMMatrixCRS.h
//  Saino
//
//  Created by Hakime Seddik on 23/03/12.
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
#import "FEMSolution.h"

@interface FEMMatrixCRS : NSObject

-(FEMMatrix * _Nonnull)createMatrixWithNumberOfRows:(int)rows totalNonZeros:(int)totalNonZeros rowNonZeros:(int * _Nonnull)rowNonZeros degreesFreedom:(int)degreesFreedom reorder:(int * _Nonnull)reorder sizeOfReorder:(int)sizeOfReorder allocateValues:(BOOL)allocateValues;
-(void)zeroRowInGlobal:(FEMSolution * _Nonnull)solution numberOfRows:(int)n;
-(void)sortGlobal:(FEMSolution * _Nonnull)solution alsoValues:(BOOL * _Nullable)alsoValues;
-(void)setElementInGlobal:(FEMSolution * _Nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)addToElementInGlobal:(FEMSolution * _Nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)glueLocalMatrix:(double * _Nonnull * _Nonnull)localMatrix inGlobal:(FEMSolution * _Nonnull)solution numberOfNodes:(int)n dofs:(int)dofs indexes:(int * _Nonnull)indexes;
-(void)setSymmetricDirichletInGlobal:(FEMSolution * _Nonnull)solution atIndex:(int)n value:(double)value;
// CRS Matrix-vector multiply
-(void)matrixVectorMultiplyInGlobal:(FEMSolution * _Nonnull)solution vector:(double * _Nonnull)u result:(double * _Nonnull)v;
-(void)complexMatrixVectorMultiplyInGlobal:(FEMSolution * _Nonnull)solution vector:(double complex * _Nonnull)u result:(double complex * _Nonnull)v;
-(void)fctlLowOrderInSolution:(FEMSolution * _Nullable)solution orMatrix:(FEMMatrix * _Nullable)matrix;

-(void)glueLocalMatrix:(double * _Nonnull * _Nonnull)localMatrix inMatrix:(FEMMatrix * _Nonnull)matrix numberOfNodes:(int)numberOfNodes dofs:(int)dofs indexes:(int * _Nonnull)indexes;
-(void)makeMatrixIndex:(FEMMatrix * _Nonnull)matrix row:(int)i col:(int)j;
-(void)zeroRowInMatrix:(FEMMatrix * _Nonnull)matrix numberOfRows:(int)n;
-(void)sortMatrix:(FEMMatrix * _Nonnull)matrix alsoValues:(BOOL * _Nullable)alsoValues;
-(void)setElementInMatrix:(FEMMatrix * _Nonnull)matrix row:(int)i col:(int)j value:(double)value;
-(void)applyProjector:(FEMMatrix * _Nonnull)pMatrix values:(double * _Nonnull)u permutation:(int * _Nullable)uperm values:(double * _Nonnull)v permutation:(int * _Nullable)vperm transpose:(BOOL * _Nullable)trans;
-(void)matrixVectorMultiply:(FEMMatrix * _Nonnull)matrix vector:(double * _Nonnull)u result:(double * _Nonnull)v;
-(void)complexMatrixVectorMultiply:(FEMMatrix * _Nonnull)matrix vector:(double complex * _Nonnull)u result:(double complex * _Nonnull)v;
-(void)zeroMatrix:(FEMMatrix * _Nonnull)matrix;


@end
