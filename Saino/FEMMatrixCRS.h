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

-(FEMMatrix * __nonnull)createMatrixWithNumberOfRows:(int)rows totalNonZeros:(int)totalNonZeros rowNonZeros:(int * __nonnull)rowNonZeros degreesFreedom:(int)degreesFreedom reorder:(int * __nonnull)reorder sizeOfReorder:(int)sizeOfReorder allocateValues:(BOOL)allocateValues;
-(void)zeroRowInGlobal:(FEMSolution * __nonnull)solution numberOfRows:(int)n;
-(void)sortGlobal:(FEMSolution * __nonnull)solution alsoValues:(BOOL * __nullable)alsoValues;
-(void)setElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)addToElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)glueLocalMatrix:(double * __nonnull * __nonnull)localMatrix inGlobal:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dofs:(int)dofs indexes:(int * __nonnull)indexes;
-(void)setSymmetricDirichletInGlobal:(FEMSolution * __nonnull)solution atIndex:(int)n value:(double)value;
// CRS Matrix-vector multiply
-(void)matrixVectorMultiplyInGlobal:(FEMSolution * __nonnull)solution vector:(double * __nonnull)u result:(double * __nonnull)v;
-(void)complexMatrixVectorMultiplyInGlobal:(FEMSolution * __nonnull)solution vector:(double complex * __nonnull)u result:(double complex * __nonnull)v;
-(void)fctlLowOrderInSolution:(FEMSolution * __nullable)solution orMatrix:(FEMMatrix * __nullable)matrix;

-(void)glueLocalMatrix:(double * __nonnull * __nonnull)localMatrix inMatrix:(FEMMatrix * __nonnull)matrix numberOfNodes:(int)numberOfNodes dofs:(int)dofs indexes:(int * __nonnull)indexes;
-(void)makeMatrixIndex:(FEMMatrix * __nonnull)matrix row:(int)i col:(int)j;
-(void)zeroRowInMatrix:(FEMMatrix * __nonnull)matrix numberOfRows:(int)n;
-(void)sortMatrix:(FEMMatrix * __nonnull)matrix alsoValues:(BOOL * __nullable)alsoValues;
-(void)setElementInMatrix:(FEMMatrix * __nonnull)matrix row:(int)i col:(int)j value:(double)value;
-(void)applyProjector:(FEMMatrix * __nonnull)pMatrix values:(double * __nonnull)u permutation:(int * __nullable)uperm values:(double * __nonnull)v permutation:(int * __nullable)vperm transpose:(BOOL * __nullable)trans;
-(void)matrixVectorMultiply:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)u result:(double * __nonnull)v;
-(void)complexMatrixVectorMultiply:(FEMMatrix * __nonnull)matrix vector:(double complex * __nonnull)u result:(double complex * __nonnull)v;
-(void)zeroMatrix:(FEMMatrix * __nonnull)matrix;


@end
