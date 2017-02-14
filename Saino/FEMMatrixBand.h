//===----------------------------------------------------------------------===//
//  FEMMatrixBand.h
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

@interface FEMMatrixBand : NSObject

-(FEMMatrix * __nonnull)createMatrixWithNumberOfRows:(int)rows subBand:(int)subBand symmetric:(BOOL)symmetric allocateValues:(BOOL)allocateValues;
-(void)zeroRowInGlobal:(FEMSolution * __nonnull)solution numberOfRows:(int)n;
-(void)setElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)addToElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)glueLocalMatrix:(double * __nonnull * __nonnull)localMatrix inGlobal:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dofs:(int)dofs indexes:(int * __nonnull)indexes;
-(void)sBand_setDirichlet:(FEMSolution * __nonnull)solution orderedNumber:(int)n value:(double)value;
// Band Matrix-vector multiply
-(void)matrixVectorMultiplyInGlobal:(FEMSolution * __nonnull)solution vector:(double * __nonnull)u result:(double * __nonnull)v;
-(void)zeroRowInMatrix:(FEMMatrix * __nonnull)matrix numberOfRows:(int)n;
-(void)setElementInMatrix:(FEMMatrix * __nonnull)matrix row:(int)i col:(int)j value:(double)value;
-(void)matrixVectorMultiply:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)u result:(double * __nonnull)v;
@end
