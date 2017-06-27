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

-(FEMMatrix * _Nonnull)createMatrixWithNumberOfRows:(int)rows subBand:(int)subBand symmetric:(BOOL)symmetric allocateValues:(BOOL)allocateValues;
-(void)zeroRowInGlobal:(FEMSolution * _Nonnull)solution numberOfRows:(int)n;
-(void)setElementInGlobal:(FEMSolution * _Nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)addToElementInGlobal:(FEMSolution * _Nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)glueLocalMatrix:(double * _Nonnull * _Nonnull)localMatrix inGlobal:(FEMSolution * _Nonnull)solution numberOfNodes:(int)n dofs:(int)dofs indexes:(int * _Nonnull)indexes;
-(void)sBand_setDirichlet:(FEMSolution * _Nonnull)solution orderedNumber:(int)n value:(double)value;
// Band Matrix-vector multiply
-(void)matrixVectorMultiplyInGlobal:(FEMSolution * _Nonnull)solution vector:(double * _Nonnull)u result:(double * _Nonnull)v;
-(void)zeroRowInMatrix:(FEMMatrix * _Nonnull)matrix numberOfRows:(int)n;
-(void)setElementInMatrix:(FEMMatrix * _Nonnull)matrix row:(int)i col:(int)j value:(double)value;
-(void)matrixVectorMultiply:(FEMMatrix * _Nonnull)matrix vector:(double * _Nonnull)u result:(double * _Nonnull)v;
@end
