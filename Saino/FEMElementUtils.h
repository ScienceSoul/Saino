//===----------------------------------------------------------------------===//
//  FEMElementUtils.h
//  Saino
//
//  Created by Seddik hakime on 28/12/12.
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
#import "FEMMesh.h"
#import "FEMModel.h"
#import "FEMSolution.h"
#import "FEMMatrix.h"

@interface FEMElementUtils : NSObject

-(FEMMatrix * _Nonnull)createMatrixInModel:(FEMModel * _Nonnull)model forSolution:(FEMSolution * _Nonnull)solution mesh:(FEMMesh * _Nonnull)mesh dofs:(int)dofs permutation:(int * _Nonnull)perm sizeOfPermutation:(int)permSize matrixFormat:(int)matrixFormat optimizeBandwidth:(BOOL)optimizeBandwidth equationName:(NSString * _Nullable)equation discontinuousGalerkinSolution:(BOOL * _Nullable)dgSolution globalBubbles:(BOOL * _Nullable)gbBubbles nodalDofsOnly:(BOOL * _Nullable)nodalDofsOnly projectorDofs:(BOOL * _Nullable)projectorDofs;
-(void)tangentDirectionsForNormal:(double * _Nonnull)normal tangent1:(double * _Nonnull)tangent1 tangent2:(double * _Nonnull)tangent2;
-(double)elementArea:(Element_t * _Nonnull)element numberOfNodes:(int)n mesh:(FEMMesh * _Nonnull)mesh nodel:(FEMModel * _Nonnull)model;

@end
