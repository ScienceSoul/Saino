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

-(FEMMatrix * __nonnull)createMatrixInModel:(FEMModel * __nonnull)model forSolution:(FEMSolution * __nonnull)solution mesh:(FEMMesh * __nonnull)mesh dofs:(int)dofs permutation:(int * __nonnull)perm sizeOfPermutation:(int)permSize matrixFormat:(int)matrixFormat optimizeBandwidth:(BOOL)optimizeBandwidth equationName:(NSString * __nullable)equation discontinuousGalerkinSolution:(BOOL * __nullable)dgSolution globalBubbles:(BOOL * __nullable)gbBubbles nodalDofsOnly:(BOOL * __nullable)nodalDofsOnly projectorDofs:(BOOL * __nullable)projectorDofs;
-(void)tangentDirectionsForNormal:(double * __nonnull)normal tangent1:(double * __nonnull)tangent1 tangent2:(double * __nonnull)tangent2;
-(double)elementArea:(Element_t * __nonnull)element numberOfNodes:(int)n mesh:(FEMMesh * __nonnull)mesh nodel:(FEMModel * __nonnull)model;

@end
