//===----------------------------------------------------------------------===//
//  FEMTimeIntegration.h
//  Saino
//
//  Created by Hakime Seddik on 14/03/12.
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
#import "FEMMesh.h"

@interface FEMTimeIntegration : NSObject


-(void)fractionalStepInSolution:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * __nonnull * __nonnull )massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution:(double * __nonnull)prevSolution rows:(int * __nonnull)rows;
-(void)bdfLocalInSolution:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * __nonnull * __nonnull)massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution:(double * __nonnull * __nonnull)prevSolution order:(int)order rows:(int * __nonnull)rows cols:(int * __nonnull)cols;
-(void)vbdfLocalInSolution:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dts:(double * __nonnull)dts massMatrix:(double * __nonnull * __nonnull)massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution:(double * __nonnull * __nonnull)prevSolution order:(int)order rows:(int * __nonnull)rows cols:(int * __nonnull)cols;
-(void)newMarkBetaInSolution:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * __nonnull * __nonnull)massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution:(double * __nonnull)prevSolution beta:(double)beta rows:(int * __nonnull)rows cols:(int * __nonnull)cols;
-(void)bossakSecondOrder:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * __nonnull * __nonnull)massMatrix dampMatrix:(double * __nonnull * __nonnull)dampMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution1:(double * __nonnull)x prevSolution2:(double * __nonnull)v prevSolution3:(double * __nonnull)a alpha:(double)alpha rows:(int * __nonnull)rows cols:(int * __nonnull)cols;

@end
