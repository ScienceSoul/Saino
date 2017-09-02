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


-(void)fractionalStepInSolution:(FEMSolution * _Nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * _Nonnull * _Nonnull )massMatrix stiffMatrix:(double * _Nonnull * _Nonnull)stiffMatrix force:(double * _Nonnull)force prevSolution:(double * _Nonnull)prevSolution rows:(int * _Nonnull)rows;
-(void)bdfLocalInSolution:(FEMSolution * _Nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * _Nonnull * _Nonnull)massMatrix stiffMatrix:(double * _Nonnull * _Nonnull)stiffMatrix force:(double * _Nonnull)force prevSolution:(double * _Nonnull * _Nonnull)prevSolution order:(int)order rows:(int * _Nonnull)rows cols:(int * _Nonnull)cols;
-(void)vbdfLocalInSolution:(FEMSolution * _Nonnull)solution numberOfNodes:(int)n dts:(double * _Nonnull)dts massMatrix:(double * _Nonnull * _Nonnull)massMatrix stiffMatrix:(double * _Nonnull * _Nonnull)stiffMatrix force:(double * _Nonnull)force prevSolution:(double * _Nonnull * _Nonnull)prevSolution order:(int)order rows:(int * _Nonnull)rows cols:(int * _Nonnull)cols;
-(void)newMarkBetaInSolution:(FEMSolution * _Nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * _Nonnull * _Nonnull)massMatrix stiffMatrix:(double * _Nonnull * _Nonnull)stiffMatrix force:(double * _Nonnull)force prevSolution:(double * _Nonnull)prevSolution beta:(double)beta rows:(int * _Nonnull)rows cols:(int * _Nonnull)cols;
-(void)bossakSecondOrder:(FEMSolution * _Nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * _Nonnull * _Nonnull)massMatrix dampMatrix:(double * _Nonnull * _Nonnull)dampMatrix stiffMatrix:(double * _Nonnull * _Nonnull)stiffMatrix force:(double * _Nonnull)force prevSolution1:(double * _Nonnull)x prevSolution2:(double * _Nonnull)v prevSolution3:(double * _Nonnull)a alpha:(double)alpha rows:(int * _Nonnull)rows cols:(int * _Nonnull)cols;

@end
