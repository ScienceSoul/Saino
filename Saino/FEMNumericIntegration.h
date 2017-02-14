//===----------------------------------------------------------------------===//
//  FEMSolver.h
//  Saino
//
//  Created by Seddik hakime on 17/06/11.
//  Copyright (c) 2011 ScienceSoul. All rights reserved.
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
#import "GaussIntegration.h"

@interface FEMNumericIntegration : NSObject {

    double _metricDeterminant;
    double *_basis;
    double **_basisFirstDerivative;
    double ***_basisSecondDerivative;
}

@property(nonatomic, assign) double metricDeterminant;
@property(nonatomic, assign) double * __nullable basis;
@property(nonatomic, assign) double * __nullable * __nullable basisFirstDerivative;
@property(nonatomic, assign) double * __nullable * __nullable * __nullable basisSecondDerivative;

-(BOOL)allocation:(FEMMesh * __nonnull)mesh;
-(void)deallocation:(FEMMesh * __nonnull)mesh;

-(BOOL)setBasisForElement:(Element_t * __nonnull)element elementNodes:(Nodes_t * __nonnull)nodes inMesh:(FEMMesh * __nonnull)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int * __nullable)degree;
-(BOOL)setBasisFirstDerivativeForElement:(Element_t * __nonnull)element elementNodes:(Nodes_t * __nonnull)nodes inMesh:(FEMMesh * __nonnull)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int * __nullable)degree;
-(BOOL)setBasisSecondDerivativeForElement:(Element_t * __nonnull)element elementNodes:(Nodes_t * __nonnull)nodes inMesh:(FEMMesh * __nonnull)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int * __nullable)degree;
-(BOOL)setMetricDeterminantForElement:(Element_t * __nonnull)element elementNodes:(Nodes_t * __nonnull)nodes inMesh:(FEMMesh * __nonnull)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w;

@end
