//
//  FEMSolver.h
//  Saino
//
//  Created by Seddik hakime on 17/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

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
@property(nonatomic, assign) double *basis;
@property(nonatomic, assign) double **basisFirstDerivative;
@property(nonatomic, assign) double ***basisSecondDerivative;

-(BOOL)allocation:(FEMMesh *)mesh;
-(void)deallocation:(FEMMesh *)mesh;

-(BOOL)setBasisForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree;
-(BOOL)setBasisFirstDerivativeForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree;
-(BOOL)setBasisSecondDerivativeForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree;
-(BOOL)setMetricDeterminantForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w;

@end
