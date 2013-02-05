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
#import "NodalBasisFunctions.h"

@interface FEMNumericIntegration : NSObject {

    GaussIntegrationPoints *_integCompound;
    double _metricDeterminant;
    double *_basis;
    double **_basisFirstDerivative;
    double ***_basisSecondDerivative;
    double **_elementMetric;
    double **_covariantMetrixTensor;
    double **_ltoGMap;
    double **_dx;                             // Partial derivatives of global coordinates with respect to local coordinates
}

@property(nonatomic, assign) GaussIntegrationPoints *integCompound;
@property(nonatomic, assign) double metricDeterminant;
@property(nonatomic, assign) double *basis;
@property(nonatomic, assign) double **basisFirstDerivative;
@property(nonatomic, assign) double ***basisSecondDerivative;
@property(nonatomic, assign) double **elementMetric;
@property(nonatomic, assign) double **covariantMetrixTensor;
@property(nonatomic, assign) double **ltoGMap;
@property(nonatomic, assign) double **dx;

-(BOOL)allocation:(FEMMesh *)mesh;
-(void)deallocation:(FEMMesh *)mesh;

-(BOOL)setBasisForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree;
-(BOOL)setBasisFirstDerivativeForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree;
-(BOOL)setBasisSecondDerivativeForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree;
-(BOOL)setMetricDeterminantForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w;
-(BOOL)setMetricForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w;

-(void)setCovariantMetrixTensorForElement:(Element_t *)element nDOFs:(int)nDOFs nodes:(Nodes_t*)nodes mesh:(FEMMesh *)mesh dLBasisdx:(double **)dLBasisdx;
-(BOOL)setLtoGMapForElement:(Element_t*)element nodes:(Nodes_t*)nodes mesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w;
-(void)setDxForElement:(Element_t *)element nDOFs:(int)nDOFs nodes:(Nodes_t*)nodes dLBasisdx:(double **)dLBasisdx;

-(double)detJForElement:(Element_t*)element nodes:(Nodes_t*)nodes mesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w;
-(void)invertMatrix3x3:(double **)GI detJ:(double)detG;
-(void)globalSecondDerivativesForElement:(Element_t*)element nodes:(Nodes_t*)nodes mesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w nodalValues:(double*)f dLBasisdx:(double**)dLBasisdx values:(double **)values;



@end
