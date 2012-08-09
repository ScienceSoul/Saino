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

//@class FEMMesh;

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

-(void)setCovariantMetrixTensor:(Element_t *)element: (int)nDOFs: (Nodes_t*)nodes: (FEMMesh *)mesh: (double **)dLBasisdx;
-(BOOL)setLtoGMap:(Element_t*)element: (Nodes_t*)nodes: (FEMMesh *)mesh: (double)u: (double)v: (double)w;
-(void)setdx:(Element_t *)element: (int)nDOFs: (Nodes_t*)nodes: (double **)dLBasisdx;

-(double)detG:(Element_t*)element: (Nodes_t*)nodes: (FEMMesh *)mesh: (double)u: (double)v: (double)w;
-(void)invertMatrix3x3:(double **)GI: (double)detG;
-(void)globalSecondDerivatives:(Element_t*)element: (Nodes_t*)nodes: (FEMMesh *)mesh: (double)u: (double)v: (double)w: (double*)f: (double**)dLBasisdx: (double **)values;



@end
