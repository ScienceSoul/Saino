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

//@private
    GaussIntegrationPoints IntegCompound;
    double *basis, **basisFirstDerivative, ***basisSecondDerivative;
    double metricDeterminant, **elementMetric, **covariantMetrixTensor, **ltoGMap;
    double **dx;         // Partial derivatives of global coordinates with respect to local coordinates
    
}

-(BOOL)allocation:(FEMMesh *)mesh;
-(void)deallocation:(FEMMesh *)mesh;

-(double)basis:(int)i;
-(double)basisFirstDerivative:(int)i: (int)j;
-(double)metricDeterminant;
-(double)elementMetric:(int)i: (int)j;
-(double)covariantMetrixTensor:(int)i: (int)j;
-(double)ltoGMap:(int)i: (int)j;
-(double)dx:(int)i: (int)j;

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

// Return pointers to data structures
-(double *)returnPointerToBasis;

@end
