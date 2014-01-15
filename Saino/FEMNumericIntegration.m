//
//  FEMSolver.m
//  Saino
//
//  Created by Seddik hakime on 17/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMNumericIntegration.h"
#import "NodalBasisFunctions.h"
#include "Utils.h"


@implementation FEMNumericIntegration

@synthesize metricDeterminant = _metricDeterminant;
@synthesize basis = _basis;
@synthesize basisFirstDerivative = _basisFirstDerivative;
@synthesize basisSecondDerivative = _basisSecondDerivative;

- (id)init
{
    self = [super init];
    if (self) {
        _metricDeterminant = 0.0;
        _basis = NULL;
        _basisFirstDerivative = NULL;
        _basisSecondDerivative = NULL;
    }
    
    return self;
}

-(BOOL)allocation:(FEMMesh *)mesh {
    
    _basis = doublevec(0, mesh.maxElementNodes-1);
    if (_basis == NULL) return NO;
    memset( _basis, 0.0, mesh.maxElementNodes*sizeof(double) );
    
    _basisFirstDerivative = doublematrix(0, mesh.maxElementNodes-1, 0, 2);
    if (_basisFirstDerivative == NULL) return NO;
    memset( *_basisFirstDerivative, 0.0, (mesh.maxElementNodes*3)*sizeof(double) );
    
    _basisSecondDerivative = d3tensor(0, mesh.maxElementNodes-1, 0, 2, 0, 2);
    if (_basisSecondDerivative == NULL) return NO;
    memset(**_basisSecondDerivative, 0.0, (mesh.maxElementNodes*3*3)*sizeof(double) );
    
    return YES;
}

-(void)deallocation:(FEMMesh *)mesh {
    
    free_dvector(_basis, 0, mesh.maxElementNodes-1);
    free_dmatrix(_basisFirstDerivative, 0, mesh.maxElementNodes-1, 0, 2);
    free_d3tensor(_basisSecondDerivative, 0, mesh.maxElementNodes-1, 0, 2, 0, 2);
}


#pragma mark Setters

-(BOOL)setBasisForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree {
/***************************************************************************************************************
    Basis function values at (u,v,w)
    
    Element_t *element ->  Element structure
    Nodes_t *nodes     ->  Element nodal coordinates
    FEMMesh *mesh      ->  Finite element mesh
    double u           ->  1st Local coordinate at which to calculate the basis function
    double v           ->  2nd local coordinate
    double w           ->  3rd local coordinate
    BOOL bubbles       ->  Are the bubbles to be avaluated
    int *degree        ->  Degree of each basis function in Basis(:) vector.
                           !! May be used with P element basis functions
 
    Retun NO if element is degenerate
 
***************************************************************************************************************/
    
    int n, dim, cdim;
    
    n = element->Type.NumberOfNodes;
    dim = element->Type.dimension;
    cdim = mesh.dimension;
    
    if (element->Type.ElementCode == 101) {
        _basis[0] = 1.0;
        return YES;
    }
    
    memset( _basis, 0.0, n*sizeof(double) );
    NodalBasisFunctions(n, _basis, element, u, v, w);
    
    return YES;
}

-(BOOL)setBasisFirstDerivativeForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree {
/********************************************************************************************
 
    Global first derivatives of basis functions at (u, v, w)
 
    Element_t *element ->  Element structure
    Nodes_t *nodes     ->  Element nodal coordinates
    FEMMesh *mesh      ->  Finite element mesh
    double u           ->  1st Local coordinate at which to calculate the basis function
    double v           ->  2nd local coordinate
    double w           ->  3rd local coordinate
    BOOL bubbles       ->  Are the bubbles to be avaluated
    int *degree        ->  Degree of each basis function in Basis(:) vector.
                           !! May be used with P element basis functions
 
    Retun NO if element is degenerate
 
********************************************************************************************/
    
    int i, j, k, n, q, dim, cdim;
    double ltoGMap[9];
    bool success;
    
    n = element->Type.NumberOfNodes;
    dim = element->Type.dimension;
    cdim = mesh.dimension;
    
    if (element->Type.ElementCode == 101) {
        for (i=0; i<2; i++) {
            self.basisFirstDerivative[0][i] = 0.0;
        }
        return YES;
    }
    
    double dLBasisdx[n*3];
    memset( dLBasisdx, 0.0, sizeof(dLBasisdx) );
    NodalFirstDerivatives(n, dLBasisdx, element, u, v, w);
    
    q = n;
    memset( *_basisFirstDerivative, 0.0, (n*3)*sizeof(double) );
    memset(ltoGMap, 0.0, sizeof(ltoGMap));
    success = localtoGlobalMap(ltoGMap, element, nodes, cdim, u, v, w, dLBasisdx);
    
    for (i=0; i<q; i++) {
        for (j=0; j<cdim; j++) {
            for (k=0; k<dim; k++) {
                _basisFirstDerivative[i][j] = _basisFirstDerivative[i][j] + dLBasisdx[3*i+k]*ltoGMap[3*j+k];
            }
        }
    }
    
    return YES;
}


/********************************************************************************************
 
    Global second derivatives of basis functions at (u, v, w)
 
    Element_t *element ->  Element structure
    Nodes_t *nodes     ->  Element nodal coordinates
    FEMMesh *mesh      ->  Finite element mesh
    double u           ->  1st Local coordinate at which to calculate the basis function
    double v           ->  2nd local coordinate
    double w           ->  3rd local coordinate
    BOOL bubbles       ->  Are the bubbles to be avaluated
    int *degree        ->  Degree of each basis function in Basis(:) vector.
    !! May be used with P element basis functions
 
    Retun NO if element is degenerate
 
********************************************************************************************/
-(BOOL)setBasisSecondDerivativeForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree {
    
    int i, j, n, q, dim, cdim;
    double dx[9], elementMetric[9], values[9];
    bool success;
    
    n = element->Type.NumberOfNodes;
    dim = element->Type.dimension;
    cdim = mesh.dimension;
    
    double NodalBasis[n];
    double dLBasisdx[n*3];
    memset( dLBasisdx, 0.0, sizeof(dLBasisdx) );
    memset(**_basisSecondDerivative, 0.0, (n*3*3)*sizeof(double) );
    NodalFirstDerivatives(n, dLBasisdx, element, u, v, w);
    
    memset(dx, 0.0, sizeof(dx));
    memset(elementMetric, 0.0, sizeof(elementMetric));
    success = contravariantMetric(elementMetric, dx, element, nodes, cdim, u, v, w, dLBasisdx);
    
    memset( NodalBasis, 0.0, sizeof(NodalBasis) );
    memset(values, 0.0, sizeof(values));
    for (q=0; q<n; q++) {
        NodalBasis[q] = 1.0;
        globalSecondDerivatives(elementMetric, element, nodes, cdim, u, v, w, NodalBasis, dLBasisdx, values);
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) {
                _basisSecondDerivative[q][i][j] = values[3*i+j];
            }
        }
        NodalBasis[q] = 0.0;
    }
    
    return YES;
}

/***************************************************************************************************************************
 
    Square root of determinant of covariant metric tensor (=sqrt(det(J^TJ))).
 
***************************************************************************************************************************/
-(BOOL)setMetricDeterminantForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w {

    int n = element->Type.NumberOfNodes;
    double covariantMetricTensor[9], dx[9];
    
    double dLBasisdx[n*3];
    memset( dLBasisdx, 0.0, sizeof(dLBasisdx) );
    NodalFirstDerivatives(n, dLBasisdx, element, u, v, w);
    
    memset(dx, 0.0, sizeof(dx));
    memset(covariantMetricTensor, 0.0, sizeof(covariantMetricTensor) );
    _metricDeterminant = detJ(covariantMetricTensor, dx, element, nodes, mesh.dimension, u, v, w, dLBasisdx);
    _metricDeterminant = sqrt(_metricDeterminant);
    
    return YES;
}

@end