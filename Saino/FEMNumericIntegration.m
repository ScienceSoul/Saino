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

@synthesize IP = _IP;
@synthesize metricDeterminant = _metricDeterminant;
@synthesize basis = _basis;
@synthesize basisFirstDerivative = _basisFirstDerivative;
@synthesize basisSecondDerivative = _basisSecondDerivative;
@synthesize elementMetric = _elementMetric;
@synthesize covariantMetrixTensor = _covariantMetrixTensor;
@synthesize ltoGMap = _ltoGMap;
@synthesize dx = _dx;

- (id)init
{
    self = [super init];
    if (self) {
        // TODO: Initialization code here.
        _metricDeterminant = 0.0;
        _basis = NULL;
        _basisFirstDerivative = NULL;
        _basisSecondDerivative = NULL;
        _elementMetric = NULL;
        _covariantMetrixTensor = NULL;
        _ltoGMap = NULL;
        _dx = NULL;
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
    
    _elementMetric = doublematrix(0, 2, 0, 2);
    if (_elementMetric == NULL) return NO;
    memset( *_elementMetric, 0.0, (3*3)*sizeof(double) );
    
    _covariantMetrixTensor = doublematrix(0, mesh.dimension-1, 0, mesh.dimension-1);
    if (_covariantMetrixTensor == NULL) return NO;
    memset( *_covariantMetrixTensor, 0.0, (mesh.dimension*mesh.dimension)*sizeof(double) );
    
    _ltoGMap = doublematrix(0, mesh.dimension-1, 0, mesh.dimension-1);
    if (_ltoGMap == NULL) return NO;
    memset( *_ltoGMap, 0.0, (mesh.dimension*mesh.dimension)*sizeof(double) );

    _dx = doublematrix(0, 2, 0, mesh.dimension-1);
    if (_dx == NULL) return NO;
    memset( *_dx, 0.0, (3*mesh.dimension)*sizeof(double) );
    
    return YES;
}

-(void)deallocation:(FEMMesh *)mesh {
    
    free_dvector(_basis, 0, mesh.maxElementNodes-1);
    free_dmatrix(_basisFirstDerivative, 0, mesh.maxElementNodes-1, 0, 2);
    free_d3tensor(_basisSecondDerivative, 0, mesh.maxElementNodes-1, 0, 2, 0, 2);
    free_dmatrix(_elementMetric, 0, 2, 0, 2);
    free_dmatrix(_covariantMetrixTensor, 0, mesh.dimension-1, 0, mesh.dimension-1);
    free_dmatrix(_ltoGMap, 0, mesh.dimension-1, 0, mesh.dimension-1);
    free_dmatrix(_dx, 0, 2, 0, mesh.dimension-1);
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
    cdim = [mesh dimension];
    
    if (element->Type.ElementCode == 101) {
        self.basis[0] = 1.0;
        return YES;
    }
    
    memset( self.basis, 0.0, n*sizeof(double) );
    NodalBasisFunctions(n, self.basis, element, u, v, w);
    
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
    double **dLBasisdx;
    
    n = element->Type.NumberOfNodes;
    dim = element->Type.dimension;
    cdim = [mesh dimension];
    
    if (element->Type.ElementCode == 101) {
        for (i=0; i<2; i++) {
            self.basisFirstDerivative[0][i] = 0.0;
        }
    }
    
    dLBasisdx = doublematrix(0, n-1, 0, 2);
    memset( *dLBasisdx, 0.0, (n*dim)*sizeof(double) );    
    NodalFirstDerivatives(n, dLBasisdx, element, u, v, w);
    
    q = n;
    memset( *_basisFirstDerivative, 0.0, (n*2)*sizeof(double) );
    [self setLtoGMapForElement:element nodes:nodes mesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
    
    for (i=0; i<q; i++) {
        for (j=0; j<cdim; j++) {
            for (k=0; k<dim; k++) {
                self.basisFirstDerivative[i][j] = self.basisFirstDerivative[i][j] * dLBasisdx[i][k]*self.ltoGMap[j][k];
            }
        }
    }
    
    free_dmatrix(dLBasisdx, 0, n-1, 0, 2);
    
    return YES;
    
}


-(BOOL)setBasisSecondDerivativeForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w withBubbles:(BOOL)bubbles basisDegree:(int *)degree {
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

    
    int i, j, n, q, dim, cdim;
    double *NodalBasis, **dLBasisdx;
    double **Values;
    
    n = element->Type.NumberOfNodes;
    dim = element->Type.dimension;
    cdim = [mesh dimension];
    
    NodalBasis = doublevec(0, n-1);
    dLBasisdx = doublematrix(0, n-1, 0, 2);
    Values = doublematrix(0, 2, 0, 2);
    
    memset( NodalBasis, 0.0, n*sizeof(double) );
    memset( *dLBasisdx, 0.0, (n*dim)*sizeof(double) );
    memset(**_basisSecondDerivative, 0.0, (n*2*2)*sizeof(double) );
    NodalFirstDerivatives(n, dLBasisdx, element, u, v, w);
    [self setMetricForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
    
    for (q=0; q<n; q++) {
        NodalBasis[q] = 1.0;
        [self globalSecondDerivativesForElement:element nodes:nodes mesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w nodalValues:NodalBasis dLBasisdx:dLBasisdx values:Values];
        for (i=0; i<2; i++) {
            for (j=0; j<2; j++) {
                self.basisSecondDerivative[q][i][j] = Values[i][j];
            }
        }
        NodalBasis[q] = 0.0;
    }
    
    free_dvector(NodalBasis, 0, n-1);
    free_dmatrix(dLBasisdx, 0, n-1, 0, 2);
    free_dmatrix(Values, 0, 2, 0, 2);
    
    return YES;
    
}

-(BOOL)setMetricDeterminantForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w {
/***************************************************************************************************************************
 
    Square root of determinant of covariant metric tensor (=sqrt(det(J^TJ))).

***************************************************************************************************************************/
    self.metricDeterminant = [self detJForElement:element nodes:nodes mesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
    
    self.metricDeterminant = sqrt(self.metricDeterminant);
    
    return YES;
}

-(BOOL)setMetricForElement:(Element_t *)element elementNodes:(Nodes_t *)nodes inMesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w {
/********************************************************************************************************
 
    Compute contravariant metric tensor (=J^TJ)^-1 of element coordinate system
 
*********************************************************************************************************/
    
    int i, j, dim, cdim;
    double detG, **GI;
    
    dim = element->Type.dimension;
    cdim = [mesh dimension];
    
    detG = [self detJForElement:element nodes:nodes mesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
    
    // Convert the metric to contravariant base
    switch (dim) {
        case 1:
            self.elementMetric[0][0] = 1.0 / detG;
            break;
        case 2:
            self.elementMetric[0][0] = self.covariantMetrixTensor[1][1] / detG;
            self.elementMetric[0][1] = -self.covariantMetrixTensor[0][1] / detG;
            self.elementMetric[1][0] = -self.covariantMetrixTensor[1][0] / detG;
            self.elementMetric[1][1] = self.covariantMetrixTensor[0][0] / detG;
            break;
        case 3:
            GI = doublematrix(0, 2, 0, 2);
            [self invertMatrix3x3:GI detJ:detG];
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    self.elementMetric[i][j] = GI[i][j];
                }
            }
            free_dmatrix(GI, 0, 2, 0, 2);
            break;
            
        default:
            errorfunct("FEMNumericIntegration:setMetricForElement", "Dimension not supported!!");
            break;
    }
    
    return YES;
    
}

-(void)setCovariantMetrixTensorForElement:(Element_t *)element nDOFs:(int)nDOFs nodes:(Nodes_t*)nodes mesh:(FEMMesh *)mesh dLBasisdx:(double **)dLBasisdx {
/******************************************************************************************************
    Compute the covariant metric tensor of the element coordinate system
 
    Arguments:
        NDOFs              -> Number of DOFs
 
        Nodes_t *nodes     -> element nodal cooridnates 
 
        double **dLBasisdx -> Derivatives of element basis function with respect to 
        local coordinates
 
******************************************************************************************************/
    
    int i, j, k, dim, cdim;
    double s;
    
    dim = element->Type.dimension;
    cdim = [mesh dimension];
    
    [self setDxForElement:element nDOFs:nDOFs nodes:nodes dLBasisdx:dLBasisdx];
    
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<cdim; k++) {
                s = s + ( self.dx[k][i] * self.dx[k][j]);
            }
            self.covariantMetrixTensor[i][j] = s;
        }
    }
}

-(BOOL)setLtoGMapForElement:(Element_t*)element nodes:(Nodes_t*)nodes mesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w {
    
    int i, j, k, dim, cdim;
    double s;
    
    dim = element->Type.dimension;
    cdim = [mesh dimension];
    
    [self setMetricForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
    
    for (i=0; i<cdim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                s = s + self.dx[i][k] * self.elementMetric[k][j];
            }
            self.ltoGMap[i][j] = s;
        }
    }
    
    return YES;
    
}

-(void)setDxForElement:(Element_t *)element nDOFs:(int)nDOFs nodes:(Nodes_t*)nodes dLBasisdx:(double **)dLBasisdx {
/***********************************************************************************
    Partial derivatives of global coordinates with respect to local coordinates
    
    Arguments:
        NDOFs              -> Number of DOFs
     
        Nodes_t *nodes     -> element nodal cooridnates 
     
        double **dLBasisdx -> Derivatives of element basis function with respect 
                              to local coordinates
     
        int el             -> Element number                                       

***********************************************************************************/
    int i, j, n, dim;
    double accum1, accum2, accum3;
    
    dim = element->Type.dimension;
    n = min(element->Type.NumberOfNodes, nDOFs);
    
    for (i=0; i<dim; i++) {
        accum1 = 0.0;
        accum2 = 0.0;
        accum3 = 0.0;
        for (j=0; j<n; j++) {
            accum1 = accum1 + (nodes->x[j] * dLBasisdx[j][i]);
            accum2 = accum2 + (nodes->y[j] * dLBasisdx[j][i]);
            accum3 = accum3 + (nodes->z[j] * dLBasisdx[j][i]);
        }
        self.dx[0][i] = accum1;
        self.dx[1][i] = accum2;
        self.dx[2][i] = accum3;
    }
}

-(double)detJForElement:(Element_t*)element nodes:(Nodes_t*)nodes mesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w {
/**********************************************************************************************
    Compute determinant of covariant metric tensor det(J^TJ) and determinant 
    of covariant metric tensor det(J^TJ).
 
    Arguments:
        Element_t *element  -> element structure
 
        Nodes_t *nodes      -> element nodal cooridnates 
 
        double u, v, w      -> Points at which evaluate the value
 
        int el              -> Element number
 
        Function return value:
        If function failure, element is degenerate
 
**********************************************************************************************/
    
    int n, dim, cdim;
    double detG, **dLBasisdx;
    
    n = element->Type.NumberOfNodes;
    dim = element->Type.dimension;
    cdim = [mesh dimension];
    
    dLBasisdx = doublematrix(0, n-1, 0, 2);
    memset( *dLBasisdx, 0.0, (n*dim)*sizeof(double) );    
    NodalFirstDerivatives(n, dLBasisdx, element, u, v, w);
    [self setCovariantMetrixTensorForElement:element nDOFs:n nodes:nodes mesh:mesh dLBasisdx:dLBasisdx];
    
    detG = 0.0;
    switch (dim) {
           
        case 1:  // Line elements
            detG = self.covariantMetrixTensor[1][1];
            if (detG <= 0.0) {
                errorfunct("FEMNumericIntegration:detJForElement", "Degenerate 1D element");
            }
            break;
            
        case 2: // Surface elements
            detG = ( self.covariantMetrixTensor[0][0]*self.covariantMetrixTensor[1][1] - self.covariantMetrixTensor[0][1]*self.covariantMetrixTensor[1][0] );
            if (detG <= 0.0) {
                if (cdim < dim) {
                    warnfunct("FEMNumericIntegration:detJForElement", "2D element in 1D coordinate systrem?");
                }
                errorfunct("FEMNumericIntegration:detJForElement", "Degenerate 2D element");
            }
            break;
            
        case 3: // Volume elements 
            detG = self.covariantMetrixTensor[0][0] * ( self.covariantMetrixTensor[1][1]*self.covariantMetrixTensor[2][2] - self.covariantMetrixTensor[1][2]*self.covariantMetrixTensor[2][1] ) + self.covariantMetrixTensor[0][1] * ( self.covariantMetrixTensor[1][2]*self.covariantMetrixTensor[2][0] - self.covariantMetrixTensor[1][0]*self.covariantMetrixTensor[2][2] ) + self.covariantMetrixTensor[0][2] * ( self.covariantMetrixTensor[1][0]*self.covariantMetrixTensor[2][1] - self.covariantMetrixTensor[1][1]*self.covariantMetrixTensor[2][0] );
            if (detG <= 0.0) {
                if (cdim < dim) {
                    warnfunct("FEMNumericIntegration:detJForElement", "2D/3D element in 1D/2D coordinate systrem?");
                }
                errorfunct("FEMNumericIntegration:detJForElement", "Degenerate 3D element");
            }
            break;
            
        default:
            errorfunct("FEMNumericIntegration:detJForElement", "Dimension not supported!!");
            break;
    }
    free_dmatrix(dLBasisdx, 0, n-1, 0, 2);
    
    return detG;
}

-(void)invertMatrix3x3:(double **)GI detJ:(double)detG {
    
    double s;
    
    s = 1.0 / detG;
    
    GI[0][0] = s * ( self.covariantMetrixTensor[1][1]*self.covariantMetrixTensor[2][2] - self.covariantMetrixTensor[2][1]*self.covariantMetrixTensor[1][2] );
    
    GI[1][0] = -s * ( self.covariantMetrixTensor[1][0]*self.covariantMetrixTensor[2][2] - self.covariantMetrixTensor[2][0]*self.covariantMetrixTensor[1][2] );
    
    GI[2][0] = s * ( self.covariantMetrixTensor[1][0]*self.covariantMetrixTensor[2][1] - self.covariantMetrixTensor[2][0]*self.covariantMetrixTensor[1][1] );
    
    GI[0][1] = -s * ( self.covariantMetrixTensor[0][1]*self.covariantMetrixTensor[2][2] - self.covariantMetrixTensor[2][1]*self.covariantMetrixTensor[0][2] );
    
    GI[1][1] = s * ( self.covariantMetrixTensor[0][0]*self.covariantMetrixTensor[2][2] - self.covariantMetrixTensor[2][0]*self.covariantMetrixTensor[0][2] );
    
    GI[2][1] = -s * ( self.covariantMetrixTensor[0][0]*self.covariantMetrixTensor[2][1] - self.covariantMetrixTensor[2][0]*self.covariantMetrixTensor[0][1] );
    
    GI[0][2] = s * ( self.covariantMetrixTensor[0][1]*self.covariantMetrixTensor[1][2] - self.covariantMetrixTensor[1][1]*self.covariantMetrixTensor[0][2] );
    
    GI[1][2] = -s * ( self.covariantMetrixTensor[0][0]*self.covariantMetrixTensor[1][2] - self.covariantMetrixTensor[1][0]*self.covariantMetrixTensor[0][2] );
    
    GI[2][2] = s * ( self.covariantMetrixTensor[0][0]*self.covariantMetrixTensor[1][1] - self.covariantMetrixTensor[1][0]*self.covariantMetrixTensor[0][1] );
}

-(void)globalSecondDerivativesForElement:(Element_t*)element nodes:(Nodes_t*)nodes mesh:(FEMMesh *)mesh firstEvaluationPoint:(double)u secondEvaluationPoint:(double)v thirdEvaluationPoint:(double)w nodalValues:(double*)f dLBasisdx:(double**)dLBasisdx values:(double **)values {
/******************************************************************************************************
    Arguments:
        Element_t* element -> structure describing the element
        Nodes_t* nodes     -> nodal coordinates
        FEMMesh *mesh      -> mesh
        double u,v,w       -> point at which to evaluate
        double *f          -> nodal values of the quantity
 
    Output: 3x3 matrix (values) of partial derivatives

******************************************************************************************************/
    
    int i, j, k, l, n, dim, cdim;
    double ***C1, ***C2, ***ddx;
    double *df;
    double **cddf, **ddf, **dxx, **bf;
    double accum1, accum2, accum3, accum4;
    double s;

    // Actually not quite correct... 
    if (element->Type.BasisFunctionDegree <= 1 ) return;
    
    n = element->Type.NumberOfNodes;
    dim = element->Type.dimension;
    cdim = [mesh dimension];
    
    C1 = d3tensor(0, 2, 0, 2, 0, 2);
    C2 = d3tensor(0, 2, 0, 2, 0, 2);
    ddx = d3tensor(0, 2, 0, 2, 0, 2);
    
    df = doublevec(0, 2);
    
    cddf = doublematrix(0, 2, 0, 2);
    ddf = doublematrix(0, 2, 0, 2);
    dxx = doublematrix(0, 2, 0, 2);

    // Partial derivatives of the basis functions are given, just
    // sum for the first partial derivatives.
    memset( df, 0.0, 3*sizeof(double) );
    memset( *dxx, 0.0, (3*3)*sizeof(double) );
    
    switch (cdim) {
        case 1:
            for (i=0; i<dim; i++) {
                accum1 = 0.0;
                accum2 = 0.0;
                for (j=0; j<n; j++) {
                    accum1 = accum1 + nodes->x[j] * dLBasisdx[j][i];
                    accum2 = accum2 + f[j] * dLBasisdx[j][i];
                }
                dxx[0][i] = accum1;
                df[i] = accum2;
            }
            break;
            
        case 2:
            for (i=0; i<dim; i++) {
                accum1 = 0.0;
                accum2 = 0.0;
                accum3 = 0.0;
                for (j=0; j<n; j++) {
                    accum1 = accum1 + nodes->x[j] * dLBasisdx[j][i];
                    accum2 = accum2 + nodes->y[j] * dLBasisdx[j][i];
                    accum3 = accum3 + f[j] * dLBasisdx[j][i];
                }
                dxx[0][i] = accum1;
                dxx[1][i] = accum2;
                df[i] = accum3;
            }
            break;
            
        case 3:
            for (i=0; i<dim; i++) {
                accum1 = 0.0;
                accum2 = 0.0;
                accum3 = 0.0;
                accum4 = 0.0;
                for (j=0; j<n; j++) {
                    accum1 = accum1 + nodes->x[j] * dLBasisdx[j][i];
                    accum2 = accum2 + nodes->y[j] * dLBasisdx[j][i];
                    accum3 = accum3 + nodes->z[j] * dLBasisdx[j][i];
                    accum4 = accum4 + f[j] * dLBasisdx[j][i];
                }
                dxx[0][i] = accum1;
                dxx[1][i] = accum2;
                dxx[2][i] = accum3;
                df[i] = accum4;
            }
            break;

        default:
            errorfunct("FEMNumericIntegration:globalSecondDerivativesForElement", "Coordinate dimension not supported!!");
            break;
    }
    
    // Get second partial derivatives with respect to local coordinates
    switch (dim) {
        case 1:
            // Line elements
            ddx[0][0][0] = SecondDerivatives1D(element, nodes->x, u);
            ddx[1][0][0] = SecondDerivatives1D(element, nodes->y, u);
            ddx[2][0][0] = SecondDerivatives1D(element, nodes->z, u);
            break;
            
        case 2:
            // Surface elements
            bf = doublematrix(0, 1, 0, 1);
            SecondDerivatives2D(bf, element, nodes->x, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddx[0][i][j] = bf[i][j];
                }
            }
            
            SecondDerivatives2D(bf, element, nodes->y, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddx[1][i][j] = bf[i][j];
                }
            }
            
            SecondDerivatives2D(bf, element, nodes->z, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddx[2][i][j] = bf[i][j];
                }
            }
            
            free_dmatrix(bf, 0, 1, 0, 1);
            break;
            
        case 3:
            // Volume elements
            bf = doublematrix(0, 2, 0, 2);
            SecondDerivatives3D(bf, element, nodes->x, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddx[0][i][j] = bf[i][j];
                }
            }
            
            SecondDerivatives3D(bf, element, nodes->y, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddx[1][i][j] = bf[i][j];
                }
            }
            
            SecondDerivatives3D(bf, element, nodes->z, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddx[2][i][j] = bf[i][j];
                }
            }

            free_dmatrix(bf, 0, 2, 0, 2);
            break;
            
        default:
            errorfunct("FEMNumericIntegration:globalSecondDerivativesForElement", "Element dimension not supported");
            break;
    }
    
    // Christoffel symbols of the second kind of the element coordinate system
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            for (k=0; k<dim; k++) {
                s = 0.0;
                for (l=0; l<cdim; l++) {
                    s = s + ddx[l][i][j]*dxx[l][k];
                }
                C2[i][j][k] = s;
            }
        }
    }
    
    // Christoffel symbols of the first kind
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            for (k=0; k<dim; k++) {
                s = 0.0;
                for (l=0; l<dim; l++) {
                    s = s + self.elementMetric[k][l]*C2[i][j][l];
                }
                C1[i][j][k] = s;   
            }
        }
    }
    
    // First add ordinary partials (change of the quantity with coordinates)...
    switch (dim) {
        case 1:
            ddf[0][0] = SecondDerivatives1D(element, f, u);
            break;
            
        case 2:
            bf = doublematrix(0, 1, 0, 1);
            SecondDerivatives2D(bf, element, f, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddf[i][j] = bf[1][2];
                }
            }
            
            free_dmatrix(bf, 0, 1, 0, 1);
            break;
            
        case 3:
            bf = doublematrix(0, 2, 0, 2);
            SecondDerivatives3D(bf, element, f, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddf[i][j] = bf[i][j];
                }
            }
            
            free_dmatrix(bf, 0, 2, 0, 2);
            break;
            
        default:
            errorfunct("FEMNumericIntegration:globalSecondDerivativesForElement", "Element dimension not supported");
            break;
    }
    
    // ... Then add change of coordinates
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                s = s + C1[i][j][k]*df[k];
            }
            ddf[i][j] = ddf[i][j] + s;
        }
    }
    
    // Convert to contravariant base
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                for (l=0; l<dim; l++) {
                    s = s + self.elementMetric[i][k]*self.elementMetric[j][l]*ddf[k][l];
                }
            }
            cddf[i][j] = s;
        }
    }
    
    // And finally transform to global coordinates
    memset( *values, 0.0, (cdim*cdim)*sizeof(double) );
    for (i=0; i<cdim; i++) {
        for (j=0; j<cdim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                for (l=0; l<dim; l++) {
                    s = s + dxx[i][k]*dxx[j][l]*cddf[k][l];
                }
            }
            values[i][j] = s;
        }
    }
    
    free_d3tensor(C2, 0, 2, 0, 2, 0, 2);
    free_d3tensor(C1, 0, 2, 0, 2, 0, 2);
    free_d3tensor(ddx, 0, 2, 0, 2, 0, 2);
    
    free_dvector(df, 0, 2);
    
    free_dmatrix(cddf, 0, 2, 0, 2);
    free_dmatrix(ddf, 0, 2, 0, 2);
    free_dmatrix(dxx, 0, 2, 0, 2);
}

@end
