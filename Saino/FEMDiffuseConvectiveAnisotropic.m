//
//  FEMDiffuseConvectiveAnisotropic.m
//  Saino
//
//  Created by Seddik hakime on 04/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMDiffuseConvectiveAnisotropic.h"

#import "FEMNumericIntegration.h"
#import "GaussIntegration.h"

@implementation FEMDiffuseConvectiveAnisotropic

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
    }
    
    return self;
}

/*********************************************************************************************************************************
    Return element local matrices and RHS vector for boundary conditions of diffusion convection equation.
    
    Arguments:
        - double **boundaryMatrix   ->  coefficient matrix if equations
        - double *boundaryVector    ->  RHS vector
        - double *loadVector        ->  coefficient of the force term
        - double *nodalAlpha        ->  coefficient for temperature dependent term
        - Element_t *element        ->  Structure describing the element (dimension, nb of nodes, interpolation degree, etc...)
        - int n                     ->  number of element nodes
        - Nodes_t *nodes            ->  Element node coordinates
*********************************************************************************************************************************/
-(void)diffuseConvectiveBoundaryMatrix:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector dimensions:(Dimensions_t)dimensions loadVector:(double *)loadVector nodalAlpha:(double *)nodalAlpha element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh; {
    
    int i, p, q, t;
    double alpha, detJ, force, s;
    BOOL stat;
    FEMNumericIntegration *integration;
    GaussIntegrationPoints *IP;
    
    memset( *boundaryMatrix, 0.0, (dimensions.mat1*dimensions.mat2)*sizeof(double) );
    memset( boundaryVector, 0.0, dimensions.vec*sizeof(double) );
    
    // Integration stuff
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) errorfunct("diffuseConvectiveBoundaryMatrix", "Allocation error in FEMNumericIntegration!");
    IP = GaussQuadrature(element, NULL, NULL);

    for (t=0; t<IP->n; t++) {
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t] withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t]];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
        force = 0.0;
        alpha = 0.0;
        for (i=0; i<n; i++) {
            force = force + (loadVector[i] * integration.basis[i]);
            alpha = alpha + (nodalAlpha[i] * integration.basis[i]);
        }
        
        for (p=0; p<n; p++) {
            for (q=0; q<n; q++) {
                boundaryMatrix[p][q] = boundaryMatrix[p][q] + s * alpha * integration.basis[q] * integration.basis[p];
            }
        }
        for (q=0; q<n; q++) {
            boundaryVector[q] = boundaryVector[q] + s * integration.basis[q] * force;
        }
    }
    GaussQuadratureDeallocation(IP);
    [integration deallocation:mesh];
}

@end
