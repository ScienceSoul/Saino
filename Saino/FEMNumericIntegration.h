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

@interface FEMNumericIntegration : FEMMesh {

//@private
    GaussIntegrationPoints IntegCompound;
    double *basis, **dBasisdx, ***ddBasisddx;
    double detJ, **elementMetric, **covariantMetrixTensor, **LtoGMap;
    double **dx; // Partial derivatives of global coordinates with respect to local coordinates
    
}

-(BOOL)setBasis: (Element_t*)element: (Nodes_t*)nodes: (double)u: (double)v: (double)w: (BOOL)Bubbles: (int)el;
-(BOOL)setdBasisdx: (Element_t*)element: (Nodes_t*)nodes: (double)u: (double)v: (double)w: (BOOL)Bubbles: (int)el;
-(double)dBasisdx: (int)i: (int)j;
-(BOOL)setddBasisddx: (Element_t*)element: (Nodes_t*)nodes: (double)u: (double)v: (double)w: (BOOL)Bubbles: (int)el;
-(void)setdx: (int)NDOFs: (Nodes_t*)nodes: (double **)dLBasisdx: (int)el;
-(double)dx: (int)i: (int)j;
-(void)setCovariantMetrixTensor: (int)NDOFs: (Nodes_t*)nodes: (double **)dLBasisdx: (int)el;
-(double)covariantMetrixTensor: (int)i: (int)j;
-(BOOL)setdetJ: (Element_t*)element: (Nodes_t*)nodes: (double)u: (double)v: (double)w: (int)el;
-(double)detJ;
-(BOOL)setElementMetric: (Element_t*)element: (Nodes_t*)nodes: (double)u: (double)v: (double)w: (int)el;
-(double)elementMetric: (int)i: (int)j;
-(double)detG: (Element_t*)element: (Nodes_t*)nodes: (double)u: (double)v: (double)w: (int)el;
-(void)InvertMatrix3x3: (double **)GI: (double)detG;
-(BOOL)setLtoGMap: (Element_t*)element: (Nodes_t*)nodes: (double)u: (double)v: (double)w: (int)el;
-(double)LtoGMap: (int)i: (int)j;
-(void)GlobalSecondDerivatives: (Element_t*)element: (Nodes_t*)nodes: (double)u: (double)v: (double)w: (double*)f: (double**)dLBasisdx: (double **)values: (int)el;

@end
