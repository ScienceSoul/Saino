//
//  FEMNavierStokesGeneral.h
//  Saino
//
//  Created by Seddik hakime on 28/05/2014.
//  Copyright (c) 2014 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMCore.h"
#import "FEMNumericIntegration.h"
#import "FEMElementDescription.h"
#import "FEMMaterialModels.h"
#import "FEMDifferentials.h"
#import "FEMElementUtils.h"

// Class computing Navier-Stokes local matrices in general coordinate system
// (i.e. not cartesian, axisymmetric or cylindrically symmetric).

@interface FEMNavierStokesGeneral : NSObject
-(void)navierStokesGeneralComposeMassMatrix:(double **) massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double **)loadVector nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz meshVelocityX:(double *)mux meshVelocityY:(double *)muy meshVelocityZ:(double *)muz isStabilize:(BOOL)stabilize isNewtonLinearization:(BOOL)newtonLinearization element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes core:(FEMCore *)core mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems materialModels:(FEMMaterialModels *)materialModels differentials:(FEMDifferentials *)differentials listUtilities:(FEMListUtilities *)listUtilities utilities:(FEMUtilities *)utilities;
-(void)navierStokesGeneralBoundary:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector loadVector:(double **)loadVector nodalAlpha:(double *)nodalAlpha nodalBeta:(double *)nodalBeta nodalExtPressure:(double *)nodalExtPressure nodalSlipCoefficient:(double **)nodalSlipCoefficient element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh  model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems;

@end
