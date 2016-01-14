//
//  FEMNavierStokes.h
//  Saino
//
//  Created by Seddik hakime on 01/05/2014.
//  Copyright (c) 2014 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMCore.h"
#import "FEMNumericIntegration.h"
#import "FEMElementDescription.h"
#import "FEMMaterialModels.h"
#import "FEMDifferentials.h"
#import "FEMElementUtils.h"

//  Class computing Navier-Stokes local matrices (cartesian coordinates)

@interface FEMNavierStokes : NSObject
-(void)navierStokesComposeMassMatrix:(double **) massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double **)loadVector nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz meshVelocityX:(double *)mux meshVelocityY:(double *)muy meshVelocityZ:(double *)muz nodalPressure:(double *)nodalPressure nodalTemperature:(double *)nodalTemperature isConvect:(BOOL)convect stabilizeFlag:(NSString *)stabilizeFlag compressibilityModel:(int)compressibilityModel isPseudoCompressible:(BOOL)pseudoCompressible nodalCompressibility:(double *)nodalCompressibility nodalGasConstant:(double *)nodalGasConstant isPorous:(BOOL)porous nodalDrag:(double **)nodalDrag isPotentialForce:(BOOL)potentialForce potentialField:(double *)potentialField potentialCoefficient:(double *)potentialCoefficient isMagneticForce:(BOOL)magneticForce isRotating:(BOOL)rotating omega:(double *)omega isDivDiscretization:(BOOL)divDiscretization isGradPDriscretization:(BOOL)gradPDriscretization isNewtonLinearization:(BOOL)newtonLinearization isTransient:(BOOL)transient element:(Element_t *)element numberOfNodes:(int)n rows:(int)rows cols:(int)cols nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution core:(FEMCore *)core mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration material:(FEMMaterial *)material elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems materialModels:(FEMMaterialModels *)materialModels differentials:(FEMDifferentials *)differentials listUtilities:(FEMListUtilities *)listUtilities utilities:(FEMUtilities *)utilities;
-(void)navierStokesBoundary:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector loadVector:(double **)loadVector nodalAlpha:(double *)nodalAlpha nodalBeta:(double *)nodalBeta nodalExtPressure:(double *)nodalExtPressure nodalSlipCoefficient:(double **)nodalSlipCoefficient isNormalTangential:(BOOL)normalTangential element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh  model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription elementUtils:(FEMElementUtils *)elementUtils;
-(void)vmsWallsBoundary:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector layerThickness:(double *)layerThickness surfaceRoughness:(double *)surfaceRoughness nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription elementUtils:(FEMElementUtils *)elementUtils;
-(void)navierStokesWallLawBoundary:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector layerThickness:(double *)layerThickness surfaceRoughness:(double *)surfaceRoughness nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription elementUtils:(FEMElementUtils *)elementUtils;
@end
