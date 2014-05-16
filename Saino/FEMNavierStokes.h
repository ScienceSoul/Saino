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
#import "FEMMaterial.h"
#import "FEMElementDescription.h"
#import "FEMMaterialModels.h"
#import "FEMDifferentials.h"

@interface FEMNavierStokes : NSObject
-(void)navierStokesComposeMassMatrix:(double **) massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double **)loadVector nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz meshVelocityX:(double *)mux meshVelocityY:(double *)muy meshVelocityZ:(double *)muz nodalPressure:(double *)nodalPressure nodalTemperature:(double *)nodalTemperature convect:(BOOL)convect stabilizeFlag:(NSString *)stabilizeFlag compressibilityModel:(int)compressibilityModel pseudoCompressible:(BOOL)pseudoCompressible nodalCompressibility:(double *)nodalCompressibility nodalGasConstant:(double *)nodalGasConstant porous:(BOOL)porous nodalDrag:(double **)nodalDrag potentialForce:(BOOL)potentialForce potentialField:(double *)potentialField potentialCoefficient:(double *)potentialCoefficient magneticForce:(BOOL)magneticForce rotating:(BOOL)rotating omega:(double *)omega divDiscretization:(BOOL)divDiscretization gradDriscretization:(BOOL)gradDriscretization newtonLinearization:(BOOL)newtonLinearization transient:(BOOL)transient element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution core:(FEMCore *)core mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration material:(FEMMaterial *)material elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems materialModels:(FEMMaterialModels*)materialModels differentials:(FEMDifferentials *)differentials listUtilities:(FEMListUtilities *)listUtilities utilities:(FEMUtilities *)utilities;
@end
