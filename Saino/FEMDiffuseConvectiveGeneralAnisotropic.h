//
//  FEMDiffuseConvectiveGeneralAnisotropic.h
//  Saino
//
//  Created by Seddik hakime on 04/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMKernel.h"

@interface FEMDiffuseConvectiveGeneralAnisotropic : NSObject

-(void)diffuseConvectiveGeneralComposeMassMatrix:(double **)massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double *)loadVector timeDerivativeTerm:(double *)nodalCT zeroDegreeTerm:(double *)nodalC0 convectionTerm:(double *)nodalC1 diffusionTerm:(double ***)nodalC2 phaseChange:(BOOL)phaseChange nodalTemperature:(double *)nodalTemperature enthalpy:(double *)enthalpy velocityX:(double *)ux velocitY:(double *)uy velocityZ:(double *)uz meshVeloX:(double *)mux meshVeloY:(double *)muy meshVeloZ:(double *)muz nodalViscosity:(double *)nodalviscosity nodaldensity:(double *)nodalDensity nodalPressure:(double *)nodalPressure nodalPressureDt:(double *)nodalPressureDt nodalPressureCoeff:(double *)nodalPressureCoeff compressible:(BOOL)compressible stabilize:(BOOL)stabilize element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution kernel:(FEMKernel *)kernel mesh:(FEMMesh *)mesh model:(FEMModel *)model listUtilities:(FEMListUtilities *)listUtilities;
-(void)diffuseConvectiveGeneralBoundaryMatrix:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector dimensions:(Dimensions_t)dimensions loadVector:(double *)loadVector nodalAlpha:(double *)nodalAlpha element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes model:(FEMModel *)model mesh:(FEMMesh *)mesh;

@end
