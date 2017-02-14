//===----------------------------------------------------------------------===//
//  FEMDiffuseConvectiveAnisotropic.h
//  Saino
//
//  Created by Seddik hakime on 04/06/13.
//  Copyright (c) 2013 ScienceSoule. All rights reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

#import <Foundation/Foundation.h>
#import "FEMCore.h"
#import "FEMNumericIntegration.h"
#import "FEMMaterialModels.h"
#import "FEMDifferentials.h"

// Class computing diffuse-convective local matrix in cartesian coordinates.

@interface FEMDiffuseConvectiveAnisotropic : NSObject

-(void)diffuseConvectiveComposeMassMatrix:(double **)massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double *)loadVector timeDerivativeTerm:(double *)nodalCT zeroDegreeTerm:(double *)nodalC0 convectionTerm:(double *)nodalC1 diffusionTerm:(double ***)nodalC2 phaseChange:(BOOL)phaseChange nodalTemperature:(double *)nodalTemperature enthalpy:(double *)enthalpy velocityX:(double *)ux velocitY:(double *)uy velocityZ:(double *)uz meshVeloX:(double *)mux meshVeloY:(double *)muy meshVeloZ:(double *)muz nodalViscosity:(double *)nodalviscosity nodaldensity:(double *)nodalDensity nodalPressure:(double *)nodalPressure nodalPressureDt:(double *)nodalPressureDt nodalPressureCoeff:(double *)nodalPressureCoeff compressible:(BOOL)compressible stabilize:(BOOL)stabilize useBubbles:(BOOL)useBubbles element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution core:(FEMCore *)core mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration materialModels:(FEMMaterialModels*)materialModels differentials:(FEMDifferentials *)differentials listUtilities:(FEMListUtilities *)listUtilities;
-(void)diffuseConvectiveBoundaryMatrix:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector dimensions:(Dimensions_t)dimensions loadVector:(double *)loadVector nodalAlpha:(double *)nodalAlpha element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh integration:(FEMNumericIntegration *)integration;

@end
