//
//  FEMDiffuseConvectiveGeneralAnisotropic.h
//  Saino
//
//  Created by Seddik hakime on 04/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"
#import "FEMMesh.h"

@interface FEMDiffuseConvectiveGeneralAnisotropic : NSObject

-(void)diffuseConvectiveGeneralBoundaryMatrix:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector dimensions:(Dimensions_t)dimensions loadVector:(double *)loadVector nodalAlpha:(double *)nodalAlpha element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes model:(FEMModel *)model mesh:(FEMMesh *)mesh;

@end
