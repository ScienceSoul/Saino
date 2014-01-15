//
//  FEMMaterialModels.h
//  Saino
//
//  Created by Seddik hakime on 13/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "FEMKernel.h"
#import "FEMModel.h"
#import "FEMMesh.h"
#import "FEMNumericIntegration.h"
#import "FEMListUtilities.h"

/***********************************************
    Class for material models of fluids mainly
***********************************************/
@interface FEMMaterialModels : NSObject

-(double)secondInvariantVelo:(double[3])velo dVelodx:(double[][3])dVelodx crtMatrix:(double[][3])crtMatrix symbols:(double[][3][3])symbols model:(FEMModel*)model;
-(double)effectiveViscosity:(double)viscosity density:(double)density velocityX:(double *)ux velocitY:(double *)uy velocityZ:(double *)uz element:(Element_t *)element nodes:(Nodes_t *)nodes numberOfNodes:(int)n numberOfPoints:(int)nd integrationU:(double)u integrationV:(double)v integrationW:(double)w muder:(double *)muder mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration;
-(double)effectiveConductivity:(double)conductivity density:(double)density element:(Element_t *)element temperature:(double *)temperature velocityX:(double *)ux velocitY:(double *)uy velocityZ:(double *)uz nodes:(Nodes_t *)nodes numberOfNodes:(int)n numberOfPoints:(int)nd integrationU:(double)u integrationV:(double)v integrationW:(double)w conductivityFlag:(NSString *)conductivityFlag kernel:(FEMKernel *)kernel mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration listUtilities:(FEMListUtilities *)listUtilities;

@end
