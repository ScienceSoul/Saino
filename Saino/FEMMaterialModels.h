//
//  FEMMaterialModels.h
//  Saino
//
//  Created by Seddik hakime on 13/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "FEMCore.h"
#import "FEMModel.h"
#import "FEMMesh.h"
#import "FEMNumericIntegration.h"
#import "FEMListUtilities.h"

/***********************************************
 
    Class for material models of fluids mainly
 
***********************************************/
@interface FEMMaterialModels : NSObject

-(double)secondInvariantVelo:(double[3])velo dVelodx:(double[][3])dVelodx crtMatrix:(double[][3])crtMatrix symbols:(double[][3][3])symbols model:(FEMModel* __nonnull)model;
-(double)effectiveViscosity:(double)viscosity density:(double)density velocityX:(double * __nonnull)ux velocitY:(double * __nonnull)uy velocityZ:(double * __nonnull)uz element:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int)n numberOfPoints:(int)nd integrationU:(double)u integrationV:(double)v integrationW:(double)w muder:(double * __nullable)muder mesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration;
-(double)effectiveConductivity:(double)conductivity density:(double)density element:(Element_t * __nonnull)element temperature:(double * __nullable)temperature velocityX:(double * __nonnull)ux velocitY:(double * __nonnull)uy velocityZ:(double * __nonnull)uz nodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int)n numberOfPoints:(int)nd integrationU:(double)u integrationV:(double)v integrationW:(double)w conductivityFlag:(NSString * __nonnull)conductivityFlag core:(FEMCore * __nonnull)core mesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration listUtilities:(FEMListUtilities * __nonnull)listUtilities;

@end
