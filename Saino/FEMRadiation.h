//
//  FEMRadiation.h
//  Saino
//
//  Created by Seddik hakime on 03/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"
#import "FEMMesh.h"

@interface FEMRadiation : NSObject

-(double)computeRadiationLoadModel:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh element:(Element_t * __nonnull)elememt temperature:(double * __nonnull)temperature reorder:(int * __nonnull)reorder emissivity:(double)emissivity angleFraction:(double * __nullable)angleFraction;
-(double)computeRadiationCoeffModel:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh element:(Element_t * __nonnull)element index:(int)k;

@end
