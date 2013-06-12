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

-(double)computeRadiationLoadModel:(FEMModel *)model mesh:(FEMMesh *)mesh element:(Element_t *)elememt temperature:(double *)temperature reorder:(int *)reorder emissivity:(double)emissivity angleFraction:(double *)angleFraction;
-(double)computeRadiationCoeffModel:(FEMModel *)model mesh:(FEMMesh *)mesh element:(Element_t *)element index:(int)k;

@end
