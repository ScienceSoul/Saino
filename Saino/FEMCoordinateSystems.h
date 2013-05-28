//
//  FEMCoordinateSystems.h
//  Saino
//
//  Created by Seddik hakime on 24/05/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"

@interface FEMCoordinateSystems : NSObject

-(double)coordinateSquareRootMetricModel:(FEMModel *)model coordX:(double)x coordY:(double)y coordZ:(double)z;

@end
