//
//  FEMCoordinateSystems.h
//  Saino
//
//  Created by Seddik hakime on 24/05/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"

@interface FEMCoordinateSystems : NSObject {
    
    int _coordinates;
}

@property(nonatomic, assign) int coordinates;

-(double)coordinateSquareRootMetricModel:(FEMModel * __nonnull)model coordX:(double)x coordY:(double)y coordZ:(double)z;
-(void)coordinateSystemInfoModel:(FEMModel * __nonnull)model metric:(double[][3])metric sqrtMetric:(double * __nonnull)sqrtMetric symbols:(double[][3][3])symbols dSymbols:(double[][3][3][3])dSymbols coordX:(double)x coordY:(double)y coordZ:(double)z;

@end
