//
//  FEMCoordinateSystems.m
//  Saino
//
//  Created by Seddik hakime on 24/05/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMCoordinateSystems.h"

@interface FEMCoordinateSystems ()
-(double)FEMCoordinateSystems_cylindricalSquareRootMetricRadial:(double)r azimut:(double)z height:(double)t;
-(double)FEMCoordinateSystems_polarSquareRootMetricModel:(FEMModel *)model radial:(double)r azimut:(double)z height:(double)t;
@end

@implementation FEMCoordinateSystems

#pragma mark Private methods

-(double)FEMCoordinateSystems_cylindricalSquareRootMetricRadial:(double)r azimut:(double)z height:(double)t {
    
    double s;
    
    return s = r;
}

-(double)FEMCoordinateSystems_polarSquareRootMetricModel:(FEMModel *)model radial:(double)r azimut:(double)z height:(double)t {
    
    double s;
    
    if (model.dimension == 2) {
        s = sqrt( pow(r, 2.0) * pow(cos(t), 2.0) );
    } else {
        s = sqrt( pow(r, 4.0) * pow(cos(t), 2.0) );
    }
    
    return s;

}

#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
    }
    
    return self;
}

-(double)coordinateSquareRootMetricModel:(FEMModel *)model coordX:(double)x coordY:(double)y coordZ:(double)z {
    
    double sqrtMetric = 0.0;
    
    if (model.coordinates == cartesian) {
        sqrtMetric = 1.0;
    } else if (model.coordinates >= cylindric && model.coordinates <= axis_symmetric) {
        sqrtMetric = [self FEMCoordinateSystems_cylindricalSquareRootMetricRadial:x azimut:y height:z];
    } else if (model.coordinates == polar) {
        sqrtMetric = [self FEMCoordinateSystems_polarSquareRootMetricModel:model radial:x azimut:y height:z];
    }
    
    return sqrtMetric;
}

@end
