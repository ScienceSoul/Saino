//===----------------------------------------------------------------------===//
//  FEMCoordinateSystems.m
//  Saino
//
//  Created by Seddik hakime on 24/05/13.
//  Copyright (c) 2013 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
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

#import "FEMCoordinateSystems.h"
#import "Utils.h"

@interface FEMCoordinateSystems ()
-(double)FEMCoordinateSystems_cylindricalSquareRootMetricRadial:(double)r azimut:(double)z height:(double)t;
-(void)FEMCoordinateSystems_cylindricalMetric:(double[][3])metric radial:(double)r azimut:(double)z height:(double)t;
-(void)FEMCoordinateSystems_cylindricalSymbols:(double[][3][3])symbols radial:(double)r azimut:(double)z height:(double)t;
-(void)FEMCoordinateSystems_cylindricalDerivSymbols:(double[][3][3][3])dSymbols radial:(double)r azimut:(double)z height:(double)t;

-(double)FEMCoordinateSystems_polarSquareRootMetricModel:(FEMModel * _Nonnull)model radial:(double)r azimut:(double)z height:(double)t;
-(void)FEMCoordinateSystems_polarModel:(FEMModel * _Nonnull)model metric:(double[][3])metric radial:(double)r azimut:(double)z height:(double)t;
-(void)FEMCoordinateSystems_polarModel:(FEMModel * _Nonnull)model symbols:(double[][3][3])symbols radial:(double)r azimut:(double)z height:(double)t;
-(void)FEMCoordinateSystems_polarModel:(FEMModel * _Nonnull)model derivSymbols:(double[][3][3][3])dSymbols radial:(double)r azimut:(double)z height:(double)t;
@end

@implementation FEMCoordinateSystems

#pragma mark Private methods

-(double)FEMCoordinateSystems_cylindricalSquareRootMetricRadial:(double)r azimut:(double)z height:(double)t {
    
    double s;
    
    return s = r;
}

-(void)FEMCoordinateSystems_cylindricalMetric:(double[][3])metric radial:(double)r azimut:(double)z height:(double)t {
    
    memset( *metric, 0.0, (3*3)*sizeof(double) );
    metric[0][0] = 1.0;
    metric[1][1] = 1.0;
    metric[2][2] = 1.0;
    if (r != 0.0) metric[2][2] = 1.0 / pow(r, 2.0);
}

-(void)FEMCoordinateSystems_cylindricalSymbols:(double[][3][3])symbols radial:(double)r azimut:(double)z height:(double)t {
    
    memset(**symbols, 0.0, (3*3*3)*sizeof(double) );
    symbols[2][2][0] = -r;
    
    if (r != 0.0) {
        symbols[0][2][2] = 1.0 / r;
        symbols[2][0][2] = 1.0 / r;
    }
}

-(void)FEMCoordinateSystems_cylindricalDerivSymbols:(double[][3][3][3])dSymbols radial:(double)r azimut:(double)z height:(double)t {
    
    memset(***dSymbols, 0.0, (3*3*3*3)*sizeof(double) );
    dSymbols[2][2][0][0] = -1.0;
    
    if (r != 0.0) {
        dSymbols[0][2][2][0] = -1.0 / pow(r, 2.0);
        dSymbols[2][0][2][0] = -1.0 / pow(r, 2.0);
    }
}

-(double)FEMCoordinateSystems_polarSquareRootMetricModel:(FEMModel * _Nonnull)model radial:(double)r azimut:(double)z height:(double)t {
    
    double s;
    
    if (model.dimension == 2) {
        s = sqrt( pow(r, 2.0) * pow(cos(t), 2.0) );
    } else {
        s = sqrt( pow(r, 4.0) * pow(cos(t), 2.0) );
    }
    
    return s;
}

-(void)FEMCoordinateSystems_polarModel:(FEMModel * _Nonnull)model metric:(double[][3])metric radial:(double)r azimut:(double)z height:(double)t {
    
    int i;
    
    memset( *metric, 0.0, (3*3)*sizeof(double) );
    for (i=0; i<3; i++) {
        metric[i][i] = 1.0;
    }
    
    if (r != 0.0) {
        metric[1][1] = 1.0 / ( pow(r, 2.0) * pow(cos(t), 2.0) );
        if (model.dimension == 3) {
            metric[2][2] = 1.0 / pow(r, 2.0);
        }
    }
}

-(void)FEMCoordinateSystems_polarModel:(FEMModel * _Nonnull)model symbols:(double[][3][3])symbols radial:(double)r azimut:(double)z height:(double)t {
    
    memset(**symbols, 0.0, (3*3*3)*sizeof(double) );
    symbols[1][1][0] = -r * pow(cos(t), 2.0);
    if (r != 0.0) {
        symbols[0][1][1] = 1.0 / r;
        symbols[1][0][1] = 1.0 / r;
    }
    
    if (model.dimension == 3) {
        symbols[2][2][0] = -r;
        symbols[1][1][2] = sin(t) * cos(t);
        
        symbols[1][2][1] = -tan(t);
        symbols[2][1][1] = -tan(t);
    
        if (r != 0.0) {
            symbols[2][0][2] = 1.0 / r;
            symbols[0][2][2] = 1.0 / r;
        }
    }
}

-(void)FEMCoordinateSystems_polarModel:(FEMModel * _Nonnull)model derivSymbols:(double[][3][3][3])dSymbols radial:(double)r azimut:(double)z height:(double)t {
    
    memset(***dSymbols, 0.0, (3*3*3*3)*sizeof(double) );
    dSymbols[1][1][0][0] = -pow(cos(t), 2.0);
    if (r != 0.0) {
        dSymbols[0][1][1][0] = -1.0 / pow(r, 2.0);
        dSymbols[1][0][1][0] = -1.0 / pow(r, 2.0);
    }
    
    if (model.dimension == 3) {
        dSymbols[1][1][0][2] = -2.0 * r * sin(t) * cos(t);
        dSymbols[2][2][0][0] = -1.0;
        dSymbols[1][1][2][2] = pow(cos(t), 2.0) - pow(sin(t), 2.0);
        
        dSymbols[1][2][1][2] = -1.0 / pow(cos(t), 2.0);
        dSymbols[2][1][1][2] = -1.0 / pow(cos(t), 2.0);
        
        if (r != 0.0) {
            dSymbols[0][2][2][0] = -1.0 / pow(r, 2.0);
            dSymbols[2][0][2][0] = -1.0 / pow(r, 2.0);
        }
    }
}

#pragma mark Public methods

@synthesize coordinates = _coordinates;

- (id)init
{
    self = [super init];
    if (self) {
        _coordinates = cartesian;
    }
    
    return self;
}

-(double)coordinateSquareRootMetricModel:(FEMModel * _Nonnull)model coordX:(double)x coordY:(double)y coordZ:(double)z {
    
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

-(void)coordinateSystemInfoModel:(FEMModel * _Nonnull)model metric:(double[_Nonnull][3])metric sqrtMetric:(double * _Nonnull)sqrtMetric symbols:(double[_Nonnull][3][3])symbols dSymbols:(double[_Nonnull][3][3][3])dSymbols coordX:(double)x coordY:(double)y coordZ:(double)z {
    
    int i;
    
    if (model.coordinates == cartesian) {
       memset( *metric, 0.0, (3*3)*sizeof(double) );         // Fixed size so we know how to initialize it
        for (i=0; i<3; i++) {
            metric[i][i] = 1.0;
        }
        *sqrtMetric = 1.0;
        memset(**symbols, 0.0, (3*3*3)*sizeof(double) );     // Fixed size so we know how to initialize it
        memset(***dSymbols, 0.0, (3*3*3*3)*sizeof(double) ); // Fixed size so we know how to initialize it
        
    } else if (model.coordinates >= cylindric && model.coordinates <= axis_symmetric) {
        *sqrtMetric = [self FEMCoordinateSystems_cylindricalSquareRootMetricRadial:x azimut:y height:z];
        [self FEMCoordinateSystems_cylindricalMetric:metric radial:x azimut:y height:z];
        [self FEMCoordinateSystems_cylindricalSymbols:symbols radial:x azimut:y height:z];
        [self FEMCoordinateSystems_cylindricalDerivSymbols:dSymbols radial:x azimut:y height:z];
    } else if (model.coordinates == polar) {
        *sqrtMetric = [self FEMCoordinateSystems_polarSquareRootMetricModel:model radial:x azimut:y height:z];
        [self FEMCoordinateSystems_polarModel:model metric:metric radial:x azimut:y height:z];
        [self FEMCoordinateSystems_polarModel:model symbols:symbols radial:x azimut:y height:z];
        [self FEMCoordinateSystems_polarModel:model derivSymbols:dSymbols radial:x azimut:y height:z];
    }
}

@end
