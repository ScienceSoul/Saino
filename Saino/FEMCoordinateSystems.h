//===----------------------------------------------------------------------===//
//  FEMCoordinateSystems.h
//  Saino
//
//  Created by Seddik hakime on 24/05/13.
//  Copyright (c) 2013 ScienceSoul. All rights reserved.
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
#import "FEMModel.h"

@interface FEMCoordinateSystems : NSObject {
    
    int _coordinates;
}

@property(nonatomic, assign) int coordinates;

-(double)coordinateSquareRootMetricModel:(FEMModel * __nonnull)model coordX:(double)x coordY:(double)y coordZ:(double)z;
-(void)coordinateSystemInfoModel:(FEMModel * __nonnull)model metric:(double[_Nonnull][3])metric sqrtMetric:(double * __nonnull)sqrtMetric symbols:(double[_Nonnull][3][3])symbols dSymbols:(double[_Nonnull][3][3][3])dSymbols coordX:(double)x coordY:(double)y coordZ:(double)z;

@end
