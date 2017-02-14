//===----------------------------------------------------------------------===//
//  FEMRadiation.h
//  Saino
//
//  Created by Seddik hakime on 03/06/13.
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
#import "FEMMesh.h"

@interface FEMRadiation : NSObject

-(double)computeRadiationLoadModel:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh element:(Element_t * __nonnull)elememt temperature:(double * __nonnull)temperature reorder:(int * __nonnull)reorder emissivity:(double)emissivity angleFraction:(double * __nullable)angleFraction;
-(double)computeRadiationCoeffModel:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh element:(Element_t * __nonnull)element index:(int)k;

@end
