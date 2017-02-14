//===----------------------------------------------------------------------===//
//  FEMDifferentials.h
//  Saino
//
//  Created by Seddik hakime on 18/06/13.
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
#import "FEMNumericIntegration.h"
#import "FEMListUtilities.h"
#import "FEMUtilities.h"
#import "FEMCoordinateSystems.h"

/*******************************************************************************
 
    This class contains some built-in material laws, and also some
    vector utilities, curl, dot, cross, etc. Some of these may be 
    of no use currently.
 
*******************************************************************************/
@interface FEMDifferentials : NSObject

-(void)lorentzForceElement:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int)n integrationU:(double)u integrationV:(double)v integrationW:(double)w lorentzForce:(double * __nonnull)lorentzForce mesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration coordinateSystems:(FEMCoordinateSystems * __nonnull)coordinateSystems listUtilities:(FEMListUtilities * __nonnull)listUtilities utilities:(FEMUtilities * __nonnull)utilities;
-(double)jouleHeatElement:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int)n integrationU:(double)u integrationV:(double)v integrationW:(double)w mesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration listUtilities:(FEMListUtilities * __nonnull)listUtilities;

@end
