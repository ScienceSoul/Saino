//===----------------------------------------------------------------------===//
//  FEMMaterialModels.h
//  Saino
//
//  Created by Seddik hakime on 13/06/13.
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
#import "FEMCore.h"
#import "FEMModel.h"
#import "FEMMesh.h"
#import "FEMNumericIntegration.h"
#import "FEMListUtilities.h"

/***********************************************
 
    Class for material models of fluids mainly
 
***********************************************/
@interface FEMMaterialModels : NSObject

-(double)secondInvariantVelo:(double[_Nonnull 3])velo dVelodx:(double[_Nonnull][3])dVelodx crtMatrix:(double[_Nullable][3])crtMatrix symbols:(double[_Nullable][3][3])symbols model:(FEMModel* _Nonnull)model;
-(double)effectiveViscosity:(double)viscosity density:(double)density velocityX:(double * _Nonnull)ux velocitY:(double * _Nonnull)uy velocityZ:(double * _Nonnull)uz element:(Element_t * _Nonnull)element nodes:(Nodes_t * _Nonnull)nodes numberOfNodes:(int)n numberOfPoints:(int)nd integrationU:(double)u integrationV:(double)v integrationW:(double)w muder:(double * _Nullable)muder mesh:(FEMMesh * _Nonnull)mesh model:(FEMModel * _Nonnull)model integration:(FEMNumericIntegration * _Nonnull)integration;
-(double)effectiveConductivity:(double)conductivity density:(double)density element:(Element_t * _Nonnull)element temperature:(double * _Nullable)temperature velocityX:(double * _Nonnull)ux velocitY:(double * _Nonnull)uy velocityZ:(double * _Nonnull)uz nodes:(Nodes_t * _Nonnull)nodes numberOfNodes:(int)n numberOfPoints:(int)nd integrationU:(double)u integrationV:(double)v integrationW:(double)w conductivityFlag:(NSString * _Nonnull)conductivityFlag core:(FEMCore * _Nonnull)core mesh:(FEMMesh * _Nonnull)mesh model:(FEMModel * _Nonnull)model integration:(FEMNumericIntegration * _Nonnull)integration listUtilities:(FEMListUtilities * _Nonnull)listUtilities;

@end
