//===----------------------------------------------------------------------===//
//  FEMBandwidthOptimize.h
//  Saino
//
//  Created by Seddik hakime on 08/01/13.
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
#import "Constructors.h"

@interface FEMBandwidthOptimize : NSObject

-(int)computeBandWidthInListMatrix:(ListMatrix_t * _Nonnull)list size:(int)n reorder:(int * _Nullable)reorder invInitialReorder:(int * _Nullable)invInitialReorder;
-(int)optimizeBandwidthInListMatrix:(ListMatrix_t * _Nonnull)listMatrix permutation:(int * _Nonnull)perm sizeOfPerm:(int)sizeOfPerm invInitialReorder:(int * _Nonnull)invInitialReorder localNodes:(int)localNodes optimize:(BOOL)optimize useOptimized:(BOOL)useOptimized equation:(NSString * _Nonnull)equation;

@end
