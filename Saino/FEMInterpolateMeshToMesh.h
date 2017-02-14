//===----------------------------------------------------------------------===//
//  FEMInterpolateMeshToMesh.h
//  Saino
//
//  Created by Seddik hakime on 12/11/2015.
//  Copyright © 2015 ScienceSoul. All rights reserved.
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
#import "FEMProjector.h"

@interface FEMInterpolateMeshToMesh : NSObject

-(void)interpolateQMesh:(FEMMesh * __nonnull)oldMesh toMesh:(FEMMesh * __nonnull)newMesh oldVariables:(NSMutableArray * __nullable)oldVar newVariables:(NSMutableArray * __nullable)newVar model:(FEMModel * __nonnull)model quadrantTree:(BOOL * __nullable)useQuandrant projector:(FEMProjector * __nullable)projector mask:(NSString * __nullable)maskName nodesPresent:(BOOL * __nullable)nodesPresent newMaskPerm:(int * __nullable)newMaskPerm;
-(void)interpolateMesh:(FEMMesh * __nonnull)oldMesh toMesh:(FEMMesh * __nonnull)newMesh oldVariables:(NSMutableArray * __nullable)oldVar newVariables:(NSMutableArray * __nullable)newVar model:(FEMModel * __nonnull)model quadrantTree:(BOOL * __nullable)useQuandrant projector:(FEMProjector * __nullable)projector mask:(NSString * __nullable)maskName unfoundNodes:(BOOL * __nullable)unfoundNodes;
-(FEMMatrix * __nonnull)weightedProjectorMesh2:(FEMMesh * __nonnull)bMesh2 mesh1:(FEMMesh * __nonnull)bMesh1 inversePermutation2:(int * __nonnull)invPerm2 sizeInversePermutation2:(int)sizeInversePermutation2 inversePermutation1:(int * __nonnull)invPerm1 sizeInversePermutation1:(int)sizeInversePermutation1 useQuadrantTree:(BOOL)useQuadrantTree repeating:(BOOL)repeating antiRepeating:(BOOL)antiRepeating periodicScale:(double)periodicScale nodalJump:(BOOL)nodalJump model:(FEMModel * __nonnull)model;

@end

