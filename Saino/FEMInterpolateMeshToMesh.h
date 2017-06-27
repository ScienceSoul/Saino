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

-(void)interpolateQMesh:(FEMMesh * _Nonnull)oldMesh toMesh:(FEMMesh * _Nonnull)newMesh oldVariables:(NSMutableArray * _Nullable)oldVar newVariables:(NSMutableArray * _Nullable)newVar model:(FEMModel * _Nonnull)model quadrantTree:(BOOL * _Nullable)useQuandrant projector:(FEMProjector * _Nullable)projector mask:(NSString * _Nullable)maskName nodesPresent:(BOOL * _Nullable)nodesPresent newMaskPerm:(int * _Nullable)newMaskPerm;
-(void)interpolateMesh:(FEMMesh * _Nonnull)oldMesh toMesh:(FEMMesh * _Nonnull)newMesh oldVariables:(NSMutableArray * _Nullable)oldVar newVariables:(NSMutableArray * _Nullable)newVar model:(FEMModel * _Nonnull)model quadrantTree:(BOOL * _Nullable)useQuandrant projector:(FEMProjector * _Nullable)projector mask:(NSString * _Nullable)maskName unfoundNodes:(BOOL * _Nullable)unfoundNodes;
-(FEMMatrix * _Nonnull)weightedProjectorMesh2:(FEMMesh * _Nonnull)bMesh2 mesh1:(FEMMesh * _Nonnull)bMesh1 inversePermutation2:(int * _Nonnull)invPerm2 sizeInversePermutation2:(int)sizeInversePermutation2 inversePermutation1:(int * _Nonnull)invPerm1 sizeInversePermutation1:(int)sizeInversePermutation1 useQuadrantTree:(BOOL)useQuadrantTree repeating:(BOOL)repeating antiRepeating:(BOOL)antiRepeating periodicScale:(double)periodicScale nodalJump:(BOOL)nodalJump model:(FEMModel * _Nonnull)model;

@end

