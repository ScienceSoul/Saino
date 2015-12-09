//
//  FEMInterpolateMeshToMesh.h
//  Saino
//
//  Created by Seddik hakime on 12/11/2015.
//  Copyright © 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMProjector.h"

@interface FEMInterpolateMeshToMesh : NSObject

-(void)interpolateQMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh oldVariables:(NSMutableArray *)oldVar newVariables:(NSMutableArray *)newVar model:(FEMModel *)model quadrantTree:(BOOL *)useQuandrant projector:(FEMProjector *)projector mask:(NSString *)maskName nodesPresent:(BOOL *)nodesPresent newMaskPerm:(int *)newMaskPerm;
-(void)interpolateMesh:(FEMMesh *)oldMesh toMesh:(FEMMesh *)newMesh oldVariables:(NSMutableArray *)oldVar newVariables:(NSMutableArray *)newVar model:(FEMModel *)model quadrantTree:(BOOL *)useQuandrant projector:(FEMProjector *)projector mask:(NSString *)maskName unfoundNodes:(BOOL *)unfoundNodes;
-(FEMMatrix *)weightedProjectorMesh2:(FEMMesh *)bMesh2 mesh1:(FEMMesh *)bMesh1 inversePermutation2:(int *)invPerm2 sizeInversePermutation2:(int)sizeInversePermutation2 inversePermutation1:(int *)invPerm1 sizeInversePermutation1:(int)sizeInversePermutation1 useQuadrantTree:(BOOL)useQuadrantTree repeating:(BOOL)repeating antiRepeating:(BOOL)antiRepeating periodicScale:(double)periodicScale nodalJump:(BOOL)nodalJump model:(FEMModel *)model;

@end

