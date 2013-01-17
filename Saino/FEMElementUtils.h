//
//  FEMElementUtils.h
//  Saino
//
//  Created by Seddik hakime on 28/12/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMesh.h"
#import "FEMModel.h"
#import "FEMSolution.h"
#import "FEMMatrix.h"
#import "FEMUtilities.h"
#import "FEMListMatrix.h"
#import "FEMBoundaryCondition.h"
#import "FEMBandwidthOptimize.h"

@interface FEMElementUtils : NSObject

-(FEMMatrix *)createMatrixInModel:(FEMModel *)model forSolution:(FEMSolution *)solution mesh:(FEMMesh *)mesh dofs:(int)dofs permutation:(int *)perm sizeOfPermutation:(int)permSize matrixFormat:(int)matrixFormat optimizeBandwidth:(BOOL)optimizeBandwidth equationName:(NSString *)equation discontinuousGalerkinSolution:(BOOL *)dgSolution globalBubbles:(BOOL *)gbBubbles;

@end
