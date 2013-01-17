//
//  FEMBandwidthOptimize.h
//  Saino
//
//  Created by Seddik hakime on 08/01/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMesh.h"
#import "FEMModel.h"
#import "FEMSolution.h"
#import "FEMMatrix.h"

@interface FEMBandwidthOptimize : NSObject

-(int)computeBandWidthInListMatrix:(ListMatrix_t *)list size:(int)n reorder:(int *)reorder invInitialReorder:(int *)invInitialReorder;
-(int)optimizeBandwidthInListMatrix:(ListMatrix_t *)listMatrix permutation:(int *)perm sizeOfPerm:(int)sizeOfPerm invInitialReorder:(int *)invInitialReorder localNodes:(int)localNodes optimize:(BOOL)optimize useOptimized:(BOOL)useOptimized equation:(NSString *)equation;

@end
