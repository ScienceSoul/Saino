//
//  FEMBandwidthOptimize.h
//  Saino
//
//  Created by Seddik hakime on 08/01/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"

@interface FEMBandwidthOptimize : NSObject

-(int)computeBandWidthInListMatrix:(ListMatrix_t * __nonnull)list size:(int)n reorder:(int * __nullable)reorder invInitialReorder:(int * __nullable)invInitialReorder;
-(int)optimizeBandwidthInListMatrix:(ListMatrix_t * __nonnull)listMatrix permutation:(int * __nonnull)perm sizeOfPerm:(int)sizeOfPerm invInitialReorder:(int * __nonnull)invInitialReorder localNodes:(int)localNodes optimize:(BOOL)optimize useOptimized:(BOOL)useOptimized equation:(NSString * __nonnull)equation;

@end
