//
//  FEMPost.h
//  Saino
//
//  Created by Seddik hakime on 29/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMesh.h"

@interface FEMPost : NSObject

-(void)writeElmerPostFile:(NSString * __nonnull)postFile resultFile:(NSString * __nonnull)resultFile model:(FEMModel * __nonnull)model timeCount:(int)timeCount append:(BOOL * __nonnull)append;

@end
