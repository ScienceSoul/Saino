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

-(void)writeElmerPostFile:(NSString *)postFile resultFile:(NSString *)resultFile model:(FEMModel *)model timeCount:(int)timeCount append:(BOOL *)append;

@end
