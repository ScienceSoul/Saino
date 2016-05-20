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

-(void)writeString:(NSString * __nonnull)string toFileHandle:(NSFileHandle * __nonnull)fileHandle;
-(void)writeInteger:(int)number toFileHandle:(NSFileHandle * __nonnull)fileHandle;
-(void)writeDouble:(double)number toFileHandle:(NSFileHandle * __nonnull)fileHandle;
-(void)writeBytes:(const void * __nonnull)bytes length:(int)length toFileHandle:(NSFileHandle * __nonnull)fileHandle;
-(void)writeElmerPostFile:(NSString * __nonnull)postFile resultFile:(NSString * __nonnull)resultFile model:(FEMModel * __nonnull)model timeCount:(int)timeCount append:(BOOL * __nonnull)append;

@end
