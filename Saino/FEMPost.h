//===----------------------------------------------------------------------===//
//  FEMPost.h
//  Saino
//
//  Created by Seddik hakime on 29/03/13.
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
#import "FEMMesh.h"

@interface FEMPost : NSObject

-(void)writeString:(NSString * __nonnull)string toFileHandle:(NSFileHandle * __nonnull)fileHandle;
-(void)writeInteger:(int)number toFileHandle:(NSFileHandle * __nonnull)fileHandle;
-(void)writeDouble:(double)number toFileHandle:(NSFileHandle * __nonnull)fileHandle;
-(void)writeBytes:(const void * __nonnull)bytes length:(int)length toFileHandle:(NSFileHandle * __nonnull)fileHandle;
-(void)writeElmerPostFile:(NSString * __nonnull)postFile resultFile:(NSString * __nonnull)resultFile model:(FEMModel * __nonnull)model timeCount:(int)timeCount append:(BOOL * __nonnull)append;

@end
