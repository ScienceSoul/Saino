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

-(void)writeString:(NSString * _Nonnull)string toFileHandle:(NSFileHandle * _Nonnull)fileHandle;
-(void)writeInteger:(int)number toFileHandle:(NSFileHandle * _Nonnull)fileHandle;
-(void)writeDouble:(double)number toFileHandle:(NSFileHandle * _Nonnull)fileHandle;
-(void)writeBytes:(const void * _Nonnull)bytes length:(int)length toFileHandle:(NSFileHandle * _Nonnull)fileHandle;
-(void)writeElmerPostFile:(NSString * _Nonnull)postFile resultFile:(NSString * _Nonnull)resultFile model:(FEMModel * _Nonnull)model timeCount:(int)timeCount append:(BOOL * _Nonnull)append;

@end
