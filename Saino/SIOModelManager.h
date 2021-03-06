//===----------------------------------------------------------------------===//
//  SIOModelManager.h
//  Saino
//
//  Created by Hakime Seddik on 20/06/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
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
#import <sys/types.h>
#import <sys/stat.h>
#import <stdlib.h>
#import <errno.h>
#import <unistd.h>
#import "sio_config.h"

@interface SIOModelManager : NSObject 

@property (nonatomic, strong, nonnull) NSString *rundir;
@property (nonatomic, strong, nullable) NSString *modeldir;
@property (nonatomic, strong, nullable) NSString *modelname;

-(int)createModel:(NSString * _Nonnull)dir;
-(int)openModel:(NSString * _Nonnull)dir;
-(int)closeModel;

-(int)openStream:(NSFileHandle * _Nullable)fstr name:(NSString * _Nonnull)name mode:(NSString * _Nonnull)mode;
-(int)closeStrem:(NSFileHandle * _Nonnull)fstr;

-(int)makeDirectory:(NSString * _Nonnull)dir;

-(void)deallocation;



@end
