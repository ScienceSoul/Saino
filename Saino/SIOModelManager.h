//
//  SIOModelManager.h
//  Saino
//
//  Created by Hakime Seddik on 20/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

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

-(int)createModel:(NSString * __nonnull)dir;
-(int)openModel:(NSString * __nonnull)dir;
-(int)closeModel;

-(int)openStream:(NSFileHandle * __nullable)fstr name:(NSString * __nonnull)name mode:(NSString * __nonnull)mode;
-(int)closeStrem:(NSFileHandle * __nonnull)fstr;

-(int)makeDirectory:(NSString * __nonnull)dir;

-(void)deallocation;



@end
