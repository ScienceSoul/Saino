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

@property (nonatomic, strong) NSString *rundir;
@property (nonatomic, strong) NSString *modeldir;
@property (nonatomic, strong) NSString *modelname;

-(int)createModel:(NSString *)dir;
-(int)openModel:(NSString *)dir;
-(int)closeModel;

-(int)openStream:(NSFileHandle *)fstr: (NSString *)name: (NSString *)mode;
-(int)closeStrem:(NSFileHandle *)fstr;

-(int)makeDirectory:(NSString *)dir;

-(void)deallocation;



@end
