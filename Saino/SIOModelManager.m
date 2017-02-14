//===----------------------------------------------------------------------===//
//  SIOModelManager.m
//  Saino
//
//  Created by Hakime Seddik on 20/06/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
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

#import "SIOModelManager.h"

@interface SIOModelManager ()

-(int)SIOModelManager_mkdir:(NSString * __nonnull)dir;
-(int)SIOModelManager_chdir:(NSString * __nonnull)dir;
-(int)SIOModelManager_checkModel:(NSString * __nonnull)model;

@end

@implementation SIOModelManager

@synthesize rundir;
@synthesize modeldir;
@synthesize modelname;

#pragma Private metods

-(int)SIOModelManager_mkdir:(NSString * __nonnull)dir {
    
    int rc;
    
    rc = mkdir([dir UTF8String], S_IRWXU | S_IRWXG);
    
    if (rc == -1) {
        switch (errno) {
            case EEXIST:
                return 0;
                break;
            default:
                fprintf(stdout, "SIOModelManager:SIOModelManager_mkdir: unexpected error at mkdir.\n");
                break;
        }
        return -1;
    }
    
    return 0;
}


-(int)SIOModelManager_chdir:(NSString * __nonnull)dir {
    
    int rc;
    
    rc = chdir([dir UTF8String]);
    
    if (rc == -1) {
        switch (errno) {
            case EACCES:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: permissions denied for directory:%s.\n", [dir UTF8String]);
                break;
            case EIO:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: i/o error in directory:%s.\n", [dir UTF8String]);
                break;
            case ENOENT:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: no such directory:%s.\n", [dir UTF8String]);
                break;
            case ENOTDIR:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: not a directory:%s.\n", [dir UTF8String]);
                break;
            default:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: unexpected error at chdir.\n");
                break;
        }
        return -1;
    }
    
    return 0;
}


-(int)SIOModelManager_checkModel:(NSString * __nonnull)model {
    
    int rc, rc_access;
    struct stat buf;
    
    // Get file status
    rc = stat([model UTF8String], &buf);
    
    if (rc == -1) {
        switch (errno) {
            case EACCES:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: permissions denied for:%s.\n", [model UTF8String]);
                break;
            case EIO:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: i/o error in:%s.\n", [model UTF8String]);
                break;
            case ENOENT:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: no such model:%s.\n", [model UTF8String]);
                break;
            case ENOTDIR:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: not a directory:%s.\n", [model UTF8String]);
                break;
            default:
                fprintf(stdout, "SIOModelManager:SIOModelManager_chdir: unexpected error at chdir.\n");
                break;
        }
        return -1;
    }
    
    // Is model a directory?
    rc = S_ISDIR(buf.st_mode);
    // Check access permissions
    rc_access = access([model UTF8String], R_OK | W_OK | X_OK);
    
    if (rc) {
        // We need read/write/exec permissions, however, since we could stat, we can search
        if (rc_access == -1) {
            fprintf(stdout, "SIOModelManager:SIOModelManager_checkModel: no permission to operate in %s.\n", [model UTF8String]);
            return -1;
        }
    } else {
        fprintf(stdout, "SIOModelManager:SIOModelManager_checkModel: %s is not a directory.\n", [model UTF8String]);
        return -1;
    }
    
    return 0;
}

#pragma mark Public methos

- (id)init
{

    NSFileManager *fileMrg;
    
    self = [super init];
    if (self) {
        // We must remember the current directory so that we can safely return
        // to it after the data base has been closed
        fileMrg = [NSFileManager defaultManager];
        self.rundir = [fileMrg currentDirectoryPath];
        self.modeldir = nil;
        self.modelname = nil;
        
        umask(S_IRWXO);
    }
    
    return self;
}

-(void)deallocation {
    
    // Return to current directory
    [self SIOModelManager_chdir:self.rundir];
    
}

-(int)createModel:(NSString * __nonnull)dir {
    
    self.modeldir = [NSString stringWithString:dir]; 
    
    if ([self SIOModelManager_chdir:self.modeldir] == -1) {
        return -1;
    }
    if ([self SIOModelManager_mkdir:self.modeldir] == -1) {
        return -1;
    }
    if ([self SIOModelManager_mkdir:self.modeldir]) {
        return -1;
    }
    
    return 0;
}

-(int)openModel:(NSString * __nonnull)dir {
    
    self.modeldir = [NSString stringWithString:dir];
    
    if ([self SIOModelManager_chdir:self.modeldir] == -1) {
        return -1;
    }
    if ([self SIOModelManager_checkModel:self.modeldir] == -1) {
        return -1;
    }
    if ([self SIOModelManager_chdir:self.modeldir] == -1) {
        return -1;
    }
    
    return 0;
}

-(int)closeModel {
    
    return 0;
}

-(int)openStream:(NSFileHandle * __nullable)fstr name:(NSString * __nonnull)name mode:(NSString * __nonnull)mode {
    
    if ([mode isEqualToString:@"read"]) {
        fstr = [NSFileHandle fileHandleForReadingAtPath:name];
    } else if ([mode isEqualToString:@"write"]) {
        fstr = [NSFileHandle fileHandleForWritingAtPath:name];
    } else {
        fprintf(stdout, "SIOModelManager:openStream: error, operation not supported.\n");
        return -1;
    }
    
    if (fstr == nil) {
        fprintf(stdout, "SIOModelManager:openStream: could not open %s.\n", [name UTF8String]);
        return -1;
    }
    
    return 0;
}

-(int)closeStrem:(NSFileHandle * __nonnull)fstr {
    
    [fstr closeFile];
    return 0;
}

-(int)makeDirectory:(NSString * __nonnull)dir {
    
    return [self SIOModelManager_mkdir:dir];
}

@end
