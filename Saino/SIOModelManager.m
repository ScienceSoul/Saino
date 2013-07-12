//
//  SIOModelManager.m
//  Saino
//
//  Created by Hakime Seddik on 20/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "SIOModelManager.h"

@interface SIOModelManager ()

-(int)SIOModelManager_mkdir:(NSString *)dir;
-(int)SIOModelManager_chdir:(NSString *)dir;
-(int)SIOModelManager_checkModel:(NSString *)model;

@end

@implementation SIOModelManager

@synthesize rundir;
@synthesize modeldir;
@synthesize modelname;

#pragma Private metods

-(int)SIOModelManager_mkdir:(NSString *)dir {
    
    int rc;
    
    rc = mkdir([dir UTF8String], S_IRWXU | S_IRWXG);
    
    if (rc == -1) {
        switch (errno) {
            case EEXIST:
                return 0;
                break;
            default:
                NSLog(@"SIOModelManager:SIOModelManager_mkdir: unexpected error at mkdir\n");
                break;
        }
        return -1;
    }
    
    return 0;
}


-(int)SIOModelManager_chdir:(NSString *)dir {
    
    int rc;
    
    rc = chdir([dir UTF8String]);
    
    if (rc == -1) {
        switch (errno) {
            case EACCES:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: permissions denied for directory:%@\n", dir);
                break;
            case EIO:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: i/o error in directory:%@\n", dir);
                break;
            case ENOENT:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: no such directory:%@\n", dir);
                break;
            case ENOTDIR:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: not a directory:%@\n", dir);
                break;
            default:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: unexpected error at chdir\n");
                break;
        }
        return -1;
    }
    
    return 0;
}


-(int)SIOModelManager_checkModel:(NSString *)model {
    
    int rc, rc_access;
    struct stat buf;
    
    // Get file status
    rc = stat([model UTF8String], &buf);
    
    if (rc == -1) {
        switch (errno) {
            case EACCES:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: permissions denied for:%@\n", model);
                break;
            case EIO:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: i/o error in:%@\n", model);
                break;
            case ENOENT:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: no such model:%@\n", model);
                break;
            case ENOTDIR:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: not a directory:%@\n", model);
                break;
            default:
                NSLog(@"SIOModelManager:SIOModelManager_chdir: unexpected error at chdir\n");
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
            NSLog(@"SIOModelManager:SIOModelManager_checkModel: no permission to operate in %@\n", model);
            return -1;
        }
    } else {
        NSLog(@"SIOModelManager:SIOModelManager_checkModel: %@ is not a directory\n", model);
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
        
        umask(S_IRWXO);
    }
    
    return self;
}

-(void)deallocation {
    
    // Return to current directory
    [self SIOModelManager_chdir:self.rundir];
    
}

-(int)createModel:(NSString *)dir {
    
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

-(int)openModel:(NSString *)dir {
    
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

-(int)openStream:(NSFileHandle *)fstr name:(NSString *)name mode:(NSString *)mode {
    
    if ([mode isEqualToString:@"read"]) {
        fstr = [NSFileHandle fileHandleForReadingAtPath:name];
    } else if ([mode isEqualToString:@"write"]) {
        fstr = [NSFileHandle fileHandleForWritingAtPath:name];
    } else {
        NSLog(@"SIOModelManager:openStream: error, operation not supported");
        return -1;
    }
    
    if (fstr == nil) {
        NSLog(@"SIOModelManager:openStream: could not open %@.", name);
        return -1;
    }
    
    return 0;
}

-(int)closeStrem:(NSFileHandle *)fstr {
    
    [fstr closeFile];
    return 0;
}

-(int)makeDirectory:(NSString *)dir {
    
    return [self SIOModelManager_mkdir:dir];
}

@end
