//===----------------------------------------------------------------------===//
//  SIOMeshIO.m
//  Saino
//
//  Created by Hakime Seddik on 28/06/12.
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

#import "SIOMeshIO.h"
#import "memory.h"

@implementation SIOMeshIO

@synthesize info = _info;
@synthesize parallelState = _parallelState;
@synthesize modelManager = _modelManager;
@synthesize meshAgent = _meshAgent;

-(id)init
{
    self = [super init];
    if (self) {
        
        _parallelState = [[SIOInfoParallel alloc] init];
        _parallelState.parallel = NO;
        _parallelState.numProc = 1;
        _parallelState.myProc = 0;
        _modelManager = [[SIOModelManager alloc] init];
        if (_modelManager == nil) {
            _info = -1;
        } else {
            _info = 0;
        }
        _meshAgent = nil;
    }
    
    return self;
}

-(id _Nonnull)initWithParallelNumberOfProcessors:(int)procs processorID:(int)me
{
    self = [super init];
    if (self) {
        
        _parallelState = [[SIOInfoParallel alloc] init];
        _parallelState.parallel = YES;
        _parallelState.numProc = procs;
        _parallelState.myProc = me;
        _modelManager = [[SIOModelManager alloc] init];
        if (_modelManager == nil) {
            _info = -1;
        } else {
            _info = 0;
        }
    }
    
    return self;
}

-(void)openMeshAtPath:(NSString * _Nonnull)directory {
    
    if (directory == nil) {
        fatal("SIOMeshIO:openMeshAtPath", "No assignment to directory name (empty object pointer).");
    }
    _meshAgent = [[SIOMeshAgent alloc] initWithManager:self.modelManager split:self.parallelState.numProc part:self.parallelState.myProc];
    
    if (self.meshAgent != nil) {
        self.info = [self.meshAgent openMesh:directory];
    } else {
        self.info = -1;
    }
}

-(void)closeMesh {
    
    [self.meshAgent closeMesh];
    self.info = 0;
}

-(void)close {
    
    self.parallelState = nil;
    self.modelManager = nil;
    self.meshAgent = nil;
    self.info = 0;
}

-(void)getMeshDescriptionNodeCount:(int * _Nonnull)nodeCount elementCount:(int * _Nonnull)elementCount boundaryElementCount:(int * _Nonnull)boundaryElementCount usedElementTypes:(int * _Nonnull)usedElementTypes elementTypeTags:(int * _Nonnull)elementTypeTags elementCountByType:(int * _Nonnull)elementCountByType {
    
    self.info = -1;
    self.info = [self.meshAgent readDescriptorNode:nodeCount element:elementCount boundaryElement:boundaryElementCount usedElementTypes:usedElementTypes usedElementTypeTags:elementTypeTags usedElementTypeCount:elementCountByType];
}

-(void)getMeshNodes:(int * _Nonnull)tags coord:(double * _Nonnull)coord {
    
    self.info = -1;
    self.info = [self.meshAgent readAllNodes:tags coord:coord];
}

-(void)getMeshElementConnection:(int * _Nonnull)tag body:(int * _Nonnull)body type:(int * _Nonnull)type pdofs:(int * _Nonnull)pdofs nodes:(int * _Nonnull)nodes colorIndex:(int * _Nullable)colorIndex parallelAssembly:(BOOL * _Nullable)parallelAssembly {
    
    int part;
    if ([self.meshAgent readNextElementConnections:tag part:&part body:body type:type pdofs:pdofs nodes:nodes colorIndex:colorIndex parallelAssembly:parallelAssembly] != -1) {
        self.info = 0;
    } else self.info = -1;
}

-(void)getMeshBoundaryElement:(int * _Nonnull)tag boundary:(int * _Nonnull)boundary leftElement:(int * _Nonnull)leftElement rightElement:(int * _Nonnull)rightElement type:(int * _Nonnull)type nodes:(int * _Nonnull)nodes coord:(double * _Nonnull)coord {
    
    int part;
    if ([self.meshAgent readNextBoundaryElement:tag part:&part boundary:boundary leftElement:leftElement rightElement:rightElement type:type nodes:nodes coord:coord] != -1) {
        self.info = 0;
    } else self.info = -1;
}

-(void)getPartDesription:(int * _Nonnull)sharedNodeCount {
    
    self.info = -1;
    self.info = [self.meshAgent readPartDescriptor:sharedNodeCount];
}

-(void)getPartNode:(int * _Nonnull)tag constraint:(int * _Nonnull)constraint coord:(double * _Nonnull)coord partCount:(int * _Nonnull)partCount parts:(int * _Nonnull)parts {
    
    if ([self.meshAgent readSharedNode:tag constraint:constraint coord:coord partCount:partCount partitions:parts] != -1) {
        self.info = 0;
    } else self.info = -1;
}

@end
