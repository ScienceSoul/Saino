//
//  SIOMeshIO.m
//  Saino
//
//  Created by Hakime Seddik on 28/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

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

-(id __nonnull)initWithParallelNumberOfProcessors:(int)procs processorID:(int)me
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

-(void)openMeshAtPath:(NSString * __nonnull)directory {
    
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

-(void)getMeshDescriptionNodeCount:(int * __nonnull)nodeCount elementCount:(int * __nonnull)elementCount boundaryElementCount:(int * __nonnull)boundaryElementCount usedElementTypes:(int * __nonnull)usedElementTypes elementTypeTags:(int * __nonnull)elementTypeTags elementCountByType:(int * __nonnull)elementCountByType {
    
    self.info = -1;
    self.info = [self.meshAgent readDescriptorNode:nodeCount element:elementCount boundaryElement:boundaryElementCount usedElementTypes:usedElementTypes usedElementTypeTags:elementTypeTags usedElementTypeCount:elementCountByType];
}

-(void)getMeshNodes:(int * __nonnull)tags coord:(double * __nonnull)coord {
    
    self.info = -1;
    self.info = [self.meshAgent readAllNodes:tags coord:coord];
}

-(void)getMeshElementConnection:(int * __nonnull)tag body:(int * __nonnull)body type:(int * __nonnull)type pdofs:(int * __nonnull)pdofs nodes:(int * __nonnull)nodes colorIndex:(int * __nullable)colorIndex parallelAssembly:(BOOL * __nullable)parallelAssembly {
    
    int part;
    if ([self.meshAgent readNextElementConnections:tag part:&part body:body type:type pdofs:pdofs nodes:nodes colorIndex:colorIndex parallelAssembly:(BOOL *)parallelAssembly] != -1) {
        self.info = 0;
    } else self.info = -1;
}

-(void)getMeshBoundaryElement:(int * __nonnull)tag boundary:(int * __nonnull)boundary leftElement:(int * __nonnull)leftElement rightElement:(int * __nonnull)rightElement type:(int * __nonnull)type nodes:(int * __nonnull)nodes coord:(double * __nonnull)coord {
    
    int part;
    if ([self.meshAgent readNextBoundaryElement:tag part:&part boundary:boundary leftElement:leftElement rightElement:rightElement type:type nodes:nodes coord:coord] != -1) {
        self.info = 0;
    } else self.info = -1;
}

-(void)getPartDesription:(int * __nonnull)sharedNodeCount {
    
    self.info = -1;
    self.info = [self.meshAgent readPartDescriptor:sharedNodeCount];
}

-(void)getPartNode:(int * __nonnull)tag constraint:(int * __nonnull)constraint coord:(double * __nonnull)coord partCount:(int * __nonnull)partCount parts:(int * __nonnull)parts {
    
    if ([self.meshAgent readSharedNode:tag constraint:constraint coord:coord partCount:partCount partitions:parts] != -1) {
        self.info = 0;
    } else self.info = -1;
}

@end
