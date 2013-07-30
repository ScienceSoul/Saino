//
//  SIOMeshIO.m
//  Saino
//
//  Created by Hakime Seddik on 28/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "SIOMeshIO.h"

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
    }
    
    return self;
}

-(id)initWithParallelNumberOfProcessors:(int)procs processorID:(int)me
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

-(void)openMeshAtPath:(NSString *)directory {
    
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
    
    self.modelManager = nil;
    self.info = 0;
}

-(void)getMeshDescriptionNodeCount:(int *)nodeCount elementCount:(int *)elementCount boundaryElementCount:(int *)boundaryElementCount usedElementTypes:(int *)usedElementTypes elementTypeTags:(int *)elementTypeTags elementCountByType:(int *)elementCountByType {
    
    self.info = -1;
    self.info = [self.meshAgent readDescriptorNode:nodeCount element:elementCount boundaryElement:boundaryElementCount usedElementTypes:usedElementTypes usedElementTypeTags:elementTypeTags usedElementTypeCount:elementCountByType];
}

-(void)getMeshNodes:(int *)tags coord:(double *)coord {
    
    self.info = -1;
    self.info = [self.meshAgent readAllNodes:tags coord:coord];
}

-(void)getMeshElementConnection:(int *)tag body:(int *)body type:(int *)type pdofs:(int *)pdofs nodes:(int *)nodes {
    
    int part;
    if ([self.meshAgent readNextElementConnections:tag part:&part body:body type:type pdofs:pdofs nodes:nodes] != -1) {
        self.info = 0;
    } else self.info = -1;
}

-(void)getMeshBoundaryElement:(int *)tag boundary:(int *)boundary leftElement:(int *)leftElement rightElement:(int *)rightElement type:(int *)type nodes:(int *)nodes coord:(double *)coord {
    
    int part;
    if ([self.meshAgent readNextBoundaryElement:tag part:&part boundary:boundary leftElement:leftElement rightElement:rightElement type:type nodes:nodes coord:coord] != -1) {
        self.info = 0;
    } else self.info = -1;
}

-(void)getPartDesription:(int *)sharedNodeCount {
    
    self.info = -1;
    self.info = [self.meshAgent readPartDescriptor:sharedNodeCount];
}

-(void)getPartNode:(int *)tag constraint:(int *)constraint coord:(double *)coord partCount:(int *)partCount parts:(int *)parts {
    
    if ([self.meshAgent readSharedNode:tag constraint:constraint coord:coord partCount:partCount partitions:parts] != -1) {
        self.info = 0;
    } else self.info = -1;
}

@end
