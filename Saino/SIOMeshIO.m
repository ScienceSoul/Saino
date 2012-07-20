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
    
    _meshAgent = [[SIOMeshAgent alloc] initWithManager:self.modelManager :self.parallelState.numProc :self.parallelState.myProc];
    
    if (self.meshAgent != nil) {
        self.info = [self.meshAgent openMesh:directory];
    } else {
        self.info = -1;
    }
}

-(void)closeMesh {
    
    [self.meshAgent closeMesh];
    self.info = 0;
    self.meshAgent = nil;
}

-(void)close {
    
    self.modelManager = nil;
    self.info = 0;
}

-(void)getMeshDescription:(int *)nodeCount :(int *)elementCount :(int *)boundaryElementCount :(int *)usedElementTypes :(int *)elementTypeTags :(int *)elementCountByType {
    
    [self.meshAgent readDescriptor:nodeCount :elementCount :boundaryElementCount :usedElementTypes :elementTypeTags :elementCountByType];
    self.info = 0;
}

-(void)getMeshNodes:(int *)tags: (double *)coord {
    
    [self.meshAgent readAllNodes:tags :coord];
    self.info = 0;
}

-(void)getMeshElementConnection:(int *)tag :(int *)body :(int *)type :(int *)pdofs :(int *)nodes {
    
    int part;
    if ([self.meshAgent readNextElementConnections:tag :&part :body :type :pdofs :nodes] != -1) {
        self.info = 0;
    } else self.info = -1;
}

-(void)getMeshBoundaryElement:(int *)tag :(int *)boundary :(int *)leftElement :(int *)rightElement :(int *)type :(int *)nodes :(double *)coord {
    
    int part;
    if ([self.meshAgent readNextBoundaryElement:tag :&part :boundary :leftElement :rightElement :type :nodes :coord] != -1) {
        self.info = 0;
    } else self.info = -1;
}

-(void)getPartDesription:(int *)sharedNodeCount {
    
    [self.meshAgent readPartDescriptor:sharedNodeCount];
    self.info = 0;
}

-(void)getPartNode:(int *)tag :(int *)constraint :(double *)coord :(int *)partCount :(int *)parts {
    
    if ([self.meshAgent readSharedNode:tag :constraint :coord :partCount :parts] != -1) {
        self.info = 0;
    } else self.info = -1;
}

@end
