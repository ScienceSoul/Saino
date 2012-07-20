//
//  SIOMeshIO.h
//  Saino
//
//  Created by Hakime Seddik on 28/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "SIOModelManager.h"
#import "SIOMeshAgent.h"
#import "SIOInfoParallel.h"

@interface SIOMeshIO : NSObject {
@private
    
    int _info;
    SIOInfoParallel *_parallelState;
    SIOModelManager *_modelManager;
    SIOMeshAgent *_meshAgent;
}

@property(nonatomic, assign) int info;
@property(nonatomic, strong) SIOInfoParallel *parallelState;
@property(nonatomic, strong) SIOModelManager *modelManager;
@property(nonatomic, strong) SIOMeshAgent *meshAgent;

-(id)initWithParallelNumberOfProcessors:(int)procs processorID:(int)me;
-(void)openMeshAtPath:(NSString *)directory;
-(void)closeMesh;
-(void)close;
-(void)getMeshDescription:(int *)nodeCount: (int *)elementCount: (int *)boundaryElementCount: (int *)usedElementTypes: (int *)elementTypeTags: (int *)elementCountByType;
-(void)getMeshNodes:(int *)tags: (double *)coord;
-(void)getMeshElementConnection:(int *)tag: (int *)body: (int *)type: (int *)pdofs: (int *)nodes;
-(void)getMeshBoundaryElement:(int *)tag: (int *)boundary: (int *)leftElement: (int *)rightElement: (int *)type: (int *)nodes: (double *)coord;
-(void)getPartDesription:(int *)sharedNodeCount;
-(void)getPartNode:(int *)tag: (int *)constraint: (double *)coord: (int *)partCount: (int *)parts;

@end
