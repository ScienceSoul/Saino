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
-(void)getMeshDescriptionNodeCount:(int *)nodeCount elementCount:(int *)elementCount boundaryElementCount:(int *)boundaryElementCount usedElementTypes:(int *)usedElementTypes elementTypeTags:(int *)elementTypeTags elementCountByType:(int *)elementCountByType;
-(void)getMeshNodes:(int *)tags coord:(double *)coord;
-(void)getMeshElementConnection:(int *)tag body:(int *)body type:(int *)type pdofs:(int *)pdofs nodes:(int *)nodes colorIndex:(int *)colorIndex parallelAssembly:(BOOL *)parallelAssembly;
-(void)getMeshBoundaryElement:(int *)tag boundary:(int *)boundary leftElement:(int *)leftElement rightElement:(int *)rightElement type:(int *)type nodes:(int *)nodes coord:(double *)coord;
-(void)getPartDesription:(int *)sharedNodeCount;
-(void)getPartNode:(int *)tag constraint:(int *)constraint coord:(double *)coord partCount:(int *)partCount parts:(int *)parts;

@end
