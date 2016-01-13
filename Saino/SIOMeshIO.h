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
@property(nonatomic, strong, nullable) SIOInfoParallel *parallelState;
@property(nonatomic, strong, nullable) SIOModelManager *modelManager;
@property(nonatomic, strong, nullable) SIOMeshAgent *meshAgent;

-(id __nonnull)initWithParallelNumberOfProcessors:(int)procs processorID:(int)me;
-(void)openMeshAtPath:(NSString * __nonnull)directory;
-(void)closeMesh;
-(void)close;
-(void)getMeshDescriptionNodeCount:(int * __nonnull)nodeCount elementCount:(int * __nonnull)elementCount boundaryElementCount:(int * __nonnull)boundaryElementCount usedElementTypes:(int * __nonnull)usedElementTypes elementTypeTags:(int * __nonnull)elementTypeTags elementCountByType:(int * __nonnull)elementCountByType;
-(void)getMeshNodes:(int * __nonnull)tags coord:(double * __nonnull)coord;
-(void)getMeshElementConnection:(int * __nonnull)tag body:(int * __nonnull)body type:(int * __nonnull)type pdofs:(int * __nonnull)pdofs nodes:(int * __nonnull)nodes colorIndex:(int * __nullable)colorIndex parallelAssembly:(BOOL * __nullable)parallelAssembly;
-(void)getMeshBoundaryElement:(int * __nonnull)tag boundary:(int * __nonnull)boundary leftElement:(int * __nonnull)leftElement rightElement:(int * __nonnull)rightElement type:(int * __nonnull)type nodes:(int * __nonnull)nodes coord:(double * __nonnull)coord;
-(void)getPartDesription:(int * __nonnull)sharedNodeCount;
-(void)getPartNode:(int * __nonnull)tag constraint:(int * __nonnull)constraint coord:(double * __nonnull)coord partCount:(int * __nonnull)partCount parts:(int * __nonnull)parts;

@end
