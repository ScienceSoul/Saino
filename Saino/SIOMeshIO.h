//===----------------------------------------------------------------------===//
//  SIOMeshIO.h
//  Saino
//
//  Created by Hakime Seddik on 28/06/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
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

-(id _Nonnull)initWithParallelNumberOfProcessors:(int)procs processorID:(int)me;
-(void)openMeshAtPath:(NSString * _Nonnull)directory;
-(void)closeMesh;
-(void)close;
-(void)getMeshDescriptionNodeCount:(int * _Nonnull)nodeCount elementCount:(int * _Nonnull)elementCount boundaryElementCount:(int * _Nonnull)boundaryElementCount usedElementTypes:(int * _Nonnull)usedElementTypes elementTypeTags:(int * _Nonnull)elementTypeTags elementCountByType:(int * _Nonnull)elementCountByType;
-(void)getMeshNodes:(int * _Nonnull)tags coord:(double * _Nonnull)coord;
-(void)getMeshElementConnection:(int * _Nonnull)tag body:(int * _Nonnull)body type:(int * _Nonnull)type pdofs:(int * _Nonnull)pdofs nodes:(int * _Nonnull)nodes colorIndex:(int * _Nullable)colorIndex parallelAssembly:(BOOL * _Nullable)parallelAssembly;
-(void)getMeshBoundaryElement:(int * _Nonnull)tag boundary:(int * _Nonnull)boundary leftElement:(int * _Nonnull)leftElement rightElement:(int * _Nonnull)rightElement type:(int * _Nonnull)type nodes:(int * _Nonnull)nodes coord:(double * _Nonnull)coord;
-(void)getPartDesription:(int * _Nonnull)sharedNodeCount;
-(void)getPartNode:(int * _Nonnull)tag constraint:(int * _Nonnull)constraint coord:(double * _Nonnull)coord partCount:(int * _Nonnull)partCount parts:(int * _Nonnull)parts;

@end
