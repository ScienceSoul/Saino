//===----------------------------------------------------------------------===//
//  SIOMeshAgent.h
//  Saino
//
//  Created by Hakime Seddik on 21/06/12.
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

#import <string.h>
#import <stdlib.h>
#import <stdio.h>
#import "FileReader.h"
#import "SIOModelManager.h"
#import "NodeCompare.h"

@interface SIOMeshAgent : NSObject {
@private
    
    SIOModelManager *_manager;
    
    // All streams. i.e., array of FileReader or NSFileHandle classes
    NSMutableArray *_meshFileStreams;
    
    int _parts;
    int _me;
    int _nodeCount;
    int _elementCount;
    int _boundaryElementCount;
    int _elementTypes;
    NSMutableArray *_elementTypeTags;    // Container for integer values
    NSMutableArray *_elementTypeCount;   // Container for integer values
    int _sharedNodeCount;
    int _borderElementCount;
    int _dim;
    int _parallel;
    int _meshFiles;
    
}

@property(nonatomic, strong, nullable) SIOModelManager *manager;
@property(nonatomic, strong, nonnull) NSMutableArray *meshFileStreams;
@property(nonatomic, strong, nonnull) NSMutableArray <NSString *> *elementTypeTags;
@property(nonatomic, strong, nonnull) NSMutableArray <NSString *> *elementTypeCount;
@property(nonatomic, assign) int parts;
@property(nonatomic, assign) int me;
@property(nonatomic, assign) int nodeCount;
@property(nonatomic, assign) int elementCount;
@property(nonatomic, assign) int boundaryElementCount;
@property(nonatomic, assign) int elementTypes;
@property(nonatomic, assign) int sharedNodeCount;
@property(nonatomic, assign) int borderElementCount;
@property(nonatomic, assign) int dim;
@property(nonatomic, assign) int parallel;
@property(nonatomic, assign) int meshFiles;

-(id _Nullable)initWithManager:(SIOModelManager * _Nonnull)mm split:(int)split part:(int)part;
-(int)createMesh:(NSString * _Nonnull)dir;
-(int)openMesh:(NSString * _Nonnull)dir;
-(int)closeMesh;

// Reading methods
-(int)readDescriptorNode:(int * _Nonnull)nodeC element:(int * _Nonnull)elementC boundaryElement:(int * _Nonnull)boundaryElementC usedElementTypes:(int * _Nonnull)usedElementTypes usedElementTypeTags:(int * _Nonnull)usedElementTypeTags usedElementTypeCount:(int * _Nonnull)usedElementTypeCount;
-(int)readNextElementConnections:(int * _Nonnull)tag part:(int * _Nonnull)part body:(int * _Nonnull)body type:(int * _Nonnull)type pdofs:(int * _Nonnull)pdofs nodes:(int * _Nonnull)nodes colorIndex:(int * _Nullable)colorIndex parallelAssembly:(BOOL * _Nullable)parallelAssembly;
-(int)readNextElementCoordinates:(int * _Nonnull)tag body:(int * _Nonnull)body type:(int * _Nonnull)type nodes:(int * _Nonnull)nodes coord:(double * _Nonnull)coord;
-(int)readNextBoundaryElement:(int * _Nonnull)tag part:(int * _Nonnull)part boundary:(int * _Nonnull)boundary leftElement:(int * _Nonnull)leftElement rightElement:(int * _Nonnull)rightElement type:(int * _Nonnull)type nodes:(int * _Nonnull)nodes coord:(double * _Nonnull)coord;
-(int)readAllNodes:(int * _Nonnull)tags coord:(double * _Nonnull)coord;

// Writing methods
//-(int)writeDescriptor:(int)nodeC: (int)elementC: (int)boundaryElementC: (int)usedElementTypes: (int *)elementTypeTags: (int *)elementCountByType;
//-(int)writeNode:(int)tag: (int)type: (double *)coord;
//-(int)writeElementConnections:(int)tag: (int)body: (int)type: (int *)nodes;
//-(int)writeBoundaryElement:(int)tag: (int)boundary: (int)leftElement: (int)rightElement: (int)type: (int *)nodes;
//-(int)writeBoundaryElement:(int)tag: (int)boundary: (int)leftElement: (int)rightElement: (int)type: (int *)nodes;

-(int)readPartDescriptor:(int * _Nonnull)shared;
-(int)readSharedNode:(int * _Nonnull)tag constraint:(int * _Nonnull)constraint coord:(double * _Nonnull)coord partCount:(int * _Nonnull)partcount partitions:(int * _Nonnull)partitions;

@end
