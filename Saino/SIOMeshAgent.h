//
//  SIOMeshAgent.h
//  Saino
//
//  Created by Hakime Seddik on 21/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

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
    
    // All streams. i.e., array of FileReader classes
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
@property(nonatomic, strong, nonnull) NSMutableArray *elementTypeTags;
@property(nonatomic, strong, nonnull) NSMutableArray *elementTypeCount;
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

-(id __nullable)initWithManager:(SIOModelManager * __nonnull)mm split:(int)split part:(int)part;
-(int)createMesh:(NSString * __nonnull)dir;
-(int)openMesh:(NSString * __nonnull)dir;
-(int)closeMesh;

// Reading methods
-(int)readDescriptorNode:(int * __nonnull)nodeC element:(int * __nonnull)elementC boundaryElement:(int * __nonnull)boundaryElementC usedElementTypes:(int * __nonnull)usedElementTypes usedElementTypeTags:(int * __nonnull)usedElementTypeTags usedElementTypeCount:(int * __nonnull)usedElementTypeCount;
-(int)readNextElementConnections:(int * __nonnull)tag part:(int * __nonnull)part body:(int * __nonnull)body type:(int * __nonnull)type pdofs:(int * __nonnull)pdofs nodes:(int * __nonnull)nodes colorIndex:(int * __nullable)colorIndex parallelAssembly:(BOOL * __nullable)parallelAssembly;
-(int)readNextElementCoordinates:(int * __nonnull)tag body:(int * __nonnull)body type:(int * __nonnull)type nodes:(int * __nonnull)nodes coord:(double * __nonnull)coord;
-(int)readNextBoundaryElement:(int * __nonnull)tag part:(int * __nonnull)part boundary:(int * __nonnull)boundary leftElement:(int * __nonnull)leftElement rightElement:(int * __nonnull)rightElement type:(int * __nonnull)type nodes:(int * __nonnull)nodes coord:(double * __nonnull)coord;
-(int)readAllNodes:(int * __nonnull)tags coord:(double * __nonnull)coord;

// Writing methods
//-(int)writeDescriptor:(int)nodeC: (int)elementC: (int)boundaryElementC: (int)usedElementTypes: (int *)elementTypeTags: (int *)elementCountByType;
//-(int)writeNode:(int)tag: (int)type: (double *)coord;
//-(int)writeElementConnections:(int)tag: (int)body: (int)type: (int *)nodes;
//-(int)writeBoundaryElement:(int)tag: (int)boundary: (int)leftElement: (int)rightElement: (int)type: (int *)nodes;
//-(int)writeBoundaryElement:(int)tag: (int)boundary: (int)leftElement: (int)rightElement: (int)type: (int *)nodes;

-(int)readPartDescriptor:(int * __nonnull)shared;
-(int)readSharedNode:(int * __nonnull)tag constraint:(int * __nonnull)constraint coord:(double * __nonnull)coord partCount:(int * __nonnull)partcount partitions:(int * __nonnull)partitions;

@end
