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

@property(nonatomic, strong) SIOModelManager *manager;
@property(nonatomic, strong) NSMutableArray *meshFileStreams;
@property(nonatomic, strong) NSMutableArray *elementTypeTags;
@property(nonatomic, strong) NSMutableArray *elementTypeCount;
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

- (id)initWithManager:(SIOModelManager *)mm: (int)split: (int)part;
-(int)createMesh:(NSString *)dir;
-(int)openMesh:(NSString *)dir;
-(int)closeMesh;

// Reading methods
-(int)readDescriptor:(int *)nodeC: (int *)elementC: (int *)boundaryElementC: (int *)usedElementTypes: (int *)usedElementTypeTags: (int *)usedElementTypeCount;
-(int)readNextElementConnections:(int *)tag: (int *)part: (int *)body: (int *)type: (int *)pdofs: (int *)nodes;
-(int)readNextElementCoordinates:(int *)tag: (int *)body: (int *)type: (int *)nodes: (double *)coord;
-(int)readNextBoundaryElement:(int *)tag: (int *)part: (int *)boundary: (int *)leftElement: (int *)rightElement: (int *)type: (int *)nodes: (double *)coord;
-(int)readAllNodes:(int *)tags: (double *)coord;

// Writing methods
//-(int)writeDescriptor:(int)nodeC: (int)elementC: (int)boundaryElementC: (int)usedElementTypes: (int *)elementTypeTags: (int *)elementCountByType;
//-(int)writeNode:(int)tag: (int)type: (double *)coord;
//-(int)writeElementConnections:(int)tag: (int)body: (int)type: (int *)nodes;
//-(int)writeBoundaryElement:(int)tag: (int)boundary: (int)leftElement: (int)rightElement: (int)type: (int *)nodes;
//-(int)writeBoundaryElement:(int)tag: (int)boundary: (int)leftElement: (int)rightElement: (int)type: (int *)nodes;

-(int)readPartDescriptor:(int *)shared;
-(int)readSharedNode:(int *)tag: (int *)constraint: (double *)coord: (int *)partcount: (int *)partitions;

@end
