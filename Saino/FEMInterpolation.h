//
//  FEMInterpolation.h
//  Saino
//
//  Created by Seddik hakime on 27/09/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"
#import "FEMMesh.h"

@interface FEMInterpolation : NSObject

-(void)findLeafElementsForPoint:(double *)point dimension:(int)dim rootQuadrant:(Quadrant_t *)root leafQuadrant:(Quadrant_t *)leaf;
-(BOOL)isPointInElement:(Element_t *)element elementNodes:(Nodes_t *)nodes point:(double *)aPoint localCoordinates:(double *)localCoords globalEpsilon:(double *)globaleps localEpsilon:(double *)localeps numericEpsilon:(double *)numericeps globalDistance:(double *)globaldist localDistance:(double *)localdist model:(FEMModel *)aModel;
-(void)putElementsInChildQuadrants:(QuadrantPointer_t *)childQuadrants motherQuadrant:(Quadrant_t *)mother mesh:(FEMMesh *)aMesh dimension:(int)dim;
-(void)createChildQuadrantFromMother:(Quadrant_t *)quadrant mesh:(FEMMesh *)aMesh dimension:(int)dim maxLeafElements:(int)maxLeafElems generation:(int *)gen;
-(void)buildQuadrantTreeForMesh:(FEMMesh *)mesh model:(FEMModel *)aModel boundingBox:(double *)box rootQuadrant:(Quadrant_t *)quadrant;

@end
