//
//  FEMPElementMaps.h
//  Saino
//
//  Created by Seddik hakime on 15/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMesh.h"
#import "Constructors.h"
#import "memory.h"
#import "Utils.h"

@interface FEMPElementMaps : NSObject

-(void)deallocation;
-(int)bubbleDofsForElement:(Element_t *)element degree:(int)p;
-(void)getTriangleEdgeMap:(int *)edge index:(int)i;
-(void)getQuadEdgeMap:(int *)edge index:(int)i;
-(void)getBrickEdgeMap:(int *)edge index:(int)i;
-(void)getBrickFaceMap:(int *)face index:(int)i;
-(int)getBrickFaceEdgeMap:(int)face localNode:(int)node;
-(void)getTetraEdgeMap:(int *)edge index:(int)i type:(int *)type;
-(void)getTetraFaceMap:(int *)face index:(int)i type:(int *)type;
-(void)getWedgeEdgeMap:(int *)edge index:(int)i;
-(void)getWedgeFaceMap:(int *)face index:(int)i;
-(void)getPyramidEdgeMap:(int *)edge index:(int)i;
-(void)getPyramidFaceMap:(int *)face index:(int)i;
-(BOOL)isPElement:(Element_t *)element;
-(void)getRefPElementNodesForElement:(Element_t *)element nodeU:(double *)u nodeV:(double *)v nodeW:(double *)w;
-(void)getFaceMapForElement:(Element_t *)element faceMap:(int **)map;
-(void)getEdgeMapForElement:(Element_t *)element edgeMap:(int **)map;
-(void)getFaceEdgeMapForElement:(Element_t *)element faceEdgeMap:(int **)map;
-(BOOL)isPPyramid:(Element_t *)element;
-(void)getBoundaryMapForElement:(Element_t *)element localNumber:(int)i resultMap:(int *)map;
-(void)getFaceEdgeMapForElement:(Element_t *)element index:(int)i resultMap:(int *)map;
-(int)getFaceDofsForElement:(Element_t *)element polyDegree:(int)p faceNumber:(int)number;
-(int)getNumberOfGaussPointsForFace:(Element_t *)face inMesh:(FEMMesh *)mesh;
-(int)getEdgePForElement:(Element_t *)element inMesh:(FEMMesh *)mesh;
-(int)getFacePForElement:(Element_t *)element inMesh:(FEMMesh *)mesh;
-(int)getNumberOfGaussPointsForElement:(Element_t *)element inMesh:(FEMMesh *)mesh;
-(int)getBubbleDofsForElement:(Element_t *)element polyDegree:(int)p;
-(int)getEdgeDofsForElement:(Element_t *)element polyDegree:(int)p;
@end