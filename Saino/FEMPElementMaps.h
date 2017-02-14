//===----------------------------------------------------------------------===//
//  FEMPElementMaps.h
//  Saino
//
//  Created by Seddik hakime on 15/08/12.
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
#import "FEMMesh.h"

@interface FEMPElementMaps : NSObject

-(void)deallocation;
-(int)bubbleDofsForElement:(Element_t * __nonnull)element degree:(int)p;
-(void)getTriangleEdgeMap:(int * __nonnull)edge index:(int)i;
-(void)getQuadEdgeMap:(int * __nonnull)edge index:(int)i;
-(void)getBrickEdgeMap:(int * __nonnull)edge index:(int)i;
-(void)getBrickFaceMap:(int * __nonnull)face index:(int)i;
-(int)getBrickFaceEdgeMap:(int)face localNode:(int)node;
-(void)getTetraEdgeMap:(int * __nonnull)edge index:(int)i type:(int * __nullable)type;
-(void)getTetraFaceMap:(int * __nonnull)face index:(int)i type:(int * __nullable)type;
-(void)getWedgeEdgeMap:(int * __nonnull)edge index:(int)i;
-(void)getWedgeFaceMap:(int * __nonnull)face index:(int)i;
-(void)getPyramidEdgeMap:(int * __nonnull)edge index:(int)i;
-(void)getPyramidFaceMap:(int * __nonnull)face index:(int)i;
-(BOOL)isPElement:(Element_t * __nonnull)element;
-(void)getRefPElementNodesForElement:(Element_t * __nonnull)element nodeU:(double * __nonnull)u nodeV:(double * __nonnull)v nodeW:(double * __nonnull)w;
-(void)getFaceMapForElement:(Element_t * __nonnull)element faceMap:(int * __nullable * __nullable)map;
-(void)getEdgeMapForElement:(Element_t * __nonnull)element edgeMap:(int * __nullable * __nullable)map;
-(void)getFaceEdgeMapForElement:(Element_t * __nonnull)element faceEdgeMap:(int * __nullable * __nullable)map;
-(BOOL)isPPyramid:(Element_t * __nonnull)element;
-(void)getBoundaryMapForElement:(Element_t * __nonnull)element localNumber:(int)i resultMap:(int * __nonnull)map;
-(void)getFaceEdgeMapForElement:(Element_t * __nonnull)element index:(int)i resultMap:(int * __nonnull)map;
-(int)getFaceDofsForElement:(Element_t * __nonnull)element polyDegree:(int)p faceNumber:(int)number;
-(int)getNumberOfGaussPointsForFace:(Element_t * __nonnull)face inMesh:(FEMMesh * __nonnull)mesh;
-(int)getEdgePForElement:(Element_t * __nonnull)element inMesh:(FEMMesh * __nonnull)mesh;
-(int)getFacePForElement:(Element_t * __nonnull)element inMesh:(FEMMesh * __nonnull)mesh;
-(int)getNumberOfGaussPointsForElement:(Element_t * __nonnull)element inMesh:(FEMMesh * __nonnull)mesh;
-(int)getBubbleDofsForElement:(Element_t * __nonnull)element polyDegree:(int)p;
-(int)getEdgeDofsForElement:(Element_t * __nonnull)element polyDegree:(int)p;
@end
