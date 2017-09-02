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
-(int)bubbleDofsForElement:(Element_t * _Nonnull)element degree:(int)p;
-(void)getTriangleEdgeMap:(int * _Nonnull)edge index:(int)i;
-(void)getQuadEdgeMap:(int * _Nonnull)edge index:(int)i;
-(void)getBrickEdgeMap:(int * _Nonnull)edge index:(int)i;
-(void)getBrickFaceMap:(int * _Nonnull)face index:(int)i;
-(int)getBrickFaceEdgeMap:(int)face localNode:(int)node;
-(void)getTetraEdgeMap:(int * _Nonnull)edge index:(int)i type:(int * _Nullable)type;
-(void)getTetraFaceMap:(int * _Nonnull)face index:(int)i type:(int * _Nullable)type;
-(void)getWedgeEdgeMap:(int * _Nonnull)edge index:(int)i;
-(void)getWedgeFaceMap:(int * _Nonnull)face index:(int)i;
-(void)getPyramidEdgeMap:(int * _Nonnull)edge index:(int)i;
-(void)getPyramidFaceMap:(int * _Nonnull)face index:(int)i;
-(BOOL)isPElement:(Element_t * _Nonnull)element;
-(void)getRefPElementNodesForElement:(Element_t * _Nonnull)element nodeU:(double * _Nonnull)u nodeV:(double * _Nonnull)v nodeW:(double * _Nonnull)w;
-(void)getFaceMapForElement:(Element_t * _Nonnull)element faceMap:(int * _Nullable * _Nullable)map;
-(void)getEdgeMapForElement:(Element_t * _Nonnull)element edgeMap:(int * _Nullable * _Nullable)map;
-(void)getFaceEdgeMapForElement:(Element_t * _Nonnull)element faceEdgeMap:(int * _Nullable * _Nullable)map;
-(BOOL)isPPyramid:(Element_t * _Nonnull)element;
-(void)getBoundaryMapForElement:(Element_t * _Nonnull)element localNumber:(int)i resultMap:(int * _Nonnull)map;
-(void)getFaceEdgeMapForElement:(Element_t * _Nonnull)element index:(int)i resultMap:(int * _Nonnull)map;
-(int)getFaceDofsForElement:(Element_t * _Nonnull)element polyDegree:(int)p faceNumber:(int)number;
-(int)getNumberOfGaussPointsForFace:(Element_t * _Nonnull)face inMesh:(FEMMesh * _Nonnull)mesh;
-(int)getEdgePForElement:(Element_t * _Nonnull)element inMesh:(FEMMesh * _Nonnull)mesh;
-(int)getFacePForElement:(Element_t * _Nonnull)element inMesh:(FEMMesh * _Nonnull)mesh;
-(int)getNumberOfGaussPointsForElement:(Element_t * _Nonnull)element inMesh:(FEMMesh * _Nonnull)mesh;
-(int)getBubbleDofsForElement:(Element_t * _Nonnull)element polyDegree:(int)p;
-(int)getEdgeDofsForElement:(Element_t * _Nonnull)element polyDegree:(int)p;
@end
