//===----------------------------------------------------------------------===//
//  FEMInterpolation.h
//  Saino
//
//  Created by Seddik hakime on 27/09/12.
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
#import "FEMModel.h"
#import "FEMMesh.h"
#import "FEMElementDescription.h"
#import "FEMPElementMaps.h"

@interface FEMInterpolation : NSObject

-(Quadrant_t * _Nullable)findLeafElementsForPoint:(double * _Nonnull)point dimension:(int)dim rootQuadrant:(Quadrant_t * _Nonnull)root;
-(BOOL)isPointInElement:(Element_t * _Nonnull)element elementNodes:(Nodes_t * _Nonnull)nodes point:(double * _Nonnull)point localCoordinates:(double * _Nonnull)localCoords globalEpsilon:(double * _Nullable)globaleps localEpsilon:(double * _Nullable)localeps numericEpsilon:(double * _Nullable)numericeps globalDistance:(double * _Nullable)globaldist localDistance:(double * _Nullable)localdist model:(FEMModel * _Nonnull)model elementDescription:(FEMElementDescription * _Nonnull)elementDescription elementMaps:(FEMPElementMaps * _Nonnull)elementMaps edgeBasis:(BOOL * _Nullable)edgeBasis;
-(void)putElementsInChildQuadrants:(QuadrantPointer_t * _Nonnull)childQuadrants motherQuadrant:(Quadrant_t * _Nonnull)mother mesh:(FEMMesh * _Nonnull)mesh dimension:(int)dim;
-(void)createChildQuadrantFromMother:(Quadrant_t * _Nonnull)quadrant mesh:(FEMMesh * _Nonnull)mesh dimension:(int)dim maxLeafElements:(int)maxLeafElems generation:(int * _Nonnull)gen;
-(void)buildQuadrantTreeForMesh:(FEMMesh * _Nonnull)mesh model:(FEMModel * _Nonnull)model boundingBox:(double * _Nonnull)box rootQuadrant:(Quadrant_t * _Nullable)quadrant;

@end
