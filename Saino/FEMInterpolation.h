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

-(Quadrant_t * __nullable)findLeafElementsForPoint:(double * __nonnull)point dimension:(int)dim rootQuadrant:(Quadrant_t * __nonnull)root;
-(BOOL)isPointInElement:(Element_t * __nonnull)element elementNodes:(Nodes_t * __nonnull)nodes point:(double * __nonnull)point localCoordinates:(double * __nonnull)localCoords globalEpsilon:(double * __nullable)globaleps localEpsilon:(double * __nullable)localeps numericEpsilon:(double * __nullable)numericeps globalDistance:(double * __nullable)globaldist localDistance:(double * __nullable)localdist model:(FEMModel * __nonnull)model elementDescription:(FEMElementDescription * __nonnull)elementDescription elementMaps:(FEMPElementMaps * __nonnull)elementMaps edgeBasis:(BOOL * __nullable)edgeBasis;
-(void)putElementsInChildQuadrants:(QuadrantPointer_t * __nonnull)childQuadrants motherQuadrant:(Quadrant_t * __nonnull)mother mesh:(FEMMesh * __nonnull)mesh dimension:(int)dim;
-(void)createChildQuadrantFromMother:(Quadrant_t * __nonnull)quadrant mesh:(FEMMesh * __nonnull)mesh dimension:(int)dim maxLeafElements:(int)maxLeafElems generation:(int * __nonnull)gen;
-(void)buildQuadrantTreeForMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model boundingBox:(double * __nonnull)box rootQuadrant:(Quadrant_t * __nullable)quadrant;

@end
