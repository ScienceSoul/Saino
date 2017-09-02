//===----------------------------------------------------------------------===//
//  FEMMeshUtils.h
//  Saino
//
//  Created by Seddik hakime on 28/08/12.
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
#import "FEMSolution.h"

@interface FEMMeshUtils : NSObject

-(void)findEdges2DInMesh:(FEMMesh * _Nonnull)mesh;
-(void)findEdges3DInMesh:(FEMMesh * _Nonnull)mesh;
-(void)findFaces3DInMesh:(FEMMesh * _Nonnull)mesh;
-(void)findEdgesForMesh:(FEMMesh * _Nonnull)mesh findEdges:(BOOL * _Nullable)present;
-(void)assignLocalNumberToEdgeElement:(Element_t * _Nonnull)edge fromElement:(Element_t * _Nonnull)element inMesh:(FEMMesh * _Nonnull)mesh;
-(FEMMatrix * _Nullable)periodicProjectorInModel:(FEMModel * _Nonnull)model forMesh:(FEMMesh * _Nonnull)mesh masterBoundary:(int)mbd targetBoundary:(int)trgt galerking:(BOOL * _Nullable)galerkin;
-(FEMMesh * _Nullable)splitMeshEqual:(FEMMesh * _Nonnull)mesh model:(FEMModel * _Nonnull)model nodal:(double * _Nullable)h sizeNodal:(int * _Nullable)sizeNodal;
-(void)setStabilizationParametersInMesh:(FEMMesh * _Nonnull)mesh model:(FEMModel * _Nonnull)model;
-(void)setCurrentMesh:(FEMMesh * _Nonnull)mesh inModel:(FEMModel * _Nonnull)model;
-(void)updateMesh:(FEMMesh * _Nonnull)mesh inSolution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model;
-(void)allocatePDefinitionsElement:(Element_t * _Nonnull)element;
-(void)setEdgeFaceDofsMesh:(FEMMesh * _Nonnull)mesh edgeDofs:(int * _Nullable)edgeDofs faceDofs:(int * _Nullable)faceDofs;
-(void)setMaximumDofsMesh:(FEMMesh * _Nonnull)mesh;
-(FEMMesh * _Nonnull)extrudeMesh:(FEMMesh * _Nonnull)mesh inLevels:(int)inLevels model:(FEMModel * _Nonnull)model;
-(FEMVariable * _Nullable)detectExtrudedStructureMesh:(FEMMesh * _Nonnull)mesh solution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model numberofNodes:(int)numberofNodes ifMask:(BOOL)ifMask mask:(variableArraysContainer * _Nullable)mask isTopActive:(BOOL)isTopActive topNodePointer:(int * _Nullable)topNodePointer isBottomActive:(BOOL)isBottomActive bottomNodePointer:(int * _Nullable)bottomNodePointer isUpActive:(BOOL)isUpActive upNodePointer:(int * _Nullable)upNodePointer isDownActive:(BOOL)isDownActive downNodePointer:(int * _Nullable)downNodePointer numberOfLayers:(int * _Nullable)numberOfLayers isMidNode:(BOOL)isMidNode midNodePointer:(int * _Nullable)midNodePointer midLayerExists:(BOOL * _Nullable)midLayerExists isNodeLayer:(BOOL)isNodeLayer nodeLayer:(int * _Nullable)nodeLayer;
-(void)preRotationalProjectorMesh1:(FEMMesh * _Nonnull)bMesh1 mesh2:(FEMMesh * _Nonnull)bMesh2 mirrorNode:(BOOL * _Nullable)mirrorNode sizeMirrorNode:(int * _Nullable)sizeMirrorNode;
-(void)postRotationalProjector:(FEMMatrix * _Nonnull)projector mirrorNode:(BOOL * _Nonnull)mirrorNode sizeMirrorNode:(int * _Nonnull)sizeMirrorNode;
-(void)saveProjector:(FEMMatrix * _Nonnull)projector saveRowSum:(BOOL)saveRowSum prefix:(NSString * _Nonnull)prefix invPerm:(int * _Nullable)invPerm;
-(void)colorMesh:(FEMMesh * _Nonnull)mesh;
-(void)saveColoredMesh:(FEMMesh * _Nonnull)mesh meshdir:(NSString * _Nonnull)dir meshName:(NSString * _Nonnull)name elementsFileName:(NSString * _Nonnull)elementsFileName saveAllElementData:(BOOL)saveAllElementData colorFileName:(NSString * _Nonnull)colorFileName;
-(void)readColoredMesh:(FEMMesh * _Nonnull)mesh name:(NSString * _Nonnull)name directory:(NSString * _Nonnull)dir readElementsFromFile:(BOOL)readElementsFromFile;
@end
