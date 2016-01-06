//
//  FEMMeshUtils.h
//  Saino
//
//  Created by Seddik hakime on 28/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMesh.h"
#import "FEMSolution.h"

@interface FEMMeshUtils : NSObject

-(void)findEdges2DInMesh:(FEMMesh * __nonnull)mesh;
-(void)findEdges3DInMesh:(FEMMesh * __nonnull)mesh;
-(void)findFaces3DInMesh:(FEMMesh * __nonnull)mesh;
-(void)findEdgesForMesh:(FEMMesh * __nonnull)mesh findEdges:(BOOL * __nullable)present;
-(void)assignLocalNumberToEdgeElement:(Element_t * __nonnull)edge fromElement:(Element_t * __nonnull)element inMesh:(FEMMesh * __nonnull)mesh;
-(FEMMatrix * __nullable)periodicProjectorInModel:(FEMModel * __nonnull)model forMesh:(FEMMesh * __nonnull)mesh masterBoundary:(int)mbd targetBoundary:(int)trgt galerking:(BOOL * __nullable)galerkin;
-(FEMMesh * __nullable)splitMeshEqual:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model nodal:(double * __nullable)h sizeNodal:(int * __nullable)sizeNodal;
-(void)setStabilizationParametersInMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model;
-(void)setCurrentMesh:(FEMMesh * __nonnull)mesh inModel:(FEMModel * __nonnull)model;
-(void)updateMesh:(FEMMesh * __nonnull)mesh inSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model;
-(void)allocatePDefinitionsElement:(Element_t * __nonnull)element;
-(void)setEdgeFaceDofsMesh:(FEMMesh * __nonnull)mesh edgeDofs:(int * __nullable)edgeDofs faceDofs:(int * __nullable)faceDofs;
-(void)setMaximumDofsMesh:(FEMMesh * __nonnull)mesh;
-(FEMMesh * __nonnull)extrudeMesh:(FEMMesh * __nonnull)mesh inLevels:(int)inLevels model:(FEMModel * __nonnull)model;
-(int)detectExtrudedStructureMesh:(FEMMesh * __nonnull)mesh solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model externVariable:(FEMVariable * __nullable)externVariable needExternVariable:(BOOL)needExternVariable isTopActive:(BOOL)isTopActive topNodePointer:(int * __nullable)topNodePointer isBottomActive:(BOOL)isBottomActive bottomNodePointer:(int * __nullable)bottomNodePointer isUpActive:(BOOL)isUpActive upNodePointer:(int * __nullable)upNodePointer isDownActive:(BOOL)isDownActive downNodePointer:(int * __nullable)downNodePointer numberOfLayers:(int * __nullable)numberOfLayers nodeLayer:(int * __nullable)nodeLayer;
-(void)preRotationalProjectorMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 mirrorNode:(BOOL * __nullable)mirrorNode sizeMirrorNode:(int * __nullable)sizeMirrorNode;
-(void)postRotationalProjector:(FEMMatrix * __nonnull)projector mirrorNode:(BOOL * __nonnull)mirrorNode sizeMirrorNode:(int * __nonnull)sizeMirrorNode;
-(void)saveProjector:(FEMMatrix * __nonnull)projector saveRowSum:(BOOL)saveRowSum prefix:(NSString * __nonnull)prefix invPerm:(int * __nullable)invPerm;
@end
