//
//  FEMMeshUtils.h
//  Saino
//
//  Created by Seddik hakime on 28/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMesh.h"
#import "FEMVariable.h"
#import "FEMProjector.h"
#import "FEMElementDescription.h"
#import "FEMPElementMaps.h"
#import "FEMBoundaryCondition.h"
#import "FEMListUtilities.h"

@interface FEMMeshUtils : NSObject

-(void)findEdges2DInMesh:(FEMMesh *)mesh;
-(void)findEdges3DInMesh:(FEMMesh *)mesh;
-(void)findFaces3DInMesh:(FEMMesh *)mesh;
-(void)findEdgesForMesh:(FEMMesh *)mesh findEdges:(BOOL *)edges;
-(void)assignLocalNumberToEdgeElement:(Element_t *)edge fromElement:(Element_t *)element inMesh:(FEMMesh *)mesh;
-(FEMMatrix *)periodicProjectorInModel:(FEMModel *)model forMesh:(FEMMesh *)mesh boundary:(int)this target:(int)trgt;

@end
