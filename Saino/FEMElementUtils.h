//
//  FEMElementUtils.h
//  Saino
//
//  Created by Seddik hakime on 28/12/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMesh.h"
#import "FEMModel.h"
#import "FEMSolution.h"
#import "FEMMatrix.h"

@interface FEMElementUtils : NSObject

-(FEMMatrix * __nonnull)createMatrixInModel:(FEMModel * __nonnull)model forSolution:(FEMSolution * __nonnull)solution mesh:(FEMMesh * __nonnull)mesh dofs:(int)dofs permutation:(int * __nonnull)perm sizeOfPermutation:(int)permSize matrixFormat:(int)matrixFormat optimizeBandwidth:(BOOL)optimizeBandwidth equationName:(NSString * __nullable)equation discontinuousGalerkinSolution:(BOOL * __nullable)dgSolution globalBubbles:(BOOL * __nullable)gbBubbles nodalDofsOnly:(BOOL * __nullable)nodalDofsOnly projectorDofs:(BOOL * __nullable)projectorDofs;
-(void)tangentDirectionsForNormal:(double * __nonnull)normal tangent1:(double * __nonnull)tangent1 tangent2:(double * __nonnull)tangent2;
-(double)elementArea:(Element_t * __nonnull)element numberOfNodes:(int)n mesh:(FEMMesh * __nonnull)mesh nodel:(FEMModel * __nonnull)model;

@end
