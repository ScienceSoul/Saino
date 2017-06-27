//===----------------------------------------------------------------------===//
//  FEMElementDescription.h
//  Saino
//
//  Created by Hakime Seddik on 05/04/12.
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
#import <Accelerate/Accelerate.h>
#import "FEMMesh.h"

@interface FEMElementDescription : NSObject {
    
}

+(id __nonnull)sharedElementDescription;
+(void)selfDestruct;

-(void)addDescriptionOfElement:(ElementType_t * __nonnull)element withBasisTerms:(int * __nonnull)terms;
-(void)initElementDescriptions;
-(int * __nonnull * __nonnull)getEdgeMap:(int)elementFamily;
-(double)elementDiameter:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes;
-(void)computeStabilizationParameterInElement:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes mesh:(FEMMesh * __nonnull)mesh numberOfNodes:(int)n mk:(double * __nonnull)mk hk:(double * __nullable)hk;
-(ElementType_t * __nullable)getElementType:(int)code inMesh:(FEMMesh * __nonnull)mesh stabilization:(BOOL * __nullable)computeStab;
-(double)firstDerivative1DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evalutationPoint:(double)u;
-(double)firstDerivativeU2DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v;
-(double)firstDerivativeV2DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v;
-(double)firstDerivativeU3DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(double)firstDerivativeV3DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(double)firstDerivativeW3DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(double)interpolateInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w withBasis:(double * __nullable)basis;
-(void)checkNormalDirectionInBDElement:(Element_t * __nonnull)boundary forNormals:(double * __nonnull)normals mesh:(FEMMesh * __nonnull)mesh x:(double)x y:(double)y z:(double)z turn:(BOOL * __nullable)turn;
-(void)normalVectorForBDElement:(Element_t * __nonnull)boundary boundaryNodes:(Nodes_t * __nonnull)nodes mesh:(FEMMesh * __nonnull)mesh paraU:(double * __nullable)u0 paraV:(double * __nullable)v0 check:(BOOL * __nullable)check normals:(double * __nonnull)normals;
-(void)globalToLocalFromElement:(Element_t * __nonnull)element elementNodes:(Nodes_t * __nonnull)nodes localU:(double * __nonnull)u localV:(double * __nonnull)v localW:(double * __nonnull)w x:(double)x y:(double)y z:(double)z model:(FEMModel * __nonnull)model;
-(double)interpolate2DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v;
-(double)interpolate3DInElement:(Element_t * __nonnull)element nodalValues:(double * __nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(void)invertMatrix3x3:(double[_Nonnull][3])G inverted:(double[_Nonnull][3])GI detG:(double)detG;
-(void)getEdgeBasisElement:(Element_t * __nonnull)element wBasis:(double * __nonnull * __nonnull)wBasis rotWBasis:(double * __nonnull * __nonnull)rotWBasis basis:(double * __nonnull)basis dBasisdx:(double * __nonnull * __nonnull)dBasisdx;

-(void)deallocation;

@end
