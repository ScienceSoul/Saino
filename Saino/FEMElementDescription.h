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

+(id _Nonnull)sharedElementDescription;
+(void)selfDestruct;

-(void)addDescriptionOfElement:(ElementType_t * _Nonnull)element withBasisTerms:(int * _Nonnull)terms;
-(void)initElementDescriptions;
-(int * _Nonnull * _Nonnull)getEdgeMap:(int)elementFamily;
-(double)elementDiameter:(Element_t * _Nonnull)element nodes:(Nodes_t * _Nonnull)nodes;
-(void)computeStabilizationParameterInElement:(Element_t * _Nonnull)element nodes:(Nodes_t * _Nonnull)nodes mesh:(FEMMesh * _Nonnull)mesh numberOfNodes:(int)n mk:(double * _Nonnull)mk hk:(double * _Nullable)hk;
-(ElementType_t * _Nullable)getElementType:(int)code inMesh:(FEMMesh * _Nonnull)mesh stabilization:(BOOL * _Nullable)computeStab;
-(double)firstDerivative1DInElement:(Element_t * _Nonnull)element nodalValues:(double * _Nonnull)x evalutationPoint:(double)u;
-(double)firstDerivativeU2DInElement:(Element_t * _Nonnull)element nodalValues:(double * _Nonnull)x evaluatedAt:(double)u andAt:(double)v;
-(double)firstDerivativeV2DInElement:(Element_t * _Nonnull)element nodalValues:(double * _Nonnull)x evaluatedAt:(double)u andAt:(double)v;
-(double)firstDerivativeU3DInElement:(Element_t * _Nonnull)element nodalValues:(double * _Nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(double)firstDerivativeV3DInElement:(Element_t * _Nonnull)element nodalValues:(double * _Nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(double)firstDerivativeW3DInElement:(Element_t * _Nonnull)element nodalValues:(double * _Nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(double)interpolateInElement:(Element_t * _Nonnull)element nodalValues:(double * _Nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w withBasis:(double * _Nullable)basis;
-(void)checkNormalDirectionInBDElement:(Element_t * _Nonnull)boundary forNormals:(double * _Nonnull)normals mesh:(FEMMesh * _Nonnull)mesh x:(double)x y:(double)y z:(double)z turn:(BOOL * _Nullable)turn;
-(void)normalVectorForBDElement:(Element_t * _Nonnull)boundary boundaryNodes:(Nodes_t * _Nonnull)nodes mesh:(FEMMesh * _Nonnull)mesh paraU:(double * _Nullable)u0 paraV:(double * _Nullable)v0 check:(BOOL * _Nullable)check normals:(double * _Nonnull)normals;
-(void)globalToLocalFromElement:(Element_t * _Nonnull)element elementNodes:(Nodes_t * _Nonnull)nodes localU:(double * _Nonnull)u localV:(double * _Nonnull)v localW:(double * _Nonnull)w x:(double)x y:(double)y z:(double)z model:(FEMModel * _Nonnull)model;
-(double)interpolate2DInElement:(Element_t * _Nonnull)element nodalValues:(double * _Nonnull)x evaluatedAt:(double)u andAt:(double)v;
-(double)interpolate3DInElement:(Element_t * _Nonnull)element nodalValues:(double * _Nonnull)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(void)invertMatrix3x3:(double[_Nonnull][3])G inverted:(double[_Nonnull][3])GI detG:(double)detG;
-(void)getEdgeBasisElement:(Element_t * _Nonnull)element wBasis:(double * _Nonnull * _Nonnull)wBasis rotWBasis:(double * _Nonnull * _Nonnull)rotWBasis basis:(double * _Nonnull)basis dBasisdx:(double * _Nonnull * _Nonnull)dBasisdx;

-(void)deallocation;

@end
