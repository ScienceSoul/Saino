//
//  FEMElementDescription.h
//  Saino
//
//  Created by Hakime Seddik on 05/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <Accelerate/Accelerate.h>

#import "FEMElementsDefinition.h"
#import "FEMSolution.h"
#import "FEMNumericIntegration.h"
#import "FEMUtilities.h"

#import "Constructors.h"
#import "memory.h"

@interface FEMElementDescription : NSObject {
    
}

-(void)deallocation;
-(void)addDescriptionOfElement:(ElementType_t)element withBasisTerms:(int *)terms;
-(void)initElementDescriptions;
-(int **)getEdgeMap:(int)elementFamily;
-(double)elementDiameter:(Element_t *)element nodes:(Nodes_t *)nodes;
-(void)computeStabilizationParameterInElement:(Element_t *)element nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh numberOfNodes:(int)n mk:(double)mk hk:(double *)hk;
-(ElementType_t *)getElementType:(int)code inMesh:(FEMMesh *)mesh stabilization:(BOOL *)computeStab;
-(double)firstDerivative1DInElement:(Element_t *)element nodalValues:(double *)x evalutationPoint:(double)u;
-(double)firstDerivativeU2DInElement:(Element_t *)element nodalValues:(double *)x evaluatedAt:(double)u andAt:(double)v;
-(double)firstDerivativeV2DInElement:(Element_t *)element nodalValues:(double *)x evaluatedAt:(double)u andAt:(double)v;
-(double)firstDerivativeU3DInElement:(Element_t *)element nodalValues:(double *)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(double)firstDerivativeV3DInElement:(Element_t *)element nodalValues:(double *)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(double)firstDerivativeW3DInElement:(Element_t *)element nodalValues:(double *)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;
-(double)interpolateInElement:(Element_t *)element nodalValues:(double *)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w withBasis:(double *)basis;
-(void)checkNormalDirectionInBDElement:(Element_t *)boundary forNormals:(double *)normals mesh:(FEMMesh *)mesh x:(double)x y:(double)y z:(double)z turn:(BOOL *)turn;
-(void)normalVectorForBDElement:(Element_t *)boundary boundaryNodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh paraU:(double *)u0 paraV:(double *)v0 check:(BOOL *)check normals:(double *)normals;
-(void)globalToLocalFromElement:(Element_t *)element elementNodes:(Nodes_t *)nodes localU:(double *)u localV:(double *)v localW:(double *)w x:(double)x y:(double)y z:(double)z model:(FEMModel *)aModel;
-(double)interpolate2DInElement:(Element_t *)element nodalValues:(double *)x evaluatedAt:(double)u andAt:(double)v;
-(double)interpolate3DInElement:(Element_t *)element nodalValues:(double *)x evaluatedAt:(double)u andAt:(double)v andAt:(double)w;

@end
