//
//  FEMUtilities.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMatrixCRS.h"
#import "FEMMatrixBand.h"
#import "FEMListUtilities.h"
#import "FEMInterpolation.h"

#import "memory.h"

@interface FEMUtilities : NSObject

-(FEMMatrix *)allocateMatrix;
-(void)zeroTheNumberOfRows:(int)n inMatrix:(FEMMatrix *)a;
-(void)setMatrixElement:(FEMMatrix *)a atIndex:(int)i andIndex:(int)j value:(double)value;
-(BOOL)checkEquationForElement:(Element_t *)element model:(FEMModel *)aModel equation:(NSString *)str;
-(int)initialPermutationInMesh:(FEMMesh *)aMesh model:(FEMModel *)aModel solution:(FEMSolution *)aSolution equation:(NSString *)str permutation:(int *)perm DGSolution:(BOOL *)dg globalBubbles:(BOOL *)gb;
-(FEMVariable *)getVariableFrom:(NSMutableArray *)anArray model:(FEMModel *)aModel name:(NSString *)name onlySearch:(BOOL *)only maskName:(NSString *)maskName info:(BOOL *)found;
-(void)addVariableTo:(NSMutableArray *)anArray mesh:(FEMMesh *)aMesh solution:(FEMSolution *)aSolution name:(NSString *)name dofs:(int )dofs container:(variableArraysContainer *)aContainer component:(BOOL)component ifOutput:(BOOL *)output ifSecondary:(BOOL *)secondary type:(int *)aType;
-(double)interpolateCurveTvalues:(double *)tValues fValues:(double *)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double *)cubicCoeff;
-(void)solveLinearSystem2x2:(double **)a afterSolve:(double *)x rightHandSide:(double *)b;
-(void)solveLinearSystem3x3:(double **)a afterSolve:(double *)x rightHandSide:(double *)b;
-(FEMMatrix *)meshProjector:(FEMMesh *)mesh1 secondmesh:(FEMMesh *)mesh2 model:(FEMModel *)aModel useQuadrantTree:(BOOL *)quadrantTree transpose:(BOOL *)trans;
-(double)cublicSplineX:(double *)x Y:(double *)y R:(double *)r T:(double)t;

// Following two methods correspond to ComponentName routines in Elmer
-(NSString *)appendNameFromString:(NSString *)string component:(int *)component;
-(NSString *)appendNameFromVariable:(FEMVariable *)variable component:(int *)component;

-(NSString *)nextFreeKeyword:(NSString *)keyword0 dictionary:(NSMutableDictionary *)dictionary;
-(void)checkOptionsInSolution:(FEMSolution *)solution;

-(void)addEquationBasicsToSolution:(FEMSolution *)solution name:(NSString *)name model:(FEMModel *)model transient:(BOOL)transient;

@end
