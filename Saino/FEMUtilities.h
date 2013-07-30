//
//  FEMUtilities.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMatrix.h"
#import "FEMVariable.h"
#import "FEMSolution.h"

@interface FEMUtilities : NSObject {
    NSString *_ext;
    NSString *_appSupportSubpath;
}

@property(nonatomic, strong) NSString *ext;
@property(nonatomic, strong) NSString *appSupportSubpath;

-(void)zeroTheNumberOfRows:(int)n inMatrix:(FEMMatrix *)a;
-(void)setMatrixElement:(FEMMatrix *)a atIndex:(int)i andIndex:(int)j value:(double)value;
-(int)initialPermutationInMesh:(FEMMesh *)aMesh model:(FEMModel *)aModel solution:(FEMSolution *)aSolution equation:(NSString *)str permutation:(int *)perm DGSolution:(BOOL *)dg globalBubbles:(BOOL *)gb;
-(FEMVariable *)getVariableFrom:(NSMutableArray *)anArray model:(FEMModel *)aModel name:(NSString *)name onlySearch:(BOOL *)only maskName:(NSString *)maskName info:(BOOL *)found;
-(void)addVariableTo:(NSMutableArray *)anArray mesh:(FEMMesh *)aMesh solution:(FEMSolution *)aSolution name:(NSString *)name dofs:(int )dofs container:(variableArraysContainer *)aContainer component:(BOOL)component ifOutput:(BOOL *)output ifSecondary:(BOOL *)secondary type:(int *)aType;
-(void)addVectorTo:(NSMutableArray *)anArray mesh:(FEMMesh *)aMesh solution:(FEMSolution *)aSolution name:(NSString *)name dofs:(int *)dofs container:(variableArraysContainer *)aContainer ifOutput:(BOOL *)output ifSecondary:(BOOL *)secondary global:(BOOL *)global initValue:(double *)initValue;
-(double)interpolateCurveTvalues:(double *)tValues fValues:(double *)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double *)cubicCoeff;
-(double)derivateCurveTvalues:(double *)tValues fValues:(double *)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double *)cubicCoeff;
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
-(void)addEquationToSolution:(FEMSolution *)solution model:(FEMModel *)model transient:(BOOL)transient;

-(BOOL)isFileNameQualified:(NSString *)file;
-(NSMutableString *)nextFreeFileName:(NSString *)fileName0 suffix:(NSString *)suffix0 lastExisting:(BOOL *)lastExisting;

// Load bundle for plug-ins support
-(NSBundle *)loadBundle:(NSString *)bundleName;

// Plug-in validation
-(BOOL)plugInClassIsValid:(Class)plugInClass;

@end
