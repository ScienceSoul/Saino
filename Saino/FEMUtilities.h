//===----------------------------------------------------------------------===//
//  FEMUtilities.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
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
#import "FEMMatrix.h"
#import "FEMSolution.h"

@interface FEMUtilities : NSObject {
    NSString * _ext;
    NSString * _appSupportSubpath;
}

@property(nonatomic, strong, nonnull) NSString *ext;
@property(nonatomic, strong, nonnull) NSString *appSupportSubpath;

-(void)zeroTheNumberOfRows:(int)n inMatrix:(FEMMatrix * __nonnull)a;
-(void)setMatrixElement:(FEMMatrix * __nonnull)a atIndex:(int)i andIndex:(int)j value:(double)value;
-(int)initialPermutationInMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model solution:(FEMSolution *__nonnull)solution equation:(NSString * __nonnull)str permutation:(int * __nonnull)perm DGSolution:(BOOL * __nullable)dg globalBubbles:(BOOL * __nullable)gb;
-(FEMVariable * __nullable)getVariableFrom:(NSMutableArray * __nonnull)anArray model:(FEMModel * __nonnull)model name:(NSString * __nonnull)name onlySearch:(BOOL * __nullable)only maskName:(NSString * __nullable)maskName info:(BOOL * __nonnull)found;
-(void)addVariableTo:(NSMutableArray * __nonnull)anArray mesh:(FEMMesh * __nullable)mesh solution:(FEMSolution * __nullable)solution name:(NSString * __nonnull)name dofs:(int )dofs container:(variableArraysContainer * __nonnull)aContainer component:(BOOL)component ifOutput:(BOOL * __nullable)output ifSecondary:(BOOL * __nullable)secondary type:(int * __nullable)aType;
-(void)addVectorTo:(NSMutableArray * __nonnull)anArray mesh:(FEMMesh * __nonnull)mesh solution:(FEMSolution * __nonnull)solution name:(NSString * __nonnull)name dofs:(int * __nullable)dofs container:(variableArraysContainer * __nonnull)aContainer ifOutput:(BOOL * __nullable)output ifSecondary:(BOOL * __nullable)secondary global:(BOOL * __nullable)global initValue:(double * __nullable)initValue;
-(double)interpolateCurveTvalues:(double * __nonnull)tValues fValues:(double * __nonnull)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double * __nullable)cubicCoeff;
-(double)derivateCurveTvalues:(double * __nonnull)tValues fValues:(double * __nonnull)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double * __nullable)cubicCoeff;
-(void)solveLinearSystem2x2:(double[][2])a afterSolve:(double * __nonnull)x rightHandSide:(double * __nonnull)b;
-(void)solveLinearSystem3x3:(double[][3])a afterSolve:(double * __nonnull)x rightHandSide:(double * __nonnull)b;
-(FEMMatrix * __nonnull)meshProjectorMesh1:(FEMMesh * __nonnull)mesh1 mesh2:(FEMMesh * __nonnull)mesh2 model:(FEMModel * __nonnull)model useQuadrantTree:(BOOL * __nullable)quadrantTree transpose:(BOOL * __nullable)transs;
-(double)cublicSplineX:(double * __nonnull)x Y:(double * __nonnull)y R:(double * __nonnull)r T:(double)t;

// Following two methods correspond to ComponentName routines in Elmer
-(NSString * __nullable)appendNameFromString:(NSString * __nonnull)string component:(int * __nullable)component;
-(NSString * __nonnull)appendNameFromVariable:(FEMVariable * __nonnull)variable component:(int * __nullable)component;

-(NSString * __nonnull)nextFreeKeyword:(NSString * __nonnull)keyword0 dictionary:(NSMutableDictionary * __nonnull)dictionary;
-(void)checkOptionsInSolution:(FEMSolution * __nonnull)solution;

-(void)addEquationBasicsToSolution:(FEMSolution * __nonnull)solution name:(NSString * __nonnull)name model:(FEMModel * __nonnull)model transient:(BOOL)transient;
-(void)addEquationToSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model transient:(BOOL)transient;

-(BOOL)isFileNameQualified:(NSString * __nonnull)file;
-(NSMutableString * __nullable)nextFreeFileName:(NSString * __nonnull)fileName0 suffix:(NSString * __nullable)suffix0 lastExisting:(BOOL * __nullable)lastExisting;

// Load bundle for plug-ins support
-(NSBundle * __nullable)loadBundle:(NSString * __nonnull)bundleName useApplicationSupportPath:(BOOL * __nullable)useApplicationSupportPath;

// Plug-in validation
-(BOOL)plugInClassIsValid:(Class __nonnull)plugInClass;

// Generate the colors to color a mesh
-(void)generateColor:(RGBColors * __nonnull)color;

@end
