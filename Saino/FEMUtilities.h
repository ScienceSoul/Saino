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

-(void)zeroTheNumberOfRows:(int)n inMatrix:(FEMMatrix * _Nonnull)a;
-(void)setMatrixElement:(FEMMatrix * _Nonnull)a atIndex:(int)i andIndex:(int)j value:(double)value;
-(int)initialPermutationInMesh:(FEMMesh * _Nonnull)mesh model:(FEMModel * _Nonnull)model solution:(FEMSolution *_Nonnull)solution equation:(NSString * _Nonnull)str permutation:(int * _Nonnull)perm DGSolution:(BOOL * _Nullable)dg globalBubbles:(BOOL * _Nullable)gb;
-(FEMVariable * _Nullable)getVariableFrom:(NSMutableArray * _Nonnull)anArray model:(FEMModel * _Nonnull)model name:(NSString * _Nonnull)name onlySearch:(BOOL * _Nullable)only maskName:(NSString * _Nullable)maskName info:(BOOL * _Nonnull)found;
-(void)addVariableTo:(NSMutableArray * _Nonnull)anArray mesh:(FEMMesh * _Nullable)mesh solution:(FEMSolution * _Nullable)solution name:(NSString * _Nonnull)name dofs:(int )dofs container:(variableArraysContainer * _Nonnull)aContainer component:(BOOL)component ifOutput:(BOOL * _Nullable)output ifSecondary:(BOOL * _Nullable)secondary type:(int * _Nullable)aType;
-(void)addVectorTo:(NSMutableArray * _Nonnull)anArray mesh:(FEMMesh * _Nonnull)mesh solution:(FEMSolution * _Nonnull)solution name:(NSString * _Nonnull)name dofs:(int * _Nullable)dofs container:(variableArraysContainer * _Nonnull)aContainer ifOutput:(BOOL * _Nullable)output ifSecondary:(BOOL * _Nullable)secondary global:(BOOL * _Nullable)global initValue:(double * _Nullable)initValue;
-(double)interpolateCurveTvalues:(double * _Nonnull)tValues fValues:(double * _Nonnull)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double * _Nullable)cubicCoeff;
-(double)derivateCurveTvalues:(double * _Nonnull)tValues fValues:(double * _Nonnull)fValues value:(double)t sizeOfTValues:(int)n cubicCoefficient:(double * _Nullable)cubicCoeff;
-(void)solveLinearSystem2x2:(double[_Nonnull][2])a afterSolve:(double * _Nonnull)x rightHandSide:(double * _Nonnull)b;
-(void)solveLinearSystem3x3:(double[_Nonnull][3])a afterSolve:(double * _Nonnull)x rightHandSide:(double * _Nonnull)b;
-(FEMMatrix * _Nonnull)meshProjectorMesh1:(FEMMesh * _Nonnull)mesh1 mesh2:(FEMMesh * _Nonnull)mesh2 model:(FEMModel * _Nonnull)model useQuadrantTree:(BOOL * _Nullable)quadrantTree transpose:(BOOL * _Nullable)transs;
-(double)cublicSplineX:(double * _Nonnull)x Y:(double * _Nonnull)y R:(double * _Nonnull)r T:(double)t;

// Following two methods correspond to ComponentName routines in Elmer
-(NSString * _Nullable)appendNameFromString:(NSString * _Nonnull)string component:(int * _Nullable)component;
-(NSString * _Nonnull)appendNameFromVariable:(FEMVariable * _Nonnull)variable component:(int * _Nullable)component;

-(NSString * _Nonnull)nextFreeKeyword:(NSString * _Nonnull)keyword0 dictionary:(NSMutableDictionary * _Nonnull)dictionary;
-(void)checkOptionsInSolution:(FEMSolution * _Nonnull)solution;

-(void)addEquationBasicsToSolution:(FEMSolution * _Nonnull)solution name:(NSString * _Nonnull)name model:(FEMModel * _Nonnull)model transient:(BOOL)transient;
-(void)addEquationToSolution:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model transient:(BOOL)transient;

-(BOOL)isFileNameQualified:(NSString * _Nonnull)file;
-(NSMutableString * _Nullable)nextFreeFileName:(NSString * _Nonnull)fileName0 suffix:(NSString * _Nullable)suffix0 lastExisting:(BOOL * _Nullable)lastExisting;

// Load bundle for plug-ins support
-(NSBundle * _Nullable)loadBundle:(NSString * _Nonnull)bundleName useApplicationSupportPath:(BOOL * _Nullable)useApplicationSupportPath;

// Plug-in validation
-(BOOL)plugInClassIsValid:(Class _Nonnull)plugInClass;

// Generate the colors to color a mesh
-(void)generateColor:(RGBColors * _Nonnull)color;

@end
