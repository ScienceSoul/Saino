//===----------------------------------------------------------------------===//
//  FEMListUtilities.h
//  Saino
//
//  Created by Hakime Seddik on 18/04/12.
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
#import "FEMModel.h"
#import "FEMValueList.h"

@interface FEMListUtilities : NSObject {
    NSMutableDictionary *_timers;
}

@property(nonatomic, strong, nonnull) NSMutableDictionary <NSString *, NSNumber *> *timers;

+(id _Nonnull)sharedListUtilities;
+(void)selfDestruct;

-(void)listParseDependencies:(NSArray * _Nonnull)dependencies index:(int)ind name:(NSString * _Nonnull)name toValues:(double * _Nonnull)t count:(int * _Nonnull)count model:(FEMModel * _Nonnull)model allGlobal:(BOOL * _Nonnull)allGlobal;

-(void)listSetNameSpace:(NSString * _Nullable)str;
-(char * _Nonnull)listGetNameSpaceForVariable:(NSString * _Nonnull)varName;
-(NSString * _Nullable)listGetString:(FEMModel * _Nullable)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName info:(BOOL * _Nonnull)found;
-(BOOL)listGetReal:(FEMModel * _Nonnull)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName numberOfNodes:(int)n indexes:(int * _Nonnull)nodeIndexes buffer:(listBuffer * _Nonnull)result minValue:(double * _Nullable)minv maxValue:(double * _Nullable)maxv;
-(double)listGetValueParameter:(FEMModel * _Nullable)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName value:(double)value info:(BOOL * _Nonnull)found minValue:(double * _Nullable)minv maxValue:(double * _Nullable)maxv;
-(BOOL)listGetRealArray:(FEMModel * _Nonnull)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName numberOfNodes:(int)n indexes:(int * _Nonnull)nodeIndexes buffer:(listBuffer * _Nonnull)result;
-(double)listGetConstReal:(FEMModel * _Nullable)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName info:(BOOL * _Nonnull)found minValue:(double * _Nullable)minv maxValue:(double * _Nullable)maxv;
-(BOOL)listGetConstRealArray:(FEMModel * _Nullable)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName buffer:(listBuffer * _Nonnull)result;

-(BOOL)listGetIntegerArray:(FEMModel * _Nullable)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName buffer:(listBuffer * _Nonnull)result;
-(int)listGetInteger:(FEMModel * _Nullable)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName info:(BOOL * _Nonnull)found minValue:(int * _Nullable)minv maxValue:(int * _Nullable)maxv;

-(BOOL)listGetLogical:(FEMModel * _Nullable)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName info:(BOOL * _Nonnull)found;

-(BOOL)listGetDerivativeValue:(FEMModel * _Nonnull)model inArray:(NSArray * _Nonnull)array forVariable:(NSString * _Nonnull)varName numberOfNodes:(int)n indexes:(int * _Nonnull)nodeIndexes buffer:(listBuffer * _Nonnull)result;

-(FEMValueList * _Nullable)listFindVariable:(NSString * _Nonnull)varName inArray:(NSArray * _Nonnull)array;
-(BOOL)listCheckPresentVariable:(NSString * _Nonnull)varName inArray:(NSArray * _Nonnull)array;
-(FEMValueList * _Nullable)listFindPrefix:(NSString * _Nonnull)prefix inArray:(NSArray * _Nonnull)array info:(BOOL * _Nonnull)found;
-(BOOL)listCheckPrefix:(NSString * _Nonnull)prefix inArray:(NSArray * _Nonnull)array;

-(void)addStringInClassList:(id _Nonnull)className theVariable:(NSString * _Nonnull)varName withValue:(NSString * _Nonnull)value;
-(void)addLogicalInClassList:(id _Nonnull)className theVariable:(NSString * _Nonnull)varName withValue:(BOOL)value;
-(void)addIntegerInClassList:(id _Nonnull)className theVariable:(NSString * _Nonnull)varName withValue:(int * _Nullable)value orUsingBlock:(double (^ _Nullable)())block;
-(void)addIntegerArrayInClassList:(id _Nonnull)className theVariable:(NSString * _Nonnull)varName withValues:(int * _Nullable)values size:(int)n orUsingBlock:(double (^ _Nullable)())block;
-(void)addConstRealInClassList:(id _Nonnull)className theVariable:(NSString * _Nonnull)varName withValue:(double * _Nullable)value orUsingBlock:(double (^ _Nullable)())block string:(NSString * _Nullable)str;
-(void)addConstRealArrayInClassList:(id _Nonnull)className theVariable:(NSString * _Nonnull)varName withValues:(double * _Nullable * _Nullable)fvalues size1:(int)m size2:(int)n orUsingBlock:(double (^ _Nullable)())block string:(NSString * _Nullable)str;
-(void)addBlockInClassList:(id _Nonnull)className theVariable:(NSString * _Nonnull)varName usingBlock:(double (^ _Nonnull)(double * _Nullable variablesValues))block dependencies:(NSArray * _Nullable)dependencies;

-(BOOL)checkElementEquation:(FEMModel * _Nonnull)model forElement:(Element_t * _Nonnull)element andEquation:(NSString * _Nonnull)equation;

-(BOOL)listCheckPresentAnyBoundaryCondition:(FEMModel * _Nonnull)model name:(NSString * _Nonnull)name;
-(BOOL)listCheckPresentAnyBodyForce:(FEMModel * _Nonnull)model name:(NSString * _Nonnull)name;
-(BOOL)listGetLogicalAnyBoundaryCondition:(FEMModel * _Nonnull)model name:(NSString * _Nonnull)name;

-(void)checkTimer:(NSString * _Nonnull)timerName deleteTimer:(BOOL * _Nullable)deleteTimer resetTimer:(BOOL * _Nullable)resetTimer model:(FEMModel * _Nonnull)model;
-(void)resetTimer:(NSString * _Nonnull)timerName model:(FEMModel * _Nonnull)model;
-(void)deletTimer:(NSString * _Nonnull)timerName;

-(void)deallocation;

@end
