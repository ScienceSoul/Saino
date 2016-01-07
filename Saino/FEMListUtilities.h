//
//  FEMListUtilities.h
//  Saino
//
//  Created by Hakime Seddik on 18/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"
#import "FEMValueList.h"

@interface FEMListUtilities : NSObject {
    NSMutableDictionary *_timers;
}

@property(nonatomic, strong, nonnull) NSMutableDictionary *timers;

-(void)listParseDependencies:(NSArray * __nonnull)dependencies index:(int)ind name:(NSString * __nonnull)name toValues:(double * __nonnull)t count:(int * __nonnull)count model:(FEMModel * __nonnull)model allGlobal:(BOOL * __nonnull)allGlobal;

-(void)listSetNameSpace:(NSString * __nonnull)str;
-(char * __nonnull)listGetNameSpaceForVariable:(char * __nonnull)varName;
-(NSString * __nullable)listGetString:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName info:(BOOL * __nonnull)found;
-(BOOL)listGetReal:(FEMModel * __nonnull)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName numberOfNodes:(int)n indexes:(int * __nonnull)nodeIndexes buffer:(listBuffer * __nonnull)result minValue:(double * __nullable)minv maxValue:(double * __nullable)maxv;
-(double)listGetValueParameter:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName value:(double)value info:(BOOL * __nonnull)found minValue:(double * __nullable)minv maxValue:(double * __nullable)maxv;
-(BOOL)listGetRealArray:(FEMModel * __nonnull)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName numberOfNodes:(int)n indexes:(int * __nonnull)nodeIndexes buffer:(listBuffer * __nonnull)result;
-(double)listGetConstReal:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName info:(BOOL * __nonnull)found minValue:(double * __nullable)minv maxValue:(double * __nullable)maxv;
-(BOOL)listGetConstRealArray:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName buffer:(listBuffer * __nonnull)result;

-(BOOL)listGetIntegerArray:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName buffer:(listBuffer * __nonnull)result;
-(int)listGetInteger:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName info:(BOOL * __nonnull)found minValue:(int * __nullable)minv maxValue:(int * __nullable)maxv;

-(BOOL)listGetLogical:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName info:(BOOL * __nonnull)found;

-(BOOL)listGetDerivativeValue:(FEMModel * __nonnull)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName numberOfNodes:(int)n indexes:(int * __nonnull)nodeIndexes buffer:(listBuffer * __nonnull)result;

-(FEMValueList * __nullable)listFindVariable:(NSString * __nonnull)varName inArray:(NSArray * __nonnull)array;
-(BOOL)listCheckPresentVariable:(NSString * __nonnull)varName inArray:(NSArray * __nonnull)array;
-(FEMValueList * __nullable)listFindPrefix:(NSString * __nonnull)prefix inArray:(NSArray * __nonnull)array info:(BOOL * __nonnull)found;
-(BOOL)listCheckPrefix:(NSString * __nonnull)prefix inArray:(NSArray * __nonnull)array;

-(void)addStringInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValue:(NSString * __nonnull)value;
-(void)addLogicalInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValue:(BOOL)value;
-(void)addIntegerInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValue:(int * __nullable)value orUsingBlock:(double (^ __nullable)())block;
-(void)addIntegerArrayInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValues:(int * __nullable)values size:(int)n orUsingBlock:(double (^ __nullable)())block;
-(void)addConstRealInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValue:(double * __nullable)value orUsingBlock:(double (^ __nullable)())block string:(NSString * __nullable)str;
-(void)addConstRealArrayInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValues:(double * __nullable * __nullable)fvalues size1:(int)m size2:(int)n orUsingBlock:(double (^ __nullable)())block string:(NSString * __nullable)str;
-(void)addBlockInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName usingBlock:(double (^ __nonnull)(double * __nullable variablesValues))block dependencies:(NSArray * __nullable)dependencies;

-(BOOL)checkElementEquation:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element andEquation:(NSString * __nonnull)equation;

-(BOOL)listCheckPresentAnyBoundaryCondition:(FEMModel * __nonnull)model name:(NSString * __nonnull)name;
-(BOOL)listCheckPresentAnyBodyForce:(FEMModel * __nonnull)model name:(NSString * __nonnull)name;

-(void)checkTimer:(NSString * __nonnull)timerName deleteTimer:(BOOL * __nullable)deleteTimer resetTimer:(BOOL * __nullable)resetTimer model:(FEMModel * __nonnull)model;
-(void)resetTimer:(NSString * __nonnull)timerName model:(FEMModel * __nonnull)model;
-(void)deletTimer:(NSString * __nonnull)timerName;

@end
