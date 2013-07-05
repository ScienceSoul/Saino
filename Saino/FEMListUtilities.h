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

@interface FEMListUtilities : NSObject

-(void)listParseStrToValues:(FEMModel *)model string:(NSString *)str index:(int)ind name:(NSString *)name values:(double *)t count:(int *)count allGlobal:(BOOL *)allGlobal;

-(void)listSetNameSpace:(NSString *)str;
-(BOOL)listGetNameSpace:(NSMutableString *)str;
-(NSString *)listGetString:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found;
-(BOOL)listGetReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes buffer:(listBuffer *)result minValue:(double *)minv maxValue:(double *)maxv;
-(BOOL)listGetRealArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes buffer:(listBuffer *)result;
-(double)listGetConstReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found minValue:(double *)minv maxValue:(double *)maxv;
-(BOOL)listGetConstRealArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName buffer:(listBuffer *)result;

-(BOOL)listGetIntegerArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName buffer:(listBuffer *)result;
-(int)listGetInteger:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found minValue:(int *)minv maxValue:(int *)maxv;

-(BOOL)listGetLogical:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found;

-(BOOL)listGetDerivativeValue:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes buffer:(listBuffer *)result;

-(FEMValueList *)listFindVariable:(NSString *)varName inArray:(NSArray *)array;
-(BOOL)listCheckPresentVariable:(NSString *)varName inArray:(NSArray *)array;

-(void)addStringInClassList:(id)className theVariable:(NSString *)varName withValue:(NSString *)value;
-(void)addLogicalInClassList:(id)className theVariable:(NSString *)varName withValue:(BOOL)value;
-(void)addIntegerInClassList:(id)className theVariable:(NSString *)varName withValue:(int)value;
-(void)addIntegerArrayInClassList:(id)className theVariable:(NSString *)varName withValues:(int *)values size:(int)n;
-(void)addConstRealInClassList:(id)className theVariable:(NSString *)varName withValue:(double)value string:(NSString *)str;
-(void)addConstRealArrayInClassList:(id)className theVariable:(NSString *)varName withValues:(double **)fvalues size1:(int)m size2:(int)n string:(NSString *)str;

-(BOOL)checkElementEquation:(FEMModel *)model forElement:(Element_t *)element andEquation:(NSString *)equation;

-(BOOL)listCheckPresentAnyBoundaryCondition:(FEMModel *)model name:(NSString *)name;
-(BOOL)listCheckPresentAnyBodyForce:(FEMModel *)model name:(NSString *)name;

@end
