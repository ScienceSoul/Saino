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
#import "FEMBodyForce.h"
#import "FEMBoundaryCondition.h"
#import "FEMSimulation.h"
#import "FEMEquation.h"
#import "FEMUtilities.h"

#import "Utils.h"

@interface FEMListUtilities : NSObject



-(void)listParseStrToValues:(FEMModel *)model: (NSString *)str: (int)ind: (NSString *)name: (double *)t: (int)count;

-(NSString *)listGetString:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found;
-(BOOL)listGetReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes buffer:(listBuffer *)result minValue:(double *)minv maxValue:(double *)maxv;
-(BOOL)listGetRealArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes buffer:(listBuffer *)result;
-(double)listGetConstReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found minValue:(double *)minv maxValue:(double *)maxv;
-(BOOL)listGetConstRealArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName buffer:(listBuffer *)result;

-(BOOL)listGetIntegerArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName buffer:(listBuffer *)result;
-(int)listGetInteger:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found minValue:(int *)minv maxValue:(int *)maxv;

-(BOOL)listGetLogical:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found;

-(FEMValueList *)listFindVariable:(NSString *)varName inArray:(NSArray *)array;
-(BOOL)listCheckPresentVariable:(NSString *)varName inArray:(NSArray *)array;

-(void)addStringInClassList:(id)className theVariable:(NSString *)varName withValue:(NSString *)value;
-(void)addLogicalInClassList:(id)className theVariable:(NSString *)varName withValue:(BOOL)value;
-(void)addIntegerInClassList:(id)className theVariable:(NSString *)varName withValue:(int)value;
-(void)addIntegerArrayInClassList:(id)className theVariable:(NSString *)varName withValues:(int *)values numberOfNodes:(int)n;
-(void)addConstRealInClassList:(id)className theVariable:(NSString *)varName withValue:(double)value string:(NSString *)str;

-(BOOL)checkElementEquation:(FEMModel *)model forElement:(Element_t *)element andEquation:(NSString *)equation;


@end
