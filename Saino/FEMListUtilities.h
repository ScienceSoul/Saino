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
#import "FEMUtilities.h"

#import "Utils.h"

@interface FEMListUtilities : NSObject 



-(void)listParseStrToValues:(FEMModel *)model: (NSString *)str: (int)ind: (NSString *)name: (double *)t: (int)count;

-(BOOL)listGetReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes resultArray:(double *)result minValue:(double *)minv maxValue:(double *)maxv;
-(BOOL)listGetRealArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes resultArray:(double ***)result;
-(double)listGetConstReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL)found minValue:(double *)minv maxValue:(double *)maxv;
-(BOOL)listGetConstRealArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName resultArray:(double **)result;

-(BOOL)listGetIntegerArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName resultArray:(int *)result;
-(int)listGetInteger:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL)found minValue:(int *)minv maxValue:(int *)maxv;

-(BOOL)listGetLogical:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL)found;

-(FEMValueList *)listFindVariable:(NSString *)varName inArray:(NSArray *)array;
-(BOOL)listCheckPresentVariable:(NSString *)varName inArray:(NSArray *)array;

-(void)addIntegerArrayInClassList:(id)className theVariable:(NSString *)varName withValues:(int *)values numberOfNodes:(int)n;


@end
