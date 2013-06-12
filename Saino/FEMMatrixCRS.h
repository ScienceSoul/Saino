//
//  FEMMatrixCRS.h
//  Saino
//
//  Created by Hakime Seddik on 23/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMSolution.h"

@interface FEMMatrixCRS : NSObject

-(FEMMatrix *)createMatrixWithNumberOfRows:(int)rows totalNonZeros:(int)totalNonZeros rowNonZeros:(int *)rowNonZeros degreesFreedom:(int)degreesFreedom reorder:(int *)reorder sizeOfReorder:(int)sizeOfReorder allocateValues:(BOOL)allocateValues;
-(void)zeroRowInGlobal:(FEMSolution *)solution numberOfRows:(int)n;
-(void)sortInGlobal:(FEMSolution *)solution alsoValues:(BOOL *)alsoValues;
-(void)setMatrixElementInGlobal:(FEMSolution *)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)addToMatrixElementInGlobal:(FEMSolution *)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)glueLocalMatrixInGlobal:(FEMSolution *)solution matrix:(double **)matrix numberOfNodes:(int)n dofs:(int)dofs indexes:(int *)indexes;
-(void)setSymmetricDirichletInGlobal:(FEMSolution *)solution atIndex:(int)n value:(double)value;
// CRS Matrix-vector multiply
-(void)matrixVectorMultiplyInGlobal:(FEMSolution *)solution multiplyVector:(double *)u resultVector:(double *)v;
-(void)complexMatrixVectorMultiplyInGlobal:(FEMSolution *)solution multiplyVector:(double complex *)u resultVector:(double complex *)v;
-(void)fctlLowOrderInSolution:(FEMSolution *)solution orMatrix:(FEMMatrix *)matrix;

-(void)glueLocalMatrixInMatrix:(FEMMatrix *)matrix localMatrix:(double **)localMatrix numberOfNodes:(int)numberOfNodes dofs:(int)dofs indexes:(int *)indexes;
-(void)makeMatrixIndex:(FEMMatrix *)a atIndex:(int)i  andIndex:(int)j;
-(void)zeroRowInMatrix:(FEMMatrix *)a numberOfRows:(int)n;
-(void)sortInMatrix:(FEMMatrix *)a alsoValues:(BOOL *)alsoValues;
-(void)setMatrixElementInMatrix:(FEMMatrix *)a atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)applyProjector:(FEMMatrix *)pMatrix values:(double *)u permutation:(int *)uperm values:(double *)v permutation:(int *)vperm transpose:(BOOL *)trans;
-(void)matrixVectorMultiplyInMatrix:(FEMMatrix *)matrix multiplyVector:(double *)u resultVector:(double *)v;
-(void)complexMatrixVectorMultiplyInMatrix:(FEMMatrix *)matrix multiplyVector:(double complex *)u resultVector:(double complex *)v;


@end
