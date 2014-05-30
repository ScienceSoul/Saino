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
-(void)sortGlobal:(FEMSolution *)solution alsoValues:(BOOL *)alsoValues;
-(void)setElementInGlobal:(FEMSolution *)solution row:(int)i col:(int)j value:(double)value;
-(void)addToElementInGlobal:(FEMSolution *)solution row:(int)i col:(int)j value:(double)value;
-(void)glueLocalMatrix:(double **)localMatrix inGlobal:(FEMSolution *)solution numberOfNodes:(int)n dofs:(int)dofs indexes:(int *)indexes;
-(void)setSymmetricDirichletInGlobal:(FEMSolution *)solution atIndex:(int)n value:(double)value;
// CRS Matrix-vector multiply
-(void)matrixVectorMultiplyInGlobal:(FEMSolution *)solution vector:(double *)u result:(double *)v;
-(void)complexMatrixVectorMultiplyInGlobal:(FEMSolution *)solution vector:(double complex *)u result:(double complex *)v;
-(void)fctlLowOrderInSolution:(FEMSolution *)solution orMatrix:(FEMMatrix *)matrix;

-(void)glueLocalMatrix:(double **)localMatrix inMatrix:(FEMMatrix *)matrix numberOfNodes:(int)numberOfNodes dofs:(int)dofs indexes:(int *)indexes;
-(void)makeMatrixIndex:(FEMMatrix *)a row:(int)i col:(int)j;
-(void)zeroRowInMatrix:(FEMMatrix *)a numberOfRows:(int)n;
-(void)sortMatrix:(FEMMatrix *)a alsoValues:(BOOL *)alsoValues;
-(void)setElementInMatrix:(FEMMatrix *)a row:(int)i col:(int)j value:(double)value;
-(void)applyProjector:(FEMMatrix *)pMatrix values:(double *)u permutation:(int *)uperm values:(double *)v permutation:(int *)vperm transpose:(BOOL *)trans;
-(void)matrixVectorMultiply:(FEMMatrix *)matrix vector:(double *)u result:(double *)v;
-(void)complexMatrixVectorMultiply:(FEMMatrix *)matrix vector:(double complex *)u result:(double complex *)v;
-(void)zeroMatrix:(FEMMatrix *)matrix;


@end
