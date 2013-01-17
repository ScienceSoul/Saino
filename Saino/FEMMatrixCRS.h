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
-(void)zeroRowInGlobal:(FEMSolution *)solution: (int)n;
-(void)sortInGlobal:(FEMSolution *)solution: (BOOL *)alsoValues;
-(void)setMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)addToMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)glueLocalMatrixInGlobal:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes;
-(void)setSymmetricDirichletInGlobal:(FEMSolution *)solution: (int)n: (double)val;

-(void)glueLocalMatrixInMatrix:(FEMMatrix *)matrix localMatrix:(double **)localMatrix numberOfNodes:(int)numberOfNodes dofs:(int)dofs indexes:(int *)indexes;
-(void)makeMatrixIndex:(FEMMatrix *)a atIndex:(int)i  andIndex:(int)j;
-(void)zeroRowInMatrix:(FEMMatrix *)a: (int)n;
-(void)sortInMatrix:(FEMMatrix *)a: (BOOL *)alsoValues;
-(void)setMatrixElementInMatrix:(FEMMatrix *)a: (int)i: (int)j: (double)value;
-(void)applyProjector:(FEMMatrix *)pMatrix: (double *)u: (int *)uperm: (double *)v: (int *)vperm: (BOOL *)trans;


@end
