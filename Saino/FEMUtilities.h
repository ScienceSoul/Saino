//
//  FEMUtilities.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMatrixCRS.h"
#import "FEMMatrixBand.h"

#import "Constructors.h"
#import "memory.h"

@interface FEMUtilities : NSObject

-(Matrix_t *)allocateMatrix;
-(void)zeroTheNumberOfRows:(int)n inMatrix:(Matrix_t *)a;
-(void)setMatrixElement:(Matrix_t *)a: (int)i: (int)j: (double)value;

-(double)interpolateCurve:(double *)tValues: (double *)fValues: (double)t: (int)n;
-(void)solveLinearSystem2x2:(double **)a: (double *)x: (double *)b;
-(void)solveLinearSystem3x3:(double **)a: (double *)x: (double *)b;

@end
