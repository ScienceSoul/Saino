//
//  FEMMatrixBand.h
//  Saino
//
//  Created by Hakime Seddik on 23/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMSolution.h"

@interface FEMMatrixBand : NSObject

-(FEMMatrix *)createMatrixWithNumberOfRows:(int)rows subBand:(int)subBand symmetric:(BOOL)symmetric allocateValues:(BOOL)allocateValues;
-(void)zeroRowInGlobal:(FEMSolution *)solution numberOfRows:(int)n;
-(void)setElementInGlobal:(FEMSolution *)solution row:(int)i col:(int)j value:(double)value;
-(void)addToElementInGlobal:(FEMSolution *)solution row:(int)i col:(int)j value:(double)value;
-(void)glueLocalMatrix:(double **)localMatrix inGlobal:(FEMSolution *)solution numberOfNodes:(int)n dofs:(int)dofs indexes:(int *)indexes;
-(void)sBand_setDirichlet:(FEMSolution *)solution orderedNumber:(int)n value:(double)value;
// Band Matrix-vector multiply
-(void)matrixVectorMultiplyInGlobal:(FEMSolution *)solution vector:(double *)u result:(double *)v;

-(void)zeroRowInMatrix:(FEMMatrix *)a numberOfRows:(int)n;
-(void)setElementInMatrix:(FEMMatrix *)a row:(int)i col:(int)j value:(double)value;
-(void)matrixVectorMultiply:(FEMMatrix *)a vector:(double *)u result:(double *)v;

@end
