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

-(FEMMatrix * __nonnull)createMatrixWithNumberOfRows:(int)rows subBand:(int)subBand symmetric:(BOOL)symmetric allocateValues:(BOOL)allocateValues;
-(void)zeroRowInGlobal:(FEMSolution * __nonnull)solution numberOfRows:(int)n;
-(void)setElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)addToElementInGlobal:(FEMSolution * __nonnull)solution row:(int)i col:(int)j value:(double)value;
-(void)glueLocalMatrix:(double * __nonnull * __nonnull)localMatrix inGlobal:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dofs:(int)dofs indexes:(int * __nonnull)indexes;
-(void)sBand_setDirichlet:(FEMSolution * __nonnull)solution orderedNumber:(int)n value:(double)value;
// Band Matrix-vector multiply
-(void)matrixVectorMultiplyInGlobal:(FEMSolution * __nonnull)solution vector:(double * __nonnull)u result:(double * __nonnull)v;
-(void)zeroRowInMatrix:(FEMMatrix * __nonnull)matrix numberOfRows:(int)n;
-(void)setElementInMatrix:(FEMMatrix * __nonnull)matrix row:(int)i col:(int)j value:(double)value;
-(void)matrixVectorMultiply:(FEMMatrix * __nonnull)matrix vector:(double * __nonnull)u result:(double * __nonnull)v;
@end
