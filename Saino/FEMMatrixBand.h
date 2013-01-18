//
//  FEMMatrixBand.h
//  Saino
//
//  Created by Hakime Seddik on 23/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMSolution.h"
#import "Constructors.h"

@interface FEMMatrixBand : NSObject

-(FEMMatrix *)createMatrixWithNumberOfRows:(int)rows subBand:(int)subBand symmetric:(BOOL)symmetric allocateValues:(BOOL)allocateValues;
-(void)zeroRowInGlobal:(FEMSolution *)solution numberOfRows:(int)n;
-(void)setMatrixElementInGlobal:(FEMSolution *)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)addToMatrixElementInGlobal:(FEMSolution *)solution atIndex:(int)i andIndex:(int)j value:(double)value;
-(void)glueLocalMatrixInGlobal:(FEMSolution *)solution matrix:(double **)matrix numberOfNodes:(int)n dofs:(int)dofs indexes:(int *)indexes;
-(void)sBand_setDirichlet:(FEMSolution *)solution orderedNumber:(int)n value:(double)value;

-(void)zeroRowInMatrix:(FEMMatrix *)a numberOfRows:(int)n;
-(void)setMatrixElementInMatrix:(FEMMatrix *)a atIndex:(int)i andIndex:(int)j value:(double)value;


@end
