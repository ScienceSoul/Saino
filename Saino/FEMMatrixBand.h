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
-(void)zeroRowInGlobal:(FEMSolution *)solution: (int)n;
-(void)setMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)addToMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)glueLocalMatrixInGlobal:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes;
-(void)sBand_setDirichlet:(FEMSolution *)solution: (int)n: (double)value;

-(void)zeroRowInMatrix:(FEMMatrix *)a: (int)n;
-(void)setMatrixElementInMatrix:(FEMMatrix *)a: (int)i: (int)j: (double)value;


@end
