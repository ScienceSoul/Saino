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

-(void)zeroRowInGlobal:(FEMSolution *)solution: (int)n;
-(void)sortInGlobal:(FEMSolution *)solution: (BOOL *)alsoValues;
-(void)setMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)addToMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)glueLocalMatrixInGlobal:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes;
-(void)setSymmetricDirichletInGlobal:(FEMSolution *)solution: (int)n: (double)val;

-(void)zeroRowInMatrix:(FEMMatrix *)a: (int)n;
-(void)sortInMatrix:(FEMMatrix *)a: (BOOL *)alsoValues;
-(void)setMatrixElementInMatrix:(FEMMatrix *)a: (int)i: (int)j: (double)value;


@end
