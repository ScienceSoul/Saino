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

-(void)zeroRowInGlobal:(FEMSolution *)solution: (int)n;
-(void)setMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)addToMatrixElementInGlobal:(FEMSolution *)solution: (int)i: (int)j: (double)value;
-(void)glueLocalMatrixInGlobal:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes;
-(void)sBand_setDirichlet:(FEMSolution *)solution: (int)n: (double)value;

-(void)zeroRowInMatrix:(Matrix_t *)a: (int)n;
-(void)setMatrixElementInMatrix:(Matrix_t *)a: (int)i: (int)j: (double)value;


@end
