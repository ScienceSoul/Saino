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

-(void)CRS_glueLocalMatrix:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes;

@end
