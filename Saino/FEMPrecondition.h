//
//  FEMPreconditioners.h
//  Saino
//
//  Created by Hakime Seddik on 27/01/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMSolution.h"

#import "memory.h"
#import "Constructors.h"
#import "Utils.h"
#import "TimeProfile.h"


@interface FEMPrecondition : NSObject

// CRS Matrix-vector multiply
-(void)CRS_MatrixVectorMultiply:(FEMSolution *)solution: (double *)u: (double *)v;
-(void)CRS_ComplexMatrixVectorMultiply:(FEMSolution *)solution: (double complex *)u: (double complex *)v;

@end

