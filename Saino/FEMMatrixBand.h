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

-(void)BAND_glueLocalMatrix:(FEMSolution *)solution: (double **)matrix: (int)n: (int)dofs: (int *)indexes;

@end
