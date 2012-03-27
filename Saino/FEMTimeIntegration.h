//
//  FEMTimeIntegration.h
//  Saino
//
//  Created by Hakime Seddik on 14/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMSolution.h"
#import "FEMMesh.h"

@interface FEMTimeIntegration : NSObject


-(void)fractionStep:(FEMSolution *)solution: (int)n: (double)dt: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double *)prevSolution: (int *)rows;
-(void)bdfLocal:(FEMSolution *)solution: (int)n: (double)dt: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double **)prevSolution: (int)order: (int *)rows: (int *)cols;
-(void)vbdfLocal:(FEMSolution *)solution: (int)n: (double *)dts: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double **)prevSolution: (int)order: (int *)rows: (int *)cols;
-(void)newMarkBeta:(FEMSolution *)solution: (int)n: (double)dt: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double *)prevSolution: (double)beta: (int *)rows;


@end