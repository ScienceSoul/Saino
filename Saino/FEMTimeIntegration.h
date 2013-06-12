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


-(void)fractionalStepInSolution:(FEMSolution *)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double **)massMatrix stiffMatrix:(double **)stiffMatrix force:(double *)force prevSolution:(double *)prevSolution rows:(int *)rows;
-(void)bdfLocalInSolution:(FEMSolution *)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double **)massMatrix stiffMatrix:(double **)stiffMatrix force:(double *)force prevSolution:(double **)prevSolution order:(int)order rows:(int *)rows cols:(int *)cols;
-(void)vbdfLocalInSolution:(FEMSolution *)solution numberOfNodes:(int)n dts:(double *)dts massMatrix:(double **)massMatrix stiffMatrix:(double **)stiffMatrix force:(double *)force prevSolution:(double **)prevSolution order:(int)order rows:(int *)rows cols:(int *)cols;
-(void)newMarkBetaInSolution:(FEMSolution *)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double **)massMatrix stiffMatrix:(double **)stiffMatrix force:(double *)force prevSolution:(double *)prevSolution beta:(double)beta rows:(int *)rows cols:(int *)cols;
-(void)bossakSecondOrder:(FEMSolution *)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double **)massMatrix dampMatrix:(double **)dampMatrix stiffMatrix:(double **)stiffMatrix force:(double *)force prevSolution1:(double *)prevSolution1 prevSolution2:(double *)prevSolution2 prevSolution3:(double *)prevSolution3 alpha:(double)alpha rows:(int *)rows cols:(int *)cols;

@end