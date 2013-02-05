//
//  HUTIter.h
//  Saino
//
//  Created by Hakime Seddik on 16/02/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <stdarg.h>
#import <string.h>

#import "FEMPrecondition.h"
#import "huti_defs.h"

@interface FEMHUTIter : NSObject {
    
    int huti_num_of_procs;
    BOOL huti_init_done;
}

-(void)hutiInit;
-(void)hutiExit;

// BI-CGSTAB
-(void)dbicgstabSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// BI-CGSTAB2
-(void)dbicgstab2SolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// TFQMR
-(void)dtfqmrSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// CG
-(void)dcgSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// CGS
-(void)dcgsSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// GNRES
-(void)dgmresSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// SGS
-(void)dsgsSolveInsolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// JACOBI
-(void)djacobiSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// BI-CGSTAB(l)
-(void)dbicgstablSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// GCR
-(void)dgcrSolveInSolution:(FEMSolution *)solution ndim:(int)ndim wrkdim:(int)wrkdim ipar:(int *)ipar dpar:(double *)dpar work:(double **)work pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;


@end

