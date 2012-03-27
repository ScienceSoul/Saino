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
-(void)dbicgstabSolve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
// BI-CGSTAB2
-(void)dbicgstab2Solve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
// TFQMR
-(void)dtfqmrSolve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
// CG
-(void)dcgSolve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
// CGS
-(void)dcgsSolve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
// GNRES
-(void)dgmresSolve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
// SGS
-(void)dsgsSolve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
// JACOBI
-(void)djacobiSolve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
// BI-CGSTAB(l)
-(void)dbicgstablSolve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;
// GCR
-(void)dgcrSolve:(FEMSolution *)solution: (int)ndim: (int)wrkdim: (int *)ipar: (double *)dpar: (double **)work: (SEL)pcondlMethod: (SEL)pcondrMethod: (SEL)matvecMethod: (SEL)mstopMethod;


@end

