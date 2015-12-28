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

#import "FEMSolution.h"
#import "FEMPrecondition.h"
#import "huti_defs.h"

@interface FEMHUTIter : NSObject {
    
    int huti_num_of_procs;
    BOOL huti_init_done;
}

-(void)hutiInit;
-(void)hutiExit;

// BI-CGSTAB
-(void)dbicgstabSolveInSolution:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double *__nonnull)b ipar:(int * __nonnull)ipar dpar:(double *__nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
// BI-CGSTAB2
-(void)dbicgstab2SolveInSolution:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double * __nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
// TFQMR
-(void)dtfqmrSolveInSolution:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double * __nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
// CG
-(void)dcgSolveInSolution:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double * __nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
// CGS
-(void)dcgsSolveInSolution:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double * __nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
// GMRES
-(void)dgmresSolveInSolution:(FEMSolution * __nonnull)solution matrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim wrkdim:(int)wrkdim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double *__nonnull)dpar work:(double * __nonnull * __nonnull)work pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;

@end

