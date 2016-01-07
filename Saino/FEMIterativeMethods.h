//
//  FEMIterativeMethods.h
//  Saino
//
//  Created by Seddik hakime on 10/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMatrix.h"

@interface FEMIterativeMethods : NSObject

// SGS
-(void)dsgsSolveMatrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double * __nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
// JACOBI
-(void)djacobiSolveMatrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double * __nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
// BI-CGSTAB(l)
-(void)dbicgstablSolveMatrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double * __nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
// GCR
-(void)dgcrSolveMatrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double * __nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;
// Richardson
-(void)drichardsonSolveMatrix:(FEMMatrix * __nonnull)matrix ndim:(int)ndim result:(double * __nonnull)x rhs:(double * __nonnull)b ipar:(int * __nonnull)ipar dpar:(double * __nonnull)dpar pcondlMethod:(SEL __nonnull)pcondlMethod pcondrMethod:(SEL __nonnull)pcondrMethod matvecMethod:(SEL __nonnull)matvecMethod mstopMethod:(SEL __nonnull)mstopMethod;


@end
