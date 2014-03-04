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
-(void)dsgsSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// JACOBI
-(void)djacobiSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// BI-CGSTAB(l)
-(void)dbicgstablSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// GCR
-(void)dgcrSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;
// Richardson
-(void)drichardsonSolveMatrix:(FEMMatrix *)matrix ndim:(int)ndim result:(double *)x rhs:(double *)b ipar:(int *)ipar dpar:(double *)dpar pcondlMethod:(SEL)pcondlMethod pcondrMethod:(SEL)pcondrMethod matvecMethod:(SEL)matvecMethod mstopMethod:(SEL)mstopMethod;


@end
