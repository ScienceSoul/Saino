//===----------------------------------------------------------------------===//
//  FEMIterativeMethods.h
//  Saino
//
//  Created by Seddik hakime on 10/06/13.
//  Copyright (c) 2013 ScienceSoul. All rights reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

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
