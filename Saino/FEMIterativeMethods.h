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
-(void)dsgsSolveMatrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double * _Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;
// JACOBI
-(void)djacobiSolveMatrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double * _Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;
// BI-CGSTAB(l)
-(void)dbicgstablSolveMatrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double * _Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;
// GCR
-(void)dgcrSolveMatrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double * _Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;
// Richardson
-(void)drichardsonSolveMatrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double * _Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;


@end
