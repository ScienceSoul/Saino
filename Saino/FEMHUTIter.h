//===----------------------------------------------------------------------===//
//  HUTIter.h
//  Saino
//
//  Created by Hakime Seddik on 16/02/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
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
-(void)dbicgstabSolveInSolution:(FEMSolution * _Nonnull)solution matrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double *_Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double *_Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;
// BI-CGSTAB2
-(void)dbicgstab2SolveInSolution:(FEMSolution * _Nonnull)solution matrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double * _Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;
// TFQMR
-(void)dtfqmrSolveInSolution:(FEMSolution * _Nonnull)solution matrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double * _Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;
// CG
-(void)dcgSolveInSolution:(FEMSolution * _Nonnull)solution matrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double * _Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;
// CGS
-(void)dcgsSolveInSolution:(FEMSolution * _Nonnull)solution matrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double * _Nonnull)dpar pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;
// GMRES
-(void)dgmresSolveInSolution:(FEMSolution * _Nonnull)solution matrix:(FEMMatrix * _Nonnull)matrix ndim:(int)ndim wrkdim:(int)wrkdim result:(double * _Nonnull)x rhs:(double * _Nonnull)b ipar:(int * _Nonnull)ipar dpar:(double *_Nonnull)dpar work:(double * _Nonnull * _Nonnull)work pcondlMethod:(SEL _Nonnull)pcondlMethod pcondrMethod:(SEL _Nonnull)pcondrMethod matvecMethod:(SEL _Nonnull)matvecMethod mstopMethod:(SEL _Nonnull)mstopMethod;

@end

