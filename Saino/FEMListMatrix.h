//===----------------------------------------------------------------------===//
//  FEMListMatrix.h
//  Saino
//
//  Created by Seddik hakime on 07/01/13.
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
#import "Constructors.h"

@interface FEMListMatrix : NSObject {
    
    int _listMatrixGrowth;
}

@property(nonatomic, assign) int listMatrixGrowth;

-(ListMatrix_t * __nonnull)allocateMatrix:(int)n;
-(void)freeMatrix:(ListMatrix_t * __nonnull)list size:(int)n;
-(ListMatrix_t * __nonnull)enlargeMatrix:(ListMatrix_t * __nonnull)matrix toSize:(int)n;
-(ListMatrixEntry_t * __nonnull)getMatrixIndexInListMatrix:(ListMatrix_t * __nullable)list atIndex:(int)k1 andIndex:(int)k2;
-(void)addToMatrixElement:(ListMatrix_t * __nonnull)list atIndex:(int)k1 andIndex:(int)k2 value:(double)value setValue:(BOOL * __nullable)setValue;
-(void)convertToCRSMatrix:(FEMMatrix * __nonnull)matrix;

@end
