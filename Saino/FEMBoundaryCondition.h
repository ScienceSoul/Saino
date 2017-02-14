//===----------------------------------------------------------------------===//
//  FEMBoundaryCondition.h
//  Saino
//
//  Created by Hakime Seddik on 27/03/12.
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
//
//------------------------------------------------------------------------------
//
// Boundary conditions
//
//------------------------------------------------------------------------------

#import <Foundation/Foundation.h>
#import "FEMMatrix.h"
#import "FEMValueList.h"

@interface FEMBoundaryCondition : NSObject {
    
    int _tag;
    NSMutableArray *_valuesList;
    FEMMatrix *_pMatrix;
}

@property(nonatomic, assign) int tag;
@property(nonatomic, strong, nonnull) NSMutableArray <FEMValueList *> *valuesList;
@property(nonatomic, strong, nullable) FEMMatrix *pMatrix;


@end
