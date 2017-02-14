//===----------------------------------------------------------------------===//
//  FEMProjector.h
//  Saino
//
//  Created by Seddik hakime on 09/08/12.
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
#import "FEMMatrix.h"
#import "FEMMesh.h"

@interface FEMProjector : NSObject {
    
    FEMMesh *_mesh;
    FEMMatrix *_matrix;
    FEMMatrix *_tMatrix;
    FEMProjector *_next;
}

@property(nonatomic, strong, nonnull) FEMMesh *mesh;
@property(nonatomic, strong, nonnull) FEMMatrix *matrix;
@property(nonatomic, strong, nonnull) FEMMatrix *tMatrix;
@property(nonatomic, strong, nullable) FEMProjector *next;

@end
