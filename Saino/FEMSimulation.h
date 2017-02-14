//===----------------------------------------------------------------------===//
//  FEMSimulation.h
//  Saino
//
//  Created by Hakime Seddik on 26/04/12.
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
#import "FEMValueList.h"

@interface FEMSimulation : NSObject {
    
    int _tag;
    NSMutableArray *_valuesList;
}

@property(nonatomic, assign) int tag;
@property(nonatomic, strong, nonnull) NSMutableArray <FEMValueList *> *valuesList;

@end
