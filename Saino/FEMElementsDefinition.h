//===----------------------------------------------------------------------===//
//  FEMElementsDefinition.h
//  Saino
//
//  Created by Hakime Seddik on 05/06/12.
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

@interface FEMElementsDefinition : NSObject 

// This class is used by FEMElementDescription for the elements definition.

@property (nonatomic, strong, nonnull) NSNumber *dimension;
@property (nonatomic, strong, nonnull) NSString *topology;
@property (nonatomic, strong, nonnull) NSNumber *code;
@property (nonatomic, strong, nonnull) NSNumber *nodes;
@property (nonatomic, strong, nonnull) NSArray <NSNumber *> *nodeU;
@property (nonatomic, strong, nonnull) NSArray <NSNumber *> *nodeV;
@property (nonatomic, strong, nonnull) NSArray <NSNumber *> *nodeW;
@property (nonatomic, strong, nonnull) NSArray <NSNumber *> *basis;
@property (nonatomic, strong, nonnull) NSArray <NSNumber *> *gaussPoints;
@property (nonatomic, strong, nonnull) NSNumber *stabilization;

@end
