//===----------------------------------------------------------------------===//
//  FEMProjector.m
//  Saino
//
//  Created by Seddik hakime on 09/08/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
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

#import "FEMProjector.h"

@implementation FEMProjector

@synthesize mesh = _mesh;
@synthesize matrix = _matrix;
@synthesize tMatrix = _tMatrix;
@synthesize next = _next;

- (id)init
{
    self = [super init];
    if (self) {
        _mesh = [[FEMMesh alloc] init];
        _matrix = [[FEMMatrix alloc] init];
        _tMatrix = [[FEMMatrix alloc] init];
        _next = nil;
    }
    
    return self;
}



@end
