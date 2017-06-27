//===----------------------------------------------------------------------===//
//  FEMValueList.m
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
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

#import "FEMValueList.h"
#import "memory.h"

@implementation FEMValueList

@synthesize model = _model;
@synthesize type = _type;
@synthesize nameLength = _nameLength;
@synthesize numberOfDependencies = _numberOfDependencies;
@synthesize lValue = _lValue;
@synthesize name = _name;
@synthesize cValue = _cValue;
@synthesize dependencies = _dependencies;
@synthesize block = _block;

- (id)init
{
    self = [super init];
    if (self) {
        _containers = (valueListArraysContainer*)malloc(sizeof(valueListArraysContainer));
        _containers->tValues = NULL;
        _containers->iValues = NULL;
        _containers->fValues = NULL;
        self.block = nil;
    }
    
    return self;
}

-(void)deallocation {
    if (_containers != NULL) {
        if (_containers->tValues != NULL) {
            free_dvector(_containers->tValues, 0, _containers->sizeTValues-1);
            _containers->tValues = NULL;
        }
        if (_containers->iValues != NULL) {
            free_ivector(_containers->iValues, 0, _containers->sizeIValues-1);
            _containers->iValues = NULL;
        }
        if (_containers->fValues != NULL) {
            free_d3tensor(_containers->fValues, 0, _containers->sizeFValues1-1, 0, _containers->sizeFValues2-1, 0, _containers->sizeFValues3-1);
            _containers->fValues = NULL;
        }
        free(_containers);
        _containers = NULL;
    }
}

-(valueListArraysContainer * _Nonnull)getContainers {
    
    return _containers;
}

@end
