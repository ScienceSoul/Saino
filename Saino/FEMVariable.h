//===----------------------------------------------------------------------===//
//  FEMVariable.h
//  Saino
//
//  Created by Seddik hakime on 27/07/12.
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
#import "Constructors.h"

@interface FEMVariable : NSObject {
    
    NSString *_name;
    id _primaryMesh;
    id _solution;
    int _nameLength;
    int _dofs;
    int _nonLinConverged;
    int _steadyConverged;
    int _nonLinIter;
    int _type;
    double _norm;
    double _prevNorm;
    double _nonLinChange;
    double _steadyChange;
    BOOL _valid;
    BOOL _output;
    BOOL _valuesChanged;
    BOOL _componentVariable;
    BOOL _secondary;
    BOOL _componentSecondaryVariable;
    
    variableArraysContainer * __nonnull _containers;
}

@property(nonatomic, strong, nullable) NSString *name;
@property(nonatomic, strong, nullable) id primaryMesh;
@property(nonatomic, strong, nullable) id solution;
@property(nonatomic, assign) int nameLength;
@property(nonatomic, assign) int dofs;
@property(nonatomic, assign) int nonLinConverged;
@property(nonatomic, assign) int steadyConverged;
@property(nonatomic, assign) int nonLinIter;
@property(nonatomic, assign) int type;
@property(nonatomic, assign) double norm;
@property(nonatomic, assign) double prevNorm;
@property(nonatomic, assign) double nonLinChange;
@property(nonatomic, assign) double steadyChange;
@property(nonatomic, assign, getter = isValid) BOOL valid;
@property(nonatomic, assign, getter = isOutput) BOOL output;
@property(nonatomic, assign, getter = isValuesChanged) BOOL valuesChanged;
@property(nonatomic, assign, getter = isComponentVariable) BOOL componentVariable;
@property(nonatomic, assign, getter =  isSecondary) BOOL secondary;
@property(nonatomic, assign, getter = isComponentSecondaryVariable) BOOL componentSecondaryVariable;

-(void)deallocation;
-(NSString * __nullable)canonicalizeName;
-(variableArraysContainer * __nonnull)getContainers;
          
@end
