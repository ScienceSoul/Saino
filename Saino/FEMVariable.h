//
//  FEMVariable.h
//  Saino
//
//  Created by Seddik hakime on 27/07/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

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
    BOOL _secondary;
    
    variableArraysContainer *_containers;
}

@property(nonatomic, strong) NSString *name;
@property(nonatomic, strong) id primaryMesh;
@property(nonatomic, strong) id solution;
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
@property(nonatomic, assign, getter =  isSecondary) BOOL secondary;

-(void)deallocation;
-(variableArraysContainer *)getContainers;
          
@end
