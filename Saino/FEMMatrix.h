//
//  FEMMatrix.h
//  Saino
//
//  Created by Hakime Seddik on 23/07/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"

@interface FEMMatrix : NSObject {
    
    FEMMatrix *_child;
    FEMMatrix *_parent;
    FEMMatrix *_constraint;
    FEMMatrix *_ematrix;
    int _numberOfRows;
    int _subband;          // >= 1
    int _format;
    int _solveCount;
    int _comm;
    int _umfPackNumeric;
    double _rhsScaling;
    BOOL _ordered;
    BOOL _lumped;
    BOOL _symmetric;
    BOOL _complexMatrix;
    BOOL _dgMatrix;
    
    matrixArraysContainer *_containers;
}

@property(nonatomic, strong) FEMMatrix *child;
@property(nonatomic, strong) FEMMatrix *parent;
@property(nonatomic, strong) FEMMatrix *constraint;
@property(nonatomic, strong) FEMMatrix *ematrix;
@property(nonatomic, assign) int numberOfRows;
@property(nonatomic, assign) int subband;
@property(nonatomic, assign) int format;
@property(nonatomic, assign) int solveCount;
@property(nonatomic, assign) int comm;
@property(nonatomic, assign) int umfPackNumeric;
@property(nonatomic, assign) double rhsScaling;
@property(nonatomic, assign, getter = isOrdered) BOOL ordered;
@property(nonatomic, assign, getter = isLumped) BOOL lumped;
@property(nonatomic, assign, getter = isSymmetric) BOOL symmetric;
@property(nonatomic, assign, getter = isComplexMatrix) BOOL complexMatrix;
@property(nonatomic, assign, getter = isDgmatrix) BOOL dgMatrix;

-(void)deallocation;
-(matrixArraysContainer *)getContainers;

@end
