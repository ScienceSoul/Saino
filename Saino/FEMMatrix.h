//
//  FEMMatrix.h
//  Saino
//
//  Created by Hakime Seddik on 23/07/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"

@interface FEMMatrix : NSObject <NSCopying> {
    
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
    int _projectorBC;
    int _projectorType;
    double _rhsScaling;
    BOOL _cholesky;
    BOOL _ordered;
    BOOL _lumped;
    BOOL _symmetric;
    BOOL _complexMatrix;
    BOOL _dgMatrix;
    
    matrixArraysContainer * __nonnull _containers;
}

@property(nonatomic, strong, nullable) FEMMatrix *child;
@property(nonatomic, strong, nullable) FEMMatrix *parent;
@property(nonatomic, strong, nullable) FEMMatrix *constraint;
@property(nonatomic, strong, nullable) FEMMatrix *ematrix;
@property(nonatomic, assign) int numberOfRows;
@property(nonatomic, assign) int subband;
@property(nonatomic, assign) int format;
@property(nonatomic, assign) int solveCount;
@property(nonatomic, assign) int comm;
@property(nonatomic, assign) int umfPackNumeric;
@property(nonatomic, assign) int projectorBC;
@property(nonatomic, assign) int projectorType;
@property(nonatomic, assign) double rhsScaling;
@property(nonatomic, assign, getter = isCholesky) BOOL cholesky;
@property(nonatomic, assign, getter = isOrdered) BOOL ordered;
@property(nonatomic, assign, getter = isLumped) BOOL lumped;
@property(nonatomic, assign, getter = isSymmetric) BOOL symmetric;
@property(nonatomic, assign, getter = isComplexMatrix) BOOL complexMatrix;
@property(nonatomic, assign, getter = isDgmatrix) BOOL dgMatrix;

-(void)deallocation;
-(matrixArraysContainer * __nonnull)getContainers;

@end
