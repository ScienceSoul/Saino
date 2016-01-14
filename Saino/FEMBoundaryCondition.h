//
//  FEMBoundaryCondition.h
//  Saino
//
//  Created by Hakime Seddik on 27/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
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
