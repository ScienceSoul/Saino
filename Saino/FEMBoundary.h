//
//  FEMBoundary.h
//  Saino
//
//  Created by Seddik hakime on 21/12/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//
//------------------------------------------------------------------------------
//
// Boundary to boundary condition mapping
//
//------------------------------------------------------------------------------

#import <Foundation/Foundation.h>

#import "FEMValueList.h"

@interface FEMBoundary : NSObject {
    
    NSMutableArray *_valuesList;
}

@property(nonatomic, strong, nonnull) NSMutableArray <FEMValueList *> *valuesList;

@end
