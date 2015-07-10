//
//  FEMBoundaryCondition.m
//  Saino
//
//  Created by Hakime Seddik on 27/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMBoundaryCondition.h"

@implementation FEMBoundaryCondition

@synthesize tag = _tag;
@synthesize valuesList = _valuesList;
@synthesize pMatrix = _pMatrix;

- (id)init
{
    self = [super init];
    if (self) {
        _valuesList = [[NSMutableArray alloc] init];
    }
    
    return self;
}


@end
