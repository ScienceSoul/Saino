//
//  FEMInitialConditions.m
//  Saino
//
//  Created by Seddik hakime on 13/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMInitialConditions.h"

@implementation FEMInitialConditions

@synthesize valuesList = _valuesList;

- (id)init
{
    self = [super init];
    if (self) {
        _valuesList = [[NSMutableArray alloc] init];
    }
    
    return self;
}

@end
