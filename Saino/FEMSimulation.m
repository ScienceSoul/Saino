//
//  FEMSimulation.m
//  Saino
//
//  Created by Hakime Seddik on 26/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMSimulation.h"

@implementation FEMSimulation

@synthesize tag = _tag;
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
