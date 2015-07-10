//
//  FEMEquation.m
//  Saino
//
//  Created by Seddik hakime on 13/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMEquation.h"

@implementation FEMEquation

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
