//
//  FEMBoundary.m
//  Saino
//
//  Created by Seddik hakime on 21/12/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMBoundary.h"

@implementation FEMBoundary

@synthesize valuesList = _valuesList;

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _valuesList = [[NSMutableArray alloc] init];
    }
    
    return self;
}

@end
