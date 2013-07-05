//
//  FEMMaterial.m
//  Saino
//
//  Created by Seddik hakime on 17/05/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMaterial.h"

@implementation FEMMaterial

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
