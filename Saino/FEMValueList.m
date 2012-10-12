//
//  FEMValueList.m
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMValueList.h"

#import "Constructors.h"
#import "memory.h"

@implementation FEMValueList

@synthesize model = _model;
@synthesize type = _type;
@synthesize nameLength = _nameLength;
@synthesize depNameLength = _depNameLength;
@synthesize method = _method;
@synthesize lValue = _lValue;
@synthesize name = _name;
@synthesize dependName = _dependName;
@synthesize cValue = _cValue;

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
        _containers = (valueListArraysContainer*)malloc(sizeof(valueListArraysContainer) * 1 );
    }
    
    return self;
}

-(void)deallocation {
    free(_containers);
}

-(valueListArraysContainer*)getContainers {
    
    return _containers;
}

@end
