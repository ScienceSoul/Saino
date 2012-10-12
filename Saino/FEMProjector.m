//
//  FEMProjector.m
//  Saino
//
//  Created by Seddik hakime on 09/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMProjector.h"

@implementation FEMProjector

@synthesize mesh = _mesh;
@synthesize matrix = _matrix;
@synthesize tMatrix = _tMatrix;

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _mesh = [[FEMMesh alloc] init];
        _matrix = [[FEMMatrix alloc] init];
        _tMatrix = [[FEMMatrix alloc] init];
    }
    
    return self;
}



@end
