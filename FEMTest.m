//
//  FEMTest.m
//  Saino
//
//  Created by Seddik hakime on 09/06/2015.
//  Copyright (c) 2015 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMTest.h"

@implementation FEMTest

@synthesize do_heatq = _do_heatq;
@synthesize do_step_stokes = _do_step_stokes;
@synthesize heatq_allDone = _heatq_allDone;
@synthesize step_stokes_allDone = _step_stokes_allDone;
@synthesize norm = _norm;

static FEMTest *sharedTest = nil;
static dispatch_once_t onceToken;

+(id)sharedTest {
    
    dispatch_once(&onceToken, ^{
        sharedTest = [[self alloc] init];
    });
    return sharedTest;
}

+(void) selfDestruct {
    sharedTest = nil;
    onceToken = 0;
}

- (id)init
{
    self = [super init];
    if (self) {
        
    }
    
    return self;
}

-(void)reset {
    _do_heatq = NO;
    _do_step_stokes = NO;
    _heatq_allDone = NO;
    _step_stokes_allDone = NO;
    _norm = 0.0;
}

@end
