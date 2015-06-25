//
//  SainoViewControler.m
//  Saino
//
//  Created by Seddik hakime on 24/06/2015.
//  Copyright (c) 2015 Institute of Low Temperature Science. All rights reserved.
//

#import "SainoViewControler.h"

#import "FEMTest.h"

@interface SainoViewControler ()

@end

@implementation SainoViewControler

@synthesize displayField;

- init {
    if (self = [super init]) {
        
    }
    return self;
}

-(IBAction)press:(id)sender {
    if ([sender isKindOfClass: [NSButton class]]) {
        int initialize = 0;
        FEMTest *test = [FEMTest sharedTest];
        [test reset];
        
        if ([[sender title] isEqualToString:@"heateq"]) {
            test.do_heatq = YES;
        } else if ([[sender title] isEqualToString:@"StepStokes"]) {
            test.do_step_stokes = YES;
        }
        _job = [[FEMJob alloc] init];
        [_job runWithInitialize:initialize];
        [_job deallocation];
        _job = nil;
        
        NSString *displayValue = [NSString stringWithFormat:@"%e",test.norm];
        [self.displayField setStringValue: displayValue];
    }
}

- (void)viewDidLoad {
    [super viewDidLoad];
    // Do view setup here.
}

@end
