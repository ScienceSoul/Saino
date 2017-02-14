//===----------------------------------------------------------------------===//
//  SainoViewControler.m
//  Saino
//
//  Created by Seddik hakime on 24/06/2015.
//  Copyright (c) 2015 ScienceSoul. All rights reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

#import "SainoViewControler.h"
#import <SainoCore/FEMTest.h>

@interface SainoViewControler ()

@end

@implementation SainoViewControler

@synthesize displayField;

- init {
    if (self = [super init]) {
        
    }
    return self;
}

-(IBAction)press:(id __nullable)sender {
    if ([sender isKindOfClass: [NSButton class]]) {
        int initialize = 0;
        FEMTest *test = [FEMTest sharedTest];
        [test reset];
        
        if ([[sender title] isEqualToString:@"heateq"]) {
            test.do_heatq = YES;
        } else if ([[sender title] isEqualToString:@"StepStokes"]) {
            test.do_step_stokes = YES;
        } else if ([[sender title] isEqualToString:@"NatConvection"]) {
            test.do_natural_convection = YES;
        } else if ([[sender title] isEqualToString:@"ISMIP-A010"]) {
            test.do_ismip_hom_A010 = YES;
        } else if ([[sender title] isEqualToString:@"ISMIP-B010"]) {
            test.do_ismip_hom_B010 = YES;
        } else if ([[sender title] isEqualToString:@"ISMIP-A010 GPU Color"]) {
            test.do_ismip_hom_A010_gpu_coloring = YES;
        } else if ([[sender title] isEqualToString:@"ISMIP-A010 GPU NZs"]) {
            test.do_ismip_hom_A010_gpu_nonzeros = YES;
        } else if ([[sender title] isEqualToString:@"ISMIP-C010"]) {
            test.do_ismip_hom_C010 = YES;
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
