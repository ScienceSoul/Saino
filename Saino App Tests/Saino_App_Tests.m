//===----------------------------------------------------------------------===//
//  Saino_App_Tests.m
//  Saino App Tests
//
//  Created by Seddik hakime on 09/06/2015.
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

#import <XCTest/XCTest.h>
#import <Cocoa/Cocoa.h>
#import <AppKit/AppKit.h>
#import <SainoCore/FEMTest.h>
#import "SainoViewControler.h"

@interface Saino_App_Tests : XCTestCase {
    SainoViewControler  * _Nonnull _saino_view_controller;
    NSView              * _Nonnull _saino_view;
    
    double _targetEps;
}

@property(nonatomic, nonnull) FEMTest *testApp;

@end

@implementation Saino_App_Tests

-(void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
    self.testApp = [FEMTest sharedTest];
    XCTAssertNotNil(self.testApp, @"Saino_App_Tests Could not create a FEMTest object.");
    
    _saino_view_controller = (SainoViewControler*)[[NSApplication sharedApplication] delegate];
    _saino_view = _saino_view_controller.view;
    
    _targetEps = 1.0e-5;
    
    NSString *dirName;
    if ([@"seddikhakime" isEqualToString:NSUserName()] || [@"hakimeseddik" isEqualToString:NSUserName()]) {
        NSString *user = [@"/Users" stringByAppendingPathComponent:NSUserName()];
        dirName = [user stringByAppendingPathComponent:@"Documents/Saino/Tests"];
    } else if ([@"_xcsbuildd" isEqualToString:NSUserName()]) {
        NSProcessInfo *processIngo = [[NSProcessInfo alloc] init];
        NSDictionary *env = [processIngo environment];
        dirName = [env[@"XCS_SOURCE_DIR"] stringByAppendingPathComponent:@"Saino/Tests"];
    } else {
        fprintf(stderr, "User name not supported for testing.\n");
    }
    self.testApp.path = [NSMutableString stringWithString:dirName];
}

-(void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

-(void)testHeateq {
    int success;
    double targetNorm = 0.768016492512e-01;
    
    [_saino_view_controller press:[_saino_view viewWithTag: 1]];  // Heateq
    
    if (self.testApp.heatq_allDone == YES) {
        if ( self.testApp.norm != -1.0 ) {
            if ( _targetEps < 0.0) {
                if ( targetNorm < self.testApp.norm )
                    success = 0;
                else
                    success = -1;
            } else  if ( 2.0 * fabs(self.testApp.norm - targetNorm) / (self.testApp.norm + targetNorm) < _targetEps ) success = 0;
            else success = -1;
        } else {
            success = 0;
        }
        if (success != 0) {
            fprintf(stdout, "Computed norm: %e.\n", self.testApp.norm);
        } else {
            fprintf(stdout, "heateq:   [Passed].\n");
        }
        XCTAssertTrue(success == 0, @"heateq:   [FAILED].\n");
    } else {
        XCTAssertTrue(self.testApp.heatq_allDone == YES, @"heateq: not reaching end of simulation:    [LOOK AT ERROR LOGS].\n");
    }
}

-(void)testStep_Stokes {
    int success;
    double targetNorm = 0.75227162e-01;
    _targetEps = 1.0e-2;
    
    [_saino_view_controller press:[_saino_view viewWithTag: 2]];  // Step_stokes

    if (self.testApp.step_stokes_allDone == YES) {
        if ( self.testApp.norm != -1.0 ) {
            if ( _targetEps < 0.0 ) {
                if ( targetNorm < self.testApp.norm )
                    success = 0;
                else
                    success = -1;
            } else  if ( 2.0 * fabs(self.testApp.norm - targetNorm) / (self.testApp.norm + targetNorm) < _targetEps ) success = 0;
            else success = -1;
        } else {
            success = 0;
        }
        if (success != 0) {
            fprintf(stdout, "Computed norm: %e.\n", self.testApp.norm);
        } else {
            fprintf(stdout, "Step_Stokes:   [Passed].\n");
        }
        XCTAssertTrue(success == 0, @"Step_Stokes:   [FAILED].\n");
    } else {
        XCTAssertTrue(self.testApp.step_stokes_allDone == YES, @"Step_Stokes: not reaching end of simulation:    [LOOK AT ERROR LOGS].\n");
    }
}

-(void)testStep_NatConvection {
    int success;
    double targetNorm = 0.275036673-03;
    
    [_saino_view_controller press:[_saino_view viewWithTag: 3]];  // Natural convection
    
    if (self.testApp.natural_convection_allDone == YES) {
        if ( self.testApp.norm != -1.0 ) {
            if ( _targetEps < 0.0 ) {
                if ( targetNorm < self.testApp.norm )
                    success = 0;
                else
                    success = -1;
            } else  if ( 2.0 * fabs(self.testApp.norm - targetNorm) / (self.testApp.norm + targetNorm) < _targetEps ) success = 0;
            else success = -1;
        } else {
            success = 0;
        }
        if (success != 0) {
            fprintf(stdout, "Computed norm: %e.\n", self.testApp.norm);
        } else {
            fprintf(stdout, "Natural convection:   [Passed].\n");
        }
        XCTAssertTrue(success == 0, @"Natural convection:   [FAILED].\n");
    } else {
        XCTAssertTrue(self.testApp.natural_convection_allDone == YES, @"Natural convection: not reaching end of simulation:    [LOOK AT ERROR LOGS].\n");
    }
}

-(void)testISMIP_HOM_A010 {
    int success;
    double targetNorm = 9.8237320560;
    
    [_saino_view_controller press:[_saino_view viewWithTag: 4]]; // ISMIP-HOM A010
    
    if (self.testApp.ismip_hom_A010_allDone == YES) {
        if (self.testApp.norm != -1.0) {
            if (_targetEps < 0.0) {
                if (targetNorm < self.testApp.norm) {
                    success = 0;
                } else
                    success = -1;
            } else if (2.0 * fabs(self.testApp.norm - targetNorm) / (self.testApp.norm + targetNorm) < _targetEps) success = 0;
            else success = -1;
        } else {
            success = 0;
        }
        if (success != 0) {
            fprintf(stdout, "Computed norm: %e.\n", self.testApp.norm);
        } else {
            fprintf(stdout, "ISMIP-HOM A010:   [Passed].\n");
        }
        XCTAssertTrue(success == 0, @"ISMIP-HOM A010:   [FAILED].\n");
    } else {
        XCTAssertTrue(self.testApp.ismip_hom_A010_allDone == YES, @"ISMIP-HOM A010: not reaching end of simulation:    [LOOK AT ERROR LOGS].\n");
    }
}

-(void)testISMIP_HOM_B010 {
    int success;
    double targetNorm = 11.5756923657;
    
    _targetEps = 1.0e-4;
    
    [_saino_view_controller press:[_saino_view viewWithTag: 5]]; // ISMIP-HOM B010
    
    if (self.testApp.ismip_hom_B010_allDone == YES) {
        if (self.testApp.norm != -1.0) {
            if (_targetEps < 0.0) {
                if (targetNorm < self.testApp.norm) {
                    success = 0;
                } else
                    success = -1;
            } else if (2.0 * fabs(self.testApp.norm - targetNorm) / (self.testApp.norm + targetNorm) < _targetEps) success = 0;
            else success = -1;
        } else {
            success = 0;
        }
        if (success != 0) {
            fprintf(stdout, "Computed norm: %e.\n", self.testApp.norm);
        } else {
            fprintf(stdout, "ISMIP-HOM B010:   [Passed].\n");
        }
        XCTAssertTrue(success == 0, @"ISMIP-HOM B010:   [FAILED].\n");
    } else {
        XCTAssertTrue(self.testApp.ismip_hom_B010_allDone == YES, @"ISMIP-HOM B010: not reaching end of simulation:    [LOOK AT ERROR LOGS].\n");
    }


}

-(void)testISMIP_HOM_C010 {
    int success;
    double targetNorm = 9.1726234962;
    
    [_saino_view_controller press:[_saino_view viewWithTag: 6]]; // ISMIP-HOM C010
    
    if (self.testApp.ismip_hom_C010_allDone == YES) {
        if (self.testApp.norm != -1.0) {
            if (_targetEps < 0.0) {
                if (targetNorm < self.testApp.norm) {
                    success = 0;
                } else
                    success = -1;
            } else if (2.0 * fabs(self.testApp.norm - targetNorm) / (self.testApp.norm + targetNorm) < _targetEps) success = 0;
            else success = -1;
        } else {
            success = 0;
        }
        if (success != 0) {
            fprintf(stdout, "Computed norm: %e.\n", self.testApp.norm);
        } else {
            fprintf(stdout, "ISMIP-HOM C010:   [Passed].\n");
        }
        XCTAssertTrue(success == 0, @"ISMIP-HOM C010:   [FAILED].\n");
    } else {
        XCTAssertTrue(self.testApp.ismip_hom_C010_allDone == YES, @"ISMIP-HOM C010: not reaching end of simulation:    [LOOK AT ERROR LOGS].\n");
    }
}

-(void)testISMIP_HOM_A010_GPU_Color {
    int success;
    double targetNorm = 9.8237320560;
    
    [_saino_view_controller press:[_saino_view viewWithTag: 7]]; // ISMIP-HOM A010 GPU Color
    
    if (self.testApp.ismip_hom_A010_gpu_color_allDone == YES) {
        if (self.testApp.norm != -1.0) {
            if (_targetEps < 0.0) {
                if (targetNorm < self.testApp.norm) {
                    success = 0;
                } else
                    success = -1;
            } else if (2.0 * fabs(self.testApp.norm - targetNorm) / (self.testApp.norm + targetNorm) < _targetEps) success = 0;
            else success = -1;
        } else {
            success = 0;
        }
        if (success != 0) {
            fprintf(stdout, "Computed norm: %e.\n", self.testApp.norm);
        } else {
            fprintf(stdout, "ISMIP-HOM A010 GPU Color:   [Passed].\n");
        }
        XCTAssertTrue(success == 0, @"ISMIP-HOM A010 GPU Color:   [FAILED].\n");
    } else {
        XCTAssertTrue(self.testApp.ismip_hom_A010_gpu_color_allDone == YES, @"ISMIP-HOM A010 GPU Color: not reaching end of simulation:    [LOOK AT ERROR LOGS].\n");
    }
}

-(void)testISMIP_HOM_A010_GPU_Nonzeros {
    int success;
    double targetNorm = 9.8237320560;
    
    [_saino_view_controller press:[_saino_view viewWithTag: 8]]; // ISMIP-HOM A010 GPU Nonzeros
    
    if (self.testApp.ismip_hom_A010_gpu_nonzeros_allDone == YES) {
        if (self.testApp.norm != -1.0) {
            if (_targetEps < 0.0) {
                if (targetNorm < self.testApp.norm) {
                    success = 0;
                } else
                    success = -1;
            } else if (2.0 * fabs(self.testApp.norm - targetNorm) / (self.testApp.norm + targetNorm) < _targetEps) success = 0;
            else success = -1;
        } else {
            success = 0;
        }
        if (success != 0) {
            fprintf(stdout, "Computed norm: %e.\n", self.testApp.norm);
        } else {
            fprintf(stdout, "ISMIP-HOM A010 GPU Nonzeros:   [Passed].\n");
        }
        XCTAssertTrue(success == 0, @"ISMIP-HOM A010 GPU Nonzeros:   [FAILED].\n");
    } else {
        XCTAssertTrue(self.testApp.ismip_hom_A010_gpu_nonzeros_allDone == YES, @"ISMIP-HOM A010 GPU Nonzeros: not reaching end of simulation:    [LOOK AT ERROR LOGS].\n");
    }
}

@end
