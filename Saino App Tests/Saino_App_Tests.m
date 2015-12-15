//
//  Saino_App_Tests.m
//  Saino App Tests
//
//  Created by Seddik hakime on 09/06/2015.
//  Copyright (c) 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <XCTest/XCTest.h>
#import <Cocoa/Cocoa.h>
#import <AppKit/AppKit.h>

#import <SainoCore/FEMTest.h>
#import "SainoViewControler.h"

@interface Saino_App_Tests : XCTestCase {
    SainoViewControler  *_saino_view_controller;
    NSView              *_saino_view;
    
    double _targetEps;
}

@property(nonatomic) FEMTest *testApp;

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
        NSLog(@"User name not supported for testing.\n");
    }
    NSLog(@"user name: %@\n", NSUserName());
    self.testApp.path = [NSMutableString stringWithString:dirName];
}

-(void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

-(void)testHeateq {
    int success;
    double targetNorm = 0.768016492512e-01;
    
    NSLog(@"I am here*********************\n");
    NSLog(@"path: %@\n", self.testApp.path);
    
    [_saino_view_controller press:[_saino_view viewWithTag: 1]];  // Heateq
    
    if (self.testApp.heatq_allDone == YES) {
        if ( self.testApp.norm != -1 ) {
            if ( _targetEps < 0 ) {
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
            NSLog(@"Computed norm: %e.\n", self.testApp.norm);
        } else {
            NSLog(@"heateq:   [Passed].\n");
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
        if ( self.testApp.norm != -1 ) {
            if ( _targetEps < 0 ) {
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
            NSLog(@"Computed norm: %e.\n", self.testApp.norm);
        } else {
            NSLog(@"Step_Stokes:   [Passed].\n");
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
        if ( self.testApp.norm != -1 ) {
            if ( _targetEps < 0 ) {
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
            NSLog(@"Computed norm: %e.\n", self.testApp.norm);
        } else {
            NSLog(@"Natural convection:   [Passed].\n");
        }
        XCTAssertTrue(success == 0, @"Natural convection:   [FAILED].\n");
    } else {
        XCTAssertTrue(self.testApp.natural_convection_allDone == YES, @"Natural convection: not reaching end of simulation:    [LOOK AT ERROR LOGS].\n");
    }
}

@end
