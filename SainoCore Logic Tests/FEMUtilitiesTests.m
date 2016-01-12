//
//  FEMUtilitiesTests.m
//  Saino
//
//  Created by Seddik hakime on 13/10/2015.
//  Copyright © 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <XCTest/XCTest.h>

#import <SainoCore/SainoCore.h>

@interface FEMUtilitiesTests : XCTestCase

@property(nonatomic, nonnull) FEMUtilities *utilities;

@end

@implementation FEMUtilitiesTests

-(void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
    self.utilities = [[FEMUtilities alloc] init];
    XCTAssertNotNil(self.utilities, @"FEMUtilitiesTests Could not create a FEMUtilities object.");
}

-(void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

-(void)testAppendNameFromString {
    NSString *testString2D = @"flow solution[velocity:2 pressure:1]";
    int dofs = 3;
    for (int i=1; i<=dofs; i++) {
        NSString *varName = [self.utilities appendNameFromString:testString2D component:&i];
        NSLog(@"%@\n", varName);
        if (i == 1) {
            XCTAssertTrue([varName isEqualToString:@"velocity 1"], @"FEMUtilitiesTests: Method appendNameFromString:component: failed to retrieve correct variable name in 2D. Should give 'velocity 1'");
        } else if (i == 2) {
            XCTAssertTrue([varName isEqualToString:@"velocity 2"], @"FEMUtilitiesTests: Method appendNameFromString:component: failed to retrieve correct variable name in 2D. Should give 'velocity 2'");
        } else if (i == 3) {
            XCTAssertTrue([varName isEqualToString:@"pressure"], @"FEMUtilitiesTests: Method appendNameFromString:component: failed to retrieve correct variable name in 2D. Should give 'pressure'");
        }
    }
    
    NSString *testString3D = @"flow solution[velocity:3 pressure:1]";
    dofs = 4;
    for (int i=1; i<=dofs; i++) {
        NSString *varName = [self.utilities appendNameFromString:testString3D component:&i];
        NSLog(@"%@\n", varName);
        if (i == 1) {
            XCTAssertTrue([varName isEqualToString:@"velocity 1"], @"FEMUtilitiesTests: Method appendNameFromString:component: failed to retrieve correct variable name in 3D. Should give 'velocity 1'");
        } else if (i == 2) {
            XCTAssertTrue([varName isEqualToString:@"velocity 2"], @"FEMUtilitiesTests: Method appendNameFromString:component: failed to retrieve correct variable name in 3D. Should give 'velocity 2'");
        } else if (i == 3) {
            XCTAssertTrue([varName isEqualToString:@"velocity 3"], @"FEMUtilitiesTests: Method appendNameFromString:component: failed to retrieve correct variable name in 3D. Should give 'velocity 3'");
        } else if (i == 4) {
            XCTAssertTrue([varName isEqualToString:@"pressure"], @"FEMUtilitiesTests: Method appendNameFromString:component: failed to retrieve correct variable name in 3D. Should give 'pressure'");
        }
    }
}

-(void)testAddVariableTo {
    
    NSMutableArray *storage = [[NSMutableArray alloc] init];
    variableArraysContainer *varContainer = (variableArraysContainer*)malloc(sizeof(variableArraysContainer));
    
    NSString *name1 = @"dummy";
    [self.utilities addVariableTo:storage mesh:nil solution:nil name:name1 dofs:1 container:varContainer component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
    
    BOOL found = NO;
    for (FEMVariable *variable in storage) {
        if ([variable.name isEqualToString:@"dummy"] == YES) {
            found = YES;
            break;
        }
    }
    XCTAssertTrue(found, @"FEMUtilitiesTests: addVariableTo: failed to properlly add the variable 'dummy'.");
}

@end
