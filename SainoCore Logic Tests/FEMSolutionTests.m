//
//  FEMSolutionTests.m
//  Saino
//
//  Created by Seddik hakime on 13/10/2015.
//  Copyright © 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <XCTest/XCTest.h>

#import <SainoCore/SainoCore.h>

@interface FEMSolutionTests : XCTestCase

@property(nonatomic) FEMSolution *solution;
@property(nonatomic) FEMUtilities *utilities;
@property(nonatomic) id<SainoSolutionsComputer> instance;

@end

@implementation FEMSolutionTests

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
    self.solution = [[FEMSolution alloc] init];
    XCTAssertNotNil(self.solution, @"FEMSolutionTests Could not create a FEMSolution object.");
    self.utilities = [[FEMUtilities alloc] init];
    XCTAssertNotNil(self.utilities, @"FEMSolutionTests Could not create a FEMUtiltiies object.");
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [self.instance deallocation:self.solution];
    [super tearDown];
}

-(void)testInstantiatePrincipalClassFromPlugIn {

    BOOL validBundle;

    BOOL useAppSupportPath = NO;
    NSString *bundleName;
    // Small hack so that testing works locally with XCode and with XCodeServer CI.
    // If the user name matches, then we are testing locally, otherwise we are in
    // XCodeServer and we need to resolve the full path
    if ([@"seddikhakime" isEqualToString:NSUserName()] || [@"hakimeseddik" isEqualToString:NSUserName()]) {
        NSString *user = [@"/Users" stringByAppendingPathComponent:NSUserName()];
        bundleName = [user stringByAppendingPathComponent:@"Documents/Saino/PlugIns/FEMStructuredMeshMapper"];
    } else bundleName = @"PlugIns/FEMStructuredMeshMapper";
    
    NSBundle *bundle = [self.utilities loadBundle:bundleName useApplicationSupportPath:&useAppSupportPath];
    if (bundle != nil) {
        Class currPrincipalClass = [bundle principalClass];
        if (currPrincipalClass) {
            validBundle = [self.utilities plugInClassIsValid:currPrincipalClass];
        }
    } else {
        NSLog(@"Can't load bundle.\n");
    }
    XCTAssertTrue(validBundle, @"FEMSolutionTests: problem loading a bundle");

    BOOL instantiateSuccess = [self.solution instantiatePrincipalClassFromPlugIn:bundle];
    XCTAssertTrue(instantiateSuccess, @"FEMSolutionTests: problem when instantiating principal class from bundle");

    self.instance = self.solution.plugInPrincipalClassInstance;
}

@end
