//
//  FEMListUtilitiesTests.m
//  Saino
//
//  Created by Seddik hakime on 11/05/2015.
//  Copyright (c) 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <XCTest/XCTest.h>

#import "FEMListUtilities.h"
#import "FEMEquation.h"
#import "memory.h"

@interface FEMListUtilitiesTests : XCTestCase {
    
}
@property (nonatomic) FEMListUtilities *listUtilities;
@property (nonatomic) FEMEquation *equation;

@end

@implementation FEMListUtilitiesTests

- (void)setUp
{
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
    self.listUtilities = [[FEMListUtilities alloc] init];
    XCTAssertNotNil(self.listUtilities, @"Could not create a list utilities object.");
    self.equation = [[FEMEquation alloc] init];
    XCTAssertNotNil(self.equation, @"Could not create an equation object.");
}

- (void)tearDown
{
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    for (FEMValueList *valueList in self.equation.valuesList) {
        [valueList deallocation];
    }
    [super tearDown];
}

- (void)testListGetString
{
    BOOL found;
    [self.listUtilities addStringInClassList:self.equation theVariable:@"test string" withValue:@"string"];
    NSString *stringTest = [self.listUtilities listGetString:nil inArray:self.equation.valuesList forVariable:@"test string" info:&found];
    XCTAssertTrue([stringTest isEqualToString:@"string"] == YES, @"Method listGetString:inArray:forVariable:info: failed to retrieve correct string.");
    //XCTFail(@"No implementation for \"%s\"", __PRETTY_FUNCTION__);
}

- (void)testListGetLogical
{
    BOOL found;
    [self.listUtilities addLogicalInClassList:self.equation theVariable:@"test logical" withValue:YES];
    BOOL logicalTest = [self.listUtilities listGetLogical:nil inArray:self.equation.valuesList forVariable:@"test logical" info:&found];
    XCTAssertTrue(logicalTest == YES, @"Method listUtilities listGetLogical:inArray:forVariable:info: failed to retrieve correct value.");
}

- (void)testListGetInteger
{
    BOOL found;
    int value = 1;
    [self.listUtilities addIntegerInClassList:self.equation theVariable:@"test integer" withValue:&value orUsingBlock:nil];
    int integerTest = [self.listUtilities listGetInteger:nil inArray:self.equation.valuesList forVariable:@"test integer" info:&found minValue:NULL maxValue:NULL];
    XCTAssertTrue(integerTest == 1, @"Method listGetInteger:inArray:forVariable:info:minValue:maxValue: failed to retrieve correct value.");
}

- (void)testListGetConstReal
{
    BOOL found;
    double value = 1.0;
    [self.listUtilities addConstRealInClassList:self.equation theVariable:@"test real" withValue:&value orUsingBlock:nil string:nil];
    int realTest = [self.listUtilities listGetConstReal:nil inArray:self.equation.valuesList forVariable:@"test real" info:&found minValue:NULL maxValue:NULL];
    XCTAssertTrue(realTest == 1.0, @"Method listGetConstReal:inArray:forVariable:info:minValue:maxValue: failed to retrieve correct value.");
}

- (void)testListGetIntegerArray
{
    BOOL found;
    listBuffer result = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    int *vector = intvec(0, 5);
    vector[0] = 1; vector[1] = 2; vector[2] = 3; vector[3] = 4; vector[4] = 5; vector[5] = 6;
    [self.listUtilities addIntegerArrayInClassList:self.equation theVariable:@"test integer array" withValues:vector size:6 orUsingBlock:nil];
    free_ivector(vector, 0, 5);
    
    found = [self.listUtilities listGetIntegerArray:nil inArray:self.equation.valuesList forVariable:@"test integer array" buffer:&result];
    XCTAssertTrue(result.ivector != NULL, @"Result vector in method listGetIntegerArray:inArray:forVariable:buffer: is null. It should have been allocated.");
    XCTAssertTrue((result.ivector[0] == 1 && result.ivector[1] == 2 && result.ivector[2] == 3 && result.ivector[3] == 4 && result.ivector[4] == 5 && result.ivector[5] == 6), @"Method listGetIntegerArray:inArray:forVariable:buffer: failed to retrieve correct value");
    free_ivector(result.ivector, 0, result.m-1);
}

- (void)testListGetConstRealArray
{
    BOOL found;
    listBuffer result = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    double **matrix = doublematrix(0, 2, 0, 2);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            matrix[i][j] = 1.0;
        }
    }
    [self.listUtilities addConstRealArrayInClassList:self.equation theVariable:@"test const real array" withValues:matrix size1:3 size2:3 orUsingBlock:nil string:nil];
    free_dmatrix(matrix, 0, 2, 0, 2);
    
    found = [self.listUtilities listGetConstRealArray:nil inArray:self.equation.valuesList forVariable:@"test const real array" buffer:&result];
    XCTAssertTrue(result.matrix != NULL, @"Result matrix in method listGetConstRealArray:inArray:forVariable:buffer: is null. It should have been allocated.");
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            XCTAssertTrue(result.matrix[i][j] == 1.0, @"Method listGetConstRealArray:inArray:forVariable:buffer: failed to retrieve correct value at matrix indexes: i: %d, j: %d\n.", i, j);
        }
    }
    free_dmatrix(result.matrix, 0, result.m-1, 0, result.n);
}

- (void)testListGetRealArray // Test only for type=LIST_TYPE_CONSTANT_TENSOR
{
    BOOL found;
    listBuffer result = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    double **matrix = doublematrix(0, 2, 0, 2);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            matrix[i][j] = 1.0;
        }
    }
    [self.listUtilities addConstRealArrayInClassList:self.equation theVariable:@"test real array" withValues:matrix size1:3 size2:3 orUsingBlock:nil string:nil];
    free_dmatrix(matrix, 0, 2, 0, 2);
    
    found = [self.listUtilities listGetRealArray:nil inArray:self.equation.valuesList forVariable:@"test real array" numberOfNodes:3 indexes:NULL buffer:&result];
    XCTAssertTrue(result.tensor != NULL, @"Result tensor in method listGetRealArray:inArray:forVariable:numberOfNodes:indexes:buffer: is null. It should have been allocated.");
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            for (int k=0; k<3; k++) {
                XCTAssertTrue(result.tensor[i][j][k] == 1.0, @"Method listGetRealArray:inArray:forVariable:numberOfNodes:indexes:buffer: failed to retrieve correct value at matrix indexes: i: %d, j: %d, k: %d\n.", i, j, k);
            }
        }
    }
    free_d3tensor(result.tensor, 0, result.m-1, 0, result.n-1, 0, result.p-1);
}

@end

