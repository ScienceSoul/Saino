//
//  FEMValueList.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"

typedef double (^inputBlock) (double *variablesValues);

@interface FEMValueList : NSObject {
    
    int _model;
    int _type;
    int _nameLength;
    int _numberOfDependencies;
    BOOL _lValue;
    NSString *_name;
    NSString *_cValue;
    NSArray *_dependencies;       // Store the names of the variable to whhich there is/are a
                                  // dependency(ies) in the exection of the block
    
    valueListArraysContainer *_containers;
    inputBlock _block;
}

@property(nonatomic, assign) int model;
@property(nonatomic, assign) int type;
@property(nonatomic, assign) int nameLength;
@property(nonatomic, assign) int numberOfDependencies;
@property(nonatomic, assign, getter = isLvalue) BOOL lValue;
@property(nonatomic, strong) NSString *name;
@property(nonatomic, strong) NSString *cValue;
@property(nonatomic, strong) NSArray *dependencies;
@property(copy) inputBlock block;

-(void)deallocation;
-(valueListArraysContainer *)getContainers;

@end
