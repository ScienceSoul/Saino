//
//  FEMValueList.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"

@interface FEMValueList : NSObject {
    
    int _model;
    int _type;
    int _nameLength;
    int _depNameLength;
    BOOL _method;
    BOOL _lValue;
    NSString *_name;
    NSString *_dependName;
    NSString *_cValue;
    
    valueListArraysContainer *_containers;
}

@property(nonatomic, assign) int model;
@property(nonatomic, assign) int type;
@property(nonatomic, assign) int nameLength;
@property(nonatomic, assign) int depNameLength;
@property(nonatomic, assign, getter = isMethod) BOOL method;
@property(nonatomic, assign, getter = isLvalue) BOOL lValue;
@property(nonatomic, strong) NSString *name;
@property(nonatomic, strong) NSString *dependName;
@property(nonatomic, strong) NSString *cValue;

-(void)deallocation;
-(valueListArraysContainer *)getContainers;

@end
