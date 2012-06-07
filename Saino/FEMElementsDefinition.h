//
//  FEMElementsDefinition.h
//  Saino
//
//  Created by Hakime Seddik on 05/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMElementsDefinition : NSObject 

// This class is used by FEMElementDescription for the elements definition.

@property (nonatomic, strong) NSNumber *dimension;
@property (nonatomic, strong) NSString *topology;
@property (nonatomic, strong) NSNumber *code;
@property (nonatomic, strong) NSNumber *nodes;
@property (nonatomic, strong) NSArray *nodeU;
@property (nonatomic, strong) NSArray *nodeV;
@property (nonatomic, strong) NSArray *nodeW;
@property (nonatomic, strong) NSArray *basis;
@property (nonatomic, strong) NSArray *gaussPoints;
@property (nonatomic, strong) NSNumber *stabilization;

@end
