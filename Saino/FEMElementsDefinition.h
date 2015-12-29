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

@property (nonatomic, strong, nonnull) NSNumber *dimension;
@property (nonatomic, strong, nonnull) NSString *topology;
@property (nonatomic, strong, nonnull) NSNumber *code;
@property (nonatomic, strong, nonnull) NSNumber *nodes;
@property (nonatomic, strong, nonnull) NSArray *nodeU;
@property (nonatomic, strong, nonnull) NSArray *nodeV;
@property (nonatomic, strong, nonnull) NSArray *nodeW;
@property (nonatomic, strong, nonnull) NSArray *basis;
@property (nonatomic, strong, nonnull) NSArray *gaussPoints;
@property (nonatomic, strong, nonnull) NSNumber *stabilization;

@end
