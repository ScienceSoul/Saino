//
//  FEMJob.h
//  Saino
//
//  Created by Seddik hakime on 09/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "FEMCore.h"
#import "FEMElementDescription.h"

@interface FEMJob : NSObject {
    
    FEMCore *_core;
    FEMElementDescription *_elementDescription;
    FEMModel *_model;
    NSString *_modelName;
}

@property(nonatomic, strong, nonnull) FEMCore *core;
@property(nonatomic, strong, nonnull) FEMElementDescription *elementDescription;
@property(nonatomic, strong, nonnull) FEMModel *model;
@property(nonatomic, strong, nonnull) NSString *modelName;

-(void)runWithInitialize:(int)initialize;
-(void)deallocation;

@end
