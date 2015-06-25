//
//  FEMTest.h
//  Saino
//
//  Created by Seddik hakime on 09/06/2015.
//  Copyright (c) 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMTest : NSObject

@property(nonatomic, assign) BOOL do_heatq;
@property(nonatomic, assign) BOOL do_step_stokes;
@property(nonatomic, assign) BOOL heatq_allDone;
@property(nonatomic, assign) BOOL step_stokes_allDone;
@property(nonatomic, assign) double norm;

// This class method retuns an instance (singleton) of FEMCore
+(id)sharedTest;
+(void)selfDestruct;

-(void)reset;

@end
