//
//  FEMTest.h
//  Saino
//
//  Created by Seddik hakime on 09/06/2015.
//  Copyright (c) 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMTest : NSObject

@property(nonatomic, strong, nonnull) NSMutableString *path;
@property(nonatomic, assign) BOOL do_heatq;
@property(nonatomic, assign) BOOL do_step_stokes;
@property(nonatomic, assign) BOOL do_natural_convection;
@property(nonatomic, assign) BOOL do_ismip_hom_A010;
@property(nonatomic, assign) BOOL do_ismip_hom_B010;
@property(nonatomic, assign) BOOL do_ismip_hom_A010_gpu_coloring;
@property(nonatomic, assign) BOOL do_ismip_hom_A010_gpu_nonzeros;
@property(nonatomic, assign) BOOL do_ismip_hom_C010;
@property(nonatomic, assign) BOOL heatq_allDone;
@property(nonatomic, assign) BOOL step_stokes_allDone;
@property(nonatomic, assign) BOOL natural_convection_allDone;
@property(nonatomic, assign) BOOL ismip_hom_A010_allDone;
@property(nonatomic, assign) BOOL ismip_hom_B010_allDone;
@property(nonatomic, assign) BOOL ismip_hom_A010_gpu_color_allDone;
@property(nonatomic, assign) BOOL ismip_hom_A010_gpu_nonzeros_allDone;
@property(nonatomic, assign) BOOL ismip_hom_C010_allDone;
@property(nonatomic, assign) double norm;

// This class method retuns an instance (singleton) of FEMTest
+(id __nonnull)sharedTest;
+(void)selfDestruct;

-(void)reset;
-(void)setUpHeateqTest:(id __nonnull)model;
-(void)setUpStepStokesTest:(id __nonnull)model;
-(void)setUpNaturalConvectionTest:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010Test:(id __nonnull)model;
-(void)setUpISMIP_HOM_B010Test:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010Test_GPU_coloring:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010Test_GPU_nonzeros:(id __nonnull)model;
-(void)setUpISMIP_HOM_C010Test:(id __nonnull)model;

@end
