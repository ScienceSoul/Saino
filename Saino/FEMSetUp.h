//
//  FEMSetUp.h
//  Saino
//
//  Created by Hakime Seddik on 08/07/2016.
//  Copyright © 2016 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMSetUp : NSObject

@property(nonatomic, strong, nonnull) NSMutableString *path;

-(void)setUpISMIP_HOM_A010:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_20x20:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_40x40:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_70x70:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_100x100:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_120x120:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_150x150:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_170x170:(id __nonnull)model;


-(void)setUpISMIP_HOM_A010_GPU_NonZeros_20x20:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros_40x40:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros_70x70:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros_100x100:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros_120x120:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros_150x150:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros_170x170:(id __nonnull)model;

-(void)setUpISMIP_HOM_B010:(id __nonnull)model;
-(void)setUpISMIP_HOM_C010:(id __nonnull)model;

@end
