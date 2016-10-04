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
-(void)setUpISMIP_HOM_A010_GPU_Coloring:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_dense1:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_dense2:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_Coloring_dense3:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros_dense1:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros_dense2:(id __nonnull)model;
-(void)setUpISMIP_HOM_A010_GPU_NonZeros_dense3:(id __nonnull)model;

@end
