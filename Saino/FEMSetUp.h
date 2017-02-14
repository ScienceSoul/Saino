//===----------------------------------------------------------------------===//
//  FEMSetUp.h
//  Saino
//
//  Created by Hakime Seddik on 08/07/2016.
//  Copyright © 2016 ScienceSoul. All rights reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

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
