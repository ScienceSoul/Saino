//===----------------------------------------------------------------------===//
//  FEMSetUp.m
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
#import "FEMSetUp.h"
#import "FEMModel.h"
#import "FEMListUtilities.h"
#import "FEMSolution.h"
#import "memory.h"

@implementation FEMSetUp

@synthesize path = _path;

-(void)setUpISMIP_HOM_A010:(id __nonnull)model {
    
    FEMModel *mod = (FEMModel *)model;
    
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    
    double L = 80.0e3;
    double slope = 0.5 * M_PI / 180.0;
    
    double yearinsec = 365.25 * 24.0 * 60.0 * 60.0;
    double rhoi = 900.0 / (1.0e6 * pow(yearinsec, 2.0));
    double gravity = -9.81 * pow(yearinsec, 2.0);
    double n = 3.0;
    double eta = pow((2.0 * 100.0), (-1.0/n));
    
    self.path = [NSMutableString stringWithString:@"/Users/hakimeseddik/Desktop/ISMIP-HOM/A010/Saino"];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"coordinate system" withValue:@"cartesian 3d"];
    
    int *vector = intvec(0, 2);
    vector[0] = 1; vector[1] = 2; vector[2] = 3;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"coordinate mapping" withValues:vector size:3 orUsingBlock:nil];
    free_ivector(vector, 0, 2);
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"simulaton type" withValue:@"steady state"];
    
    int ivalue = 16;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"extruded mesh levels" withValue:&ivalue orUsingBlock:nil];
    
    int *intervals = intvec(0, 0);
    intervals[0] = 1;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"output intervals" withValues:intervals size:1 orUsingBlock:nil];
    free_ivector(intervals, 0, 0);
    
    ivalue = 1;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"steady state maximum iterations" withValue:&ivalue orUsingBlock:nil];
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"steady state minimum iterations" withValue:&ivalue orUsingBlock:nil];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"output file" withValue:@"ismip_A010_.result"];
    [listUtilities addStringInClassList:mod.simulation theVariable:@"post file" withValue:@"ismip_A010_.ep"];
    
    ivalue = 4;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"maximum output level" withValue:&ivalue orUsingBlock:nil];
    
    mod.bodies = @[ @{@"equation" : @1,
                      @"body force" : @1,
                      @"material" : @1,
                      @"initial condition" : @1}];
    mod.numberOfBodies = 1;
    
    FEMInitialConditions *initialConditions = [[FEMInitialConditions alloc] init];
    double value = 0.0;
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"pressure" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 3" withValue:&value orUsingBlock:nil string:nil];
    mod.initialConditions = @[initialConditions];
    mod.numberOfInitialConditions = 1;
    
    FEMBodyForce *bodyForce = [[FEMBodyForce alloc] init];
    [listUtilities addConstRealInClassList:bodyForce theVariable:@"flow bodyforce 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:bodyForce theVariable:@"flow bodyforce 2" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:bodyForce theVariable:@"flow bodyforce 3" withValue:&gravity orUsingBlock:nil string:nil];
    mod.bodyForces = @[bodyForce];
    mod.numberOfBodyForces = 1;
    
    FEMMaterial *material = [[FEMMaterial alloc] init];
    [listUtilities addConstRealInClassList:material theVariable:@"density" withValue:&rhoi orUsingBlock:nil string:nil];
    [listUtilities addStringInClassList:material theVariable:@"viscosity model" withValue:@"power law"];
    [listUtilities addConstRealInClassList:material theVariable:@"viscosity" withValue:&eta orUsingBlock:nil string:nil];
    value = 1.0 / n;
    [listUtilities addConstRealInClassList:material theVariable:@"viscosity exponent" withValue:&value orUsingBlock:nil string:nil];
    value = 1.0e-10;
    [listUtilities addConstRealInClassList:material theVariable:@"critical shear rate" withValue:&value orUsingBlock:nil string:nil];
    
    [listUtilities addLogicalInClassList:material theVariable:@"cauchy" withValue:NO];
    mod.materials = @[material];
    mod.numberOfMaterials = 1;
    
    FEMSolution *solution1 = [[FEMSolution alloc] init];
    NSString *user = [@"/Users" stringByAppendingPathComponent:NSUserName()];
    solution1.plugInName = [user stringByAppendingPathComponent:@"Documents/Saino/PlugIns/FEMStructuredMeshMapper"];
    [solution1.solutionInfo setObject:@"mapcoordinate" forKey:@"equation"];
    [solution1.solutionInfo setObject:@3 forKey:@"active coordinate"];
    //    [solution1.solutionInfo setObject:@"dsdt" forKey:@"mesh velocity variable"];
    //    [solution1.solutionInfo setObject:@"ds" forKey:@"mesh update variable"];
    [solution1.solutionInfo setObject:@YES forKey:@"mesh velocity first zero"];
    
    FEMSolution *solution2 = [[FEMSolution alloc] init];
    [solution2.solutionInfo setObject:@"navier-stokes" forKey:@"equation"];
    [solution2.solutionInfo setObject:@"stabilized" forKey:@"stabilization method"];
    [solution2.solutionInfo setObject:@"stokes" forKey:@"flow model"];
    //    [solution2.solutionInfo setObject:@"-dofs 1 dsdt" forKey:@"exported variable 1"];
    //    [solution2.solutionInfo setObject:@"-dofs 1 ds" forKey:@"exported variable 1"];
    [solution2.solutionInfo setObject:@"iterative" forKey:@"linear system solver"];
    [solution2.solutionInfo setObject:@"bi-cgstab" forKey:@"linear system iterative method"];
    [solution2.solutionInfo setObject:@"ilu0" forKey:@"linear system preconditioning"];
    [solution2.solutionInfo setObject:@1.0e-8 forKey:@"linear system convergence tolerance"];
    [solution2.solutionInfo setObject:@500 forKey:@"linear system maximum iterations"];
    [solution2.solutionInfo setObject:@500 forKey:@"linear system residual output"];
    [solution2.solutionInfo setValue:@YES forKey:@"linear system abort not converged"];
    [solution2.solutionInfo setValue:@50 forKey:@"nonlinear system maximum iterations"];
    [solution2.solutionInfo setValue:@1.0e-5 forKey:@"nonlinear system convergence tolerance"];
    [solution2.solutionInfo setValue:@5 forKey:@"nonlinear system newton after iterations"];
    [solution2.solutionInfo setValue:@1.0e-2 forKey:@"nonlinear system newton after tolerance"];
    [solution2.solutionInfo setValue:@1.0 forKey:@"nonlinear system relaxation factor"];
    [solution2.solutionInfo setValue:@1.0e-3 forKey:@"steady state convergence tolerance"];
    
    mod.solutions = @[solution1, solution2];
    mod.numberOfSolutions = 2;
    
    FEMEquation *equation = [[FEMEquation alloc] init];
    vector = intvec(0, 1);
    vector[0] = 1;
    vector[1] = 2;
    [listUtilities addIntegerArrayInClassList:equation theVariable:@"active solvers" withValues:vector size:2 orUsingBlock:nil];
    free_ivector(vector, 0, 1);
    mod.equations = @[equation];
    mod.numberOfEquations = 1;
    
    mod.boundaryID = @[@1, @2, @3, @4, @5, @6];
    
    // BC y=y0
    FEMBoundaryCondition *boundaryCondition1 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 1;
    [listUtilities addIntegerArrayInClassList:boundaryCondition1 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
    double **translate = doublematrix(0, 2, 0, 0);
    
    // BC x=xmax
    FEMBoundaryCondition *boundaryCondition2 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 2;
    [listUtilities addIntegerArrayInClassList:boundaryCondition2 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    ivalue = 4;
    [listUtilities addIntegerInClassList:boundaryCondition2 theVariable:@"periodic bc" withValue:&ivalue orUsingBlock:nil];
    translate[0][0] = L;
    translate[1][0] = 0.0;
    translate[2][0] = 0.0;
    [listUtilities addConstRealArrayInClassList:boundaryCondition2 theVariable:@"periodic bc translate" withValues:translate size1:3 size2:1 orUsingBlock:nil string:nil];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc velocity 1" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc velocity 2" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc velocity 3" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc pressure" withValue:YES];
    
    // BC y=ymax
    FEMBoundaryCondition *boundaryCondition3 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 3;
    [listUtilities addIntegerArrayInClassList:boundaryCondition3 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    ivalue = 1;
    [listUtilities addIntegerInClassList:boundaryCondition3 theVariable:@"periodic bc" withValue:&ivalue orUsingBlock:nil];
    translate[0][0] = 0.0;
    translate[1][0] = L;
    translate[2][0] = 0.0;
    [listUtilities addConstRealArrayInClassList:boundaryCondition3 theVariable:@"periodic bc translate" withValues:translate size1:3 size2:1 orUsingBlock:nil string:nil];
    [listUtilities addLogicalInClassList:boundaryCondition3 theVariable:@"periodic bc velocity 1" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition3 theVariable:@"periodic bc velocity 2" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition3 theVariable:@"periodic bc velocity 3" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition3 theVariable:@"periodic bc pressure" withValue:YES];
    
    free_dmatrix(translate, 0, 2, 0, 0);
    
    // BC x=x0
    FEMBoundaryCondition *boundaryCondition4 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 4;
    [listUtilities addIntegerArrayInClassList:boundaryCondition4 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
    // BC bedrock
    FEMBoundaryCondition *boundaryCondition5 = [[FEMBoundaryCondition alloc] init];
    value = 0.0;
    [listUtilities addConstRealInClassList:boundaryCondition5 theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:boundaryCondition5 theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:boundaryCondition5 theVariable:@"velocity 3" withValue:&value orUsingBlock:nil string:nil];
    NSArray *dependencies1 = @[@"coordinate 1", @"coordinate 2"];
    double (^block1) (double *) = ^(double *t) {
        return -t[0] * tan(slope) - 1000.0 + 500.0 * sin(2.0*M_PI*t[0]/L) * sin(2.0*M_PI*t[1]/L);
    };
    [listUtilities addBlockInClassList:boundaryCondition5 theVariable:@"bottom surface" usingBlock:block1 dependencies:dependencies1];
    
    // BC free surface
    FEMBoundaryCondition *boundaryCondition6 = [[FEMBoundaryCondition alloc] init];
    NSArray *dependencies2 = @[@"coordinate 1"];
    double (^block2) (double *) = ^(double *t) {
        return -t[0] * tan(slope);
    };
    [listUtilities addBlockInClassList:boundaryCondition6 theVariable:@"top surface" usingBlock:block2 dependencies:dependencies2];
    
    boundaryCondition1.tag = 1;
    boundaryCondition2.tag = 2;
    boundaryCondition3.tag = 3;
    boundaryCondition4.tag = 4;
    boundaryCondition5.tag = 5;
    boundaryCondition6.tag = 6;
    mod.boundaryConditions = @[boundaryCondition1, boundaryCondition2, boundaryCondition3, boundaryCondition4, boundaryCondition5, boundaryCondition6];
    mod.numberOfBoundaryConditions = 6;
    
    mod.meshDir = [NSMutableString stringWithString:self.path];
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80"];
}

-(void)setUpISMIP_HOM_A010_GPU_Coloring_20x20:(id __nonnull)model {
    
    [self setUpISMIP_HOM_A010:model];
    
    FEMModel *mod = (FEMModel *)model;
    
    double yearinsec = 365.25 * 24.0 * 60.0 * 60.0;
    double rhoi = 900.0 / (1.0e6 * pow(yearinsec, 2.0));
    double gravity = -9.81 * pow(yearinsec, 2.0);
    double n = 3.0;
    double eta = pow((2.0 * 100.0), (-1.0/n));
    
    NSString *kernelPath = @"/Users/hakimeseddik/Documents/Saino/Saino/NavierStokesAssemblyKernel_opt";
    
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setValue:@100 forKey:@"nonlinear system newton after iterations"];
            [solution.solutionInfo setValue:@1.0e-10 forKey:@"nonlinear system newton after tolerance"];
            
            [solution.solutionInfo setObject:@NO forKey:@"optimize bandwidth"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly"];
            [solution.solutionInfo setObject:@"element coloring" forKey:@"parallel assembly method"];
            [solution.solutionInfo setObject:@YES forKey:@"color mesh"];
            [solution.solutionInfo setObject:kernelPath forKey:@"gpu kernel source file"];
            [solution.solutionInfo setObject:@"ice flow" forKey:@"gpu flow type"];
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@(rhoi) forKey:@"gpu ice density"];
            [solution.solutionInfo setObject:@(eta) forKey:@"gpu ice viscosity"];
            [solution.solutionInfo setObject:@(gravity) forKey:@"gpu ice gravity"];
            [solution.solutionInfo setObject:@NO forKeyedSubscript:@"enable newton linearization"];
            [solution.solutionInfo setObject:@YES forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@YES forKey:@"use gpu local memory"];
            
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@16 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@YES forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            [solution.solutionInfo setObject:@YES forKey:@"precompute nonzero indices"];
            [solution.solutionInfo setObject:@YES forKey:@"use global stiff and force"];
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_colored"];
}

-(void)setUpISMIP_HOM_A010_GPU_Coloring_40x40:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_Coloring_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@YES forKey:@"precompute nonzero indices"];
            [solution.solutionInfo setObject:@NO forKey:@"use global stiff and force"];
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_40x40_colored"];
}

-(void)setUpISMIP_HOM_A010_GPU_Coloring_70x70:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_Coloring_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];

            [solution.solutionInfo setObject:@YES forKey:@"precompute nonzero indices"];
            [solution.solutionInfo setObject:@NO forKey:@"use global stiff and force"];
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_dense1_colored"];
}

-(void)setUpISMIP_HOM_A010_GPU_Coloring_100x100:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_Coloring_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@YES forKey:@"precompute nonzero indices"];
            [solution.solutionInfo setObject:@NO forKey:@"use global stiff and force"];
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_100x100_colored"];
}

-(void)setUpISMIP_HOM_A010_GPU_Coloring_120x120:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_Coloring_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@YES forKey:@"precompute nonzero indices"];
            [solution.solutionInfo setObject:@NO forKey:@"use global stiff and force"];
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_dense2_colored"];
}

-(void)setUpISMIP_HOM_A010_GPU_Coloring_150x150:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_Coloring_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@YES forKey:@"precompute nonzero indices"];
            [solution.solutionInfo setObject:@NO forKey:@"use global stiff and force"];
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_dense3_colored"];
}

-(void)setUpISMIP_HOM_A010_GPU_Coloring_170x170:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_Coloring_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@YES forKey:@"precompute nonzero indices"];
            [solution.solutionInfo setObject:@NO forKey:@"use global stiff and force"];
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_dense4_colored"];
}

-(void)setUpISMIP_HOM_A010_GPU_NonZeros_20x20:(id __nonnull)model {
    
    [self setUpISMIP_HOM_A010:model];
    
    FEMModel *mod = (FEMModel *)model;
    
    double yearinsec = 365.25 * 24.0 * 60.0 * 60.0;
    double rhoi = 900.0 / (1.0e6 * pow(yearinsec, 2.0));
    double gravity = -9.81 * pow(yearinsec, 2.0);
    double n = 3.0;
    double eta = pow((2.0 * 100.0), (-1.0/n));
    
    NSString *kernelPath = @"/Users/hakimeseddik/Documents/Saino/Saino/NavierStokesAssemblyKernel_opt";
    
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@NO forKey:@"optimize bandwidth"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly"];
            [solution.solutionInfo setObject:@"nonzero entries" forKey:@"parallel assembly method"];
            [solution.solutionInfo setObject:kernelPath forKey:@"gpu kernel source file"];
            [solution.solutionInfo setObject:@"ice flow" forKey:@"gpu flow type"];
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@(rhoi) forKey:@"gpu ice density"];
            [solution.solutionInfo setObject:@(eta) forKey:@"gpu ice viscosity"];
            [solution.solutionInfo setObject:@(gravity) forKey:@"gpu ice gravity"];
            [solution.solutionInfo setObject:@NO forKeyedSubscript:@"enable newton linearization"];
            [solution.solutionInfo setObject:@YES forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            [solution.solutionInfo setObject:@YES forKey:@"use global basis functions coefficients"];
            
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros assembly work-group size"];
            [solution.solutionInfo setObject:@64 forKey:@"nonzeros assembly thread block size"];
            [solution.solutionInfo setObject:@4 forKey:@"nonzeros assembly nonzeros per thread"];
            [solution.solutionInfo setObject:@1024 forKey:@"maximum allowed allocation size for reduction matrix buffer (MB)"];
            
            [solution.solutionInfo setObject:@YES forKey:@"enable GPU debug mode"];
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_nzs"];
}

-(void)setUpISMIP_HOM_A010_GPU_NonZeros_40x40:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_NonZeros_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros assembly work-group size"];
            [solution.solutionInfo setObject:@64 forKey:@"nonzeros assembly thread block size"];
            [solution.solutionInfo setObject:@4 forKey:@"nonzeros assembly nonzeros per thread"];
            
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_40x40_nzs"];
}

-(void)setUpISMIP_HOM_A010_GPU_NonZeros_70x70:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_NonZeros_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];

            [solution.solutionInfo setObject:@256 forKey:@"nonzeros adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros assembly work-group size"];
            [solution.solutionInfo setObject:@64 forKey:@"nonzeros assembly thread block size"];
            [solution.solutionInfo setObject:@4 forKey:@"nonzeros assembly nonzeros per thread"];

            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_dense1_nzs"];
}

-(void)setUpISMIP_HOM_A010_GPU_NonZeros_100x100:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_NonZeros_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros assembly work-group size"];
            [solution.solutionInfo setObject:@64 forKey:@"nonzeros assembly thread block size"];
            [solution.solutionInfo setObject:@4 forKey:@"nonzeros assembly nonzeros per thread"];
            
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_100x100_nzs"];
}

-(void)setUpISMIP_HOM_A010_GPU_NonZeros_120x120:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_NonZeros_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"double" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros assembly work-group size"];
            [solution.solutionInfo setObject:@64 forKey:@"nonzeros assembly thread block size"];
            [solution.solutionInfo setObject:@4 forKey:@"nonzeros assembly nonzeros per thread"];
            
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_dense2_nzs"];
}

-(void)setUpISMIP_HOM_A010_GPU_NonZeros_150x150:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_NonZeros_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"single" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros assembly work-group size"];
            [solution.solutionInfo setObject:@64 forKey:@"nonzeros assembly thread block size"];
            [solution.solutionInfo setObject:@4 forKey:@"nonzeros assembly nonzeros per thread"];
            
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_dense3_nzs"];
}

-(void)setUpISMIP_HOM_A010_GPU_NonZeros_170x170:(id __nonnull)model {
    [self setUpISMIP_HOM_A010_GPU_NonZeros_20x20:model];
    FEMModel *mod = (FEMModel *)model;
    for (FEMSolution *solution in mod.solutions) {
        if ([solution.solutionInfo[@"equation"] isEqualToString:@"navier-stokes"] == YES) {
            [solution.solutionInfo setObject:@"single" forKey:@"gpu floating-point precision"];
            [solution.solutionInfo setObject:@NO forKey:@"compute basis and basis derivatives in separate kernel"];
            [solution.solutionInfo setObject:@NO forKey:@"use gpu local memory"];
            [solution.solutionInfo setObject:@YES forKey:@"parallel assembly enable work-groups"];
            // If single precision and use of local memory or if work-groups are enabled:
            [solution.solutionInfo setObject:@64 forKey:@"adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@64 forKey:@"coloring assembly/locals compute work-group size"];
            
            [solution.solutionInfo setObject:@NO forKey:@"use global basis functions coefficients"];
            [solution.solutionInfo setObject:@YES forKey:@"use element nodal data"];
            
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros adjust global work size to be a multiple of"];
            [solution.solutionInfo setObject:@256 forKey:@"nonzeros assembly work-group size"];
            [solution.solutionInfo setObject:@64 forKey:@"nonzeros assembly thread block size"];
            [solution.solutionInfo setObject:@4 forKey:@"nonzeros assembly nonzeros per thread"];
            
            
            [solution.solutionInfo setObject:@YES forKey:@"enable gpu multiply-and-add operations"];
        }
    }
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80_dense4_nzs"];
}

-(void)setUpISMIP_HOM_B010:(id __nonnull)model {
    
    FEMModel *mod = (FEMModel *)model;
    
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    
    double L = 80.0e3;
    double slope = 0.5 * M_PI / 180.0;
    
    double yearinsec = 365.25 * 24.0 * 60.0 * 60.0;
    double rhoi = 900.0 / (1.0e6 * pow(yearinsec, 2.0));
    double gravity = -9.81 * pow(yearinsec, 2.0);
    double n = 3.0;
    double eta = pow((2.0 * 100.0), (-1.0/n));
    
    self.path = [NSMutableString stringWithString:@"/Users/hakimeseddik/Desktop/ISMIP-HOM/B010"];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"coordinate system" withValue:@"cartesian 2d"];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"simulaton type" withValue:@"steady state"];
    
    int ivalue = 1;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"steady state maximum iterations" withValue:&ivalue orUsingBlock:nil];
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"steady state minimum iterations" withValue:&ivalue orUsingBlock:nil];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"output file" withValue:@"ismip_B010_.result"];
    [listUtilities addStringInClassList:mod.simulation theVariable:@"post file" withValue:@"ismip_B010_.ep"];
    
    ivalue = 3;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"maximum output level" withValue:&ivalue orUsingBlock:nil];
    
    mod.bodies = @[ @{@"equation" : @1,
                      @"body force" : @1,
                      @"material" : @1,
                      @"initial condition" : @1}];
    mod.numberOfBodies = 1;
    
    FEMInitialConditions *initialConditions = [[FEMInitialConditions alloc] init];
    double value = 0.0;
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"pressure" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    mod.initialConditions = @[initialConditions];
    mod.numberOfInitialConditions = 1;
    
    FEMBodyForce *bodyForce = [[FEMBodyForce alloc] init];
    [listUtilities addConstRealInClassList:bodyForce theVariable:@"flow bodyforce 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:bodyForce theVariable:@"flow bodyforce 2" withValue:&gravity orUsingBlock:nil string:nil];
    mod.bodyForces = @[bodyForce];
    mod.numberOfBodyForces = 1;
    
    FEMMaterial *material = [[FEMMaterial alloc] init];
    [listUtilities addConstRealInClassList:material theVariable:@"density" withValue:&rhoi orUsingBlock:nil string:nil];
    [listUtilities addStringInClassList:material theVariable:@"viscosity model" withValue:@"power law"];
    [listUtilities addConstRealInClassList:material theVariable:@"viscosity" withValue:&eta orUsingBlock:nil string:nil];
    value = 1.0 / n;
    [listUtilities addConstRealInClassList:material theVariable:@"viscosity exponent" withValue:&value orUsingBlock:nil string:nil];
    value = 1.0e-10;
    [listUtilities addConstRealInClassList:material theVariable:@"critical shear rate" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addLogicalInClassList:material theVariable:@"cauchy" withValue:NO];
    mod.materials = @[material];
    mod.numberOfMaterials = 1;
    
    FEMSolution *solution1 = [[FEMSolution alloc] init];
    NSString *user = [@"/Users" stringByAppendingPathComponent:NSUserName()];
    solution1.plugInName = [user stringByAppendingPathComponent:@"Documents/Saino/PlugIns/FEMStructuredMeshMapper"];
    [solution1.solutionInfo setObject:@"mapcoordinate" forKey:@"equation"];
    [solution1.solutionInfo setObject:@2 forKey:@"active coordinate"];
    //    [solution1.solutionInfo setObject:@"dsdt" forKey:@"mesh velocity variable"];
    //    [solution1.solutionInfo setObject:@"ds" forKey:@"mesh update variable"];
    [solution1.solutionInfo setObject:@YES forKey:@"mesh velocity first zero"];
    
    FEMSolution *solution2 = [[FEMSolution alloc] init];
    [solution2.solutionInfo setObject:@"navier-stokes" forKey:@"equation"];
    [solution2.solutionInfo setObject:@"stabilized" forKey:@"stabilization method"];
    [solution2.solutionInfo setObject:@"stokes" forKey:@"flow model"];
    //    [solution2.solutionInfo setObject:@"-dofs 1 dsdt" forKey:@"exported variable 1"];
    //    [solution2.solutionInfo setObject:@"-dofs 1 ds" forKey:@"exported variable 1"];
    [solution2.solutionInfo setObject:@"iterative" forKey:@"linear system solver"];
    [solution2.solutionInfo setObject:@"bi-cgstab" forKey:@"linear system iterative method"];
    [solution2.solutionInfo setObject:@"ilu1" forKey:@"linear system preconditioning"];
    [solution2.solutionInfo setObject:@1.0e-8 forKey:@"linear system convergence tolerance"];
    [solution2.solutionInfo setObject:@500 forKey:@"linear system maximum iterations"];
    [solution2.solutionInfo setObject:@500 forKey:@"linear system residual output"];
    [solution2.solutionInfo setValue:@YES forKey:@"linear system abort not converged"];
    [solution2.solutionInfo setValue:@100 forKey:@"nonlinear system maximum iterations"];
    [solution2.solutionInfo setValue:@1.0e-5 forKey:@"nonlinear system convergence tolerance"];
    [solution2.solutionInfo setValue:@5 forKey:@"nonlinear system newton after iterations"];
    [solution2.solutionInfo setValue:@1.0e-2 forKey:@"nonlinear system newton after tolerance"];
    [solution2.solutionInfo setValue:@1.0 forKey:@"nonlinear system relaxation factor"];
    [solution2.solutionInfo setValue:@1.0e-3 forKey:@"steady state convergence tolerance"];
    [solution2.solutionInfo setValue:@YES forKeyPath:@"nonlinear system reset newton"];
    
    mod.solutions = @[solution1, solution2];
    mod.numberOfSolutions = 2;
    
    FEMEquation *equation = [[FEMEquation alloc] init];
    int *vector = intvec(0, 1);
    vector[0] = 1;
    vector[1] = 2;
    [listUtilities addIntegerArrayInClassList:equation theVariable:@"active solvers" withValues:vector size:2 orUsingBlock:nil];
    free_ivector(vector, 0, 1);
    mod.equations = @[equation];
    mod.numberOfEquations = 1;
    
    mod.boundaryID = @[@1, @2, @3, @4];
    
    // BC bedrock
    FEMBoundaryCondition *boundaryCondition1 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 1;
    [listUtilities addIntegerArrayInClassList:boundaryCondition1 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    value = 0.0;
    [listUtilities addConstRealInClassList:boundaryCondition1 theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:boundaryCondition1 theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    NSArray *dependencies1 = @[@"coordinate 1"];
    double (^block1) (double *) = ^(double *t) {
        return -t[0] * tan(slope) - 1000.0 + 500.0 * sin(2.0*M_PI*t[0]/L);
    };
    [listUtilities addBlockInClassList:boundaryCondition1 theVariable:@"bottom surface" usingBlock:block1 dependencies:dependencies1];
    
    
    // BC periodic right
    FEMBoundaryCondition *boundaryCondition2 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 2;
    [listUtilities addIntegerArrayInClassList:boundaryCondition2 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    ivalue = 4;
    [listUtilities addIntegerInClassList:boundaryCondition2 theVariable:@"periodic bc" withValue:&ivalue orUsingBlock:nil];
    double **translate = doublematrix(0, 1, 0, 0);
    translate[0][0] = L;
    translate[1][0] = 0.0;
    [listUtilities addConstRealArrayInClassList:boundaryCondition2 theVariable:@"periodic bc translate" withValues:translate size1:2 size2:1 orUsingBlock:nil string:nil];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc velocity 1" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc velocity 2" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc pressure" withValue:YES];
    free_dmatrix(translate, 0, 1, 0, 0);
    
    // BC free surface
    FEMBoundaryCondition *boundaryCondition3 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 3;
    [listUtilities addIntegerArrayInClassList:boundaryCondition3 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    NSArray *dependencies2 = @[@"coordinate 1"];
    double (^block2) (double *) = ^(double *t) {
        return -t[0] * tan(slope);
    };
    [listUtilities addBlockInClassList:boundaryCondition3 theVariable:@"top surface" usingBlock:block2 dependencies:dependencies2];
    
    // BC periodic left
    FEMBoundaryCondition *boundaryCondition4 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 4;
    [listUtilities addIntegerArrayInClassList:boundaryCondition4 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
    boundaryCondition1.tag = 1;
    boundaryCondition2.tag = 2;
    boundaryCondition3.tag = 3;
    boundaryCondition4.tag = 4;
    mod.boundaryConditions = @[boundaryCondition1, boundaryCondition2, boundaryCondition3, boundaryCondition4];
    mod.numberOfBoundaryConditions = 4;
    
    mod.meshDir = [NSMutableString stringWithString:self.path];
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80"];
}

-(void)setUpISMIP_HOM_C010:(id __nonnull)model {
    
    FEMModel *mod = (FEMModel *)model;
    
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    
    double L = 80.0e3;
    double slope = 0.1 * M_PI / 180.0;
    
    double yearinsec = 365.25 * 24.0 * 60.0 * 60.0;
    double rhoi = 900.0 / (1.0e6 * pow(yearinsec, 2.0));
    double gravity = -9.81 * pow(yearinsec, 2.0);
    double n = 3.0;
    double eta = pow((2.0 * 100.0), (-1.0/n));
    
    self.path = [NSMutableString stringWithString:@"/Users/hakimeseddik/Desktop/ISMIP-HOM/C010"];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"coordinate system" withValue:@"cartesian 3d"];
    
    int *vector = intvec(0, 2);
    vector[0] = 1; vector[1] = 2; vector[2] = 3;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"coordinate mapping" withValues:vector size:3 orUsingBlock:nil];
    free_ivector(vector, 0, 2);
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"simulaton type" withValue:@"steady state"];
    
    int ivalue = 16;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"extruded mesh levels" withValue:&ivalue orUsingBlock:nil];
    
    int *intervals = intvec(0, 0);
    intervals[0] = 1;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"output intervals" withValues:intervals size:1 orUsingBlock:nil];
    free_ivector(intervals, 0, 0);
    
    ivalue = 1;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"steady state maximum iterations" withValue:&ivalue orUsingBlock:nil];
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"steady state minimum iterations" withValue:&ivalue orUsingBlock:nil];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"output file" withValue:@"ismip_C010_.result"];
    [listUtilities addStringInClassList:mod.simulation theVariable:@"post file" withValue:@"ismip_C010_.ep"];
    
    ivalue = 20;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"maximum output level" withValue:&ivalue orUsingBlock:nil];
    
    mod.bodies = @[ @{@"equation" : @1,
                      @"body force" : @1,
                      @"material" : @1,
                      @"initial condition" : @1}];
    mod.numberOfBodies = 1;
    
    FEMInitialConditions *initialConditions = [[FEMInitialConditions alloc] init];
    double value = 0.0;
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"pressure" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 3" withValue:&value orUsingBlock:nil string:nil];
    value = 10.0;
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    mod.initialConditions = @[initialConditions];
    mod.numberOfInitialConditions = 1;
    
    FEMBodyForce *bodyForce = [[FEMBodyForce alloc] init];
    [listUtilities addConstRealInClassList:bodyForce theVariable:@"flow bodyforce 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:bodyForce theVariable:@"flow bodyforce 2" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:bodyForce theVariable:@"flow bodyforce 3" withValue:&gravity orUsingBlock:nil string:nil];
    mod.bodyForces = @[bodyForce];
    mod.numberOfBodyForces = 1;
    
    FEMMaterial *material = [[FEMMaterial alloc] init];
    [listUtilities addConstRealInClassList:material theVariable:@"density" withValue:&rhoi orUsingBlock:nil string:nil];
    [listUtilities addStringInClassList:material theVariable:@"viscosity model" withValue:@"power law"];
    [listUtilities addConstRealInClassList:material theVariable:@"viscosity" withValue:&eta orUsingBlock:nil string:nil];
    value = 1.0 / n;
    [listUtilities addConstRealInClassList:material theVariable:@"viscosity exponent" withValue:&value orUsingBlock:nil string:nil];
    value = 1.0e-10;
    [listUtilities addConstRealInClassList:material theVariable:@"critical shear rate" withValue:&value orUsingBlock:nil string:nil];
    
    [listUtilities addLogicalInClassList:material theVariable:@"cauchy" withValue:NO];
    mod.materials = @[material];
    mod.numberOfMaterials = 1;
    
    FEMSolution *solution1 = [[FEMSolution alloc] init];
    NSString *user = [@"/Users" stringByAppendingPathComponent:NSUserName()];
    solution1.plugInName = [user stringByAppendingPathComponent:@"Documents/Saino/PlugIns/FEMStructuredMeshMapper"];
    [solution1.solutionInfo setObject:@"mapcoordinate" forKey:@"equation"];
    [solution1.solutionInfo setObject:@3 forKey:@"active coordinate"];
    //    [solution1.solutionInfo setObject:@"dsdt" forKey:@"mesh velocity variable"];
    //    [solution1.solutionInfo setObject:@"ds" forKey:@"mesh update variable"];
    [solution1.solutionInfo setObject:@YES forKey:@"mesh velocity first zero"];
    
    FEMSolution *solution2 = [[FEMSolution alloc] init];
    [solution2.solutionInfo setObject:@"navier-stokes" forKey:@"equation"];
    [solution2.solutionInfo setObject:@"stabilized" forKey:@"stabilization method"];
    [solution2.solutionInfo setObject:@"stokes" forKey:@"flow model"];
    //    [solution2.solutionInfo setObject:@"-dofs 1 dsdt" forKey:@"exported variable 1"];
    //    [solution2.solutionInfo setObject:@"-dofs 1 ds" forKey:@"exported variable 1"];
    [solution2.solutionInfo setObject:@"iterative" forKey:@"linear system solver"];
    [solution2.solutionInfo setObject:@"bi-cgstab" forKey:@"linear system iterative method"];
    [solution2.solutionInfo setObject:@"ilu0" forKey:@"linear system preconditioning"];
    [solution2.solutionInfo setObject:@1.0e-6 forKey:@"linear system convergence tolerance"];
    [solution2.solutionInfo setObject:@1000 forKey:@"linear system maximum iterations"];
    [solution2.solutionInfo setValue:@YES forKey:@"linear system abort not converged"];
    [solution2.solutionInfo setValue:@50 forKey:@"nonlinear system maximum iterations"];
    [solution2.solutionInfo setValue:@1.0e-5 forKey:@"nonlinear system convergence tolerance"];
    [solution2.solutionInfo setValue:@5 forKey:@"nonlinear system newton after iterations"];
    [solution2.solutionInfo setValue:@1.0e-2 forKey:@"nonlinear system newton after tolerance"];
    [solution2.solutionInfo setValue:@1.0 forKey:@"nonlinear system relaxation factor"];
    [solution2.solutionInfo setValue:@1.0e-3 forKey:@"steady state convergence tolerance"];
    //[solution2.solutionInfo setObject:@"never" forKey:@"invoke solution computer"];
    
    mod.solutions = @[solution1, solution2];
    mod.numberOfSolutions = 2;
    
    FEMEquation *equation = [[FEMEquation alloc] init];
    vector = intvec(0, 1);
    vector[0] = 1;
    vector[1] = 2;
    [listUtilities addIntegerArrayInClassList:equation theVariable:@"active solvers" withValues:vector size:2 orUsingBlock:nil];
    free_ivector(vector, 0, 1);
    mod.equations = @[equation];
    mod.numberOfEquations = 1;
    
    mod.boundaryID = @[@1, @2, @3, @4, @5, @6];
    
    // BC y=y0
    FEMBoundaryCondition *boundaryCondition1 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 1;
    [listUtilities addIntegerArrayInClassList:boundaryCondition1 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
    double **translate = doublematrix(0, 2, 0, 0);
    
    // BC x=xmax
    FEMBoundaryCondition *boundaryCondition2 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 2;
    [listUtilities addIntegerArrayInClassList:boundaryCondition2 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    ivalue = 4;
    [listUtilities addIntegerInClassList:boundaryCondition2 theVariable:@"periodic bc" withValue:&ivalue orUsingBlock:nil];
    translate[0][0] = L;
    translate[1][0] = 0.0;
    translate[2][0] = 0.0;
    [listUtilities addConstRealArrayInClassList:boundaryCondition2 theVariable:@"periodic bc translate" withValues:translate size1:3 size2:1 orUsingBlock:nil string:nil];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc velocity 1" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc velocity 2" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc velocity 3" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition2 theVariable:@"periodic bc pressure" withValue:YES];
    
    // BC y=ymax
    FEMBoundaryCondition *boundaryCondition3 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 3;
    [listUtilities addIntegerArrayInClassList:boundaryCondition3 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    ivalue = 1;
    [listUtilities addIntegerInClassList:boundaryCondition3 theVariable:@"periodic bc" withValue:&ivalue orUsingBlock:nil];
    translate[0][0] = 0.0;
    translate[1][0] = L;
    translate[2][0] = 0.0;
    [listUtilities addConstRealArrayInClassList:boundaryCondition3 theVariable:@"periodic bc translate" withValues:translate size1:3 size2:1 orUsingBlock:nil string:nil];
    [listUtilities addLogicalInClassList:boundaryCondition3 theVariable:@"periodic bc velocity 1" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition3 theVariable:@"periodic bc velocity 2" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition3 theVariable:@"periodic bc velocity 3" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition3 theVariable:@"periodic bc pressure" withValue:YES];
    
    free_dmatrix(translate, 0, 2, 0, 0);
    
    // BC x=x0
    FEMBoundaryCondition *boundaryCondition4 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 4;
    [listUtilities addIntegerArrayInClassList:boundaryCondition4 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
    // BC bedrock
    FEMBoundaryCondition *boundaryCondition5 = [[FEMBoundaryCondition alloc] init];
    [listUtilities addLogicalInClassList:boundaryCondition5 theVariable:@"flow force bc" withValue:YES];
    [listUtilities addLogicalInClassList:boundaryCondition5 theVariable:@"normal-tangential velocity" withValue:YES];
    
    value = 0.0;
    [listUtilities addConstRealInClassList:boundaryCondition5 theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    
    NSArray *dependencies1 = @[@"coordinate 1", @"coordinate 2"];
    double (^block0) (double *) = ^(double *t) {
        return 1.0e-3*(1.0 + sin(2.0*M_PI*t[0]/L) * sin(2.0*M_PI*t[1]/L));
    };
    [listUtilities addBlockInClassList:boundaryCondition5 theVariable:@"slip coefficient 2" usingBlock:block0 dependencies:dependencies1];
    [listUtilities addBlockInClassList:boundaryCondition5 theVariable:@"slip coefficient 3" usingBlock:block0 dependencies:dependencies1];
    
    NSArray *dependencies2 = @[@"coordinate 1"];
    double (^block1) (double *) = ^(double *t) {
        return -t[0] * tan(slope) - 1000.0;
    };
    [listUtilities addBlockInClassList:boundaryCondition5 theVariable:@"bottom surface" usingBlock:block1 dependencies:dependencies2];
    
    // BC free surface
    FEMBoundaryCondition *boundaryCondition6 = [[FEMBoundaryCondition alloc] init];
    double (^block2) (double *) = ^(double *t) {
        return -t[0] * tan(slope);
    };
    [listUtilities addBlockInClassList:boundaryCondition6 theVariable:@"top surface" usingBlock:block2 dependencies:dependencies2];
    
    boundaryCondition1.tag = 1;
    boundaryCondition2.tag = 2;
    boundaryCondition3.tag = 3;
    boundaryCondition4.tag = 4;
    boundaryCondition5.tag = 5;
    boundaryCondition6.tag = 6;
    mod.boundaryConditions = @[boundaryCondition1, boundaryCondition2, boundaryCondition3, boundaryCondition4, boundaryCondition5, boundaryCondition6];
    mod.numberOfBoundaryConditions = 6;
    
    mod.meshDir = [NSMutableString stringWithString:self.path];
    mod.meshName = [NSMutableString stringWithString:@"rectangle_L80"];
}



@end
