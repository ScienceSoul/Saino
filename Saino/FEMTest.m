//
//  FEMTest.m
//  Saino
//
//  Created by Seddik hakime on 09/06/2015.
//  Copyright (c) 2015 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMTest.h"

#import "FEMModel.h"
#import "FEMListUtilities.h"
#import "FEMSolution.h"
#import "memory.h"

@implementation FEMTest

@synthesize path = _path;
@synthesize do_heatq = _do_heatq;
@synthesize do_step_stokes = _do_step_stokes;
@synthesize do_natural_convection = _do_natural_convection;
@synthesize do_ismip_hom_A010 = _do_ismip_hom_A010;
@synthesize heatq_allDone = _heatq_allDone;
@synthesize step_stokes_allDone = _step_stokes_allDone;
@synthesize natural_convection_allDone = _natural_convection_allDone;
@synthesize ismip_hom_allDone = _ismip_hom_allDone;
@synthesize norm = _norm;

static FEMTest * __nullable sharedTest = nil;
static dispatch_once_t onceToken;

+(id __nonnull)sharedTest {
    
    dispatch_once(&onceToken, ^{
        sharedTest = [[self alloc] init];
    });
    return sharedTest;
}

+(void) selfDestruct {
    sharedTest = nil;
    onceToken = 0;
}

- (id)init
{
    self = [super init];
    if (self) {
        
    }
    
    return self;
}

-(void)reset {
    _do_heatq = NO;
    _do_step_stokes = NO;
    _do_natural_convection = NO;
    _do_ismip_hom_A010 = NO;
    _heatq_allDone = NO;
    _step_stokes_allDone = NO;
    _natural_convection_allDone = NO;
    _ismip_hom_allDone = NO;
    _norm = 0.0;
}

-(void)setUpHeateqTest:(id __nonnull)model {
    
    FEMModel *mod = (FEMModel *)model;
    
    FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"coordinate system" withValue:@"cartesian 2d"];
    
    int *vector = intvec(0, 2);
    vector[0] = 1; vector[1] = 2; vector[2] = 3;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"coordinate mapping" withValues:vector size:3 orUsingBlock:nil];
    free_ivector(vector, 0, 2);
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"simulation type" withValue:@"steady state"];
    //[listUtilities addStringInClassList:mod.simulation theVariable:@"simulation type" withValue:@"transient"];
    [listUtilities addStringInClassList:mod.simulation theVariable:@"time stepping method" withValue:@"bdf"];
    
    int ivalue = 1;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"bdf order" withValue:&ivalue orUsingBlock:nil];
    ivalue = 2;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"time step intervals" withValue:&ivalue orUsingBlock:nil];
    
    double value = 1.0;
    [listUtilities addConstRealInClassList:mod.simulation theVariable:@"time step size" withValue:&value orUsingBlock:nil string:nil];
    
    ivalue = 1;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"steady state maximum iterations" withValue:&ivalue orUsingBlock:nil];
    
    int *intervals = intvec(0, 0);
    intervals[0] = 1;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"output intervals" withValues:intervals size:1 orUsingBlock:nil];
    free_ivector(intervals, 0, 0);
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"output file" withValue:@"TempDist.result"];
    [listUtilities addStringInClassList:mod.simulation theVariable:@"post file" withValue:@"TempDist.ep"];
    
    double **gravity = doublematrix(0, 3, 0, 0);
    gravity[0][0] = 0.0;
    gravity[1][0] = -1.0;
    gravity[2][0] = 0.0;
    gravity[3][0] = 9.82;
    [listUtilities addConstRealArrayInClassList:mod.constants theVariable:@"gravity" withValues:gravity size1:4 size2:1 orUsingBlock:nil string:nil];
    
    free_dmatrix(gravity, 0, 3, 0, 0);
    value = 5.67e-08;
    [listUtilities addConstRealInClassList:mod.constants theVariable:@"stefan boltzmann" withValue:&value orUsingBlock:nil string:nil];
    
    mod.bodies = @[ @{ @"name" : @"body1",
                        @"body force" : @1,
                        @"equation" : @1,
                        @"material" : @1 }];
    mod.numberOfBodies = 1;
    
    
    FEMEquation *equation = [[FEMEquation alloc] init];
    [listUtilities addStringInClassList:equation theVariable:@"name" withValue:@"equation1"];
    [listUtilities addLogicalInClassList:equation theVariable:@"heat equation" withValue:YES];
    mod.equations = @[equation];
    mod.numberOfEquations = 1;
    
    FEMSolution *solution = [[FEMSolution alloc] init];
    [solution.solutionInfo setObject:@"always" forKey:@"invoke solution computer"];
    [solution.solutionInfo setObject:@"heat equation" forKey:@"equation"];
    [solution.solutionInfo setObject:@"temperature" forKey:@"variable"];
    [solution.solutionInfo setObject:@1 forKey:@"variable dofs"];
    [solution.solutionInfo setObject:@"iterative" forKey:@"linear system solver"];
    [solution.solutionInfo setObject:@"bi-cgstab" forKey:@"linear system iterative method"];
    [solution.solutionInfo setObject:@1000 forKey:@"linear system maximum iterations"];
    [solution.solutionInfo setObject:@1.0e-08 forKey:@"linear system convergence tolerance"];
    [solution.solutionInfo setObject:@YES forKey:@"linear system abort not converged"];
    [solution.solutionInfo setObject:@"ilu0" forKey:@"linear system preconditioning"];
    [solution.solutionInfo setObject:@1 forKey:@"linear system residual output"];
    [solution.solutionInfo setObject:@1.0e-05 forKey:@"steady state convergence tolerance"];
    [solution.solutionInfo setObject:@YES forKey:@"stabilize"];
    [solution.solutionInfo setObject:@1.0e-05 forKey:@"nonlinear system convergence tolerance"];
    [solution.solutionInfo setObject:@1 forKey:@"nonlinear system maximum iterations"];
    [solution.solutionInfo setObject:@3 forKey:@"nonlinear system newton after iterations"];
    [solution.solutionInfo setObject:@1.0e-02 forKey:@"nonlinear system newton after tolerance"];
    [solution.solutionInfo setObject:@1.0 forKey:@"nonlinear system relaxation factor"];
    [solution.solutionInfo setObject:@NO forKey:@"parallel assembly"];
    mod.solutions = @[solution];
    mod.numberOfSolutions = 1;
    
    FEMMaterial *material = [[FEMMaterial alloc] init];
    [listUtilities addStringInClassList:material theVariable:@"name" withValue:@"material1"];
    value = 1.0;
    [listUtilities addConstRealInClassList:material theVariable:@"density" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:material theVariable:@"heat conductivity" withValue:&value orUsingBlock:nil string:nil];
    
    mod.materials = @[material];
    mod.numberOfMaterials = 1;
    
    FEMBodyForce *bodyForce = [[FEMBodyForce alloc] init];
    [listUtilities addStringInClassList:bodyForce theVariable:@"name" withValue:@"bodyforce1"];
    [listUtilities addConstRealInClassList:bodyForce theVariable:@"heat source" withValue:&value orUsingBlock:nil string:nil];
    mod.bodyForces = @[bodyForce];
    mod.numberOfBodyForces = 1;
    
    mod.boundaryID = @[@1];
    
    FEMBoundaryCondition *boundaryCondition = [[FEMBoundaryCondition alloc] init];
    [listUtilities addStringInClassList:boundaryCondition theVariable:@"name" withValue:@"constraint1"];
    vector = intvec(0, 5);
    vector[0] = 1; vector[1] = 2; vector[2] = 3; vector[3] = 4; vector[4] = 5; vector[5] = 6;
    [listUtilities addIntegerArrayInClassList:boundaryCondition theVariable:@"target boundaries" withValues:vector size:6 orUsingBlock:nil];
    free_ivector(vector, 0, 5);
    value = 0.0;
    [listUtilities addConstRealInClassList:boundaryCondition theVariable:@"temperature" withValue:&value orUsingBlock:nil string:nil];
    
    boundaryCondition.tag = 1;
    mod.boundaryConditions = @[boundaryCondition];
    mod.numberOfBoundaryConditions = 1;
    
    mod.meshDir = [NSMutableString stringWithString:[self.path stringByAppendingPathComponent:@"heateq"]];
    mod.meshName = [NSMutableString stringWithString:[self.path stringByAppendingPathComponent:@"heateq/Mesh"]];
}

-(void)setUpStepStokesTest:(id __nonnull)model {
    
    FEMModel *mod = (FEMModel *)model;
    
    FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"coordinate system" withValue:@"cartesian 2d"];
    int ivalue = 9;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"maximum output level" withValue:&ivalue orUsingBlock:nil];
    
    int *vector = intvec(0, 2);
    vector[0] = 1; vector[1] = 2; vector[2] = 3;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"coordinate mapping" withValues:vector size:3 orUsingBlock:nil];
    free_ivector(vector, 0, 2);
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"simulation type" withValue:@"steady state"];
    ivalue = 1;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"steady state maximum iterations" withValue:&ivalue orUsingBlock:nil];
    
    int *intervals = intvec(0, 0);
    intervals[0] = 1;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"output intervals" withValues:intervals size:1 orUsingBlock:nil];
    free_ivector(intervals, 0, 0);
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"output file" withValue:@"Step.result"];
    [listUtilities addStringInClassList:mod.simulation theVariable:@"post file" withValue:@"Step.ep"];
    
    double **gravity = doublematrix(0, 3, 0, 0);
    gravity[0][0] = 0.0;
    gravity[1][0] = -1.0;
    gravity[2][0] = 0.0;
    gravity[3][0] = 9.82;
    [listUtilities addConstRealArrayInClassList:mod.constants theVariable:@"gravity" withValues:gravity size1:4 size2:1 orUsingBlock:nil string:nil];
    free_dmatrix(gravity, 0, 3, 0, 0);
    
    double value = 5.67e-08;
    [listUtilities addConstRealInClassList:mod.constants theVariable:@"stefan boltzmann" withValue:&value orUsingBlock:nil string:nil];
    
    mod.bodies = @[ @{@"name": @"body1",
                       @"equation" : @1,
                       @"material" : @1,
                       @"initial condition" : @1}];
    mod.numberOfBodies = 1;
    
    FEMInitialConditions *initialConditions = [[FEMInitialConditions alloc] init];
    value = 0.0;
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    mod.initialConditions = @[initialConditions];
    mod.numberOfInitialConditions = 1;
    
    FEMEquation *equation = [[FEMEquation alloc] init];
    [listUtilities addStringInClassList:equation theVariable:@"name" withValue:@"equation1"];
    vector = intvec(0, 0);
    vector[0] = 1;
    [listUtilities addIntegerArrayInClassList:equation theVariable:@"active solvers" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    [listUtilities addLogicalInClassList:equation theVariable:@"ns convect" withValue:NO];
    mod.equations = @[equation];
    mod.numberOfEquations = 1;
    
    FEMSolution *solution = [[FEMSolution alloc] init];
    [solution.solutionInfo setObject:@"navier-stokes" forKey:@"equation"];
    [solution.solutionInfo setObject:@"iterative" forKey:@"linear system solver"];
    [solution.solutionInfo setObject:@"bi-cgstab(l)" forKey:@"linear system iterative method"];
    [solution.solutionInfo setObject:@500 forKey:@"linear system maximum iterations"];
    [solution.solutionInfo setObject:@"stabilized" forKey:@"stabilization method"];
    [solution.solutionInfo setObject:@1.0e-08 forKey:@"linear system convergence tolerance"];
    [solution.solutionInfo setObject:@"ilu0" forKey:@"linear system preconditioning"];
    [solution.solutionInfo setObject:@1 forKey:@"linear system residual output"];
    [solution.solutionInfo setObject:@1.0e-05 forKey:@"steady state convergence tolerance"];
    [solution.solutionInfo setObject:@1 forKey:@"nonlinear system maximum iterations"];
    [solution.solutionInfo setObject:@1.0e-03 forKey:@"nonlinear system convergence tolerance"];
    [solution.solutionInfo setObject:@1.0e-03 forKey:@"nonlinear system newton after tolerance"];
    [solution.solutionInfo setObject:@3 forKey:@"nonlinear system newton after iterations"];
    mod.solutions = @[solution];
    mod.numberOfSolutions = 1;
    
    FEMMaterial *material = [[FEMMaterial alloc] init];
    [listUtilities addStringInClassList:material theVariable:@"name" withValue:@"material1"];
    value = 1.0;
    [listUtilities addConstRealInClassList:material theVariable:@"density" withValue:&value orUsingBlock:nil string:nil];
    value = 1.0e20;
    [listUtilities addConstRealInClassList:material theVariable:@"viscosity" withValue:&value orUsingBlock:nil string:nil];
    mod.materials = @[material];
    mod.numberOfMaterials = 1;
    
    mod.boundaryID = @[@1, @2, @3];
    
    FEMBoundaryCondition *boundaryCondition1 = [[FEMBoundaryCondition alloc] init];
    [listUtilities addStringInClassList:boundaryCondition1 theVariable:@"name" withValue:@"inflow"];
    vector = intvec(0, 0);
    vector[0] = 1;
    [listUtilities addIntegerArrayInClassList:boundaryCondition1 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
    NSArray *dependencies = @[@"coordinate 2"];
    double (^block) (double *) = ^(double *t) {
        return (2.0-t[0]) * (t[0]-1.0);
    };
    [listUtilities addBlockInClassList:boundaryCondition1 theVariable:@"velocity 1" usingBlock:block dependencies:dependencies];
    
    value = 0.0;
    [listUtilities addConstRealInClassList:boundaryCondition1 theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    
    FEMBoundaryCondition *boundaryCondition2 = [[FEMBoundaryCondition alloc] init];
    [listUtilities addStringInClassList:boundaryCondition2 theVariable:@"name" withValue:@"outflow"];
    vector = intvec(0, 0);
    vector[0] = 2;
    [listUtilities addIntegerArrayInClassList:boundaryCondition2 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    [listUtilities addConstRealInClassList:boundaryCondition2 theVariable:@"pressure" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:boundaryCondition2 theVariable:@"velocity 2"  withValue:&value orUsingBlock:nil string:nil];
    
    FEMBoundaryCondition *boundaryCondition3 = [[FEMBoundaryCondition alloc] init];
    [listUtilities addStringInClassList:boundaryCondition3 theVariable:@"name" withValue:@"wall"];
    vector = intvec(0, 0);
    vector[0] = 3;
    [listUtilities addIntegerArrayInClassList:boundaryCondition3 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    [listUtilities addConstRealInClassList:boundaryCondition3 theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:boundaryCondition3 theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    
    boundaryCondition1.tag = 1;
    boundaryCondition2.tag = 2;
    boundaryCondition3.tag = 3;
    mod.boundaryConditions = @[boundaryCondition1, boundaryCondition2, boundaryCondition3];
    mod.numberOfBoundaryConditions = 3;
    
    mod.meshDir = [NSMutableString stringWithString:[self.path stringByAppendingPathComponent:@"Step_stokes"]];
    mod.meshName = [NSMutableString stringWithString:[self.path stringByAppendingPathComponent:@"Step_stokes/Step"]];
}

-(void)setUpNaturalConvectionTest:(id __nonnull)model {
    
    FEMModel *mod = (FEMModel *)model;
    
    FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"coordinate system" withValue:@"cartesian"];
    int ivalue = 4;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"maximum output level" withValue:&ivalue orUsingBlock:nil];
    
    int *vector = intvec(0, 2);
    vector[0] = 1; vector[1] = 2; vector[2] = 3;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"coordinate mapping" withValues:vector size:3 orUsingBlock:nil];
    free_ivector(vector, 0, 2);
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"simulation type" withValue:@"transient"];
    ivalue = 20;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"steady state maximum iterations" withValue:&ivalue orUsingBlock:nil];
    [listUtilities addStringInClassList:mod.simulation theVariable:@"time stepping method" withValue:@"bdf"];
    ivalue = 2;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"bdf order" withValue:&ivalue orUsingBlock:nil];
    ivalue = 3;
    [listUtilities addIntegerInClassList:mod.simulation theVariable:@"time step intervals" withValue:&ivalue orUsingBlock:nil];
    double value = 100.0;
    [listUtilities addConstRealInClassList:mod.simulation theVariable:@"time step size" withValue:&value orUsingBlock:nil string:nil];
    
    int *intervals = intvec(0, 0);
    intervals[0] = 0;
    [listUtilities addIntegerArrayInClassList:mod.simulation theVariable:@"output intervals" withValues:intervals size:1 orUsingBlock:nil];
    free_ivector(intervals, 0, 0);
    
    [listUtilities addStringInClassList:mod.simulation theVariable:@"output file" withValue:@"case.result"];
    [listUtilities addStringInClassList:mod.simulation theVariable:@"post file" withValue:@"case.ep"];
    
    double **gravity = doublematrix(0, 3, 0, 0);
    gravity[0][0] = 0.0;
    gravity[1][0] = -1.0;
    gravity[2][0] = 0.0;
    gravity[3][0] = 9.82;
    [listUtilities addConstRealArrayInClassList:mod.constants theVariable:@"gravity" withValues:gravity size1:4 size2:1 orUsingBlock:nil string:nil];
    free_dmatrix(gravity, 0, 3, 0, 0);
    
    value = 5.67e-08;
    [listUtilities addConstRealInClassList:mod.constants theVariable:@"stefan boltzmann" withValue:&value orUsingBlock:nil string:nil];
    
    NSArray *bList = @[@1];
    mod.bodies = @[ @{@"target bodies" : bList,
                      @"name": @"body1",
                      @"equation" : @1,
                      @"material" : @1,
                      @"body force" : @1,
                      @"initial condition" : @1}];
    mod.numberOfBodies = 1;
    
    FEMSolution *solution1 = [[FEMSolution alloc] init];
    [solution1.solutionInfo setObject:@"heat equation" forKey:@"equation"];
    [solution1.solutionInfo setObject:@"-dofs 1 temperature" forKey:@"variable"];
    [solution1.solutionInfo setObject:@YES forKey:@"stabilize"];
    [solution1.solutionInfo setObject:@NO forKey:@"bubbles"];
    [solution1.solutionInfo setObject:@YES forKey:@"optimize bandwidth"];
    [solution1.solutionInfo setObject:@1.0e-05 forKey:@"steady state convergence tolerance"];
    [solution1.solutionInfo setObject:@YES forKey:@"linear system symmetric"];
    [solution1.solutionInfo setObject:@"iterative" forKey:@"linear system solver"];
    [solution1.solutionInfo setObject:@"gcr" forKey:@"linear system iterative method"];
    [solution1.solutionInfo setObject:@500 forKey:@"linear system maximum iterations"];
    [solution1.solutionInfo setObject:@1.0e-08 forKey:@"linear system convergence tolerance"];
    [solution1.solutionInfo setObject:@NO forKey:@"linear system abort not converged"];
    [solution1.solutionInfo setObject:@"ilu0" forKey:@"linear system preconditioning"];
    [solution1.solutionInfo setObject:@10 forKey:@"linear system residual output"];
    [solution1.solutionInfo setObject:@1.0e-08 forKey:@"nonlinear system convergence tolerance"];
    [solution1.solutionInfo setObject:@1 forKey:@"nonlinear system maximum iterations"];
    [solution1.solutionInfo setObject:@3 forKey:@"nonlinear system newton after iterations"];
    [solution1.solutionInfo setObject:@1.0e-03 forKey:@"nonlinear system newton after tolerance"];
    [solution1.solutionInfo setObject:@0.7 forKey:@"nonlinear system relaxation factor"];
    [solution1.solutionInfo setObject:@NO forKey:@"parallel assembly"];
    
    FEMSolution *solution2 = [[FEMSolution alloc] init];
    [solution2.solutionInfo setObject:@"navier-stokes" forKey:@"equation"];
    [solution2.solutionInfo setObject:@YES forKey:@"stabilize"];
    [solution2.solutionInfo setObject:@NO forKey:@"bubbles"];
    [solution2.solutionInfo setObject:@YES forKey:@"optimize bandwidth"];
    [solution2.solutionInfo setObject:@1.0e-05 forKey:@"steady state convergence tolerance"];
    [solution2.solutionInfo setObject:@YES forKey:@"linear system symmetric"];
    [solution2.solutionInfo setObject:@"iterative" forKey:@"linear system solver"];
    [solution2.solutionInfo setObject:@"bi-cgstab(l)" forKey:@"linear system iterative method"];
    [solution2.solutionInfo setObject:@500 forKey:@"linear system maximum iterations"];
    [solution2.solutionInfo setObject:@1.0e-10 forKey:@"linear system convergence tolerance"];
    [solution2.solutionInfo setObject:@"ilu0" forKey:@"linear system preconditioning"];
    [solution2.solutionInfo setObject:@NO forKey:@"linear system abort not converged"];
    [solution2.solutionInfo setObject:@10 forKey:@"linear system residual output"];
    [solution2.solutionInfo setObject:@1.0e-08 forKey:@"nonlinear system convergence tolerance"];
    [solution2.solutionInfo setObject:@1 forKey:@"nonlinear system maximum iterations"];
    [solution2.solutionInfo setObject:@3 forKey:@"nonlinear system newton after iterations"];
    [solution2.solutionInfo setObject:@1.0e-03 forKey:@"nonlinear system newton after tolerance"];
    [solution2.solutionInfo setObject:@0.7 forKey:@"nonlinear system relaxation factor"];

    mod.solutions = @[solution1, solution2];
    mod.numberOfSolutions = 2;
    
    FEMEquation *equation = [[FEMEquation alloc] init];
    [listUtilities addStringInClassList:equation theVariable:@"name" withValue:@"natural convection"];
    vector = intvec(0, 1);
    vector[0] = 1;
    vector[1] = 2;
    [listUtilities addIntegerArrayInClassList:equation theVariable:@"active solvers" withValues:vector size:2 orUsingBlock:nil];
    free_ivector(vector, 0, 1);
    [listUtilities addStringInClassList:equation theVariable:@"convection" withValue:@"computed"];
    mod.equations = @[equation];
    mod.numberOfEquations = 1;
    
    FEMMaterial *material = [[FEMMaterial alloc] init];
    [listUtilities addStringInClassList:material theVariable:@"name" withValue:@"water (room temperature)"];
    value = 998.3;
    [listUtilities addConstRealInClassList:material theVariable:@"density" withValue:&value orUsingBlock:nil string:nil];
    value = 1.002e-3;
    [listUtilities addConstRealInClassList:material theVariable:@"viscosity" withValue:&value orUsingBlock:nil string:nil];
    value = 0.58;
    [listUtilities addConstRealInClassList:material theVariable:@"heat conductivity" withValue:&value orUsingBlock:nil string:nil];
    value = 293.0;
    [listUtilities addConstRealInClassList:material theVariable:@"reference temperature" withValue:&value orUsingBlock:nil string:nil];
    value = 4183.0;
    [listUtilities addConstRealInClassList:material theVariable:@"heat capacity" withValue:&value orUsingBlock:nil string:nil];
    value = 0.207e-3;
    [listUtilities addConstRealInClassList:material theVariable:@"heat expansion coefficient" withValue:&value orUsingBlock:nil string:nil];
    mod.materials = @[material];
    mod.numberOfMaterials = 1;
    
    FEMBodyForce *bodyForce = [[FEMBodyForce alloc] init];
    [listUtilities addStringInClassList:bodyForce theVariable:@"name" withValue:@"bouancy"];
    [listUtilities addLogicalInClassList:bodyForce theVariable:@"boussinesq" withValue:YES];
    mod.bodyForces = @[bodyForce];
    mod.numberOfBodyForces = 1;
    
    FEMInitialConditions *initialConditions = [[FEMInitialConditions alloc] init];
    [listUtilities addStringInClassList:initialConditions theVariable:@"name" withValue:@"initial guess"];
    value = 1.0e-9;
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    value = 0.0;
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    value = 295.0;
    [listUtilities addConstRealInClassList:initialConditions theVariable:@"temperature" withValue:&value orUsingBlock:nil string:nil];
    mod.initialConditions = @[initialConditions];
    mod.numberOfInitialConditions = 1;
    
    mod.boundaryID = @[@1, @2, @3, @4];
    
    FEMBoundaryCondition *boundaryCondition1 = [[FEMBoundaryCondition alloc] init];
    [listUtilities addStringInClassList:boundaryCondition1 theVariable:@"name" withValue:@"bottom"];
    vector = intvec(0, 0);
    vector[0] = 1;
    [listUtilities addIntegerArrayInClassList:boundaryCondition1 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    value = 0.0;
    [listUtilities addConstRealInClassList:boundaryCondition1 theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:boundaryCondition1 theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    
    FEMBoundaryCondition *boundaryCondition2 = [[FEMBoundaryCondition alloc] init];
    [listUtilities addStringInClassList:boundaryCondition2 theVariable:@"name" withValue:@"right"];
    vector = intvec(0, 0);
    vector[0] = 2;
    [listUtilities addIntegerArrayInClassList:boundaryCondition2 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    value = 0.0;
    [listUtilities addConstRealInClassList:boundaryCondition2 theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:boundaryCondition2 theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    value = 297.0;
    [listUtilities addConstRealInClassList:boundaryCondition2 theVariable:@"temperature" withValue:&value orUsingBlock:nil string:nil];
    
    FEMBoundaryCondition *boundaryCondition3 = [[FEMBoundaryCondition alloc] init];
    [listUtilities addStringInClassList:boundaryCondition3 theVariable:@"name" withValue:@"top"];
    vector = intvec(0, 0);
    vector[0] = 3;
    [listUtilities addIntegerArrayInClassList:boundaryCondition3 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    value = 0.0;
    [listUtilities addConstRealInClassList:boundaryCondition3 theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    
    FEMBoundaryCondition *boundaryCondition4 = [[FEMBoundaryCondition alloc] init];
    [listUtilities addStringInClassList:boundaryCondition4 theVariable:@"name" withValue:@"left"];
    vector = intvec(0, 0);
    vector[0] = 4;
    [listUtilities addIntegerArrayInClassList:boundaryCondition4 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    value = 0.0;
    [listUtilities addConstRealInClassList:boundaryCondition4 theVariable:@"velocity 1" withValue:&value orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:boundaryCondition4 theVariable:@"velocity 2" withValue:&value orUsingBlock:nil string:nil];
    value = 293.0;
    [listUtilities addConstRealInClassList:boundaryCondition4 theVariable:@"temperature" withValue:&value orUsingBlock:nil string:nil];

    boundaryCondition1.tag = 1;
    boundaryCondition2.tag = 2;
    boundaryCondition3.tag = 3;
    boundaryCondition4.tag = 4;
    mod.boundaryConditions = @[boundaryCondition1, boundaryCondition2, boundaryCondition3, boundaryCondition4];
    mod.numberOfBoundaryConditions = 4;
    
    mod.meshDir = [NSMutableString stringWithString:[self.path stringByAppendingPathComponent:@"NaturalConvection"]];
    mod.meshName = [NSMutableString stringWithString:[self.path stringByAppendingPathComponent:@"NaturalConvection/square"]];
}

-(void)setUpISMIP_HOM_A010Test:(id __nonnull)model {
    
    FEMModel *mod = (FEMModel *)model;
    
    FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
    
    double L = 10.0e3;
    double slope = 0.5 * M_PI / 180.0;
    
    double yearinsec = 365.25 * 24.0 * 24.0 * 60.0 * 60.0;
    double rhoi = 900.0 / (1.0e6 * pow(yearinsec, 2.0));
    double gravity = -9.81 * pow(yearinsec, 2.0);
    double n = 3.0;
    double eta = pow((2.0 * 100.0), (-1.0/n));
    
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
    [solution1.solutionInfo setObject:@"MapCoordinate" forKey:@"equation"];
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
    [solution2.solutionInfo setObject:@"never" forKey:@"invoke solution computer"];
    
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
    
    FEMBoundaryCondition *boundaryCondition1 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 1;
    [listUtilities addIntegerArrayInClassList:boundaryCondition1 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
    FEMBoundaryCondition *boundaryCondition2 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 2;
    [listUtilities addIntegerArrayInClassList:boundaryCondition2 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
    FEMBoundaryCondition *boundaryCondition3 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 3;
    [listUtilities addIntegerArrayInClassList:boundaryCondition3 theVariable:@"target boudnaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
    FEMBoundaryCondition *boundaryCondition4 = [[FEMBoundaryCondition alloc] init];
    vector = intvec(0, 0);
    vector[0] = 4;
    [listUtilities addIntegerArrayInClassList:boundaryCondition4 theVariable:@"target boundaries" withValues:vector size:1 orUsingBlock:nil];
    free_ivector(vector, 0, 0);
    
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
    
    mod.meshDir = [NSMutableString stringWithString:[self.path stringByAppendingPathComponent:@"ISMIP-HOM-A010"]];
    mod.meshName = [NSMutableString stringWithString:[self.path stringByAppendingPathComponent:@"ISMIP-HOM-A010/rectangle"]];
}

@end
