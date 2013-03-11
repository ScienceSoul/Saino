//
//  FEMJob.m
//  Saino
//
//  Created by Seddik hakime on 09/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMJob.h"
#import "FEMSolution.h"
#import "FEMListUtilities.h"
#import "FEMEquation.h"
#import "FEMUtilities.h"
#import "memory.h"

@interface FEMJob ()
-(void)FEMJob_addSolutionsModel:(FEMModel *)model;
-(void)FEMJob_addMeshCoordinatesAndTimeModel:(FEMModel *)model;
@end

@implementation FEMJob {
    
    BOOL _transient;
    double *_sTime;
    double *_sStep;
    double *_sInterval;
    double *_sSize;
    double **_sPrevSizes;
    double *_steadyIt;
    double *_nonLinIt;
    int _sizeSTime;
    int _sizeSStep;
    int _sizeSInterval;
    int _sizeSSize;
    int _sizeSteadyIt;
    int _sizeNonLinIt;
    int _size1SPrevSizes;
    int _size22SPrevSizes;
}

@synthesize model = _model;

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _transient = NO;
    }
    
    return self;
}

-(void)deallocation {

}

/***********************************
    Add flags for active sokutions
***********************************/
-(void)FEMJob_addSolutionsModel:(FEMModel *)model {
    
    int i, j;
    BOOL initSolution, found;
    NSString *eq;
    NSNumber *aNumber;
    FEMListUtilities *listUtilities;
    FEMUtilities *utilities;
    listBuffer activeSolvers = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    listUtilities = [[FEMListUtilities alloc] init];
    
    i = 1;
    for (FEMSolution *solution in model.solutions) {
        
        // This is a hack that sets Equation flags true for the active solvers.
        // The equation flag is the legacy way of setting a solution active and
        // is still used internally.
        if ((solution.solutionInfo)[@"equation"] != nil) {
            eq = (solution.solutionInfo)[@"equation"];
            for (FEMEquation *equation in model.equations) {
                found = [listUtilities listGetIntegerArray:model inArray:equation.valuesList forVariable:@"active solvers" buffer:&activeSolvers];
                if (found == YES) {
                    for (j=0; j<activeSolvers.m; j++) {
                        if (activeSolvers.ivector[j] == i) {
                            [listUtilities addLogicalInClassList:equation.valuesList theVariable:eq withValue:YES];
                            break;
                        }
                    }
                }
                free_ivector(activeSolvers.ivector, 0, activeSolvers.m-1);
            }
        }
        i++;
    }
    
    utilities = [[FEMUtilities alloc] init];
    initSolution = NO;
    for (FEMSolution *solution in model.solutions) {
        if ((solution.solutionInfo)[@"equation"] != nil) {
            eq = (solution.solutionInfo)[@"equation"];
        }
        if ((solution.solutionInfo)[@"initialize"] != nil) {
            if ([(solution.solutionInfo)[@"initialize"] boolValue] == YES) {
                [solution.matrix deallocation];
                solution.matrix = nil;
                aNumber = @YES;
                [solution.solutionInfo setObject:aNumber forKey:@"initialize"];
            }
        }
        
        if (solution.plugInPrincipalClassInstance == nil || initSolution == YES) {
            //TODO: Make sure that this is alsways correct. Here if the solution has no mesh
            // we assigned it to the first mesh listed in model.meshes
            if (solution.mesh == nil) solution.mesh = model.meshes[0];
        }
        [utilities addEquationBasicsToSolution:solution name:eq model:model transient:_transient];
        [utilities addEquationToSolution:solution model:model transient:_transient];
    }
}

/*******************************************************************
    Add coordinate and time variables to the meshe(s) in the model
*******************************************************************/
-(void)FEMJob_addMeshCoordinatesAndTimeModel:(FEMModel *)model {
    
    FEMVariable *dtVar;
    FEMSolution *solution;
    FEMUtilities *utilities;
    Nodes_t *nodes;
    variableArraysContainer *varContainers = NULL, *dtVarContainers = NULL;
    BOOL found;
    
    utilities = [[FEMUtilities alloc] init];
    
    varContainers = allocateVariableContainer();
    solution = nil;
    for (FEMMesh *mesh in model.meshes) {
        nodes = mesh.getNodes;
        
        varContainers->Values = nodes->x;
        varContainers->sizeValues = nodes->numberOfNodes;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"coordinate 1" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        
        varContainers->Values = nodes->y;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"coordinate 2" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        
        varContainers->Values = nodes->z;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"coordinate 3" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        
        varContainers->Values = _sTime;
        varContainers->sizeValues = _sizeSTime;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"time" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        
        varContainers->Values = _sStep;
        varContainers->sizeValues = _sizeSStep;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"time step" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        
        varContainers->Values = _sSize;
        varContainers->sizeValues = _sizeSSize;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"time step size" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        
        varContainers->Values = _sInterval;
        varContainers->sizeValues = _sizeSInterval;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"time step interval" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        
        // Save some previous time steps for variable time step multistep methods
        dtVar = [utilities getVariableFrom:mesh.variables model:model name:@"time step size" onlySearch:NULL maskName:NULL info:&found];
        dtVarContainers = dtVar.getContainers;
        dtVarContainers->PrevValues = _sPrevSizes;
        dtVarContainers->size1PrevValues = _size1SPrevSizes;
        dtVarContainers->size2PrevValues = _size22SPrevSizes;
        
        varContainers->Values = _nonLinIt;
        varContainers->sizeValues = _sizeNonLinIt;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"nonlin iter" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        
        varContainers->Values = _steadyIt;
        varContainers->sizeValues = _sizeSteadyIt;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"coupled iter" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
    }
    free(varContainers);
}

@end
