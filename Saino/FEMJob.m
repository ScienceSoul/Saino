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
-(void)FEMJob_addSolutions:(FEMModel *)model;
@end

@implementation FEMJob {
    
    BOOL _transient;
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

-(void)FEMJob_addSolutions:(FEMModel *)model {
    
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

@end
