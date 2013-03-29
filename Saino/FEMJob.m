//
//  FEMJob.m
//  Saino
//
//  Created by Seddik hakime on 09/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMJob.h"
#import "FEMKernel.h"
#import "FEMSolution.h"
#import "FEMListUtilities.h"
#import "FEMEquation.h"
#import "FEMUtilities.h"
#import "FEMMeshUtils.h"
#import "FEMElementDescription.h"
#import "FEMElementUtils.h"
#import "FEMInitialConditions.h"
#import "Utils.h"
#import "TimeProfile.h"

@interface FEMJob ()
-(void)FEMJob_addSolutionsModel:(FEMModel *)model;
-(void)FEMJob_addMeshCoordinatesAndTimeModel:(FEMModel *)model;
-(void)FEMJob_setInitialConditionsModel:(FEMModel *)model;
-(void)FEMJob_initCondModel:(FEMModel *)model;
-(void)FEMJob_restartModel:(FEMModel *)model;
-(void)FEMJob_runSimulation:(FEMModel *)model timeIntervals:(int)timeIntervals coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter outputIntervals:(int* )outputIntervals transient:(BOOL)transient scanning:(BOOL)scanning;
-(void)FEMJob_saveToPostModel:(FEMModel *)model currentStep:(int)currentStep;
-(void)FEMJob_saveCurrent:(FEMModel *)model currentStep:(int)currentStep;
@end

@implementation FEMJob {
    
    BOOL _transient;
    BOOL _initDirichlet;
    BOOL _lastSaved;
    int _totalTimeSteps;
    int *_timeSteps;
    double _s;
    double *_sTime;
    double *_sStep;
    double *_sInterval;
    double *_sSize;
    double **_sPrevSizes;
    double *_steadyIt;
    double *_nonLinIt;
    double **_timeStepSizes;
    int _sizeSTime;
    int _sizeSStep;
    int _sizeSInterval;
    int _sizeSSize;
    int _sizeSteadyIt;
    int _sizeNonLinIt;
    int _size1SPrevSizes;
    int _size2SPrevSizes;
    
    NSMutableString *_outputName;
    NSMutableString *_postName;
    NSString *_when;
}

@synthesize model = _model;

#pragma mark Private methods

/***********************************
    Add flags for active solutions
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
        
        if (solution.plugInPrincipalClassInstance == nil || solution.selector == NULL || initSolution == YES) {
            //TODO: Make sure that this is alsways correct. Here if the solution has no mesh
            // we assigned it to the first mesh listed in model.meshes
            if (solution.mesh == nil) solution.mesh = model.meshes[0];
            model.solution = solution;
            [utilities addEquationBasicsToSolution:solution name:eq model:model transient:_transient];
            [utilities addEquationToSolution:solution model:model transient:_transient];
        }
    }
}

/*******************************************************************
    Add coordinate and time variables to the meshe(s) in the model
*******************************************************************/
-(void)FEMJob_addMeshCoordinatesAndTimeModel:(FEMModel *)model {
    
    FEMVariable *dtVar = nil;
    FEMSolution *solution = nil;
    FEMUtilities *utilities;
    Nodes_t *nodes;
    variableArraysContainer *varContainers = NULL, *dtVarContainers = NULL;
    BOOL found;
    
    utilities = [[FEMUtilities alloc] init];
    
    varContainers = allocateVariableContainer();
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
        dtVarContainers->size2PrevValues = _size2SPrevSizes;
        
        varContainers->Values = _nonLinIt;
        varContainers->sizeValues = _sizeNonLinIt;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"nonlin iter" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        
        varContainers->Values = _steadyIt;
        varContainers->sizeValues = _sizeSteadyIt;
        [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"coupled iter" dofs:1 container:varContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
    }
    free(varContainers);
}

-(void)FEMJob_initCondModel:(FEMModel *)model {
    
    int i, j, k, k1, l, n, t, dofs;
    int *indexes;
    double integral;
    BOOL found;
    NSString *str = nil;
    NSMutableString *varName;
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0};
    FEMKernel *kernel;
    FEMSolution *solution;
    FEMInitialConditions *initialCondition;
    FEMListUtilities *listUtilities;
    FEMMeshUtils *meshUtilities;
    Element_t *elements = NULL, *edges = NULL;
    variableArraysContainer *variableContainers = NULL;
    
    kernel = [FEMKernel sharedKernel];
    listUtilities = [[FEMListUtilities alloc] init];
    meshUtilities = [[FEMMeshUtils alloc] init];
    
    for (i=1; i<=model.numberOfBodies; i++) {
        for (FEMMesh *mesh in model.meshes) {
            indexes = intvec(0, mesh.maxElementDofs-1);
            [meshUtilities setCurrentMesh:mesh inModel:model];
            
            elements = mesh.getElements;
            for (t=0; t<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; t++) {
                if (elements[t].BodyID == i) {
                    if ((model.bodies)[i-1][@"initial condition"] != nil) {
                        j = [(model.bodies)[i-1][@"initial condition"] intValue];
                        if (j < 1) j = 1;
                        if (j > model.numberOfInitialConditions) j = model.numberOfInitialConditions;
                        
                        n = [kernel getNumberOfNodesForElement:&elements[t]];
                        for (FEMVariable *variable in mesh.variables) {
                            variableContainers = variable.getContainers;
                            if (variable.solution != nil) {
                                dofs = [kernel getElementDofsSolution:solution model:model forElement:&elements[t] atIndexes:indexes];
                            } else {
                                dofs = [kernel getElementDofsSolution:nil model:model forElement:&elements[t] atIndexes:indexes];
                            }
                            
                            solution = (FEMSolution *)variable.solution;
                            if (solution == nil) solution = (FEMSolution *)model.solution;
                            
                            if ((solution.solutionInfo)[@"namespace"] != nil) {
                                str = (solution.solutionInfo)[@"namespace"];
                                [listUtilities listSetNameSpace:str];
                            }
                            
                            initialCondition = (model.initialConditions)[j-1];
                            if (variableContainers->sizeValues == variable.dofs) {
                                for (k=0; k<variableContainers->sizeValues; k++) {
                                    variableContainers->Values[k] = [listUtilities listGetConstReal:model inArray:initialCondition.valuesList forVariable:variable.name info:&found minValue:NULL maxValue:NULL];
                                }
                            } else if (variable.dofs <= 1) {
                                [listUtilities getReal:model inArray:initialCondition.valuesList forVariable:variable.name element:&elements[t] buffer:&work info:&found];
                                if (found == YES) {
                                    for (k=0; k<n; k++) {
                                        k1 = indexes[k];
                                        if (variableContainers->Perm !=NULL) k1 = variableContainers->Perm[k1];
                                        if (k1 >= 0) variableContainers->Values[k1] = work.vector[k];
                                    }
                                    free_dvector(work.vector, 0, work.m-1);
                                }
                                
                                if (_transient == YES && solution.timeOrder == 2) {
                                    varName = [NSMutableString stringWithString:variable.name];
                                    [varName appendString:@" velocity"];
                                    [listUtilities getReal:model inArray:initialCondition.valuesList forVariable:varName element:&elements[t] buffer:&work info:&found];
                                    if (found == YES) {
                                        for (k=0; k<n; k++) {
                                            k1 = indexes[k];
                                            if (variableContainers->Perm != NULL) k1 = variableContainers->Perm[k1];
                                            if (k1 >= 0) variableContainers->PrevValues[k1][0] = work.vector[k];
                                        }
                                        free_dvector(work.vector, 0, work.m-1);
                                    }
                                    varName = [NSMutableString stringWithString:variable.name];
                                    [varName appendString:@" acceleration"];
                                    [listUtilities getReal:model inArray:initialCondition.valuesList forVariable:varName element:&elements[t] buffer:&work info:&found];
                                    if (found == YES) {
                                        for (k=0; k<n; k++) {
                                            k1 = indexes[k];
                                            if (variableContainers->Perm != NULL) k1 = variableContainers->Perm[k1];
                                            if (k1 >= 0) variableContainers->PrevValues[k1][1] = work.vector[k];
                                        }
                                         free_dvector(work.vector, 0, work.m-1);
                                    }
                                }
                                
                                edges = mesh.getEdges;
                                if (edges != NULL) {
                                    if (i <= mesh.numberOfBulkElements) {
                                        varName = [NSMutableString stringWithString:variable.name];
                                        [varName appendString:@" {e}"];
                                        found = [listUtilities listCheckPresentVariable:varName inArray:initialCondition.valuesList];
                                        if (found == YES) {
                                            for (k=0; k<elements[t].Type.NumberOfEdges; k++) {
                                                l = variableContainers->Perm[elements[t].EdgeIndexes[k]+mesh.numberOfNodes];
                                                if (l >= 0) {
                                                    [kernel localBoundaryIntegral:model inSolution:solution atBoundary:initialCondition.valuesList forElement:&edges[elements[t].EdgeIndexes[k]] withNumberOfNodes:edges[elements[t].EdgeIndexes[k]].Type.NumberOfNodes andParent:&elements[t] withNumberOfNodes:n boundaryName:varName functionIntegral:&integral];
                                                    variableContainers->Values[l] = integral;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            else {
                                found = [listUtilities listGetRealArray:model inArray:initialCondition.valuesList forVariable:variable.name numberOfNodes:n indexes:elements[t].NodeIndexes buffer:&work];
                                if (found == YES) {
                                    for (k=0; k<n; k++) {
                                        k1 = indexes[k];
                                        for (l=0; l<min(work.m, variable.dofs); l++) {
                                            if (variableContainers->Perm != NULL) k1 = variableContainers->Perm[k1];
                                            if (k1 >= 0) variableContainers->Values[variable.dofs*k1+l] = work.tensor[l][0][k];
                                        }
                                    }
                                    free_d3tensor(work.tensor, 0, work.m-1, 0, work.n-1, 0, work.p-1);
                                }
                            }
                            [listUtilities listSetNameSpace:@""];
                        }
                    }
                }
            }
            free_ivector(indexes, 0, mesh.maxElementDofs-1);
        }
    }
}

/*************************************************************
    Check if we are restarting, and if yes read field values
*************************************************************/
-(void)FEMJob_restartModel:(FEMModel *)model {
    
    int k, minVal;
    double startTime;
    BOOL found;
    NSString *restartFile;
    FEMVariable *var = nil;
    FEMListUtilities *listUtilities;
    FEMMeshUtils *meshUtilities;
    FEMUtilities *utilities;
    variableArraysContainer *varContainers = NULL;
    
    listUtilities = [[FEMListUtilities alloc] init];
    
    restartFile = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"restart file" info:&found];
    
    if (found == YES) {
        minVal = 0;
        k = [listUtilities listGetInteger:model inArray:model.simulation.valuesList forVariable:@"restart position" info:&found minValue:&minVal maxValue:NULL];
        
        for (FEMMesh *mesh in model.meshes) {
            
            if (mesh.name.length > 0) {
                [_outputName setString:mesh.name];
                [_outputName appendString:@"/"];
                [_outputName appendString:restartFile];
            } else {
                [_outputName setString:restartFile];
            }
            //TODO: add support for parallel run
            
            meshUtilities = [[FEMMeshUtils alloc] init];
            [meshUtilities setCurrentMesh:mesh inModel:model];
            //TODO: this method is incomplete (and so doesn't work), we need to load the restart file here...
            
            startTime = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"restart time" info:&found minValue:NULL maxValue:NULL];
            if (found == YES) {
                utilities = [[FEMUtilities alloc] init];
                var = [utilities getVariableFrom:mesh.variables model:model name:@"time" onlySearch:NULL maskName:NULL info:&found];
                if (var != nil) {
                    varContainers = var.getContainers;
                    varContainers->Values[0] = startTime;
                }
            }
        }
    }
}

/******************************************
    Set initial conditions for the fields
******************************************/
-(void)FEMJob_setInitialConditionsModel:(FEMModel *)model {
    
    int i, j, k, l, m, n, t, dim, vectDof, realDof;
    double udot, parU, parV, *nrm, *t1, *t2, *vec, *tmp;
    BOOL found, ntBoundary, pointed, check;
    listBuffer work = { NULL, NULL, NULL, NULL, 0, 0, 0};
    NSArray *bc;
    NSString *str = nil;
    NSMutableString *varName;
    FEMVariable *vectVariable = nil;
    FEMSolution *solution = nil;
    FEMKernel *kernel;
    FEMListUtilities *listUtilities;
    FEMMeshUtils *meshUtilities;
    FEMElementDescription *elementDescription;
    FEMElementUtils *elementUtils;
    Element_t *elements = NULL;
    Nodes_t *nodes = NULL;
    variableArraysContainer *variableContainers = NULL, *vectVarContainers = NULL;
    
    dim = model.dimension;
    
    kernel = [FEMKernel sharedKernel];

    listUtilities = [[FEMListUtilities alloc] init];
    if ([listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"restart before initial conditions" info:&found] == YES) {
        [self FEMJob_restartModel:model];
        [self FEMJob_initCondModel:model];
    } else {
        [self FEMJob_initCondModel:model];
        [self FEMJob_restartModel:model];
    }
    
    // Make sure that initial values at boundaries are set correctly
    // NOTE: This overrides the inital condition settings for field variables!!!
    _initDirichlet = [listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"initialize dirichlet conditions" info:&found];
    if (found == NO) _initDirichlet = YES;
    
    elementDescription = [[FEMElementDescription alloc] init];
    nrm = doublevec(0, 2);
    t1 = doublevec(0, 2);
    t2 = doublevec(0, 2);
    vec = doublevec(0, 2);
    tmp = doublevec(0, 2);
    
    if (_initDirichlet == YES) {
        meshUtilities = [[FEMMeshUtils alloc] init];
        for (FEMMesh *mesh in model.meshes) {
            [meshUtilities setCurrentMesh:mesh inModel:model];
            
            elements = mesh.getElements;
            for (t=mesh.numberOfBulkElements; t<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; t++) {
                n = elements[t].Type.NumberOfNodes;
                bc = [kernel getBoundaryCondition:model forElement:&elements[t]];
                
                for (FEMVariable *variable in mesh.variables) {
                    variableContainers = variable.getContainers;
                    solution = (FEMSolution *)variable.solution;
                    if (solution == nil) solution = (FEMSolution *)model.solution;
                    
                    if ((solution.solutionInfo)[@"namespace"] != nil) {
                        str = (solution.solutionInfo)[@"namespace"];
                        [listUtilities listSetNameSpace:str];
                    }
                    
                    if (variable.dofs <= 1) {
                        [listUtilities getReal:model inArray:bc forVariable:variable.name element:&elements[t] buffer:&work info:&found];
                        if (found == YES) {
                            ntBoundary = NO;
                            if ([kernel getElementFamily:&elements[t]] != 1) {
                                k = (int)(variable.name.length);
                                vectDof = [variable.name characterAtIndex:k-1] - '0';
                                if (vectDof >= 1 && vectDof <= 3) {
                                    varName = [NSMutableString stringWithString:@"normal-tangential "];
                                    [varName appendString:[variable.name substringToIndex:k-1]];
                                    ntBoundary = [listUtilities listGetLogical:model inArray:bc forVariable:varName info:&found];
                                    
                                    if (ntBoundary == YES) {
                                        ntBoundary = NO;
                                        for (FEMVariable *vectVar in mesh.variables) {
                                            vectVarContainers = vectVar.getContainers;
                                            if (vectVar.dofs >= dim) {
                                                for (realDof=0; realDof<vectVar.dofs; realDof++) {
                                                    k = 0;
                                                    for (j=realDof; j<vectVarContainers->sizeValues; j+=vectVar.dofs) {
                                                        if (variableContainers->ComponentValues[k] == &vectVarContainers->Values[j]) {
                                                            pointed = YES;
                                                        } else pointed = NO;
                                                         k++;
                                                    }
                                                    if (pointed == YES) {
                                                        ntBoundary = YES;
                                                        break;
                                                    }
                                                }
                                                if (ntBoundary == YES) {
                                                    vectVariable = vectVar;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                    
                                    if (ntBoundary == YES) {
                                        nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
                                        initNodes(nodes);
                                        [kernel getNodes:solution model:model inElement:&elements[t] resultNodes:nodes numberOfNodes:NULL];
                                        parU = 0.0;
                                        parV = 0.0;
                                        check = YES;
                                        [elementDescription normalVectorForBDElement:&elements[t] boundaryNodes:nodes mesh:mesh paraU:&parU paraV:&parV check:&check normals:nrm];
                                        switch (elements[t].Type.dimension) {
                                            case 1:
                                                t1[0] = nrm[1];
                                                t1[1] = -nrm[0];
                                                break;
                                            case 2:
                                                elementUtils = [[FEMElementUtils alloc] init];
                                                [elementUtils tangentDirectionsForNormal:nrm tangent1:t1 tangent2:t2];
                                                break;
                                        }
                                        
                                        switch (vectDof) {
                                            case 1:
                                                memcpy(vec, nrm, 3);
                                                break;
                                            case 2:
                                                memcpy(vec, t1, 3);
                                                break;
                                            case 3:
                                                memcpy(vec, t2, 3);
                                                break;
                                        }
                                        free_dvector(nodes->x, 0, nodes->numberOfNodes-1);
                                        free_dvector(nodes->y, 0, nodes->numberOfNodes-1);
                                        free_dvector(nodes->z, 0, nodes->numberOfNodes-1);
                                        free(nodes);
                                    }
                                }
                            }
                            
                            for (j=0; j<n; j++) {
                                k = elements[t].NodeIndexes[j];
                                if (variableContainers->Perm != NULL) k = variableContainers->Perm[k];
                                if (k >= 0) {
                                    if (ntBoundary == YES) {
                                        for (l=0; l<dim; l++) {
                                            m = l + realDof - (--vectDof);
                                            tmp[l] = vectVarContainers->Values[vectVariable.dofs*k+m];
                                        }
                                        udot = 0.0;
                                        for (i=0; i<dim; i++) {
                                            udot = udot + (vec[i]*tmp[i]);
                                        }
                                        for (i=0; i<dim; i++) {
                                            tmp[i] = tmp[i]+(work.vector[j]-udot)*vec[i];
                                        }
                                        for (l=0; l<dim; l++) {
                                            m = l + realDof - (--vectDof);
                                            vectVarContainers->Values[vectVariable.dofs*k+m] = tmp[l];
                                        }
                                    } else {
                                        variableContainers->Values[k] = work.vector[j];
                                    }
                                }
                            }
                            free_dvector(work.vector, 0, work.m-1);
                        }
                        
                        if (_transient == YES && solution.timeOrder == 2) {
                            varName = [NSMutableString stringWithString:variable.name];
                            [varName appendString:@" velocity"];
                            [listUtilities getReal:model inArray:bc forVariable:varName element:&elements[t] buffer:&work info:&found];
                            if (found == YES) {
                                for (j=0; j<n; j++) {
                                    k = elements[t].NodeIndexes[j];
                                    if (variableContainers->Perm != NULL) k = variableContainers->Perm[k];
                                    if (k >= 0) variableContainers->PrevValues[k][0] = work.vector[j];
                                }
                                free_dvector(work.vector, 0, work.m-1);
                            }
                            
                            varName = [NSMutableString stringWithString:variable.name];
                            [varName appendString:@" acceleration"];
                            [listUtilities getReal:model inArray:bc forVariable:varName element:&elements[t] buffer:&work info:&found];
                            if (found == YES) {
                                for (j=0; j<n; j++) {
                                    k = elements[t].NodeIndexes[j];
                                    if (variableContainers->Perm != NULL) k = variableContainers->Perm[k];
                                    if (k >= 0) variableContainers->PrevValues[k][1] = work.vector[j];
                                }
                                free_dvector(work.vector, 0, work.m-1);
                            }
                        }
                    } else {
                        found = [listUtilities listGetRealArray:model inArray:bc forVariable:variable.name numberOfNodes:n indexes:elements[t].NodeIndexes buffer:&work];
                        if (found == YES) {
                            for (j=0; j<n; j++) {
                                k = elements[t].NodeIndexes[j];
                                for (l=0; l<min(work.m, variable.dofs); l++) {
                                    if (variableContainers->Perm != NULL) k = variableContainers->Perm[k];
                                    if (k >= 0) variableContainers->Values[variable.dofs*k+l] = work.tensor[l][0][j];
                                }
                            }
                            free_d3tensor(work.tensor, 0, work.m-1, 0, work.n-1, 0, work.p-1);
                        }
                    }
                    
                    [listUtilities listSetNameSpace:@""];
                }
            }
        }
    }
    
    [elementDescription deallocation];
    free_dvector(nrm, 0, 2);
    free_dvector(t1, 0, 2);
    free_dvector(t2, 0, 2);
    free_dvector(vec, 0, 2);
    free_dvector(tmp, 0, 2);
}

-(void)FEMJob_runSimulation:(FEMModel *)model timeIntervals:(int)timeIntervals coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter outputIntervals:(int* )outputIntervals transient:(BOOL)transient scanning:(BOOL)scanning{
    
    int i, j, k, jj, kk, n, interval, timeStep, stepCount = 0, cumTimeStep, realTimeStep, timeLeft, adaptiveKeepSmallest, minVal, smallestCount,
        stepControl=-1;
    double dt = 0, ddt, dtFunc, newTime, maxTime, prevTime=0, adaptiveLimit, adaptiveMaxTimeStep, adaptiveMinTimeStep, cumTime, maxErr,
           exitCond, arg;
    double **xx, *xxnrm, *yynrm, ***prevxx;
    BOOL found, adaptiveTime, steadyStateReached=NO, execThis;
    FEMKernel *kernel;
    FEMListUtilities *listUtilities;
    variableArraysContainer *varContainers = NULL;
    
    kernel = [FEMKernel sharedKernel];
    listUtilities = [[FEMListUtilities alloc] init];
    
    for (FEMSolution *solution in model.solutions) {
        if (solution.selector == NULL || solution.plugInPrincipalClassInstance == nil) continue;
        if (solution.solutionSolveWhen == SOLUTION_SOLVE_AHEAD_ALL) {
            [kernel activateSolution:solution model:model timeStep:dt transientSimulation:transient];
        }
    }
    
    for (interval=0; interval<timeIntervals; interval++) {
        stepCount = stepCount + _timeSteps[interval];
    }
    
    cumTimeStep = 0;
    ddt = 0.0;
    
    for (interval=1; interval<=timeIntervals; interval++) {
        
        if (transient == YES || scanning == YES) {
            dt = _timeStepSizes[interval-1][0];
        } else {
            dt = 1;
        }
        
        // Go through number of time steps within an interval
        realTimeStep = 1;
        for (timeStep=1; timeStep<=_timeSteps[interval-1]; timeStep++) {
            
            cumTimeStep++;;
            _sStep[0] = cumTimeStep;
            
            dtFunc = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"time step size" info:&found minValue:NULL maxValue:NULL];
            
            if (found == YES) dt = dtFunc;
            
            _sTime[0] = _sTime[0] + dt;
            
            // Move the old time steps one step down the ladder
            if (timeStep > 1 || interval > 1) {
                for (i=_size2SPrevSizes-1; i>=1; i--) {
                    _sPrevSizes[0][i] = _sPrevSizes[0][i-1];
                }
                _sPrevSizes[0][0] = _sSize[0];
            }
            _sSize[0] = dt;
            
            _sInterval[0] = interval;
            if (transient == NO) _steadyIt[0]++;
            
            // TODO: When we will support parallel run, we need to be sure that the code below
            // is only executed by only one processor
            NSLog(@"JOB: \n");
            NSLog(@"JOB: ----------------------------------------------------------------------\n");
            
            if (transient == YES || scanning == YES) {
                NSLog(@"JOB: Time: %d / %d %f\n", cumTimeStep, stepCount, _sTime[0]);
                
                newTime = realtime();
                
                if (cumTimeStep > 1) {
                    maxTime = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"real time max" info:&found minValue:NULL maxValue:NULL];
                    if (found == YES) {
                        NSLog(@"JOB: Fraction of real time left: %f\n", 1.0-realtime()/maxTime);
                    } else {
                        timeLeft = round((stepCount-(cumTimeStep-1))*(newTime-prevTime)/60.0);
                        if (timeLeft > 120) {
                            NSLog(@"JOB: Estimated time left: %d hours.\n", timeLeft/60);
                        } else if (timeLeft > 60) {
                            NSLog(@"JOB: Estimated time left: 1 hour %d minutes.\n", timeLeft % 60);
                        } else if (timeLeft >= 1) {
                            NSLog(@"JOB: Estimated time left: %d minutes.\n", timeLeft);
                        } else {
                            NSLog(@"JOB: Estimated time left: less than a minute.\n");
                        }
                    }
                }
                prevTime = newTime;
            } else {
                NSLog(@"JOB: Steady state iteration: %d\n", cumTimeStep);
            }
            
            NSLog(@"JOB: ----------------------------------------------------------------------\n");
            NSLog(@"JOB: \n"); // End of code that should only be executed by one processor if parallel run
            
            // Solve any and all governing equations in the system
            adaptiveTime = [listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"adaptive time stepping" info:&found];
            
            if (transient == YES && adaptiveTime == YES) { // Adaptive time stepping
                adaptiveLimit = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"adaptive time error" info:&found minValue:NULL maxValue:NULL];
                
                if (found == NO) {
                    errorfunct("FEMJob_runSimulation", "Adaptive time limit must be be given for adaptive stepping scheme.");
                }
                
                adaptiveMaxTimeStep = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"adaptive max time step" info:&found minValue:NULL maxValue:NULL];
                if (found == NO) adaptiveMaxTimeStep = dt;
                adaptiveMaxTimeStep = min(adaptiveMaxTimeStep, dt);
                
                adaptiveMinTimeStep = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"adaptive min time step" info:&found minValue:NULL maxValue:NULL];
                
                minVal = 0;
                adaptiveKeepSmallest = [listUtilities listGetInteger:model inArray:model.simulation.valuesList forVariable:@"adaptive keep smallest" info:&found minValue:&minVal maxValue:NULL];
                
                jj = 0;
                kk = 0;
                for (FEMSolution *solution in model.solutions) {
                    varContainers = solution.variable.getContainers;
                    if (varContainers->Values != NULL) {
                        if (varContainers->PrevValues != NULL) {
                            jj = max(jj, varContainers->size2PrevValues);
                        }
                        kk = max(kk, varContainers->sizeValues);
                    }
                }
                xx = doublematrix(0, model.numberOfSolutions-1, 0, kk-1);
                yynrm = doublevec(0, model.numberOfSolutions-1);
                xxnrm = doublevec(0, model.numberOfSolutions-1);
                prevxx = d3tensor(0, model.numberOfSolutions-1, 0, kk-1, 0, jj-1);
                
                cumTime = 0.0;
                if (ddt == 0.0 || ddt > adaptiveMaxTimeStep) ddt = adaptiveMaxTimeStep;
                
                _s = _sTime[0] = dt;
                smallestCount = 0;
                while (cumTime < dt-1.0e-12) {
                    ddt = min(dt-cumTime, ddt);
                    i = 0;
                    for (FEMSolution *solution in model.solutions) {
                        varContainers = solution.variable.getContainers;
                        if (varContainers->Values != NULL) {
                            n = varContainers->sizeValues;
                            for (j=0; j<n; j++) {
                                xx[i][j] = varContainers->Values[j];
                            }
                            xxnrm[i] = solution.variable.norm;
                            if (varContainers->PrevValues != NULL) {
                                for (j=0; j<varContainers->size2PrevValues; j++) {
                                    for (k=0; k<n; k++) {
                                        prevxx[i][k][j] = varContainers->PrevValues[k][j];
                                    }
                                }
                            }
                        }
                        i++;
                    }
                    
                    _sTime[0] = _s + cumTime + ddt;
                    _sSize[0] = ddt;
                    [kernel solveEquationsModel:model timeStep:&ddt transientSimulation:transient coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter steadyStateReached:&steadyStateReached realTimeStep:&realTimeStep];
                    
                    maxErr = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"adaptive error measure" info:&found minValue:NULL maxValue:NULL];
                    
                    i = 0;
                    for (FEMSolution *solution in model.solutions) {
                        varContainers = solution.variable.getContainers;
                        if (varContainers->Values != NULL) {
                            n = varContainers->sizeValues;
                            yynrm[i] = solution.variable.norm;
                            for (j=0; j<n; j++) {
                                varContainers->Values[j] = xx[i][j];
                            }
                            if (varContainers->PrevValues != NULL) {
                                for (j=0; j<varContainers->size2PrevValues; j++) {
                                    for (k=0; k<n; k++) {
                                        varContainers->PrevValues[k][j] = prevxx[i][k][j];
                                    }
                                }
                            }
                        }
                        i++;
                    }
                    
                    _sStep[0] = ddt / 2;
                    _sTime[0] = _s + cumTime + ddt/2;
                    arg = ddt / 2;
                    [kernel solveEquationsModel:model timeStep:&arg transientSimulation:transient coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter steadyStateReached:&steadyStateReached realTimeStep:&realTimeStep];
                    _sTime[0] = _s + cumTime + ddt;
                    [kernel solveEquationsModel:model timeStep:&arg transientSimulation:transient coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter steadyStateReached:&steadyStateReached realTimeStep:&realTimeStep];
                    
                    maxErr = fabs(maxErr - [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"adaptive error measure" info:&found minValue:NULL maxValue:NULL]);
                    
                    if (found == NO) {
                        maxErr = 0.0;
                        i = 0;
                        for (FEMSolution *solution in model.solutions) {
                            varContainers = solution.variable.getContainers;
                            if (varContainers->Values != NULL) {
                                if (yynrm[i] != solution.variable.norm) {
                                    maxErr = max(maxErr, fabs(yynrm[i]-solution.variable.norm)/yynrm[i]);
                                }
                            }
                            i++;
                        }
                    }
                    
                    if (maxErr < adaptiveLimit || ddt <= adaptiveMinTimeStep) {
                        cumTime = cumTime + ddt;
                        realTimeStep++;
                        if (smallestCount >= adaptiveKeepSmallest || stepControl > 0) {
                            ddt = min(2.0*ddt, adaptiveMaxTimeStep);
                            stepControl = 1;
                            smallestCount = 0;
                        } else {
                            stepControl = 0;
                            smallestCount++;
                        }
                    } else {
                        i = 0;
                        for (FEMSolution *solution in model.solutions) {
                            varContainers = solution.variable.getContainers;
                            if (varContainers->Values != NULL) {
                                n = varContainers->sizeValues;
                                solution.variable.norm = xxnrm[i];
                                for (j=0; j<n; j++) {
                                    varContainers->Values[j] = xx[i][j];
                                }
                                if (varContainers->PrevValues != NULL) {
                                    for (j=0; j<varContainers->size2PrevValues; j++) {
                                        for (k=0; k<n; k++) {
                                            varContainers->PrevValues[k][j] = prevxx[i][k][j];
                                        }
                                    }
                                }
                            }
                            i++;
                        }
                        ddt = ddt / 2.0;
                        stepControl = -1;
                    }
                    NSLog(@"Adaptive(cum, ddt, err): %f, %f, %f\n", cumTime, ddt, maxErr);
                }
                
                _sSize[0] = dt;
                _sTime[0] = _s + dt;
                
                free_dmatrix(xx, 0, model.numberOfSolutions-1, 0, kk-1);
                free_dvector(yynrm, 0, model.numberOfSolutions-1);
                free_dvector(xxnrm, 0, model.numberOfSolutions-1);
                free_d3tensor(prevxx, 0, model.numberOfSolutions-1, 0, kk-1, 0, jj-1);
            } else {
                [kernel solveEquationsModel:model timeStep:&dt transientSimulation:transient coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter steadyStateReached:&steadyStateReached realTimeStep:&realTimeStep];
                realTimeStep++;
            }
            
            // Save results to disk if requested
            _lastSaved = NO;
            if (outputIntervals[interval-1] != 0) {
                
                [self FEMJob_saveToPostModel:model currentStep:0];
                k = (timeStep-1) % outputIntervals[interval-1];
                if (k == 0 || steadyStateReached == YES) {
                    for (FEMSolution *solution in model.solutions) {
                        if (solution.selector == NULL || solution.plugInPrincipalClassInstance == nil) continue;
                        execThis = (solution.solutionSolveWhen == SOLUTION_SOLVE_AHEAD_SAVE) ? YES : NO;
                        if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                             _when = (solution.solutionInfo)[@"invoke solution computer"];
                            execThis = ([_when isEqualToString:@"before saving"] == YES) ? YES : NO;
                            if (execThis) [kernel activateSolution:solution model:model timeStep:dt transientSimulation:transient];
                        }
                    }
                    
                    // TODO: Call save current here
                    _lastSaved = YES;
                    
                    for (FEMSolution *solution in model.solutions) {
                        if (solution.selector == NULL || solution.plugInPrincipalClassInstance == nil) continue;
                        execThis = (solution.solutionSolveWhen == SOLUTION_SOLVE_AFTER_SAVE) ? YES : NO;
                        if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                            _when = (solution.solutionInfo)[@"invoke solution computer"];
                            execThis = ([_when isEqualToString:@"after saving"] == YES) ? YES : NO;
                        }
                        if (execThis == YES) [kernel activateSolution:solution model:model timeStep:dt transientSimulation:transient];
                    }
                }
            }
            
            maxTime = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"real time max" info:&found minValue:NULL maxValue:NULL];
            if (found == YES && realtime() > maxTime) {
                NSLog(@"JOB: Reached allowed maximum real time, exiting...");
                goto jump;
            }
            
            exitCond = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"exit condition" info:&found minValue:NULL maxValue:NULL];
            if (found == YES && exitCond > 0.0) {
                NSLog(@"JOB: Found a positive exit condition, exiting...");
                goto jump;
            }
            
            if (steadyStateReached == YES && (transient == NO || scanning == NO)) {
                if (timeStep >= coupledMinIter) break;
            }
        } // Time step within an interval
    } // Time step intervals, i.e., the simulation
    
jump:
    for (FEMSolution *solution in model.solutions) {
        if (solution.selector == NULL || solution.plugInPrincipalClassInstance == nil) continue;
        if ( (solution.solutionInfo)[@"invoke solution computer"] != nil) {
            _when = (solution.solutionInfo)[@"invoke solution computer"];
            if ([_when isEqualToString:@"after simulation"] || [_when isEqualToString:@"after all"]) {
                [kernel activateSolution:solution model:model timeStep:dt transientSimulation:transient];
                _lastSaved = NO;
            }
        } else {
            if (solution.solutionSolveWhen == SOLUTION_SOLVE_AFTER_ALL) {
                [kernel activateSolution:solution model:model timeStep:dt transientSimulation:transient];
                _lastSaved = NO;
            }
        }
    }
    
    if (_lastSaved == NO) {
        for (FEMSolution *solution in model.solutions) {
            if (solution.selector == NULL || solution.plugInPrincipalClassInstance == nil) continue;
            execThis = (solution.solutionSolveWhen == SOLUTION_SOLVE_AHEAD_SAVE) ? YES : NO;
            if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                _when = (solution.solutionInfo)[@"invoke solution computer"];
                execThis = ([_when isEqualToString:@"before saving"] == YES) ? YES : NO;
            }
            if (execThis == YES) [kernel activateSolution:solution model:model timeStep:dt transientSimulation:transient];
        }
    }
}

/***********************************************************************************
    Saves results file to post proecessing file of ElmerPost format if requested
***********************************************************************************/
-(void)FEMJob_saveToPostModel:(FEMModel *)model currentStep:(int)currentStep {
    
    int j, k, timeSteps, savedEigenValues;
    BOOL found, lastExisting, eigenAnal=NO;
    NSString *outputFile, *postFile, *saveWhich;
    FEMListUtilities *listUtilities;
    FEMUtilities *utilities;
    FEMMeshUtils *meshUtilities;
    variableArraysContainer *varContainers = NULL;
    
    listUtilities = [[FEMListUtilities alloc] init];
    
    outputFile = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"output file" info:&found];
    if (found == YES) {
        // TODO: add support for parallel run
    }
    
    postFile = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"post file" info:&found];
    if (found == NO) return;
    
    // TODO: add support for parallel run
    
    utilities = [[FEMUtilities alloc] init];
    meshUtilities = [[FEMMeshUtils alloc] init];
    
    // Loop over all meshes
    for (FEMMesh *mesh in model.meshes) {
        if (currentStep == 0 && mesh.savesDone > 0) continue;
        
        // Check whether this mesh is active for saving
        if (mesh.outputActive == YES) {
            if ((int)[mesh.name length] == 0 || [utilities isFileNameQualified:outputFile]) {
                [_outputName setString:outputFile];
            } else {
                [_outputName setString:mesh.name];
                [_outputName appendString:@"/"];
                [_outputName appendString:outputFile];
            }
            
            if ((int)[mesh.name length] == 0 || [utilities isFileNameQualified:postFile]) {
                [_postName setString:postFile];
            } else {
                [_postName setString:mesh.name];
                [_postName appendString:@"/"];
                [_postName appendString:postFile];
            }
            
            if ([listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"filename numbering" info:&found] == YES) {
                if (currentStep == 0) {
                    [_postName setString:[utilities nextFreeFileName:_postName suffix:nil lastExisting:NULL]];
                } else {
                    lastExisting = YES;
                    [_postName setString:[utilities nextFreeFileName:_postName suffix:nil lastExisting:&lastExisting]];
                }
            }
            
            // Set the current mesh pointer in the model
            [meshUtilities setCurrentMesh:mesh inModel:model];
            
            // Use number of time steps or number of eigenmodes
            timeSteps = _totalTimeSteps;
            for (FEMSolution *solution in model.solutions) {
                if (solution.mesh == mesh) {
                    if ((solution.solutionInfo)[@"eigen analysis"] != nil) {
                        eigenAnal = [(solution.solutionInfo)[@"eigen analysis"] boolValue];
                    }
                    eigenAnal = (eigenAnal == YES || [(solution.solutionInfo)[@"harmonic analysis"] boolValue] == YES) ? YES : NO;
                    if (eigenAnal == YES) timeSteps = max(timeSteps, solution.nOfEigenValues);
                }
            }
            
            for (FEMSolution *solution in model.solutions) {
                if (solution.mesh == mesh) {
                    eigenAnal = [(solution.solutionInfo)[@"eigen analysis"] boolValue];
                    eigenAnal = (eigenAnal == YES || [(solution.solutionInfo)[@"harmonic analysis"] boolValue] == YES) ? YES : NO;
                    if (eigenAnal == YES) {
                        if ((solution.solutionInfo)[@"eigen and harmonic solution output"] != nil) {
                            saveWhich = (solution.solutionInfo)[@"eigen and harmonic solution output"];
                        }
                        savedEigenValues = solution.nOfEigenValues;
                        for (j=0; j<savedEigenValues; j++) {
                            for (FEMVariable *variable in mesh.variables) {
                                varContainers = variable.getContainers;
                                if (varContainers->EigenValues == NULL) continue;
                                if (solution.matrix.complexMatrix == YES) {
                                    for (k=0; k<varContainers->sizeValues/2; k++) {
                                        varContainers->Values[2*k] = creal(varContainers->EigenVectors[j][k]);
                                        varContainers->Values[2*k+1] = cimag(varContainers->EigenVectors[j][k]);
                                    }
                                } else {
                                    if ([saveWhich isEqualToString:@"real part"] == YES) {
                                        for (k=0; k<varContainers->sizeValues; k++) {
                                            // TODO: is that really correct, we dont take the real part?
                                            varContainers->Values[k] = varContainers->EigenVectors[j][k];
                                        }
                                    } else if ([saveWhich isEqualToString:@"imag part"] == YES) {
                                        for (k=0; k<varContainers->sizeValues; k++) {
                                            varContainers->Values[k] = cimag(varContainers->EigenVectors[j][k]);
                                        }
                                    } else if ([saveWhich isEqualToString:@"abs value"] == YES) {
                                        for (k=0; k<varContainers->sizeValues; k++) {
                                            varContainers->Values[k] = cabs(varContainers->EigenVectors[j][k]);
                                        }
                                    } else if ([saveWhich isEqualToString:@"phase angle"] == YES) {
                                        for (k=0; k<varContainers->sizeValues; k++) {
                                            varContainers->Values[k] = atan2(cimag(varContainers->EigenVectors[j][k]), creal(varContainers->EigenVectors[j][k]));
                                        }
                                    } else {
                                        varContainers->CValues = malloc ( varContainers->size2EigenVectors * sizeof ( double * ));
                                        for (k=0; k<varContainers->size2EigenVectors; k++) {
                                            varContainers->CValues[k] = &varContainers->EigenVectors[j][k];
                                        }
                                        varContainers->sizeCValues = varContainers->size2EigenVectors;
                                    }
                                }
                            }
                            
                            if (currentStep > 0) {
                                if (mesh.savesDone != 0) {
                                    if (_totalTimeSteps == 1) {
                                        mesh.savesDone = j+1;
                                    } else {
                                        mesh.savesDone = currentStep;
                                    }
                                }
                                // TODO: call write to post file here
                            }
                        }
                        break;
                    }
                }
            }
            
            // If this mesh has not geen saved, then do so
            if (eigenAnal == NO || currentStep == 0) {
                // TODO: call write to post file here
            }
            
            for (FEMVariable *variable in mesh.variables) {
                varContainers = variable.getContainers;
                if (varContainers->EigenValues != NULL) {
                    memset( varContainers->Values, 0.0, (varContainers->sizeValues*sizeof(varContainers->Values)) );
                    varContainers->CValues = NULL;
                }
            }
        }
    }
}

// Saves current time step to external files
-(void)FEMJob_saveCurrent:(FEMModel *)model currentStep:(int)currentStep {
    
    int j, k;
    BOOL found, binaryOutput, saveAll, eigenAnal=NO;
    NSString *outputFile;
    FEMListUtilities *listUtilities;
    variableArraysContainer *varContainers = NULL;
    
    listUtilities = [[FEMListUtilities alloc] init];
    
    outputFile = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"output file" info:&found];
    if (found == YES) {
        // TODO: add support for parallel run
        
        binaryOutput = [listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"binary output" info:&found];
        if (found == NO) binaryOutput = NO;
        
        saveAll = ([listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"omit unchanged variables in output" info:&found] == YES) ? NO : YES;
        if (found == NO) saveAll = YES;
        
        for (FEMMesh *mesh in model.meshes) {
            if (mesh.outputActive == YES) {
                if ((int)[mesh.name length] > 0) {
                    [_outputName setString:mesh.name];
                    [_outputName appendString:@"/"];
                    [_outputName appendString:outputFile];
                } else {
                    [_outputName setString:outputFile];
                }
                
                for (FEMSolution *solution in model.solution) {
                    if (solution.mesh == mesh) {
                        eigenAnal = [(solution.solutionInfo)[@"eigen analysis"] boolValue];
                        eigenAnal = (eigenAnal == YES|| [(solution.solutionInfo)[@"harmonic analysis"] boolValue] == YES) ? YES : NO;
                        
                        if (eigenAnal == YES) {
                            varContainers = solution.variable.getContainers;
                            if (varContainers->EigenValues != NULL) {
                                if (_totalTimeSteps == 1) {
                                    for (j=0; j<solution.nOfEigenValues; j++) {
                                        if (solution.matrix.complexMatrix == YES) {
                                            for (k=0; k<varContainers->sizeValues/2; k++) {
                                                varContainers->Values[2*k] = creal(varContainers->EigenVectors[j][k]);
                                                varContainers->Values[2*k+1] = cimag(varContainers->EigenVectors[j][k]);
                                            }
                                        } else {
                                            for (k=0; k<varContainers->sizeValues; k++) {
                                                varContainers->Values[k] = creal(varContainers->EigenVectors[j][k]);
                                            }
                                        }
                                        // TODO: call save results here
                                    }
                                } else {
                                    j = min(currentStep, varContainers->size1EigenVectors);
                                    if (solution.matrix.complexMatrix == YES) {
                                        for (k=0; k<varContainers->sizeValues/2; k++) {
                                            varContainers->Values[2*k] = creal(varContainers->EigenVectors[j][k]);
                                            varContainers->Values[2*k+1] = cimag(varContainers->EigenVectors[j][k]);
                                        }
                                    } else {
                                        for (k=0; k<varContainers->sizeValues; k++) {
                                            varContainers->Values[k] = creal(varContainers->EigenVectors[j][k]);
                                        }
                                    }
                                    // TODO: call save results here
                                }
                                memset( varContainers->Values, 0.0, (varContainers->sizeValues*sizeof(varContainers->Values)) );
                            }
                        }
                    }
                }
                
                if (eigenAnal == NO) {
                    // TODO: call save results here
                }
            }
        }
    } else {
        for (FEMMesh *mesh in model.meshes) {
            if (mesh.outputActive == YES) mesh.savesDone++;
        }
    }
    [self FEMJob_saveToPostModel:model currentStep:currentStep];
}
  
#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _transient = NO;
        _outputName = [NSMutableString stringWithString:@""];
        _postName = [NSMutableString stringWithString:@""];
    }
    
    return self;
}

-(void)deallocation {
    
}

@end
