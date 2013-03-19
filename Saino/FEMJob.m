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

@interface FEMJob ()
-(void)FEMJob_addSolutionsModel:(FEMModel *)model;
-(void)FEMJob_addMeshCoordinatesAndTimeModel:(FEMModel *)model;
-(void)FEMJob_setInitialConditionsModel:(FEMModel *)model;
-(void)FEMJob_initCondModel:(FEMModel *)model;
-(void)FEMJob_restartModel:(FEMModel *)model;
-(void)FEMJob_runModel:(FEMModel *)model timeIntervals:(int)timeIntervals coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter outputIntervals:(int* )outputIntervals transient:(BOOL)transient scanning:(BOOL)scanning;
@end

@implementation FEMJob {
    
    BOOL _transient;
    BOOL _initDirichlet;
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
    
    NSMutableString *_outputName;
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

-(void)FEMJob_runModel:(FEMModel *)model timeIntervals:(int)timeIntervals coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter outputIntervals:(int* )outputIntervals transient:(BOOL)transient scanning:(BOOL)scanning{
    
    for (FEMSolution *solution in model.solutions) {
        
    }
}

#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _transient = NO;
        _outputName = [NSMutableString stringWithString:@""];
    }
    
    return self;
}

-(void)deallocation {
    
}

@end
