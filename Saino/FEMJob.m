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
#import "FEMMeshUtils.h"
#import "FEMElementDescription.h"
#import "FEMElementUtils.h"
#import "FEMPost.h"
#import "Utils.h"
#import "TimeProfile.h"
#import "GaussIntegration.h"

#ifdef TEST
    #import "FEMTest.h"
#endif

@interface FEMJob ()
-(void)FEMJob_addSolutionsModel:(FEMModel * __nonnull)model;
-(void)FEMJob_addMeshCoordinatesAndTimeModel:(FEMModel * __nonnull)model;
-(void)FEMJob_setInitialConditionsModel:(FEMModel * __nonnull)model;
-(void)FEMJob_initCondModel:(FEMModel * __nonnull)model;
-(void)FEMJob_restartModel:(FEMModel * __nonnull)model;
-(void)FEMJob_runSimulation:(FEMModel * __nonnull)model timeIntervals:(int)timeIntervals coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter outputIntervals:(int * __nonnull)outputIntervals transient:(BOOL)transient scanning:(BOOL)scanning post:(FEMPost * __nonnull)post;
-(void)FEMJob_saveToPostModel:(FEMModel * __nonnull)model currentStep:(int)currentStep post:(FEMPost * __nonnull)post;
-(int)FEMJob_saveResult:(NSString * __nonnull)fileName model:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh time:(int)time simulationTime:(double)simulationTime binary:(BOOL)binary saveAll:(BOOL)saveAll freeSurface:(BOOL * __nullable)freeSurface post:(FEMPost * __nonnull)post;
-(void)FEMJob_saveCurrent:(FEMModel * __nonnull)model currentStep:(int)currentStep post:(FEMPost * __nonnull)post;
@end

@implementation FEMJob {
    
    BOOL _scanning;
    BOOL _transient;
    BOOL _initDirichlet;
    BOOL _lastSaved;
    BOOL _firstLoad;
    BOOL _firstTime;
    BOOL _silent;
    BOOL _version;
    BOOL _gotModelName;
    int _timeIntervals;
    int _totalTimeSteps;
    int _savedSteps;
    int _coupledMaxIter;
    int _coupledMinIter;
    int * __nullable _timeSteps;
    int * __nullable _outputIntervals;
    double _s;
    double _dt;
    double * __nullable _sTime;
    double * __nullable _sStep;
    double * __nullable _sInterval;
    double * __nullable _sSize;
    double * __nullable * __nullable _sPrevSizes;
    double * __nullable _steadyIt;
    double * __nullable _nonLinIt;
    double * __nullable * __nullable _timeStepSizes;
    int _sizeTimeSteps;
    int _sizeOutputIntervals;
    int _sizeSTime;
    int _sizeSStep;
    int _sizeSInterval;
    int _sizeSSize;
    int _sizeSteadyIt;
    int _sizeNonLinIt;
    int _size1SPrevSizes;
    int _size2SPrevSizes;
    int _size1TimeStepSizes;
    int _size2TimeStepSizes;
    
    NSMutableString * __nonnull _outputName;
    NSMutableString * __nonnull _postName;
    NSString * __nullable _when;
}

@synthesize core = _core;
@synthesize elementDescription = _elementDescription;
@synthesize listUtilities = _listUtilities;
@synthesize model = _model;
@synthesize modelName = _modelName;

#pragma mark Private methods

/**********************************************************************
 
    Add flags for active solutions
    
    Method corresponds partially to Elmer from git on October 27 2015

**********************************************************************/
-(void)FEMJob_addSolutionsModel:(FEMModel * __nonnull)model {
    
    int i, j;
    BOOL initSolution, found;
    NSString *eq;
    FEMListUtilities *listUtilities;
    FEMUtilities *utilities;
    listBuffer activeSolvers = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    
    fprintf(stdout, "FEMJob:FEMJob_addSolutionsModel: setting up %d solvers.\n", model.numberOfSolutions);
    
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
                            [listUtilities addLogicalInClassList:equation theVariable:eq withValue:YES];
                            break;
                        }
                    }
                }
            }
        }
        i++;
    }
    
    if (activeSolvers.ivector != NULL) {
        free_ivector(activeSolvers.ivector, 0, activeSolvers.m-1);
    }
    
    utilities = [[FEMUtilities alloc] init];
    initSolution = NO;
    i = 1;
    for (FEMSolution *solution in model.solutions) {
        if ((solution.solutionInfo)[@"equation"] != nil) {
            eq = (solution.solutionInfo)[@"equation"];
            fprintf(stdout, "FEMJob:FEMJob_addSolutionsModel: setting up solver %d : %s.\n", i, [eq UTF8String]);
        }
        if ((solution.solutionInfo)[@"initialize"] != nil) {
            if ([(solution.solutionInfo)[@"initialize"] boolValue] == YES) {
                [solution.matrix deallocation];
                solution.matrix = nil;
                [solution.solutionInfo setObject:@YES forKey:@"initialize"];
            }
        }
        
        if ((solution.plugInPrincipalClassInstance == nil && solution.hasBuiltInSolution == NO) || initSolution == YES) {
            //TODO: Make sure that this is always correct. Here if the solution has no mesh
            // we assigned it to the first mesh listed in model.meshes
            if (solution.mesh == nil) solution.mesh = model.meshes[0];
            model.solution = solution;
            [utilities addEquationBasicsToSolution:solution name:eq model:model transient:_transient];
            [utilities addEquationToSolution:solution model:model transient:_transient];
        }
        i++;
    }
}

/*******************************************************************
    Add coordinate and time variables to the meshe(s) in the model
*******************************************************************/
-(void)FEMJob_addMeshCoordinatesAndTimeModel:(FEMModel * __nonnull)model {
    
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

-(void)FEMJob_initCondModel:(FEMModel * __nonnull)model {
    
    int i, j, k, k1, l, n, t, dofs;
    int *indexes;
    double integral;
    BOOL found;
    NSString *str = nil;
    NSMutableString *varName;
    listBuffer vector = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer tensor = { NULL, NULL, NULL, NULL, 0, 0, 0};
    FEMSolution *solution;
    FEMInitialConditions *initialCondition;
    FEMListUtilities *listUtilities;
    FEMMeshUtils *meshUtilities;
    Element_t *elements = NULL, *edges = NULL;
    variableArraysContainer *variableContainers = NULL;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
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
                        
                        n = [self.core getNumberOfNodesForElement:&elements[t]];
                        for (FEMVariable *variable in mesh.variables) {
                            variableContainers = variable.getContainers;
                            if (variable.solution != nil) {
                                dofs = [self.core getElementDofsSolution:solution model:model forElement:&elements[t] atIndexes:indexes disableDiscontinuousGalerkin:NULL];
                            } else {
                                dofs = [self.core getElementDofsSolution:nil model:model forElement:&elements[t] atIndexes:indexes disableDiscontinuousGalerkin:NULL];
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
                                found = [self.core getReal:model forElement:&elements[t] inArray:initialCondition.valuesList variableName:variable.name buffer:&vector listUtilities:listUtilities];
                                if (found == YES) {
                                    for (k=0; k<n; k++) {
                                        k1 = indexes[k];
                                        if (variableContainers->Perm != NULL) k1 = variableContainers->Perm[k1];
                                        if (k1 >= 0) {
                                             if (variable.isComponentVariable == YES) {
                                                 *(variableContainers->ComponentValues[k1]) = vector.vector[k];
                                             } else {
                                                 variableContainers->Values[k1] = vector.vector[k];
                                             }
                                        }
                                    }
                                }
                                
                                if (_transient == YES && solution.timeOrder == 2) {
                                    varName = [NSMutableString stringWithString:variable.name];
                                    [varName appendString:@" velocity"];
                                    found = [self.core getReal:model forElement:&elements[t] inArray:initialCondition.valuesList variableName:varName buffer:&vector listUtilities:listUtilities];
                                    if (found == YES) {
                                        for (k=0; k<n; k++) {
                                            k1 = indexes[k];
                                            if (variableContainers->Perm != NULL) k1 = variableContainers->Perm[k1];
                                            if (k1 >= 0) variableContainers->PrevValues[k1][0] = vector.vector[k];
                                        }
                                    }
                                    varName = [NSMutableString stringWithString:variable.name];
                                    [varName appendString:@" acceleration"];
                                    found = [self.core getReal:model forElement:&elements[t] inArray:initialCondition.valuesList variableName:varName buffer:&vector listUtilities:listUtilities];
                                    if (found == YES) {
                                        for (k=0; k<n; k++) {
                                            k1 = indexes[k];
                                            if (variableContainers->Perm != NULL) k1 = variableContainers->Perm[k1];
                                            if (k1 >= 0) variableContainers->PrevValues[k1][1] = vector.vector[k];
                                        }
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
                                                    [self.core localBoundaryIntegral:model inSolution:solution atBoundary:initialCondition.valuesList forElement:&edges[elements[t].EdgeIndexes[k]] withNumberOfNodes:edges[elements[t].EdgeIndexes[k]].Type.NumberOfNodes andParent:&elements[t] withNumberOfNodes:n boundaryName:varName functionIntegral:&integral];
                                                    variableContainers->Values[l] = integral;
                                                }
                                            }
                                        }
                                    }
                                }
                            } else {
                                found = [listUtilities listGetRealArray:model inArray:initialCondition.valuesList forVariable:variable.name numberOfNodes:n indexes:elements[t].NodeIndexes buffer:&tensor];
                                if (found == YES) {
                                    for (k=0; k<n; k++) {
                                        k1 = indexes[k];
                                        for (l=0; l<min(tensor.m, variable.dofs); l++) {
                                            if (variableContainers->Perm != NULL) k1 = variableContainers->Perm[k1];
                                            if (k1 >= 0) variableContainers->Values[variable.dofs*k1+l] = tensor.tensor[l][0][k];
                                        }
                                    }
                                }
                            }
                            [listUtilities listSetNameSpace:nil];
                        }
                    }
                }
            }
            free_ivector(indexes, 0, mesh.maxElementDofs-1);
        }
    }
    if (vector.vector != NULL) free_dvector(vector.vector, 0, vector.m-1);
    if (tensor.tensor != NULL) free_d3tensor(tensor.tensor, 0, tensor.m-1, 0, tensor.n-1, 0, tensor.p-1);
}

/*************************************************************
    Check if we are restarting, and if yes read field values
*************************************************************/
-(void)FEMJob_restartModel:(FEMModel * __nonnull)model {
    
    int k, minVal;
    double startTime;
    BOOL found;
    NSString *restartFile;
    FEMVariable *var = nil;
    FEMListUtilities *listUtilities;
    FEMMeshUtils *meshUtilities;
    FEMUtilities *utilities;
    variableArraysContainer *varContainers = NULL;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    
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
-(void)FEMJob_setInitialConditionsModel:(FEMModel * __nonnull)model {
    
    int i, j, k, l, m, n, t, dim, vectDof=0, realDof=-1;
    double udot, parU, parV, *nrm, *t1, *t2, *vec, *tmp;
    BOOL found, ntBoundary, pointed=NO, check;
    listBuffer vector = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer tensor = { NULL, NULL, NULL, NULL, 0, 0, 0};
    NSArray *bc;
    NSString *str = nil;
    NSMutableString *varName;
    FEMVariable *vectVariable = nil;
    FEMSolution *solution = nil;
    FEMListUtilities *listUtilities;
    FEMMeshUtils *meshUtilities;
    FEMElementUtils *elementUtils;
    Element_t *elements = NULL;
    Nodes_t *nodes = NULL;
    variableArraysContainer *variableContainers = NULL, *vectVarContainers = NULL;
    
    dim = model.dimension;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
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
    
    if (_initDirichlet == YES) {
        FEMElementDescription * elementDescription = [FEMElementDescription sharedElementDescription];
        nrm = doublevec(0, 2);
        t1 = doublevec(0, 2);
        t2 = doublevec(0, 2);
        vec = doublevec(0, 2);
        tmp = doublevec(0, 2);

        meshUtilities = [[FEMMeshUtils alloc] init];
        for (FEMMesh *mesh in model.meshes) {
            [meshUtilities setCurrentMesh:mesh inModel:model];
            
            elements = mesh.getElements;
            for (t=mesh.numberOfBulkElements; t<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; t++) {
                n = elements[t].Type.NumberOfNodes;
                bc = [self.core getBoundaryCondition:model forElement:&elements[t]];
                
                for (FEMVariable *variable in mesh.variables) {
                    variableContainers = variable.getContainers;
                    solution = (FEMSolution *)variable.solution;
                    if (solution == nil) solution = (FEMSolution *)model.solution;
                    
                    if ((solution.solutionInfo)[@"namespace"] != nil) {
                        str = (solution.solutionInfo)[@"namespace"];
                        [listUtilities listSetNameSpace:str];
                    }
                    
                    if (variable.dofs <= 1) {
                        found = [self.core getReal:model forElement:&elements[t] inArray:bc variableName:variable.name buffer:&vector listUtilities:listUtilities];
                        if (found == YES) {
                            ntBoundary = NO;
                            if ([self.core getElementFamily:&elements[t]] != 1) {
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
                                        [self.core getNodes:solution model:model inElement:&elements[t] resultNodes:nodes numberOfNodes:NULL mesh:nil];
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
                                                memcpy(vec, nrm, 3*sizeof(double));
                                                break;
                                            case 2:
                                                memcpy(vec, t1, 3*sizeof(double));
                                                break;
                                            case 3:
                                                memcpy(vec, t2, 3*sizeof(double));
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
                                        udot = cblas_ddot(dim, vec, 1, tmp, 1);
                                        for (i=0; i<dim; i++) {
                                            tmp[i] = tmp[i]+(vector.vector[j]-udot)*vec[i];
                                        }
                                        for (l=0; l<dim; l++) {
                                            m = l + realDof - (--vectDof);
                                            vectVarContainers->Values[vectVariable.dofs*k+m] = tmp[l];
                                        }
                                    } else {
                                        if (variable.isComponentVariable == YES) {
                                            *(variableContainers->ComponentValues[k]) = vector.vector[j];
                                        } else {
                                            variableContainers->Values[k] = vector.vector[j];
                                        }
                                    }
                                }
                            }
                        }
                        
                        if (_transient == YES && solution.timeOrder == 2) {
                            varName = [NSMutableString stringWithString:variable.name];
                            [varName appendString:@" velocity"];
                            found = [self.core getReal:model forElement:&elements[t] inArray:bc variableName:varName buffer:&vector listUtilities:listUtilities];
                            if (found == YES) {
                                for (j=0; j<n; j++) {
                                    k = elements[t].NodeIndexes[j];
                                    if (variableContainers->Perm != NULL) k = variableContainers->Perm[k];
                                    if (k >= 0) variableContainers->PrevValues[k][0] = vector.vector[j];
                                }
                            }
                            
                            varName = [NSMutableString stringWithString:variable.name];
                            [varName appendString:@" acceleration"];
                            found = [self.core getReal:model forElement:&elements[t] inArray:bc variableName:varName buffer:&vector listUtilities:listUtilities];
                            if (found == YES) {
                                for (j=0; j<n; j++) {
                                    k = elements[t].NodeIndexes[j];
                                    if (variableContainers->Perm != NULL) k = variableContainers->Perm[k];
                                    if (k >= 0) variableContainers->PrevValues[k][1] = vector.vector[j];
                                }
                            }
                        }
                    } else {
                        found = [listUtilities listGetRealArray:model inArray:bc forVariable:variable.name numberOfNodes:n indexes:elements[t].NodeIndexes buffer:&tensor];
                        if (found == YES) {
                            for (j=0; j<n; j++) {
                                k = elements[t].NodeIndexes[j];
                                for (l=0; l<min(tensor.m, variable.dofs); l++) {
                                    if (variableContainers->Perm != NULL) k = variableContainers->Perm[k];
                                    if (k >= 0) variableContainers->Values[variable.dofs*k+l] = tensor.tensor[l][0][j];
                                }
                            }
                        }
                    }
                    
                    [listUtilities listSetNameSpace:nil];
                }
            }
        }
        free_dvector(nrm, 0, 2);
        free_dvector(t1, 0, 2);
        free_dvector(t2, 0, 2);
        free_dvector(vec, 0, 2);
        free_dvector(tmp, 0, 2);
    }
    
    if (vector.vector != NULL) free_dvector(vector.vector, 0, vector.m-1);
    if (tensor.tensor != NULL) free_d3tensor(tensor.tensor, 0, tensor.m-1, 0, tensor.n-1, 0, tensor.p-1);
}

-(void)FEMJob_runSimulation:(FEMModel * __nonnull)model timeIntervals:(int)timeIntervals coupledMinIteration:(int)coupledMinIter coupleMaxIteration:(int)coupleMaxIter outputIntervals:(int * __nonnull)outputIntervals transient:(BOOL)transient scanning:(BOOL)scanning post:(FEMPost * __nonnull)post {
    
    int i, j, k, jj, kk, n, interval, timeStep, stepCount = 0, cumTimeStep, realTimeStep, timeLeft, adaptiveKeepSmallest, minVal, smallestCount,
        stepControl=-1;
    double dt = 0, ddt, dtFunc, newTime, maxTime, prevTime=0, adaptiveLimit, adaptiveMaxTimeStep, adaptiveMinTimeStep, cumTime, maxErr,
           exitCond, arg;
    double **xx, *xxnrm, *yynrm, ***prevxx;
    BOOL found, adaptiveTime, steadyStateReached=NO, execThis;
    FEMListUtilities *listUtilities;
    variableArraysContainer *varContainers = NULL;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    
    for (FEMSolution *solution in model.solutions) {
        if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
        if (solution.solutionSolveWhen == SOLUTION_SOLVE_AHEAD_ALL) {
            [self.core activateSolution:solution model:model timeStep:dt transientSimulation:transient];
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
            fprintf(stdout, "JOB: \n");
            fprintf(stdout, "JOB: ----------------------------------------------------------------------\n");
            
            if (transient == YES || scanning == YES) {
                fprintf(stdout, "JOB: Time: %d / %d %f.\n", cumTimeStep, stepCount, _sTime[0]);
                
                newTime = realtime();
                
                if (cumTimeStep > 1) {
                    maxTime = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"real time max" info:&found minValue:NULL maxValue:NULL];
                    if (found == YES) {
                        fprintf(stdout, "JOB: Fraction of real time left: %f.\n", 1.0-realtime()/maxTime);
                    } else {
                        timeLeft = round((stepCount-(cumTimeStep-1))*(newTime-prevTime)/60.0);
                        if (timeLeft > 120) {
                            fprintf(stdout, "JOB: Estimated time left: %d hours.\n", timeLeft/60);
                        } else if (timeLeft > 60) {
                            fprintf(stdout, "JOB: Estimated time left: 1 hour %d minutes.\n", timeLeft % 60);
                        } else if (timeLeft >= 1) {
                            fprintf(stdout, "JOB: Estimated time left: %d minutes.\n", timeLeft);
                        } else {
                            fprintf(stdout, "JOB: Estimated time left: less than a minute.\n");
                        }
                    }
                }
                prevTime = newTime;
            } else {
                fprintf(stdout, "JOB: Steady state iteration: %d.\n", cumTimeStep);
            }
            
            fprintf(stdout, "JOB: ----------------------------------------------------------------------\n");
            fprintf(stdout, "JOB: \n"); // End of code that should only be executed by one processor if parallel run
            
            // Solve any and all governing equations in the system
            adaptiveTime = [listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"adaptive time stepping" info:&found];
            
            if (transient == YES && adaptiveTime == YES) { // Adaptive time stepping
                adaptiveLimit = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"adaptive time error" info:&found minValue:NULL maxValue:NULL];
                
                if (found == NO) {
                    fatal("FEMJob:FEMJob_runSimulation", "Adaptive time limit must be be given for adaptive stepping scheme.");
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
                    [self.core solveEquationsModel:model timeStep:&ddt transientSimulation:transient coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter steadyStateReached:&steadyStateReached realTimeStep:&realTimeStep];
                    
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
                    [self.core solveEquationsModel:model timeStep:&arg transientSimulation:transient coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter steadyStateReached:&steadyStateReached realTimeStep:&realTimeStep];
                    _sTime[0] = _s + cumTime + ddt;
                    [self.core solveEquationsModel:model timeStep:&arg transientSimulation:transient coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter steadyStateReached:&steadyStateReached realTimeStep:&realTimeStep];
                    
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
                    fprintf(stdout, "FEMJob:FEMJob_runSimulation: adaptive(cum, ddt, err): %f, %f, %f.\n", cumTime, ddt, maxErr);
                }
                
                _sSize[0] = dt;
                _sTime[0] = _s + dt;
                
                free_dmatrix(xx, 0, model.numberOfSolutions-1, 0, kk-1);
                free_dvector(yynrm, 0, model.numberOfSolutions-1);
                free_dvector(xxnrm, 0, model.numberOfSolutions-1);
                free_d3tensor(prevxx, 0, model.numberOfSolutions-1, 0, kk-1, 0, jj-1);
            } else {
                [self.core solveEquationsModel:model timeStep:&dt transientSimulation:transient coupledMinIteration:coupledMinIter coupleMaxIteration:coupleMaxIter steadyStateReached:&steadyStateReached realTimeStep:&realTimeStep];
                realTimeStep++;
            }
            
            // Save results to disk if requested
            _lastSaved = NO;
            if (outputIntervals[interval-1] != 0) {
                
                [self FEMJob_saveToPostModel:model currentStep:0 post:post];
                k = (timeStep-1) % outputIntervals[interval-1];
                if (k == 0 || steadyStateReached == YES) {
                    for (FEMSolution *solution in model.solutions) {
                        if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
                        execThis = (solution.solutionSolveWhen == SOLUTION_SOLVE_AHEAD_SAVE) ? YES : NO;
                        if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                             _when = (solution.solutionInfo)[@"invoke solution computer"];
                            execThis = ([_when isEqualToString:@"before saving"] == YES) ? YES : NO;
                        }
                        if (execThis) [self.core activateSolution:solution model:model timeStep:dt transientSimulation:transient];
                    }
                    
                    // TODO: Call save current here
                    [self FEMJob_saveCurrent:model currentStep:timeStep post:post];
                    _lastSaved = YES;
                    
                    for (FEMSolution *solution in model.solutions) {
                        if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
                        execThis = (solution.solutionSolveWhen == SOLUTION_SOLVE_AFTER_SAVE) ? YES : NO;
                        if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                            _when = (solution.solutionInfo)[@"invoke solution computer"];
                            execThis = ([_when isEqualToString:@"after saving"] == YES) ? YES : NO;
                        }
                        if (execThis == YES) [self.core activateSolution:solution model:model timeStep:dt transientSimulation:transient];
                    }
                }
            }
            
            maxTime = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"real time max" info:&found minValue:NULL maxValue:NULL];
            if (found == YES && realtime() > maxTime) {
                fprintf(stdout, "JOB: Reached allowed maximum real time, exiting...");
                goto jump;
            }
            
            exitCond = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"exit condition" info:&found minValue:NULL maxValue:NULL];
            if (found == YES && exitCond > 0.0) {
                fprintf(stdout, "JOB: Found a positive exit condition, exiting...");
                goto jump;
            }
            
            if (steadyStateReached == YES && !(transient == YES || scanning == YES)) {
                if (timeStep >= coupledMinIter) break;
            }
        } // Time step within an interval
    } // Time step intervals, i.e., the simulation
    
jump:
    for (FEMSolution *solution in model.solutions) {
        if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
        if ( (solution.solutionInfo)[@"invoke solution computer"] != nil) {
            _when = (solution.solutionInfo)[@"invoke solution computer"];
            if ([_when isEqualToString:@"after simulation"] == YES || [_when isEqualToString:@"after all"] == YES) {
                [self.core activateSolution:solution model:model timeStep:dt transientSimulation:transient];
                _lastSaved = NO;
            }
        } else {
            if (solution.solutionSolveWhen == SOLUTION_SOLVE_AFTER_ALL) {
                [self.core activateSolution:solution model:model timeStep:dt transientSimulation:transient];
                _lastSaved = NO;
            }
        }
    }
    
    if (_lastSaved == NO) {
        for (FEMSolution *solution in model.solutions) {
            if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
            execThis = (solution.solutionSolveWhen == SOLUTION_SOLVE_AHEAD_SAVE) ? YES : NO;
            if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                _when = (solution.solutionInfo)[@"invoke solution computer"];
                execThis = ([_when isEqualToString:@"before saving"] == YES) ? YES : NO;
            }
            if (execThis == YES) [self.core activateSolution:solution model:model timeStep:dt transientSimulation:transient];
        }
    }
}

/***********************************************************************************
    Saves results file to post proecessing file of ElmerPost format if requested
***********************************************************************************/
-(void)FEMJob_saveToPostModel:(FEMModel * __nonnull)model currentStep:(int)currentStep post:(FEMPost * __nonnull)post {
    
    int j, k, timeSteps, savedEigenValues;
    BOOL found, lastExisting, eigenAnal=NO, append;
    NSString *outputFile, *postFile, *saveWhich;
    variableArraysContainer *varContainers = NULL;
    
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    
    outputFile = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"output file" info:&found];
    if (found == YES) {
        // TODO: add support for parallel run
    }
    
    postFile = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"post file" info:&found];
    if (found == NO) return;
    
    // TODO: add support for parallel run
    
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    FEMMeshUtils *meshUtilities = [[FEMMeshUtils alloc] init];
    
    // Loop over all meshes
    for (FEMMesh *mesh in model.meshes) {
        if (currentStep == 0 && mesh.savesDone > 0) continue;
        
        // Check whether this mesh is active for saving
        if (mesh.outputActive == YES) {
            if ([mesh.name length] == 0 || [utilities isFileNameQualified:outputFile]) {
                [_outputName setString:outputFile];
            } else {
                [_outputName setString:mesh.name];
                [_outputName appendString:@"/"];
                [_outputName appendString:outputFile];
            }
            
            if ([mesh.name length] == 0 || [utilities isFileNameQualified:postFile]) {
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
            
            // Use number of time steps or number of eigen modes
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
                                        if (varContainers->CValues == NULL) {
                                            varContainers->CValues = malloc ( varContainers->size2EigenVectors * sizeof ( double complex * ));
                                            varContainers->sizeCValues = varContainers->size2EigenVectors;
                                        }
                                        for (k=0; k<varContainers->size2EigenVectors; k++) {
                                            varContainers->CValues[k] = &varContainers->EigenVectors[j][k];
                                        }
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
                                append = YES;
                                [post writeElmerPostFile:_postName resultFile:_outputName model:model timeCount:solution.nOfEigenValues append:&append];
                            }
                        }
                        break;
                    }
                }
            }
            
            // If this mesh has not been saved, then do so
            if (eigenAnal == NO || currentStep == 0) {
                append = YES;
                [post writeElmerPostFile:_postName resultFile:_outputName model:model timeCount:timeSteps append:&append];
            }
            
            for (FEMVariable *variable in mesh.variables) {
                varContainers = variable.getContainers;
                if (varContainers->EigenValues != NULL) {
                    memset( varContainers->Values, 0.0, varContainers->sizeValues*sizeof(double) );
                    free(varContainers->CValues);
                    varContainers->CValues = NULL;
                }
            }
        }
    }
}

/***********************************************************************************
    Save fields in a file that may be used for restarting a simulation.
    The data is saved in ascii format (no binary format support yet).
***********************************************************************************/
-(int)FEMJob_saveResult:(NSString * __nonnull)fileName model:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh time:(int)time simulationTime:(double)simulationTime binary:(BOOL)binary saveAll:(BOOL)saveAll freeSurface:(BOOL * __nullable)freeSurface post:(FEMPost * __nonnull)post {
    
    int i, k, dofs, n, saveCount;
    BOOL found, freeSurfaceFlag, moveBoundary, sameAsPrev, all;
    NSMutableString *fName;
    FEMListUtilities *listUtilities;
    FEMUtilities *utilities;
    NSFileManager *fileManager;
    NSFileHandle *outputFileHandle;
    variableArraysContainer *varContainers = NULL, *prev = NULL;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    utilities = [[FEMUtilities alloc] init];
    
    // If first time here, count the number of variables
    if ([utilities isFileNameQualified:fileName] == NO) {
        if ([model.outputPath length] > 0) {
            fName = [NSMutableString stringWithString:model.outputPath];
            [fName appendString:@"/"];
            [fName appendString:fileName];
        }
    }
    
    if (freeSurface != NULL) {
        freeSurfaceFlag = *freeSurface;
    } else {
        freeSurfaceFlag = NO;
        moveBoundary = NO;
        for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
            freeSurfaceFlag = (freeSurfaceFlag == YES || [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"free surface" info:&found] == YES) ? YES : NO;
            if (freeSurfaceFlag == YES) {
                moveBoundary = [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"internal move boundary" info:&found];
                if (found == NO) moveBoundary = YES;
                freeSurfaceFlag = (freeSurfaceFlag == YES && moveBoundary == YES) ? YES : NO;
            }
            if (freeSurfaceFlag == YES) break;
        }
    }
    
    fileManager = [NSFileManager defaultManager];
    
    char space = ' ';
    char newLine = '\n';
    NSData *spaceBuff = [NSMutableData dataWithBytes:&space length:sizeof(space)];
    NSData *newLineBuff = [NSMutableData dataWithBytes:&newLine length:sizeof(newLine)];

    if (mesh.savesDone == 0) {
        if ([fileManager fileExistsAtPath:fName] == NO) {
            if ([fileManager createFileAtPath:fName contents:nil attributes:nil] == YES) {
                outputFileHandle = [NSFileHandle fileHandleForWritingAtPath:fName];
            } else {
                fatal("FEMJob:FEMJob_saveResult", "Can't create result file.");
            }
        } else { // File already exists, erase it and start from scratch
            [fileManager removeItemAtPath:fName error:nil];
            if ([fileManager createFileAtPath:fName contents:nil attributes:nil] == YES) {
                outputFileHandle = [NSFileHandle fileHandleForWritingAtPath:fName];
            } else {
                fatal("FEMJob:FEMJob_saveResult", "Can't create result file.");
            }
        }
        // The first time, we start by writing the header
        [post writeString:@"ACSII 1" toFileHandle:outputFileHandle]; [outputFileHandle writeData:newLineBuff];
        [post writeString:@"!File started at: " toFileHandle:outputFileHandle];
        NSString *dateString = [NSString stringWithCString:dateAndTime() encoding:NSASCIIStringEncoding];
        [post writeString:dateString toFileHandle:outputFileHandle];
        [outputFileHandle writeData:newLineBuff];
        
        [post writeString:@"Degrees of freedom:" toFileHandle:outputFileHandle];
        [outputFileHandle writeData:newLineBuff];
        dofs = 0;
        for (FEMVariable *variable in mesh.variables) {
            if (variable.output == YES) {
                varContainers = variable.getContainers;
                if (variable.dofs > 1 && varContainers->sizeValues > 1) {
                    if ([variable.name length] < 10 || [[variable.name substringToIndex:10] isEqualToString:@"coordinate"] == NO || freeSurfaceFlag == YES) {
                        [post writeString:variable.name toFileHandle:outputFileHandle]; [outputFileHandle writeData:spaceBuff];
                        [post writeInteger:variable.dofs toFileHandle:outputFileHandle]; [outputFileHandle writeData:spaceBuff];
                        [post writeString:@" :fs" toFileHandle:outputFileHandle];
                        [outputFileHandle writeData:newLineBuff];
                    }
                }
            }
        }
        
        for (FEMVariable *variable in mesh.variables) {
            if (variable.output == YES) {
                varContainers = variable.getContainers;
                if (variable.dofs == 1 && varContainers->sizeValues > 1) {
                    if ([variable.name length] < 10 || [[variable.name substringToIndex:10] isEqualToString:@"coordinate"] == NO || freeSurfaceFlag == YES) {
                        [post writeString:variable.name toFileHandle:outputFileHandle]; [outputFileHandle writeData:spaceBuff];
                        [post writeInteger:variable.dofs toFileHandle:outputFileHandle]; [outputFileHandle writeData:spaceBuff];
                        [post writeString:@" :fs" toFileHandle:outputFileHandle];
                        if (variable.dofs == 1) dofs++;
                        [outputFileHandle writeData:newLineBuff];
                    }
                }
            }
        }
        [post writeString:@"Total Dofs: " toFileHandle:outputFileHandle];
        [post writeInteger:dofs toFileHandle:outputFileHandle];
        [outputFileHandle writeData:newLineBuff];
        [post writeString:@"Number of Nodes: " toFileHandle:outputFileHandle];
        [post writeInteger:mesh.numberOfNodes toFileHandle:outputFileHandle];
        [outputFileHandle writeData:newLineBuff];
    } else {
        outputFileHandle = [NSFileHandle fileHandleForWritingAtPath:fName];
        [outputFileHandle seekToEndOfFile];
    }
    
    [post writeString:@"Time: " toFileHandle:outputFileHandle];
    [post writeInteger:mesh.savesDone+1 toFileHandle:outputFileHandle]; [outputFileHandle writeData:spaceBuff];
    [post writeInteger:time toFileHandle:outputFileHandle]; [outputFileHandle writeData:spaceBuff];
    [post writeDouble:simulationTime toFileHandle:outputFileHandle];
    [outputFileHandle writeData:newLineBuff];
    
    // Write data to disk
    prev = allocateVariableContainer();
    for (FEMVariable *variable in mesh.variables) {
        varContainers = variable.getContainers;
        if (variable.output == YES && variable.dofs == 1 && varContainers->sizeValues > 1) {
            if ([variable.name length] < 10 || [[variable.name substringToIndex:10] isEqualToString:@"coordinate"] == NO || freeSurfaceFlag == YES) {
                if (saveAll == YES || variable.valuesChanged == YES) {
                    [post writeString:variable.name toFileHandle:outputFileHandle];
                    [outputFileHandle writeData:newLineBuff];
                    // Permutations...
                    if (varContainers->Perm == NULL) {
                        [post writeString:@"Perm: NULL" toFileHandle:outputFileHandle];
                        [outputFileHandle writeData:newLineBuff];
                    } else {
                        sameAsPrev = NO;
                        if (varContainers->Perm == prev->Perm) {
                            sameAsPrev = YES;
                        } else if (prev->Perm != NULL) {
                            if (varContainers->sizePerm == prev->sizePerm) {
                                all = YES;
                                for (i=0; i<varContainers->sizePerm; i++) {
                                    if (varContainers->Perm[i] != prev->Perm[i]) {
                                        all = NO;
                                        break;
                                    }
                                }
                                if (all == YES) sameAsPrev = YES;
                            }
                        }
                        if (sameAsPrev == YES) {
                            [post writeString:@"Perm: use previous" toFileHandle:outputFileHandle];
                            [outputFileHandle writeData:newLineBuff];
                        } else {
                            prev->Perm = varContainers->Perm;
                            n = 0;
                            for (i=0; i<varContainers->sizePerm; i++) {
                                if (varContainers->Perm[i] >= 0) n++;
                            }
                            [post writeString:@"Perm: " toFileHandle:outputFileHandle];
                            [post writeInteger:varContainers->sizePerm toFileHandle:outputFileHandle]; [outputFileHandle writeData:spaceBuff];
                            [post writeInteger:n toFileHandle:outputFileHandle];
                            [outputFileHandle writeData:newLineBuff];
                            for (i=0; i<varContainers->sizePerm; i++) {
                                if (varContainers->Perm[i] >= 0) {
                                    [post writeInteger:i+1 toFileHandle:outputFileHandle]; [outputFileHandle writeData:spaceBuff];
                                    [post writeInteger:varContainers->Perm[i]+1 toFileHandle:outputFileHandle];
                                    [outputFileHandle writeData:newLineBuff];
                                }
                            }
                        }
                    }
                    
                    if (varContainers->Perm != NULL) {
                        n = varContainers->sizePerm;
                    } else {
                        n = mesh.numberOfNodes;
                    }
                    
                    // and the values...
                    for (i=0; i<n; i++) {
                        k = 0;
                        if (varContainers->Perm != NULL) k = varContainers->Perm[i];
                        if (k >= 0) {
                            [post writeDouble:varContainers->Values[k] toFileHandle:outputFileHandle];
                            [outputFileHandle writeData:newLineBuff];
                        }
                    }
                    variable.valuesChanged = NO;
                }
            }
        }
    }
    
    [outputFileHandle closeFile];
    mesh.savesDone++;
    saveCount = mesh.savesDone;
    free(prev);
    
    return saveCount;
}


/***********************************************************************************
    Saves current time step to external files
***********************************************************************************/
-(void)FEMJob_saveCurrent:(FEMModel * __nonnull)model currentStep:(int)currentStep post:(FEMPost * __nonnull)post {
    
    int j, k;
    BOOL found, binaryOutput, saveAll, eigenAnal=NO;
    NSString *outputFile;
    FEMListUtilities *listUtilities;
    variableArraysContainer *varContainers = NULL;
    
    listUtilities = [FEMListUtilities sharedListUtilities];
    
    outputFile = [listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"output file" info:&found];
    if (found == YES) {
        // TODO: add support for parallel run
        
        binaryOutput = [listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"binary output" info:&found];
        if (found == NO) binaryOutput = NO;
        
        saveAll = ([listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"omit unchanged variables in output" info:&found] == YES) ? NO : YES;
        if (found == NO) saveAll = YES;
        
        for (FEMMesh *mesh in model.meshes) {
            if (mesh.outputActive == YES) {
                if ([mesh.name length] > 0) {
                    [_outputName setString:mesh.name];
                    [_outputName appendString:@"/"];
                    [_outputName appendString:outputFile];
                } else {
                    [_outputName setString:outputFile];
                }
                
                for (FEMSolution *solution in model.solutions) {
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
                                        _savedSteps = [self FEMJob_saveResult:_outputName model:model mesh:mesh time:j+1 simulationTime:_sTime[0] binary:binaryOutput saveAll:saveAll freeSurface:NULL post:post];
                                    }
                                } else {
                                    j = min(currentStep, varContainers->size1EigenVectors);
                                    if (solution.matrix.complexMatrix == YES) {
                                        for (k=0; k<varContainers->sizeValues/2; k++) {
                                            varContainers->Values[2*k] = creal(varContainers->EigenVectors[j-1][k]);
                                            varContainers->Values[2*k+1] = cimag(varContainers->EigenVectors[j-1][k]);
                                        }
                                    } else {
                                        for (k=0; k<varContainers->sizeValues; k++) {
                                            varContainers->Values[k] = creal(varContainers->EigenVectors[j-1][k]);
                                        }
                                    }
                                    _savedSteps = [self FEMJob_saveResult:_outputName model:model mesh:mesh time:currentStep simulationTime:_sTime[0] binary:binaryOutput saveAll:saveAll freeSurface:NULL post:post];
                                }
                                memset( varContainers->Values, 0.0, varContainers->sizeValues*sizeof(double) );
                            }
                        }
                    }
                }
                
                if (eigenAnal == NO) {
                    _savedSteps = [self FEMJob_saveResult:_outputName model:model mesh:mesh time:round(_sTime[0]) simulationTime:_sTime[0] binary:binaryOutput saveAll:saveAll freeSurface:NULL post:post];
                }
            }
        }
    } else {
        for (FEMMesh *mesh in model.meshes) {
            if (mesh.outputActive == YES) mesh.savesDone++;
        }
    }
    [self FEMJob_saveToPostModel:model currentStep:currentStep post:post];
}
  
#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {

        // Instanciate a kernel class for this job
        _core = [FEMCore sharedCore];
        
        // Instanciate an element description class for this job
        _elementDescription = [FEMElementDescription sharedElementDescription];
        
        // Instanciate a list utilities class for this job
        _listUtilities = [FEMListUtilities sharedListUtilities];

        _transient = NO;
        _outputName = [NSMutableString stringWithString:@""];
        _postName = [NSMutableString stringWithString:@""];
        
        _firstLoad = YES;
        _firstTime = YES;
        _silent = NO;
        _version = NO;
        _gotModelName = NO;
        
        _timeSteps = NULL;
        _outputIntervals = NULL;
        _sTime = NULL;
        _sStep = NULL;
        _sInterval = NULL;
        _sSize = NULL;
        _sPrevSizes = NULL;
        _steadyIt = NULL;
        _nonLinIt = NULL;
        _timeStepSizes = NULL;
    }
    
    return self;
}

-(void)deallocation {
    
    [_core deallocation];
    [FEMCore selfDestruct];
    
    [_elementDescription deallocation];
    [FEMElementDescription selfDestruct];
    
    [_listUtilities deallocation];
    [FEMListUtilities selfDestruct];
    
    [_model deallocation];
    
    if (_timeStepSizes != NULL) {
        free_ivector(_timeSteps, 0, _sizeTimeSteps-1);
        _timeSteps = NULL;
    }
    if (_outputIntervals != NULL) {
        free_ivector(_outputIntervals, 0, _sizeOutputIntervals-1);
        _outputIntervals = NULL;
    }
    if (_sTime != NULL) {
        free_dvector(_sTime, 0, 0);
        _sTime = NULL;
    }
    if (_sStep != NULL) {
        free_dvector(_sStep, 0, 0);
        _sStep = NULL;
    }
    if (_sInterval != NULL) {
        free_dvector(_sInterval, 0, 0);
        _sInterval = NULL;
    }
    if (_sSize != NULL) {
        free_dvector(_sSize, 0, 0);
        _sSize = NULL;
    }
    if (_steadyIt != NULL) {
        free_dvector(_steadyIt, 0, 0);
        _steadyIt = NULL;
    }
    if (_nonLinIt != NULL) {
        free_dvector(_nonLinIt, 0, 0);
        _nonLinIt = NULL;
    }
    if (_sPrevSizes != NULL) {
        free_dmatrix(_sPrevSizes, 0, 0, 0, 4);
        _sPrevSizes = NULL;
    }
    if (_timeStepSizes == NULL) {
        free_dmatrix(_timeStepSizes, 0, _size1TimeStepSizes-1, 0, _size2TimeStepSizes-1);
        _timeStepSizes = NULL;
    }
    
    // Finally deallocate the Gauss quadrature points
    GaussQuadratureDeallocation();
}

-(void)runWithInitialize:(int)initialize {
    
    int i, j, k, extrudeLevels, interval, minVal, timeStep=0;
    NSString *eq, *when;
    FEMMesh *extrudedMesh;
    listBuffer listBuffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    BOOL colorMesh=NO, execThis, found, parallelAssembly=NO;
    
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    FEMPost *post = [[FEMPost alloc] init];
    
    // TODO: Add support for parallel run
    
    NSArray *args = [[NSProcessInfo processInfo] arguments];
    if (_firstTime == YES) {
        if ([args count]-3 > 0) {
            for (NSString *argument in args) {
                _silent = (_silent == YES || ([argument isEqualToString:@"-s"] == YES || [argument isEqualToString:@"--silent"] == YES)) ? YES : NO;
                _version = (_version == YES || ([argument isEqualToString:@"-v"] == YES || [argument isEqualToString:@"--version"] == YES)) ? YES : NO;
            }
        }
        
        if (_silent == NO) {
            fprintf(stdout, "JOB: \n");
            fprintf(stdout, "JOB: =====================================================================\n");
            fprintf(stdout, "JOB: Saino finite element software, Welcome.\n");
            fprintf(stdout, "JOB: Saino is an object oriented, GPU based implementation of ElmerSolver.\n");
            fprintf(stdout, "JOB: This program is free software under (L)GPL.\n");
            fprintf(stdout, "JOB: Copyright 15 April 2011 - ScienceSoul Hakime Seddik.\n");
            fprintf(stdout, "JOB: Copyright 1st April 1995 - CSC IT Center for Science Ltd.\n");
            // TODO: Add support for parallel info output
            fprintf(stdout, "JOB: =====================================================================\n");
        }
        
        if (_version == YES) return;
        _firstTime = NO;
    }
   
    if ([args count]-3 > 0) { // TODO: Add support for parallel run
        self.modelName = [NSString stringWithString:args[1]];
        if ([self.modelName characterAtIndex:0] != '-') {
            _gotModelName = YES;
        }
    }
    
    if (_gotModelName == NO) {
        NSFileHandle *startFile = [NSFileHandle fileHandleForReadingAtPath:@"SAINO_STARTINFO"];
        if (startFile == nil) {
            fatal("FEMJob:runWithInitialize", "SAINO_STARTINFO file not found.");
        }
        NSData *buffer = [startFile readDataToEndOfFile];
        self.modelName = [[NSString alloc] initWithData:buffer encoding:NSUTF8StringEncoding];
        [startFile closeFile];
    }
    
    if (initialize == 1) {
        [self.model deallocation];
        _firstLoad = YES;
    }
    
    while (1) {
        if (initialize == 2) goto jump;
        
        if (_firstLoad == YES) {
            if (_silent == NO) {
                fprintf(stdout, "JOB: \n");
                fprintf(stdout, "JOB: \n");
                fprintf(stdout, "JOB: ---------------------------------------------------------------------\n");
                fprintf(stdout, "JOB: Reading model: %s.\n", [self.modelName UTF8String]);
            }
            
            self.model = [[FEMModel alloc] init];
            self.model.mdf = [[FileReader alloc] initWithFilePath:self.modelName];
            if (self.model.mdf == nil) {
                fprintf(stderr, "FEMJob:runWithInitialize: Unable to find model description file [' %s '].\n", [self.modelName UTF8String]);
                fatal("FEMJob:runWithInitialize");
            }
            [self.model loadModelName:self.modelName boundariesOnly:NO dummy:NULL dummy:NULL];
            
            // Check whether we need to color the mesh if parallel assembly is required
            
            for (FEMSolution *solution in self.model.solutions) {
                if ([(solution.solutionInfo)[@"parallel assembly"] boolValue] == YES) {
                    parallelAssembly = YES;
                    colorMesh = [(solution.solutionInfo)[@"color mesh"] boolValue];
                    break;
                }
            }
            
            // Optionally perform simple extrusion to increase the dimension of the mesh
            extrudeLevels = [listUtilities listGetInteger:self.model inArray:self.model.simulation.valuesList forVariable:@"extruded mesh levels" info:&found minValue:NULL maxValue:NULL];
            if (extrudeLevels > 1) {
                if (parallelAssembly == YES && colorMesh == NO) {
                    fprintf(stderr, "FEMJob:runWithInitialize: if the < color mesh> option is set to < NO >, the extruded mesh won't be colored.\n");
                    fprintf(stderr, "FEMJob:runWithInitialize: this means that the non-extruded mesh should already be colored which does not make sense.\n");
                    fprintf(stderr, "FEMJob:runWithInitialize: if the mesh is extruded, it should not be colored before extrusion.\n");
                    fatal("FEMJob:runWithInitialize");
                }
                FEMMeshUtils *meshUtils = [[FEMMeshUtils alloc] init];
                extrudedMesh = [meshUtils extrudeMesh:self.model.meshes[0] inLevels:extrudeLevels-2 model:self.model];
                for (FEMSolution *solution in self.model.solutions) {
                    if (solution.mesh == self.model.meshes[0]) {
                        solution.mesh = extrudedMesh;
                    }
                }
                // Deallocate and replace the 2D mesh we extruded with the resulting 3D mesh
                // since presumably it won't be used anymore.
                FEMMesh *mesh = self.model.meshes[0];
                [mesh deallocation];
                [self.model.meshes replaceObjectAtIndex:0 withObject:extrudedMesh];
                
                // If periodic BC given, compute the boundary mesh projector
                i = 0;
                for (FEMBoundaryCondition *boundary in self.model.boundaryConditions) {
                    k = [listUtilities listGetInteger:self.model inArray:boundary.valuesList forVariable:@"periodic bc" info:&found minValue:NULL maxValue:NULL];
                    if (found == YES) {
                        boundary.pMatrix = [meshUtils periodicProjectorInModel:self.model forMesh:extrudedMesh masterBoundary:i targetBoundary:k-1 galerking:NULL];
                    }
                    i++;
                }
            }
            // If parallel assembly and if required, color the mesh here
            // Store the result in the mesh directory so that we do the coloring only once
            if (parallelAssembly == YES && colorMesh == YES) {
                FEMMeshUtils *meshUtils = [[FEMMeshUtils alloc] init];
                FEMMesh *mesh = self.model.meshes[0];
                // Fist check whether we already have done the coloring before for this mesh
                BOOL alreadyColored = NO;
                NSFileManager *fileManager = [NSFileManager defaultManager];
                NSString *file1 = [[self.model.meshDir stringByAppendingPathComponent:self.model.meshName] stringByAppendingPathComponent:@"mesh.colors"];
                NSString *file2 = [[self.model.meshDir stringByAppendingPathComponent:self.model.meshName] stringByAppendingPathComponent:@"mesh.colored_elements"];
                if ([fileManager fileExistsAtPath:file1] == YES && [fileManager fileExistsAtPath:file2] == YES) alreadyColored = YES;
                if (alreadyColored == NO) {
                    [meshUtils colorMesh:mesh];
                    // Save the coloring
                    [meshUtils saveColoredMesh:mesh meshdir:self.model.meshDir meshName:self.model.meshName elementsFileName:@"mesh.colored_elements" saveAllElementData:NO colorFileName:@"mesh.colors"];
                }
                [meshUtils readColoredMesh:mesh name:self.model.meshName directory:self.model.meshDir readElementsFromFile:YES];
            }
            if (_silent == NO) {
                fprintf(stdout, "JOB: ---------------------------------------------------------------------\n");
            }
        } else {
            if (initialize == 3) {
                if (self.model.mdf != nil) {
                    [self.model.mdf closeHandle];
                    self.model.mdf = nil;
                }
                self.model.mdf = [[FileReader alloc] initWithFilePath:self.modelName];
            }
            
            // TODO: Add support for reloading MDF if we really need that
            // For now, we always get out of the loop
            if(/* DISABLES CODE */ (YES)) break;
            
            for (FEMMesh *mesh in self.model.meshes) {
                mesh.savesDone = 0;
            }
        }
        
    jump:
        [listUtilities addLogicalInClassList:self.model.simulation theVariable:@"initialization phase" withValue:YES];
        
        // Check for transient case
        eq = [listUtilities listGetString:self.model inArray:self.model.simulation.valuesList forVariable:@"simulation type" info:&found];
        _scanning = ([eq isEqualToString:@"scanning"] == YES) ? YES : NO;
        _transient = ([eq isEqualToString:@"transient"] == YES) ? YES : NO;
        
        // Figure out what (flow, heat, stress...) shoud be computed and get memory for the dofs
        [self FEMJob_addSolutionsModel:self.model];
        
         // Time integration and/or steady state steps
        if (_transient == YES || _scanning == YES) {
            found = [listUtilities listGetIntegerArray:self.model inArray:self.model.simulation.valuesList forVariable:@"time step intervals" buffer:&listBuffer];
            if (found == NO) fatal("FEMJob:runWithInitialize", "Keyword > time step intervals < must be defined for transient and scanning simulations.");
            
            _timeSteps = intvec(0, listBuffer.m-1);
            _sizeTimeSteps = listBuffer.m;
            memcpy(_timeSteps, listBuffer.ivector, listBuffer.m*sizeof(int));
            free_ivector(listBuffer.ivector, 0, listBuffer.m-1);
            listBuffer.ivector = NULL;
            
            found = [listUtilities listGetConstRealArray:self.model inArray:self.model.simulation.valuesList forVariable:@"time step sizes" buffer:&listBuffer];
            if (found == NO) {
                if (_scanning == YES || [listUtilities listCheckPresentVariable:@"time step size" inArray:self.model.simulation.valuesList] == YES) {
                    _timeStepSizes = doublematrix(0, _sizeTimeSteps-1, 0, 0);
                    _size1TimeStepSizes = _sizeTimeSteps;
                    _size2TimeStepSizes = 1;
                    for (i=0; i<_size1TimeStepSizes; i++) {
                        for (j=0; j<_size2TimeStepSizes; j++) {
                            _timeStepSizes[i][j] = 1.0;
                        }
                    }
                }
            } else {
                _timeStepSizes = doublematrix(0, listBuffer.m-1, 0, listBuffer.n-1);
                _size1TimeStepSizes = listBuffer.m;
                _size2TimeStepSizes = listBuffer.n;
                for (i=0; i<listBuffer.m; i++) {
                    for (j=0; j<listBuffer.n; j++) {
                        _timeStepSizes[i][j] = listBuffer.matrix[i][j];
                    }
                }
                free_dmatrix(listBuffer.matrix, 0, listBuffer.m-1, 0, listBuffer.n-1);
                listBuffer.matrix = NULL;
            }
            
            _timeIntervals = _sizeTimeSteps;
            
            minVal = 1;
            _coupledMaxIter = [listUtilities listGetInteger:self.model inArray:self.model.simulation.valuesList forVariable:@"steady state maximum iterations" info:&found minValue:&minVal maxValue:NULL];
            if (found == NO) _coupledMaxIter = 1;
        } else { // Steady state
            _timeSteps = intvec(0, 0);
            _sizeTimeSteps = 1;
            
            minVal = 1;
            _timeSteps[0] = [listUtilities listGetInteger:self.model inArray:self.model.simulation.valuesList forVariable:@"steady state maximum iterations" info:&found minValue:&minVal maxValue:NULL];
            if (found == NO) _timeSteps[0] = 1;
            
            _timeStepSizes = doublematrix(0, 0, 0, 0);
            _size1TimeStepSizes = 1;
            _size2TimeStepSizes = 1;
            _timeStepSizes[0][0] = 1.0;
            
            _coupledMaxIter = 1;
            _timeIntervals = 1;
        }
        
        if (_firstLoad == YES) {
            _sTime = doublevec(0, 0);
            _sizeSTime = 1;
            
            _sStep = doublevec(0, 0);
            _sizeSStep = 1;
            
            _sInterval = doublevec(0, 0);
            _sizeSInterval = 1;
            
            _sSize = doublevec(0, 0);
            _sizeSSize = 1;
            
            _steadyIt = doublevec(0, 0);
            _sizeSteadyIt = 1;
            
            _nonLinIt = doublevec(0, 0);
            _sizeNonLinIt = 1;
            
            _sPrevSizes = doublematrix(0, 0, 0, 4);
            _size1SPrevSizes = 1;
            _size2SPrevSizes = 5;
        }
        
        _dt = 0.0;
        
        memset( _sTime, 0.0, 1*sizeof(double) );
        memset( _sStep, 0.0, 1*sizeof(double) );
        memset( _sSize, _dt, 1*sizeof(double) );
        memset( _sInterval, 0.0, 1*sizeof(double) );
        memset( _steadyIt, 0.0, 1*sizeof(double) );
        memset( _nonLinIt, 0.0, 1*sizeof(double) );
        memset( *_sPrevSizes, 0.0, (1*5)*sizeof(double) );
        
        _coupledMinIter = [listUtilities listGetInteger:self.model inArray:self.model.simulation.valuesList forVariable:@"steady state min iterations" info:&found minValue:NULL maxValue:NULL];
        
        // Add coordinates and simulation time to list of variables so that coordinate dependent
        // parameter computing methods can ask for them
        if (_firstLoad == YES) [self FEMJob_addMeshCoordinatesAndTimeModel:self.model];
        
        // Get output file options
        found = [listUtilities listGetIntegerArray:self.model inArray:self.model.simulation.valuesList forVariable:@"output intervals" buffer:&listBuffer];
        if (found == NO) {
            _outputIntervals = intvec(0, _sizeTimeSteps-1);
            _sizeOutputIntervals = _sizeTimeSteps;
            for (i=0; i<_sizeOutputIntervals; i++) {
                _outputIntervals[i] = 1;
            }
        } else {
            _outputIntervals = intvec(0, listBuffer.m-1);
            _sizeOutputIntervals = listBuffer.m;
            memcpy(_outputIntervals, listBuffer.ivector, listBuffer.m*sizeof(int));
            free_ivector(listBuffer.ivector, 0, listBuffer.m-1);
            listBuffer.ivector = NULL;
        }
        
        // Initial conditions
        if (_firstLoad == YES) [self FEMJob_setInitialConditionsModel:self.model];
        
        // Compute the total number of steps that will be saved to the file.
        // Particularly look if the last step will be saved or if it has to be saved separately.
        _totalTimeSteps = 0;
        _lastSaved = YES;
        for (interval=0; interval<_timeIntervals; interval++) {
            for (timeStep=1; timeStep<=_timeSteps[interval]; timeStep++) {
                if (_outputIntervals[interval] == 0) continue;
                _lastSaved = NO;
                if ((timeStep-1) % _outputIntervals[interval] == 0) {
                    _lastSaved = YES;
                    _totalTimeSteps++;
                }
            }
        }
        
        for (FEMSolution *solution in self.model.solutions) {
            if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                when = (solution.solutionInfo)[@"invoke solution computer"];
                if ([when isEqualToString:@"after simulation"] == YES || [when isEqualToString:@"after all"] == YES) {
                    _lastSaved = NO;
                }
            } else {
                if (solution.solutionSolveWhen == SOLUTION_SOLVE_AFTER_ALL) {
                    _lastSaved = NO;
                }
            }
        }
        
        if (_lastSaved == NO) _totalTimeSteps++;
        if (_totalTimeSteps == 0) _totalTimeSteps = 1;
        
        [listUtilities addLogicalInClassList:self.model.simulation theVariable:@"initialization phase" withValue:NO];
        
        _firstLoad = NO;
        if (initialize == 1) break;
        
        // Here we actually start the simulation....
        // First go through time intervals
        [self FEMJob_runSimulation:self.model timeIntervals:_timeIntervals coupledMinIteration:_coupledMinIter coupleMaxIteration:_coupledMaxIter outputIntervals:_outputIntervals transient:_transient scanning:_scanning post:post];
        
        // Always save the last step to output
        if (_lastSaved == NO) {
            for (FEMSolution *solution in self.model.solutions) {
                if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
                execThis = (solution.solutionSolveWhen == SOLUTION_SOLVE_AHEAD_SAVE) ? YES : NO;
                if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                    when = (solution.solutionInfo)[@"invoke solution computer"];
                    execThis = ([when isEqualToString:@"before saving"] == YES) ? YES : NO;
                }
                if (execThis == YES) [self.core activateSolution:solution model:self.model timeStep:_dt transientSimulation:_transient];
            }
            
            [self FEMJob_saveToPostModel:self.model currentStep:0 post:post];
            [self FEMJob_saveCurrent:self.model currentStep:timeStep post:post];
            
            for (FEMSolution *solution in self.model.solutions) {
                if (solution.hasBuiltInSolution == NO && solution.plugInPrincipalClassInstance == nil) continue;
                 execThis = (solution.solutionSolveWhen == SOLUTION_SOLVE_AFTER_SAVE) ? YES : NO;
                if ((solution.solutionInfo)[@"invoke solution computer"] != nil) {
                    when = (solution.solutionInfo)[@"invoke solution computer"];
                    execThis = ([when isEqualToString:@"after saving"] == YES) ? YES : NO;
                }
                if (execThis == YES) [self.core activateSolution:solution model:self.model timeStep:_dt transientSimulation:_transient];
            }
        }
        
        if (initialize >= 2) break;
    }
    
    // -----------------------------------------------------------------------
    // THIS IS THE END LITTLE APPRENTICE, SAY GOOD BYE TO THE DARK SAMOURAI...
    // -----------------------------------------------------------------------
    if (initialize != 1) fprintf(stdout, "JOB: *** ALL DONE ***\n");
    
    // TODO: add support for parallel runs
    fprintf(stdout, "JOB: The end.\n");
    
#ifdef TEST
    FEMTest *test = [FEMTest sharedTest];
    if (test.do_heatq == YES) {
        test.heatq_allDone = YES;
    } else if (test.do_step_stokes == YES) {
        test.step_stokes_allDone = YES;
    } else if (test.do_natural_convection == YES) {
        test.natural_convection_allDone = YES;
    } else if (test.do_ismip_hom_A010 == YES) {
        test.ismip_hom_A010_allDone = YES;
    } else if (test.do_ismip_hom_A010_gpu == YES) {
        test.ismip_hom_A010_gpu_allDone = YES;
    }
#endif
}

@end
