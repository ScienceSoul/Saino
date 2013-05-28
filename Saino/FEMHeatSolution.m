//
//  FEMHeatEquation.m
//  Saino
//
//  Created by Seddik hakime on 19/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMHeatSolution.h"
#import "FEMKernel.h"
#import "FEMBoundaryCondition.h"
#import "FEMListUtilities.h"
#import "FEMUtilities.h"
#import "FEMElementUtils.h"
#import "FEMBodyForce.h"
#import "FEMMaterial.h"
#import "FEMEquation.h"
#import "Utils.h"

@implementation FEMHeatSolution {
    
    BOOL _allocationDone;
    int _doneTime;
    int *_indexes;
    int *_saveIndexes;
    double _powerScaling;
    double _prevPowerScaling;
    double *_u, *_v, *_w;
    double **_mu;
    double *_pressure;
    double *_dPressureDt;
    double *_pressureCoeff;
    double *_density;
    double *_work;
    double *_latentHeat;
    double **_phaseVelocity;
    double *_electricConductivity;
    double *_permeability;
    double *_viscosity;
    double *_c0;
    double *_heatTransferCoeff;
    double *_heatExpansionCoeff;
    double *_referenceTemperature;
    double **_mass;
    double *_localTemperature;
    double *_heatCapacity;
    double *_enthalpy;
    double *_nodalEmissivity;
    double *_gasConstant;
    double *_aText;
    double ***_heatConductivity;
    double **_stiff;
    double *_load;
    double *_force;
    double *_timeForce;
    double *_perfusionRate;
    double *_perfusionDensity;
    double *_perfusionHeatCapacity;
    double *_perfusionRefTemperature;
    double *_heaterArea;
    double *_heaterDensity;
    double *_heaterSource;
    double *_heaterScaling;
    double *_heaterTarget;
    double *_xx;
    double *_yy;
    double *_forceHeater;
    bool *_smarterHeaters;
    bool *_integralHeaters;
    Nodes_t *_elementNodes;
}

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _allocationDone = NO;
        _doneTime = 0;
        _powerScaling = 1.0;
        _prevPowerScaling = 1.0;
    }
    
    return self;
}

-(void)fieldSolutionComputer:(FEMSolution *)solution model:(FEMModel *)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
    int i, j, k, l, n, t, bf_id, body_id, compressibilityModel, eq_id, iter, localNodes, mat_id, nsdofs, nonLinearIter, newtonIter, smartHeaterNode, smartHeaterBC;
    int *tempPerm, *flowPerm;
    double dist, dt, dt0, jx, jy, jz, cumulativeTime, newtonTol, meltPoint, minDist, nonLinearTol, norm, powerTimeScale, relax, s, saveRelax,
           smartTol, stefanBoltzmann, sum;
    double controlPoint[3], *temperature, *forceVector, *flowSolution, *prevSolution, *tSolution=NULL, *tSolution1;
    BOOL checkLatentHeatRelease=NO, constantBulk=NO, found, gotIt, isRadiation, firstTime, gotMeltPoint, integralHeaterControl, heaterControlLocal,
         newtonLinearization, phaseChange=NO, saveBulk, smartHeaterControl, smartHeaterAverage, smartTolReached, stabilize, transientAssembly,
         transientHeaterControl, useBubbles;
    NSString *convectionField, *stabilizeFlag, *convectionFlag, *compressibilityFlag, *phaseModel;
    Element_t *elements=NULL, *element = NULL;
    Nodes_t *meshNodes=NULL;
    FEMKernel *kernel;
    FEMMesh *mesh;
    FEMVariable *flowSol, *densitySol;
    FEMListUtilities *listUtilities;
    FEMUtilities *utilities;
    FEMElementUtils *elementUtils;
    FEMBodyForce *bodyForceAtID = nil;
    FEMBoundaryCondition *boundaryConditionAtID = nil;
    FEMMaterial *materialAtID = nil;
    FEMEquation *equationAtID = nil;
    matrixArraysContainer *matContainers=NULL;
    variableArraysContainer *tempContainers=NULL, *flowSolContainers=NULL;
    listBuffer realWork = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    kernel = [FEMKernel sharedKernel];
    
    listUtilities = [[FEMListUtilities alloc] init];
    utilities = [[FEMUtilities alloc] init];
    
    mesh = (FEMMesh *)model.mesh;
    elements = mesh.getElements;
    meshNodes = mesh.getNodes;
    
    isRadiation = NO;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        isRadiation = ([listUtilities listCheckPresentVariable:@"radiation" inArray:boundaryCondition.valuesList] == YES) ? YES : NO;
        if (isRadiation == YES) break;
    }
    
    if (isRadiation == YES) {
        //TODO: implement the radiation factors computation
    }
    
    // Get variables needed for solution
    
    if (solution.matrix == nil) return;
    
    matContainers = solution.matrix.getContainers;
    forceVector = matContainers->RHS;
    
    tempContainers = solution.variable.getContainers;
    tempPerm = tempContainers->Perm;
    temperature = tempContainers->Values;
    
    localNodes = 0;
    for (i=0; i<tempContainers->sizePerm; i++) {
        if (tempPerm[i] >= 0) localNodes++;
    }
    if (localNodes <= 0) return;
    
    if ((solution.solutionInfo)[@"temperature convection field"] != nil) {
        convectionField = (solution.solutionInfo)[@"temperature convection field"];
        flowSol = [utilities getVariableFrom:solution.mesh.variables model:model name:convectionField onlySearch:NULL maskName:nil info:&found];
    } else {
        flowSol = [utilities getVariableFrom:solution.mesh.variables model:model name:@"flow solution" onlySearch:NULL maskName:nil info:&found];
    }
    
    if (flowSol != nil) {
        flowSolContainers = flowSol.getContainers;
        flowPerm = flowSolContainers->Perm;
        nsdofs = flowSol.dofs;
        flowSolution = flowSolContainers->Values;
    }
    
    densitySol = [utilities getVariableFrom:solution.mesh.variables model:model name:@"density" onlySearch:NULL maskName:nil info:&found];
    
    // Check whether we have some heater controls. This will affect initialization stuff.
    smartHeaterControl = ([listUtilities listCheckPresentAnyBodyForce:model name:@"smart heater control"] == YES) ? YES : NO;
    integralHeaterControl = ([listUtilities listCheckPresentAnyBodyForce:model name:@"integral heat source"] == YES) ? YES : NO;
    
    // Allocate some permanent storage, this is done first time only
    if (_allocationDone == NO || solution.mesh.changed == YES) {
        n = solution.mesh.maxElementDofs;
        if (_allocationDone == YES) {
            free_ivector(_indexes, 0, n-1);
            free_ivector(_saveIndexes, 0, n-1);
            free_dvector(_u, 0, n-1);
            free_dvector(_v, 0, n-1);
            free_dvector(_w, 0, n-1);
            free_dmatrix(_mu, 0, 2, 0, n-1);
            free_dvector(_pressure, 0, n-1);
            free_dvector(_dPressureDt, 0, n-1);
            free_dvector(_pressureCoeff, 0, n-1);
            free_dvector(_density, 0, n-1);
            free_dvector(_work, 0, n-1);
            free_dvector(_latentHeat, 0, n-1);
            free_dmatrix(_phaseVelocity, 0, 2, 0, n-1);
            free_dvector(_electricConductivity, 0, n-1);
            free_dvector(_permeability, 0, n-1);
            free_dvector(_viscosity, 0, n-1);
            free_dvector(_c0, 0, n-1);
            free_dvector(_heatTransferCoeff, 0, n-1);
            free_dvector(_heatExpansionCoeff, 0, n-1);
            free_dvector(_referenceTemperature, 0, n-1);
            free_dmatrix(_mass, 0, (2*n)-1, 0, (2*n)-1);
            free_dvector(_localTemperature, 0, n-1);
            free_dvector(_heatCapacity, 0, n-1);
            free_dvector(_enthalpy, 0, n-1);
            free_dvector(_nodalEmissivity, 0, n-1);
            free_dvector(_gasConstant, 0, n-1);
            free_d3tensor(_heatConductivity, 0, 2, 0, 2, 0, n-1);
            free_dmatrix(_stiff, 0, (2*n)-1, 0, (2*n)-1);
            free_dvector(_load, 0, n-1);
            free_dvector(_force, 0, (2*n)-1);
            free_dvector(_timeForce, 0, (2*n)-1);
            free_dvector(_perfusionRate, 0, n-1);
            free_dvector(_perfusionDensity, 0, n-1);
            free_dvector(_perfusionHeatCapacity, 0, n-1);
            free_dvector(_perfusionRefTemperature, 0, n-1);
            free_dvector(_elementNodes->x, 0, n-1);
            free_dvector(_elementNodes->y, 0, n-1);
            free_dvector(_elementNodes->z, 0, n-1);
            free(_elementNodes);
        }
        _indexes = intvec(0, n-1);
        _saveIndexes = intvec(0, n-1);
        _u = doublevec(0, n-1);
        _v = doublevec(0, n-1);
        _w = doublevec(0, n-1);
        _mu = doublematrix(0, 2, 0, n-1);
        _pressure = doublevec(0, n-1);
        _dPressureDt = doublevec(0, n-1);
        _pressureCoeff = doublevec(0, n-1);
        _density = doublevec(0, n-1);
        _work = doublevec(0, n-1);
        _latentHeat = doublevec(0, n-1);
        _phaseVelocity = doublematrix(0, 2, 0, n-1);
        _electricConductivity = doublevec(0, n-1);
        _permeability = doublevec(0, n-1);
        _viscosity = doublevec(0, n-1);
        _c0 = doublevec(0, n-1);
        _heatTransferCoeff = doublevec(0, n-1);
        _heatExpansionCoeff = doublevec(0, n-1);
        _referenceTemperature = doublevec(0, n-1);
        _mass = doublematrix(0, (2*n)-1, 0, (2*n)-1);
        _localTemperature = doublevec(0, n-1);
        _heatCapacity = doublevec(0, n-1);
        _enthalpy = doublevec(0, n-1);
        _nodalEmissivity = doublevec(0, n-1);
        _gasConstant = doublevec(0, n-1);
        _heatConductivity = d3tensor(0, 2, 0, 2, 0, n-1);
        _stiff = doublematrix(0, (2*n)-1, 0, (2*n)-1);
        _load = doublevec(0, n-1);
        _force = doublevec(0, (2*n)-1);
        _timeForce = doublevec(0, (2*n)-1);
        _perfusionRate = doublevec(0, n-1);
        _perfusionDensity = doublevec(0, n-1);
        _perfusionHeatCapacity = doublevec(0, n-1);
        _perfusionRefTemperature = doublevec(0, n-1);
        _elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
        _elementNodes->x = doublevec(0, n-1);
        _elementNodes->y = doublevec(0, n-1);
        _elementNodes->z = doublevec(0, n-1);
        if (_indexes == NULL || _saveIndexes == NULL || _u == NULL || _v == NULL || _w == NULL || _mu == NULL || _pressure == NULL ||
            _dPressureDt == NULL || _pressureCoeff == NULL || _density == NULL || _work == NULL || _latentHeat == NULL || _phaseVelocity == NULL ||
            _electricConductivity == NULL || _permeability == NULL || _viscosity == NULL || _c0 == NULL || _heatTransferCoeff == NULL ||
            _heatExpansionCoeff == NULL || _referenceTemperature == NULL || _mass == NULL || _localTemperature == NULL || _heatCapacity == NULL ||
            _enthalpy == NULL || _nodalEmissivity == NULL || _gasConstant == NULL || _heatConductivity == NULL || _stiff == NULL || _load == NULL ||
            _force == NULL || _timeForce == NULL || _perfusionRate == NULL || _perfusionDensity == NULL || _perfusionHeatCapacity == NULL ||
            _perfusionRefTemperature == NULL || _elementNodes->x == NULL || _elementNodes->y == NULL || _elementNodes->z == NULL) {
            errorfunct("FEMHeatSolution", "Memory allocation error");
        }
        
        if (smartHeaterControl == YES || integralHeaterControl == YES) {
            n = model.numberOfBodyForces;
            if (_allocationDone == YES) {
                free_dvector(_heaterArea, 0, n-1);
                free_dvector(_heaterDensity, 0, n-1);
                free_dvector(_heaterSource, 0, n-1);
                free_dvector(_heaterScaling, 0, n-1);
                free_dvector(_heaterTarget, 0, n-1);
                free_bvector(_smarterHeaters, 0, n-1);
                free_bvector(_integralHeaters, 0, n-1);
            }
            _heaterArea = doublevec(0, n-1);
            _heaterDensity = doublevec(0, n-1);
            _heaterSource = doublevec(0, n-1);
            _heaterScaling = doublevec(0, n-1);
            _heaterTarget = doublevec(0, n-1);
            _smarterHeaters = boolvec(0, n-1);
            _integralHeaters = boolvec(0, n-1);
            if (_heaterArea == NULL || _heaterDensity == NULL || _heaterSource == NULL || _heaterScaling == NULL || _heaterTarget == NULL ||
                _smarterHeaters == NULL || _integralHeaters == NULL) {
                 errorfunct("FEMHeatSolution", "Memory allocation error");
            }
            memset( _smarterHeaters, false, n*sizeof(bool) );
            memset( _integralHeaters, false, n*sizeof(bool) );
        }
        
        if (smartHeaterControl == YES) {
            if (_allocationDone == YES) {
                free_dvector(_xx, 0, tempContainers->sizeValues-1);
                free_dvector(_yy, 0, tempContainers->sizeValues-1);
                free_dvector(_forceHeater, 0, tempContainers->sizeValues-1);
            }
            
            _xx = doublevec(0, tempContainers->sizeValues-1);
            _yy = doublevec(0, tempContainers->sizeValues-1);
            _forceHeater = doublevec(0, tempContainers->sizeValues-1);
            if (_xx == NULL || _yy == NULL || _forceHeater == NULL) {
                errorfunct("FEMHeatSolution", "Memory allocation error");
            }
            memset( _xx, 0.0, tempContainers->sizeValues*sizeof(double) );
            memset( _yy, 0.0, tempContainers->sizeValues*sizeof(double) );
            memset( _forceHeater, 0.0, tempContainers->sizeValues*sizeof(double) );
        }
        
        _allocationDone = YES;
    }
    
    // Do some additional initialization and go for it
    dt = timeStep;
    if (isRadiation == YES) {
        stefanBoltzmann = [listUtilities listGetConstReal:model inArray:model.constants.valuesList forVariable:@"stefan boltzmann" info:&found minValue:NULL maxValue:NULL];
    }
    
    if ((solution.solutionInfo)[@"stabilize"] != nil) {
        stabilize = [(solution.solutionInfo)[@"stabilize"] boolValue];
    }
    
    if ((solution.solutionInfo)[@"bubbles"] != nil) {
        useBubbles = [(solution.solutionInfo)[@"bubbles"] boolValue];
    } else useBubbles = YES;
    
    if ((solution.solutionInfo)[@"stabilization method"] != nil) {
        stabilizeFlag = (solution.solutionInfo)[@"stabilization method"];
    }
    if ([stabilizeFlag isEqualToString:@"vms"] == YES) {
        stabilize = NO;
        useBubbles = NO;
    } else if ([stabilizeFlag isEqualToString:@"stabilized"] == YES) {
        stabilize = YES;
        useBubbles = NO;
    } else if ([stabilizeFlag isEqualToString:@"bubbles"] == YES) {
        stabilize = NO;
        useBubbles = YES;
    }
    
    if ((solution.solutionInfo)[@"nonlinear system max iterations"] != nil) {
        nonLinearIter = [(solution.solutionInfo)[@"nonlinear system max iterations"] intValue];
    } else nonLinearIter = 1;
    
    if ((solution.solutionInfo)[@"nonlinear system convergence tolerance"] != nil) {
        nonLinearTol = [(solution.solutionInfo)[@"nonlinear system convergence tolerance"] doubleValue];
    }
    
    if (isRadiation == YES) {
        if ((solution.solutionInfo)[@"nonlinear system newton after tolerance"] != nil) {
            newtonTol = [(solution.solutionInfo)[@"nonlinear system newton after tolerance"] doubleValue];
        } else newtonTol = 1.0;
        if ((solution.solutionInfo)[@"nonlinear system newton after iterations"] != nil) {
            newtonIter = [(solution.solutionInfo)[@"nonlinear system newton after iterations"] intValue];
        } else newtonIter = 0;
    } else {
        newtonTol = 1.0;
        newtonIter = 0;
    }
    if (newtonIter == 0) newtonLinearization = YES;
    
    if ((solution.solutionInfo)[@"nonlinear system relaxation factor"] != nil) {
        relax = [(solution.solutionInfo)[@"nonlinear system relaxation factor"] doubleValue];
    } else relax = 1.0;
    
    transientAssembly = transient;
    found = NO;
    if ((solution.solutionInfo)[@"steady state transition time step"] != nil) {
        dt0 = [(solution.solutionInfo)[@"steady state transition time step"] doubleValue];
        found = YES;
    } else {
        if ((solution.solutionInfo)[@"smart heater time scale"] != nil) {
            dt0 = [(solution.solutionInfo)[@"smart heater time scale"] doubleValue];
            found = YES;
        }
    }
    if (found == YES && dt > dt0) transientAssembly = NO;
    
    transientHeaterControl = NO;
    if (smartHeaterControl == YES) {
        
        // Mark the smart heaters
        for (i=0; i<model.numberOfBodyForces; i++) {
            _smarterHeaters[i] = false;
        }
        bf_id = 0;
        i = 0;
        for (FEMBodyForce *bodyForce in model.bodyForces) {
            if ([listUtilities listGetLogical:model inArray:bodyForce.valuesList forVariable:@"smart heater control" info:&found] == YES) {
                _smarterHeaters[i] = true;
                bf_id = i;
            }
            i++;
        }
        
        // Find the BC that controls the heater
        // If not found assume that smart heater is related to phase change
        bodyForceAtID = (model.bodyForces)[bf_id];
        meltPoint = [listUtilities listGetConstReal:model inArray:bodyForceAtID.valuesList forVariable:@"smart heater temperature" info:&gotMeltPoint minValue:NULL maxValue:NULL];
        
        smartHeaterAverage = NO;
        minDist = 0.0;
        smartHeaterNode = -1;
        smartHeaterNode = [listUtilities listGetInteger:model inArray:bodyForceAtID.valuesList forVariable:@"smart heater control node" info:&found minValue:NULL maxValue:NULL];
        if (found == NO) {
            found = [listUtilities listGetConstRealArray:model inArray:bodyForceAtID.valuesList forVariable:@"smart heater control point" buffer:&realWork];
            if (found == YES) {
                for (i=0; i<3; i++) {
                    controlPoint[i] = realWork.matrix[i][0];
                }
                minDist = HUGE_VAL;
                for (l=0; l<model.numberOfNodes; l++) {
                    if (tempPerm[l] < 0) continue;
                    
                    jx = meshNodes->x[l];
                    jy = meshNodes->y[l];
                    jz = meshNodes->z[l];
                    
                    dist = pow((controlPoint[0]-jx), 2.0) + pow((controlPoint[1]-jy), 2.0) + pow((controlPoint[2]-jz), 2.0);
                    if (dist < minDist) {
                        minDist = dist;
                        smartHeaterNode = l;
                    }
                }
            }
            NSLog(@"FEMHeatSolution: found control point at distance: %f\n", sqrt(minDist));
            NSLog(@"FEMHeatSolution: control point index: %d\n", smartHeaterNode);
        }
        
        if (gotMeltPoint == NO || smartHeaterNode < 0) {
            gotIt = NO;
            found = NO;
            smartHeaterBC = -1;
            
            i = 0;
            for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                gotIt = [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"smart heater boundary" info:&found];
                if (gotIt) {
                    smartHeaterBC = i;
                    break;
                }
                i++;
            }
            if (gotIt == NO) {
                i = 0;
                for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
                    gotIt = [listUtilities listGetLogical:model inArray:boundaryCondition.valuesList forVariable:@"phase change" info:&found];
                    if (gotIt == YES) {
                        smartHeaterBC = i;
                        break;
                    }
                    i++;
                }
            }
            if (smartHeaterBC < 0) {
                errorfunct("FEMHeatSolution", "Smart heater boundary / Phase change is undefined");
            }
            
            boundaryConditionAtID = (model.boundaryConditions)[smartHeaterBC];
            meltPoint = [listUtilities listGetConstReal:model inArray:boundaryConditionAtID.valuesList forVariable:@"smart heater temperature" info:&found minValue:NULL maxValue:NULL];
            if (found == NO) {
                for (FEMMaterial *material in model.materials) {
                    meltPoint = [listUtilities listGetConstReal:model inArray:material.valuesList forVariable:@"melting point" info:&found minValue:NULL maxValue:NULL];
                    if (found == YES) break;
                }
                if (found == NO) errorfunct("FEMHeatSolution", "Smart heater temperature / melting point is undefined.");
            }
            
            // Find the node related to temperature control
            if ((solution.solutionInfo)[@"smart heater average"] != nil) {
                smartHeaterAverage = [(solution.solutionInfo)[@"smart heater average"] boolValue];
            }
            if (smartHeaterAverage == NO) {
                jx = -HUGE_VALL;
                for (k=mesh.numberOfBulkElements; k<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; k++) {
                    if (elements[k].BoundaryInfo->Constraint == smartHeaterBC+1) {
                        for (l=0; l<elements[k].Type.NumberOfNodes; l++) {
                            if (meshNodes->x[elements[k].NodeIndexes[l]] >= jx) {
                                j = elements[k].NodeIndexes[l];
                                jx = meshNodes->x[elements[k].NodeIndexes[l]];
                            }
                        }
                    }
                }
                smartHeaterNode = j;
            }
        }
        
        if ((solution.solutionInfo)[@"smart heater control after tolerance"] != nil) {
            smartTol = [(solution.solutionInfo)[@"smart heater control after tolerance"] doubleValue];
        } else {
            smartTolReached = YES;
            smartTol = 1.0;
        }
    
        powerTimeScale = 0.0;
        if ((solution.solutionInfo)[@"smart heater time scale"] != nil) {
            powerTimeScale = [(solution.solutionInfo)[@"smart heater time scale"] doubleValue];
        }
    
        if (transient == YES && dt < powerTimeScale) {
            transientHeaterControl = YES;
            NSLog(@"FEMHeatSolution: using transient heater control.\n");
        } else {
            transientHeaterControl = NO;
            NSLog(@"FEMHeatSolution: using steady-state heater control.\n");
        }
    
        if (solution.doneTime != _doneTime) {
            _prevPowerScaling = _powerScaling;
            _doneTime = solution.doneTime;
        }
    }
    
    if (integralHeaterControl == YES) {
        NSLog(@"FEMHeatSolution: using integral heater control");
        for (i=0; i<model.numberOfBodyForces; i++) {
            _integralHeaters[i] = false;
        }
        i = 0;
        for (FEMBodyForce *bodyForce in model.bodyForces) {
            if ([listUtilities listCheckPresentVariable:@"integral heat source" inArray:bodyForce.valuesList] == YES) {
                _integralHeaters[i] = true;
            }
            i++;
        }
    }
    
    if ((solution.solutionInfo)[@"constant bulk system"] != nil) {
        constantBulk = [(solution.solutionInfo)[@"constant bulk system"] boolValue];
    }
    saveBulk = (constantBulk == YES || [(solution.solutionInfo)[@"save bulk system"] boolValue] == YES) ? YES : NO;
    saveBulk = (constantBulk == YES || [(solution.solutionInfo)[@"calculate loads"] boolValue] == YES) ? YES : NO;
    
    saveRelax = relax;
    cumulativeTime = 0.0;
    
    firstTime = YES;
    prevSolution = doublevec(0, localNodes-1);
    
    elementUtils = [[FEMElementUtils alloc] init];
    
    while (cumulativeTime < timeStep-1.0e-12 || transient == NO) {
        // The first time around this has been done by the caller...
        if (transient == YES && firstTime == NO) [kernel initializeTimeStepInSolution:solution model:model];
        firstTime = NO;
        
        // Save current solution
        for (i=0; i<localNodes; i++) {
            prevSolution[i] = temperature[i];
        }
        if (transient == YES) {
            if (tSolution == NULL) {
                tSolution = doublevec(0, localNodes-1);
                tSolution1 = doublevec(0, localNodes-1);
                for (i=0; i<tempContainers->size1PrevValues; i++) {
                    tSolution[i] = tempContainers->PrevValues[i][0];
                }
            }
        }
        
        norm = solution.variable.norm;
        
        for (iter=1; iter<=nonLinearIter; iter++) {
            NSLog(@"FEMHeatSolution\n");
            NSLog(@"FEMHeatSolution\n");
            NSLog(@"FEMHeatSolution: --------------------------------------------------------\n");
            NSLog(@"FEMHeatSolution: TEMPERATURE ITERATION %d\n", iter);
            NSLog(@"FEMHeatSolution: --------------------------------------------------------\n");
            NSLog(@"FEMHeatSolution\n");
            NSLog(@"FEMHeatSolution: Starting Assembly...\n");
            
            if (constantBulk == YES && matContainers->BulkValues != NULL) {
                memcpy(matContainers->Values, matContainers->BulkValues, matContainers->sizeBulkValues*sizeof(double));
                memcpy(matContainers->RHS, matContainers->BulkRHS, matContainers->sizeBulkRHS*sizeof(double));
                goto jump;
            }
            
            [kernel defaultInitializeSolution:solution model:model];
            
            if (smartHeaterControl == YES || integralHeaterControl == YES) {
                if (smartHeaterControl == YES) memset( _forceHeater, 0.0, tempContainers->sizeValues*sizeof(double) );
                memset( _heaterArea, 0.0, model.numberOfBodyForces*sizeof(double) );
                memset( _heaterSource, 0.0, model.numberOfBodyForces*sizeof(double) );
                memset( _heaterScaling, 0.0, model.numberOfBodyForces*sizeof(double) );
                memset( _heaterTarget, 0.0, model.numberOfBodyForces*sizeof(double) );
                heaterControlLocal = NO;
                
                for (t=0; t<solution.numberOfActiveElements; t++) {
                    element = [kernel getActiveElement:t solution:solution model:model];
                    bf_id = [kernel getBodyForceIDForElement:element model:model];
                    bodyForceAtID = (model.bodyForces)[bf_id-1];
                    if (bodyForceAtID == nil) continue;
                    if (!_smarterHeaters[bf_id-1] || !_integralHeaters[bf_id-1]) continue;
                    
                    mat_id = [kernel getMaterialIDForElement:element model:model];
                    materialAtID = (model.materials)[mat_id-1];
                    
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&buffer];
                    memcpy(_density, buffer.vector, buffer.m*sizeof(double));
                    
                    found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"heat source" buffer:&buffer];
                    memcpy(_load, buffer.vector, buffer.m*sizeof(double));
                    
                    s = [elementUtils elementArea:element numberOfNodes:n mesh:solution.mesh nodel:model];
                    
                    if (model.coordinates == axis_symmetric || model.coordinates == cylindric_symmetric) s = 2.0 * pi * s;
                    
                    sum = 0.0;
                    for (i=0; i<n; i++) {
                        sum = sum + (_density[i] * _load[i]);
                    }
                    _heaterSource[bf_id-1] = _heaterSource[bf_id-1] + s * sum / n;
                    _heaterArea[bf_id-1] = _heaterArea[bf_id-1] + s;
                    sum = 0.0;
                    for (i=0; i<n; i++) {
                        sum = sum + _density[i];
                    }
                    _heaterDensity[bf_id-1] = _heaterDensity[bf_id-1] + s * sum / n;
                }
                i = 0;
                for (FEMBodyForce *bodyForce in model.bodyForces) {
                    if (_integralHeaters[i] || _smarterHeaters[i]) {
                        _heaterDensity[i] = _heaterDensity[i] / _heaterDensity[i];
                    }
                    if (_integralHeaters[i]) {
                        _heaterTarget[i] = [listUtilities listGetConstReal:model inArray:bodyForce.valuesList forVariable:@"integral heat source" info:&found minValue:NULL maxValue:NULL];
                        _heaterScaling[i] = _heaterTarget[i] / _heaterSource[i];
                    }
                    i++;
                }
            }
            
            body_id = -1;
            materialAtID = nil;
            // Bulk elements
            for (t=0; t<solution.numberOfActiveElements; t++) {
                
                // Check if this element belongs to a body where temperature
                // should be calculated
                element = [kernel getActiveElement:t solution:solution model:model];
                if (element->BodyID != body_id) {
                    eq_id = [kernel getEquationIDForElement:element model:model];
                    equationAtID = (model.equations)[eq_id-1];
                    convectionFlag = [listUtilities listGetString:model inArray:equationAtID.valuesList forVariable:@"convection" info:&found];
                    
                    mat_id = [kernel getMaterialIDForElement:element model:model];
                    materialAtID = (model.materials)[mat_id-1];
                    
                    compressibilityFlag = [listUtilities listGetString:model inArray:equationAtID.valuesList forVariable:@"compressibility model" info:&found];
                    if (found == NO) compressibilityModel = imcrompressible;
                    if ([compressibilityFlag isEqualToString:@"incompressible"] == YES) {
                        compressibilityModel = imcrompressible;
                    } else if ([compressibilityFlag isEqualToString:@"user defined"] == YES) {
                        compressibilityModel = user_defined1;
                    } else if ([compressibilityFlag isEqualToString:@"perfect gas"] == YES || [compressibilityFlag isEqualToString:@"perfect gas equation 1"] == YES) {
                        compressibilityModel = perfect_gas1;
                    } else if ([compressibilityFlag isEqualToString:@"thermal"] == YES) {
                        compressibilityModel = thermal;
                    } else {
                        compressibilityModel = imcrompressible;
                    }
                    
                    phaseModel = [listUtilities listGetString:model inArray:equationAtID.valuesList forVariable:@"phase change model" info:&found];
                    if (found == NO) phaseModel = [listUtilities listGetString:model inArray:materialAtID.valuesList forVariable:@"phase change model" info:&found];
                    phaseChange = (found == YES && [phaseModel isEqualToString:@"none"] == NO) ? YES : NO;
                    if (phaseChange == YES) {
                        checkLatentHeatRelease = [listUtilities listGetLogical:model inArray:equationAtID.valuesList forVariable:@"check latent heat release" info:&found];
                    }
                }
                
                n = [kernel getNumberOfNodesForElement:element];
                [kernel getNodes:solution model:model inElement:element resultNodes:_elementNodes numberOfNodes:NULL];
                [kernel getScalarLocalField:_localTemperature sizeField:solution.mesh.maxElementDofs name:nil element:element solution:solution model:model timeStep:NULL];
                
                // Get element material parameters
                found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"heat capacity" buffer:&buffer];
                 memcpy(_heatCapacity, buffer.vector, buffer.m*sizeof(double));
                
                found = [listUtilities listGetRealArray:model inArray:materialAtID.valuesList forVariable:@"heat conductivity" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer];
                memset( **_heatConductivity, 0.0, (3*3*solution.mesh.maxElementDofs)*sizeof(double) );
            }
            
            
            
            
        jump:
            NSLog(@"");
            
        }
        if (buffer.vector != NULL) free_dvector(buffer.vector, 0, buffer.m-1);
        if (buffer.tensor != NULL) free_d3tensor(buffer.tensor, 0, buffer.m-1, 0, buffer.n-1, 0, buffer.p-1);
        
        
        
    }
    
    
}

@end
