//===----------------------------------------------------------------------===//
//  FEMHeatSolution_OpenCL.m
//  Saino
//
//  Created by Seddik hakime on 26/09/13.
//  Copyright (c) 2013 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
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

#import "FEMHeatSolution_OpenCL.h"
#import "FEMCore.h"
#import "FEMListUtilities.h"
#import "FEMElementUtils.h"
#import "FEMRadiation.h"
#import "FEMNumericIntegration.h"
#import "FEMCoordinateSystems.h"
#import "FEMElementDescription.h"
#import "FEMDiffuseConvectiveAnisotropic.h"
#import "FEMDiffuseConvectiveGeneralAnisotropic.h"
#import "GaussIntegration.h"
#import "Utils.h"
#import "TimeProfile.h"
#import "OpenCLUtils.h"

static int k1 = 0, n1 = 0;
static double * _Nullable * _Nullable stiff = NULL, * _Nullable * _Nullable mass = NULL, * _Nullable * _Nullable x = NULL;

enum {
    PHASE_SPATIAL_1 = 1,
    PHASE_SPATIAL_2,
    PHASE_TEMPORAL
};

@interface FEMHeatSolution_OpenCL ()
-(void)FEMHeatSolution_nullify;
@end

@implementation FEMHeatSolution_OpenCL {
    
    BOOL _allocationDone;
    BOOL _constantBulk;
    BOOL _heatGapBC;
    BOOL _newtonLinearization;
    BOOL _phaseSpatial;
    BOOL _transientAssembly;
    int _doneTime;
    int _phaseChangeModel;
    int _localNodes;
    int * _Nullable _indexes;
    int * _Nullable _saveIndexes;
    int * _Nullable _tempPerm;
    double _dt;
    double _emissivity;
    double _powerScaling;
    double _prevPowerScaling;
    double _s;
    double _stefanBoltzmann;
    double _text;
    double _visibleFraction;
    double * _Nullable _u, * _Nullable _v, * _Nullable _w;
    double * _Nullable * _Nullable _mu;
    double * _Nullable _pressure;
    double * _Nullable _dPressureDt;
    double * _Nullable _pressureCoeff;
    double * _Nullable _density;
    double * _Nullable _work;
    double * _Nullable _latentHeat;
    double * _Nullable * _Nullable _phaseVelocity;
    double * _Nullable _electricConductivity;
    double _normal[3];
    double * _Nullable _permeability;
    double * _Nullable _viscosity;
    double * _Nullable _c0;
    double * _Nullable _heatTransferCoeff;
    double * _Nullable _heatExpansionCoeff;
    double * _Nullable _referenceTemperature;
    double * _Nullable * _Nullable _mass;
    double * _Nullable _localTemperature;
    double * _Nullable _heatCapacity;
    double * _Nullable _enthalpy;
    double * _Nullable _nodalEmissivity;
    double * _Nullable _gasConstant;
    double * _Nullable _aText;
    double * _Nullable * _Nullable * _Nullable _heatConductivity;
    double * _Nullable * _Nullable _stiff;
    double * _Nullable _load;
    double * _Nullable _force;
    double * _Nullable _timeForce;
    double * _Nullable _perfusionRate;
    double * _Nullable _perfusionDensity;
    double * _Nullable _perfusionHeatCapacity;
    double * _Nullable _perfusionRefTemperature;
    double * _Nullable _heaterArea;
    double * _Nullable _heaterDensity;
    double * _Nullable _heaterSource;
    double * _Nullable _heaterScaling;
    double * _Nullable _heaterTarget;
    double * _Nullable _xx;
    double * _Nullable _yy;
    double * _Nullable _forceHeater;
    double * _Nullable _prevSolution;
    double * _Nullable _temperature;
    double * _Nullable _tSolution;
    double * _Nullable _tSolution1;
    bool * _Nullable _smarterHeaters;
    bool * _Nullable _integralHeaters;
    Nodes_t * _Nullable _elementNodes;
    NSString * _Nullable _phaseModel;
    NSString * _Nullable _radiationFlag;
}

-(void)FEMHeatSolution_nullify {
    
    _indexes = NULL;
    _saveIndexes = NULL;
    _u = NULL;
    _v = NULL;
    _w = NULL;
    _mu = NULL;
    _pressure = NULL;
    _dPressureDt = NULL;
    _pressureCoeff = NULL;
    _density = NULL;
    _work = NULL;
    _latentHeat = NULL;
    _phaseVelocity = NULL;
    _electricConductivity = NULL;
    _permeability = NULL;
    _viscosity = NULL;
    _c0 = NULL;
    _heatTransferCoeff = NULL;
    _heatExpansionCoeff = NULL;
    _referenceTemperature = NULL;
    _mass = NULL;
    _localTemperature = NULL;
    _heatCapacity = NULL;
    _enthalpy = NULL;
    _nodalEmissivity = NULL;
    _gasConstant = NULL;
    _aText = NULL;
    _heatConductivity = NULL;
    _stiff = NULL;
    _load = NULL;
    _force = NULL;
    _timeForce = NULL;
    _perfusionRate = NULL;
    _perfusionDensity = NULL;
    _perfusionHeatCapacity = NULL;
    _perfusionRefTemperature = NULL;
    if (_elementNodes != NULL) {
        _elementNodes->x = NULL;
        _elementNodes->y = NULL;
        _elementNodes->z = NULL;
    }
    _elementNodes = NULL;
}

- (id)init
{
    self = [super init];
    if (self) {
        _allocationDone = NO;
        _constantBulk = NO;
        _newtonLinearization = NO;
        _doneTime = 0;
        _powerScaling = 1.0;
        _prevPowerScaling = 1.0;
        _tempPerm = NULL;
        _temperature = NULL;
        
        [self FEMHeatSolution_nullify];
        
        _heaterArea = NULL;
        _heaterDensity = NULL;
        _heaterSource = NULL;
        _heaterScaling = NULL;
        _heaterTarget = NULL;
        _smarterHeaters = NULL;
        _integralHeaters = NULL;
        _xx = NULL;
        _yy = NULL;
        _forceHeater = NULL;
        
        _prevSolution = NULL;
        _tSolution = NULL;
        _tSolution1 = NULL;
        
        _phaseModel = nil;
        _radiationFlag = nil;
    }
    
    return self;
}

-(void)deallocation:(FEMSolution * _Nonnull)solution {
    
    int n = solution.mesh.maxElementDofs;
    variableArraysContainer *tempContainers = solution.variable.getContainers;
    
    if (stiff != NULL) {
        free_dmatrix(stiff, 0, 0, 0, n1-1);
    }
    if (mass != NULL) {
        free_dmatrix(mass, 0, 0, 0, n1-1);
    }
    if (x != NULL) {
        free_dmatrix(x, 0, n1-1, 0, k1-1);
    }
    
    if (_indexes != NULL) free_ivector(_indexes, 0, n-1);
    if (_saveIndexes != NULL) free_ivector(_saveIndexes, 0, n-1);
    if (_u != NULL) free_dvector(_u, 0, n-1);
    if (_v != NULL) free_dvector(_v, 0, n-1);
    if (_w != NULL) free_dvector(_w, 0, n-1);
    if (_mu != NULL) free_dmatrix(_mu, 0, 2, 0, n-1);
    if (_pressure != NULL) free_dvector(_pressure, 0, n-1);
    if (_dPressureDt != NULL) free_dvector(_dPressureDt, 0, n-1);
    if (_pressureCoeff != NULL) free_dvector(_pressureCoeff, 0, n-1);
    if (_density != NULL) free_dvector(_density, 0, n-1);
    if (_work != NULL) free_dvector(_work, 0, n-1);
    if (_latentHeat != NULL) free_dvector(_latentHeat, 0, n-1);
    if (_phaseVelocity != NULL) free_dmatrix(_phaseVelocity, 0, 2, 0, n-1);
    if (_electricConductivity != NULL) free_dvector(_electricConductivity, 0, n-1);
    if (_permeability != NULL) free_dvector(_permeability, 0, n-1);
    if (_viscosity != NULL) free_dvector(_viscosity, 0, n-1);
    if (_c0 != NULL) free_dvector(_c0, 0, n-1);
    if (_heatTransferCoeff != NULL) free_dvector(_heatTransferCoeff, 0, n-1);
    if (_heatExpansionCoeff != NULL) free_dvector(_heatExpansionCoeff, 0, n-1);
    if (_referenceTemperature != NULL) free_dvector(_referenceTemperature, 0, n-1);
    if (_mass != NULL) free_dmatrix(_mass, 0, (2*n)-1, 0, (2*n)-1);
    if (_localTemperature != NULL) free_dvector(_localTemperature, 0, n-1);
    if (_heatCapacity != NULL) free_dvector(_heatCapacity, 0, n-1);
    if (_enthalpy != NULL) free_dvector(_enthalpy, 0, n-1);
    if (_nodalEmissivity != NULL) free_dvector(_nodalEmissivity, 0, n-1);
    if (_gasConstant != NULL) free_dvector(_gasConstant, 0, n-1);
    if (_aText != NULL) free_dvector(_aText, 0, n-1);
    if (_heatConductivity != NULL) free_d3tensor(_heatConductivity, 0, 2, 0, 2, 0, n-1);
    if (_stiff != NULL) free_dmatrix(_stiff, 0, (2*n)-1, 0, (2*n)-1);
    if (_load != NULL) free_dvector(_load, 0, n-1);
    if (_force != NULL) free_dvector(_force, 0, (2*n)-1);
    if (_timeForce != NULL) free_dvector(_timeForce, 0, (2*n)-1);
    if (_perfusionRate != NULL) free_dvector(_perfusionRate, 0, n-1);
    if (_perfusionDensity != NULL) free_dvector(_perfusionDensity, 0, n-1);
    if (_perfusionHeatCapacity != NULL) free_dvector(_perfusionHeatCapacity, 0, n-1);
    if (_perfusionRefTemperature != NULL) free_dvector(_perfusionRefTemperature, 0, n-1);
    if (_elementNodes->x != NULL) free_dvector(_elementNodes->x, 0, n-1);
    if (_elementNodes->y != NULL) free_dvector(_elementNodes->y, 0, n-1);
    if (_elementNodes->z != NULL) free_dvector(_elementNodes->z, 0, n-1);
    if (_elementNodes != NULL) free(_elementNodes);
    
    if (_heaterArea != NULL) free_dvector(_heaterArea, 0, n-1);
    if (_heaterDensity != NULL) free_dvector(_heaterDensity, 0, n-1);
    if (_heaterSource != NULL) free_dvector(_heaterSource, 0, n-1);
    if (_heaterScaling != NULL) free_dvector(_heaterScaling, 0, n-1);
    if (_heaterTarget != NULL) free_dvector(_heaterTarget, 0, n-1);
    if (_smarterHeaters != NULL) free_bvector(_smarterHeaters, 0, n-1);
    if (_integralHeaters != NULL) free_bvector(_integralHeaters, 0, n-1);
    
    if (_xx != NULL) free_dvector(_xx, 0, tempContainers->sizeValues-1);
    if (_yy != NULL) free_dvector(_yy, 0, tempContainers->sizeValues-1);
    if (_forceHeater != NULL) free_dvector(_forceHeater, 0, tempContainers->sizeValues-1);
    
    if (_tSolution != NULL) free_dvector(_tSolution, 0, _localNodes-1);
    if (_tSolution1 != NULL) free_dvector(_tSolution1, 0, _localNodes-1);
}

-(void)solutionComputer:(FEMSolution * _Nonnull)solution model:(FEMModel * _Nonnull)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
    int i, n, nb, t, body_id, cols, iter, nBasis, nonLinearIter, newtonIter, position, returnValue, rows;
    static int firstTimeCL = 0;
    int *colorMapping = NULL, *elementNodeIndexesStore = NULL;
    char *program_source;
    double at, at0, ct, mt, mmt, dt0, cumulativeTime, newtonTol, nonLinearTol, norm, prevNorm, relax, relativeChange, saveRelax, st, totat, totst;
    double *forceVector;
    BOOL found, heatFluxBC, firstTime, stabilize = YES, useBubbles;
    NSString *stabilizeFlag;
    NSArray *bc;
    cl_context         context;
	cl_command_queue   cmd_queue;
	cl_device_id       devices;
	cl_int             err;
	size_t src_len;
    Element_t *elements = NULL, *element = NULL;
    Nodes_t *meshNodes = NULL;
    FEMMesh *mesh;
    FEMBodyForce *bodyForceAtID = nil;
    FEMMaterial *materialAtID = nil;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *tempContainers = NULL;
    FEMCore *core = [FEMCore sharedCore];
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    FEMTimeIntegration *timeIntegration;
    
    if (transient == YES) {
        timeIntegration = [[FEMTimeIntegration alloc] init];
    }
    
    mesh = (FEMMesh *)model.mesh;
    elements = mesh.getElements;
    meshNodes = mesh.getNodes;
    colorMapping = mesh.getColorMapping;
    elementNodeIndexesStore = mesh.getElementNodeIndexesStore;

    if (solution.matrix == nil) return;
    
    // Get variables needed for solution
    
    matContainers = solution.matrix.getContainers;
    forceVector = matContainers->RHS;
    
    tempContainers = solution.variable.getContainers;
    _tempPerm = tempContainers->Perm;
    _temperature = tempContainers->Values;
    
    _localNodes = 0;
    for (i=0; i<tempContainers->sizePerm; i++) {
        if (_tempPerm[i] >= 0) _localNodes++;
    }
    if (_localNodes <= 0) return;
    
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
            free_dvector(_aText, 0, n-1);
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
            [self FEMHeatSolution_nullify];
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
        _aText = doublevec(0, n-1);
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
        initNodes(_elementNodes);
        _elementNodes->x = doublevec(0, n-1);
        _elementNodes->y = doublevec(0, n-1);
        _elementNodes->z = doublevec(0, n-1);
        if (_indexes == NULL || _saveIndexes == NULL || _u == NULL || _v == NULL || _w == NULL || _mu == NULL || _pressure == NULL ||
            _dPressureDt == NULL || _pressureCoeff == NULL || _density == NULL || _work == NULL || _latentHeat == NULL || _phaseVelocity == NULL ||
            _electricConductivity == NULL || _permeability == NULL || _viscosity == NULL || _c0 == NULL || _heatTransferCoeff == NULL ||
            _heatExpansionCoeff == NULL || _referenceTemperature == NULL || _mass == NULL || _localTemperature == NULL || _heatCapacity == NULL ||
            _enthalpy == NULL || _nodalEmissivity == NULL || _gasConstant == NULL || _aText == NULL || _heatConductivity == NULL || _stiff == NULL |
            _load == NULL || _force == NULL || _timeForce == NULL || _perfusionRate == NULL || _perfusionDensity == NULL ||
            _perfusionHeatCapacity == NULL || _perfusionRefTemperature == NULL || _elementNodes->x == NULL || _elementNodes->y == NULL ||
            _elementNodes->z == NULL) {
            fatal("FEMHeatSolution:solutionComputer", "Memory allocation error.");
        }
        
        cols = 2*n;
        rows = 2*n;
        
        _allocationDone = YES;
    }
    
    // Do some additional initialization and go for it
    _dt = timeStep;
    if (solution.solutionInfo[@"stabilize"] != nil) {
        stabilize = [solution.solutionInfo[@"stabilize"] boolValue];
    }
    
    if (solution.solutionInfo[@"bubbles"] != nil) {
        useBubbles = [solution.solutionInfo[@"bubbles"] boolValue];
    } else useBubbles = YES;
    
    if (solution.solutionInfo[@"stabilization method"] != nil) {
        stabilizeFlag = solution.solutionInfo[@"stabilization method"];
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
    
    if (solution.solutionInfo[@"nonlinear system maximum iterations"] != nil) {
        nonLinearIter = [solution.solutionInfo[@"nonlinear system maximum iterations"] intValue];
    } else nonLinearIter = 1;
    
    nonLinearTol = 0.0;
    if (solution.solutionInfo[@"nonlinear system convergence tolerance"] != nil) {
        nonLinearTol = [solution.solutionInfo[@"nonlinear system convergence tolerance"] doubleValue];
    }
    
    newtonTol = 1.0;
    newtonIter = 0;
    
    if (newtonIter == 0) _newtonLinearization = YES;
    
    if (solution.solutionInfo[@"nonlinear system relaxation factor"] != nil) {
        relax = [solution.solutionInfo[@"nonlinear system relaxation factor"] doubleValue];
    } else relax = 1.0;
    
    _transientAssembly = transient;
    found = NO;
    dt0 = 0.0;
    if (solution.solutionInfo[@"steady state transition time step"] != nil) {
        dt0 = [solution.solutionInfo[@"steady state transition time step"] doubleValue];
        found = YES;
    }
    if (found == YES && _dt > dt0) _transientAssembly = NO;
    
    saveRelax = relax;
    cumulativeTime = 0.0;
    
    firstTime = YES;
    _prevSolution = doublevec(0, _localNodes-1);
    
    FEMNumericIntegration *integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) fatal("FEMHeatSolution_OpenCL:solutionComputer", "Allocation error in FEMNumericIntegration.");
    
	// Connect to a compute devise
	err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, 1, &devices, NULL);

    size_t returned_size = 0;
	cl_char vendor_name[1024] = {0};
	cl_char device_name[1024] = {0};
	err = clGetDeviceInfo(devices, CL_DEVICE_VENDOR, sizeof(vendor_name), vendor_name, &returned_size);
	err = clGetDeviceInfo(devices, CL_DEVICE_NAME, sizeof(device_name), device_name, &returned_size);
    
    if (firstTimeCL == 0) {
        fprintf(stdout, "FEMHeatSolution_OpenCL:solutionComputer: Connecting to %s %s...\n\n", vendor_name, device_name);
		device_stats(devices);
	}
    
    returnValue = LoadFileIntoString("/Users/seddikhakime/Documents/Saino/Saino/heatSolutiomAssemblyKernel.cl", &program_source, &src_len);
	if (returnValue) {
        fprintf(stderr, "FEMHeatSolution_OpenCL:solutionComputer: Error: Can't load kernel source.\n");
        exit(-1);
	}
    
    // Create the context of the command queue
	context = clCreateContext(0, 1, &devices, NULL, NULL, &err);
	cmd_queue = clCreateCommandQueue(context, devices, 0, NULL);
	
	// Allocate memory for program and kernels
    // Create the program .cl file
	cl_program program = clCreateProgramWithSource(context, 1, (const char**)&program_source, NULL, &err);
	if (err) {
        fprintf(stderr, "FEMHeatSolution_OpenCL:solutionComputer: Can't create program. Error was: %d.\n", err);
		exit(-1);
	}
    
    // Build the program (compile it)
	err = clBuildProgram(program, 0, NULL, NULL, NULL, &err);
	char build[2048];
	clGetProgramBuildInfo(program, devices, CL_PROGRAM_BUILD_LOG, 2048, build, NULL);
	if (err) {
        fprintf(stderr, "FEMHeatSolution_OpenCL:solutionComputer: Can't build program. Error was: %d.\n", err);
        fprintf(stderr, "FEMHeatSolution_OpenCL:solutionComputer: Build Log:\n%s.\n", build);
		exit(-1);
	}

    // Create the kernel
	cl_kernel CLkernel = clCreateKernel(program, "heatSolutionAssembly", &err);
	if (err) {
        fprintf(stderr, "FEMHeatSolution_OpenCL:solutionComputer:heatAssembly: Can't create kernel. Error was: %d.\n", err);
		exit(-1);
	}
    
    size_t thread_size;
	clGetKernelWorkGroupInfo(CLkernel, devices, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &thread_size, NULL);
    
    if (firstTimeCL == 0) {
        fprintf(stdout, "FEMHeatSolution_OpenCL:solutionComputer:heatAssembly: Recommended Work Group Size: %lu.\n", thread_size);
	}
    
    // Allocate memory and queue it to be written to the device
	cl_mem matDiag = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_int)*matContainers->sizeDiag, NULL, NULL);
	cl_mem matRows = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_int)*matContainers->sizeRows, NULL, NULL);
	cl_mem matCols = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_int)*matContainers->sizeCols, NULL, NULL);
	cl_mem matValues = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_double)*matContainers->sizeValues, NULL, NULL);
	cl_mem matRhs = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_double)*matContainers->sizeRHS, NULL, NULL);
    cl_mem colorMappingIn = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_int)*mesh.numberOfBulkElements, NULL, NULL);
    cl_mem elementNodeIndexesStoreIn = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_int)*(mesh.numberOfBulkElements*mesh.maxElementDofs), NULL, NULL);
    cl_mem permutationIn = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_int)*(tempContainers->sizePerm), NULL, NULL);
    cl_mem nodesX = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_double)*(mesh.numberOfNodes), NULL, NULL);
    cl_mem nodesY = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_double)*(mesh.numberOfNodes), NULL, NULL);
    cl_mem nodesZ = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_double)*(mesh.numberOfNodes), NULL, NULL);

    firstTimeCL = 1;
    
    double allocationSize = sizeof(cl_int)*matContainers->sizeDiag + sizeof(cl_int)*matContainers->sizeRows + sizeof(cl_int)*matContainers->sizeCols + sizeof(cl_double)*matContainers->sizeValues +
                          + sizeof(cl_double)*matContainers->sizeRHS
                          + sizeof(cl_int)*mesh.numberOfBulkElements
                          + sizeof(cl_int)*(mesh.numberOfBulkElements*mesh.maxElementDofs)
                          + sizeof(cl_int)*(tempContainers->sizePerm)
                          + sizeof(cl_double)*(mesh.numberOfNodes) + sizeof(cl_double)*(mesh.numberOfNodes) + sizeof(cl_double)*(mesh.numberOfNodes);
    
    fprintf(stdout, "FEMHeatSolution:solutionComputer: Allocation for the CRS matrix (MB): %f.\n",
          (sizeof(cl_double)*matContainers->sizeValues+sizeof(cl_double)*matContainers->sizeRHS+sizeof(cl_int)*matContainers->sizeCols+sizeof(cl_int)*matContainers->sizeRows+sizeof(cl_int)*matContainers->sizeDiag)/(1024.0*1024.0));
    fprintf(stdout, "FEMHeatSolution:solutionComputer: Total allocation to the device (MB): %f.\n", allocationSize/(1024.0*1024.0));

    while (cumulativeTime < timeStep-1.0e-12 || transient == NO) {
        // The first time around this has been done by the caller...
        if (transient == YES && firstTime == NO) [core initializeTimeStepInSolution:solution model:model];
        firstTime = NO;
        
        // Save current solution
        memcpy(_prevSolution, _temperature, _localNodes*sizeof(double));
        if (transient == YES) {
            if (_tSolution == NULL) {
                _tSolution = doublevec(0, _localNodes-1);
                _tSolution1 = doublevec(0, _localNodes-1);
                for (i=0; i<tempContainers->size1PrevValues; i++) {
                    _tSolution[i] = tempContainers->PrevValues[i][0];
                }
            }
        }
        
        totat = 0.0;
        totst = 0.0;
        
        norm = solution.variable.norm;
        
        for (iter=1; iter<=nonLinearIter; iter++) {
            at0 = realtime();
            
            fprintf(stdout, "FEMHeatSolution:solutionComputer:\n");
            fprintf(stdout, "FEMHeatSolution:solutionComputer:\n");
            fprintf(stdout, "FEMHeatSolution:solutionComputer: -----------------------------------------------------------\n");
            fprintf(stdout, "FEMHeatSolution:solutionComputer: TEMPERATURE ITERATION %d.\n", iter);
            fprintf(stdout, "FEMHeatSolution:solutionComputer: -----------------------------------------------------------\n");
            fprintf(stdout, "FEMHeatSolution:solutionComputer:\n");
            fprintf(stdout, "FEMHeatSolution:solutionComputer: Starting Assembly...\n");
            
            [core defaultInitializeSolution:solution model:model];
            
            body_id = -1;
            materialAtID = nil;
            bodyForceAtID = nil;
            nb = 0;
            
            at = cputime();
            
            mt = cputime();
            
            // Queue memory to be written to the device
            err = clEnqueueWriteBuffer(cmd_queue, matDiag, CL_TRUE, 0, sizeof(cl_int)*matContainers->sizeDiag, (void*)matContainers->Diag, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, matRows, CL_TRUE, 0, sizeof(cl_int)*matContainers->sizeRows, (void*)matContainers->Rows, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, matCols, CL_TRUE, 0, sizeof(cl_int)*matContainers->sizeCols, (void*)matContainers->Cols, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, matValues, CL_TRUE, 0, sizeof(cl_double)*matContainers->sizeValues, (void*)matContainers->Values, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, matRhs, CL_TRUE, 0, sizeof(cl_double)*matContainers->sizeRHS, (void*)matContainers->RHS, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, colorMappingIn, CL_TRUE, 0, sizeof(cl_int)*mesh.numberOfBulkElements, (void*)colorMapping, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, elementNodeIndexesStoreIn, CL_TRUE, 0, sizeof(cl_int)*(mesh.numberOfBulkElements*mesh.maxElementDofs), (void*)elementNodeIndexesStore, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, permutationIn, CL_TRUE, 0, sizeof(cl_int)*(tempContainers->sizePerm), (void*)tempContainers->Perm, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, nodesX, CL_TRUE, 0, sizeof(cl_double)*(mesh.numberOfNodes), (void*)meshNodes->x, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, nodesY, CL_TRUE, 0, sizeof(cl_double)*(mesh.numberOfNodes), (void*)meshNodes->y, 0, NULL, NULL);
            err = clEnqueueWriteBuffer(cmd_queue, nodesZ, CL_TRUE, 0, sizeof(cl_double)*(mesh.numberOfNodes), (void*)meshNodes->z, 0, NULL, NULL);
            
            // Push the data out to the device
            clFinish(cmd_queue);
            
            mt = cputime() -  mt;
            
            n = 3;
            nBasis = 3;
            int dimension = model.dimension;
            int varDofs = solution.variable.dofs;
            
            // Set kernel arguments
            err  = clSetKernelArg(CLkernel, 0, sizeof(cl_mem), &matValues);
            err |= clSetKernelArg(CLkernel, 1, sizeof(cl_mem), &matRhs);
            err |= clSetKernelArg(CLkernel, 2, sizeof(cl_mem), &matDiag);
            err |= clSetKernelArg(CLkernel, 3, sizeof(cl_mem), &matRows);
            err |= clSetKernelArg(CLkernel, 4, sizeof(cl_mem), &matCols);
            err |= clSetKernelArg(CLkernel, 5, sizeof(cl_mem), &colorMappingIn);
            err |= clSetKernelArg(CLkernel, 6, sizeof(cl_mem), &elementNodeIndexesStoreIn);
            err |= clSetKernelArg(CLkernel, 7, sizeof(cl_mem), &permutationIn);
            err |= clSetKernelArg(CLkernel, 8, sizeof(cl_mem), &nodesX);
            err |= clSetKernelArg(CLkernel, 9, sizeof(cl_mem), &nodesY);
            err |= clSetKernelArg(CLkernel, 10, sizeof(cl_mem), &nodesZ);
            err |= clSetKernelArg(CLkernel, 12, sizeof(int), &dimension);
            err |= clSetKernelArg(CLkernel, 13, sizeof(int), &n);
            err |= clSetKernelArg(CLkernel, 14, sizeof(int), &n);
            err |= clSetKernelArg(CLkernel, 15, sizeof(int), &nBasis);
            err |= clSetKernelArg(CLkernel, 16, sizeof(int), &varDofs);
            
            ct = cputime();
                
            for (NSMutableArray *color in mesh.colors) {
                
                position = 0;
                for (i=0; i<[color[1] intValue]; i++) {
                    NSMutableArray *cc = mesh.colors[i];
                    position = position + [cc[0] intValue];
                }
                
                size_t global_work_size = [color[0] intValue];
                //size_t local_work_size = 64;
                err |= clSetKernelArg(CLkernel, 11, sizeof(int), &position);
                
                //Queue up the kernels
                err = clEnqueueNDRangeKernel(cmd_queue, CLkernel, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
                
                // Finish the calculation
                clFinish(cmd_queue);
            }
            
            ct = cputime() -  ct;
            mmt = cputime();
            
            // Read results data from the device
            err = clEnqueueReadBuffer(cmd_queue, matValues, CL_TRUE, 0, sizeof(cl_double)*matContainers->sizeValues, matContainers->Values, 0, NULL, NULL);
            err = clEnqueueReadBuffer(cmd_queue, matRhs, CL_TRUE, 0, sizeof(cl_double)*matContainers->sizeRHS, matContainers->RHS, 0, NULL, NULL);
            clFinish(cmd_queue);

            mmt = cputime() -  mmt;
            
            at = cputime() -  at;

            // Neumann & Newton boundary conditions
            for (t=0; t<solution.mesh.numberOfBoundaryElements; t++) {
                element = [core getBoundaryElement:solution atIndex:t];
                if ([core isActiveBoundaryElement:element inSolution:solution model:model] == NO) continue;
                
                n = element->Type.NumberOfNodes;
                if ([core getElementFamily:element] == 1) continue;
                
                bc = [core getBoundaryCondition:model forElement:element];
                if (bc == nil) continue;
                
                heatFluxBC = [listUtilities listGetLogical:model inArray:bc forVariable:@"heat fluc bc" info:&found];
                if (found == YES && heatFluxBC == NO) continue;
            } // Neumann & Newton BCs
            
            [core defaultFinishAssemblySolution:solution model:model timeIntegration:timeIntegration utilities:utilities];
            fprintf(stdout, "FEMHeatSolution:solutionComputer: Assembly done.\n");
            
            [core dirichletBoundaryConditions:model inSolution:solution usingOffset:NULL offDiaginalMatrix:NULL];
            
            // Solve the system and check for convergence
            st = cputime();
            
            prevNorm = norm;
            norm = [core findSolution:solution model:model backRorateNT:NULL];
            
            st = cputime() - st;
            totat = totat + at;
            totst = totst + st;
            fprintf(stdout, "FEMHeatSolution:solutionComputer: iter: %d, Assembly (compute, mem up, mem down, all, tot) (s): %f %f %f %f %f.\n", iter, ct, mt, mmt, at, totat);
            fprintf(stdout, "FEMHeatSolution:solutionComputer: iter: %d, Solve (s): %f %f.\n", iter, st, totst);
            
            relativeChange = solution.variable.nonLinChange;
            fprintf(stdout, "FEMHeatSolution:solutionComputer: result norm: %e.\n", norm);
            fprintf(stdout, "FEMHeatSolution:solutionComputer: relative change: %e.\n", relativeChange);
            
            if (relativeChange < newtonTol || iter >= newtonIter) _newtonLinearization = YES;
            if (relativeChange < nonLinearTol) break;
            
        } // Non linear iteration loop
        
        // Compute cumulative time done by now and time remaining
        if (transient == NO) break;
        cumulativeTime = cumulativeTime + _dt;
        _dt = timeStep - cumulativeTime;
    } // time interval
    
    // Release kernel, program and memory objects
	clReleaseKernel(CLkernel);
	clReleaseProgram(program);
	clReleaseCommandQueue(cmd_queue);
	clReleaseContext(context);
	
	clReleaseMemObject(matDiag);
	clReleaseMemObject(matRows);
	clReleaseMemObject(matCols);
	clReleaseMemObject(matValues);
	clReleaseMemObject(matRhs);
    clReleaseMemObject(colorMappingIn);
    clReleaseMemObject(elementNodeIndexesStoreIn);
    clReleaseMemObject(permutationIn);
    clReleaseMemObject(nodesX);
    clReleaseMemObject(nodesY);
    clReleaseMemObject(nodesZ);
    
    [integration deallocation:mesh];
    
    solution.dt = timeStep;
    [solution.solutionInfo setValue:@(saveRelax) forKey:@"nonlinear system relaxation factor"];
    free_dvector(_prevSolution, 0, _localNodes-1);
}

@end

