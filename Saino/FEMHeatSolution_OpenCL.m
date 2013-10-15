//
//  FEMHeatSolution_OpenCL.m
//  Saino
//
//  Created by Seddik hakime on 26/09/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMHeatSolution_OpenCL.h"
#import "FEMKernel.h"
#import "FEMBoundaryCondition.h"
#import "FEMListUtilities.h"
#import "FEMUtilities.h"
#import "FEMElementUtils.h"
#import "FEMBodyForce.h"
#import "FEMMaterial.h"
#import "FEMEquation.h"
#import "FEMRadiation.h"
#import "FEMNumericIntegration.h"
#import "FEMCoordinateSystems.h"
#import "FEMElementDescription.h"
#import "FEMDiffuseConvectiveAnisotropic.h"
#import "FEMDiffuseConvectiveGeneralAnisotropic.h"
#import "FEMTimeIntegration.h"
#import "GaussIntegration.h"
#import "Utils.h"
#import "TimeProfile.h"
#import "heatSolutiomAssemblyKernel.cl.h"

static int k1 = 0, n1 = 0;
static double **stiff = NULL, **mass = NULL, **x = NULL;

enum {
    PHASE_SPATIAL_1 = 1,
    PHASE_SPATIAL_2,
    PHASE_TEMPORAL
};

@interface FEMHeatSolution_OpenCL ()

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
    int *_indexes;
    int *_saveIndexes;
    int *_tempPerm;
    double _dt;
    double _emissivity;
    double _powerScaling;
    double _prevPowerScaling;
    double _s;
    double _stefanBoltzmann;
    double _text;
    double _visibleFraction;
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
    double _normal[3];
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
    double *_prevSolution;
    double *_temperature;
    double *_tSolution;
    double *_tSolution1;
    bool *_smarterHeaters;
    bool *_integralHeaters;
    Nodes_t *_elementNodes;
    NSString *_phaseModel;
    NSString *_radiationFlag;
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
        
        _indexes = NULL;
        _saveIndexes = NULL;
        _tempPerm = NULL;
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
        _heaterArea = NULL;
        _heaterDensity = NULL;
        _heaterSource = NULL;
        _heaterScaling = NULL;
        _heaterTarget = NULL;
        _xx = NULL;
        _yy = NULL;
        _forceHeater = NULL;
        _prevSolution = NULL;
        _temperature = NULL;
        _tSolution = NULL;
        _tSolution1 = NULL;
        _smarterHeaters = NULL;
        _integralHeaters = NULL;
        _elementNodes = NULL;
        
        _tSolution = NULL;
        _tSolution1 = NULL;
        
        _phaseModel = nil;
        _radiationFlag = nil;
    }
    
    return self;
}

-(void)deallocation:(FEMSolution *)solution {
    
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

-(void)fieldSolutionComputer:(FEMSolution *)solution model:(FEMModel *)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
    int i, j, k, n, nb, t, bf_id, body_id, cols, compressibilityModel, eq_id, indx1, indx2, indx3, iter, kernelDof, mat_id, nBasis, nn, nnn, nonLinearIter, newtonIter, rows;
    int *elementPermutationStore, *indexes = NULL, indexStore[511];
    double at, at0, C1, dt0, cumulativeTime, detJ, newtonTol, nonLinearTol, norm,
    prevNorm, referencePressure, relax, relativeChange, saveRelax, specificHeatRatio, s, st, totat, totst, u, v, w;
    double *elementBasis, *elementdBasis, *forceVector, *nodesInfo;
    BOOL all, bubbles, found, heatFluxBC, firstTime, stabilize = YES, stat, useBubbles;
    NSString *stabilizeFlag, *convectionFlag, *compressibilityFlag;
    NSArray *bc;
    Element_t *elements = NULL, *element = NULL;
    Nodes_t *meshNodes = NULL;
    FEMMesh *mesh;
    FEMVariable *densitySol;
    FEMBodyForce *bodyForceAtID = nil;
    FEMMaterial *materialAtID = nil;
    FEMEquation *equationAtID = nil;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *tempContainers = NULL;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    GaussIntegrationPoints *IP;
    
    FEMKernel *kernel = [FEMKernel sharedKernel];
    FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
    
    mesh = (FEMMesh *)model.mesh;
    elements = mesh.getElements;
    meshNodes = mesh.getNodes;
    
    
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
            errorfunct("FEMHeatSolution:fieldSolutionComputer", "Memory allocation error");
        }
        
        cols = 2*n;
        rows = 2*n;
        
        _allocationDone = YES;
    }
    
    // Do some additional initialization and go for it
    _dt = timeStep;
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
    
    if ((solution.solutionInfo)[@"nonlinear system maximum iterations"] != nil) {
        nonLinearIter = [(solution.solutionInfo)[@"nonlinear system maximum iterations"] intValue];
    } else nonLinearIter = 1;
    
    if ((solution.solutionInfo)[@"nonlinear system convergence tolerance"] != nil) {
        nonLinearTol = [(solution.solutionInfo)[@"nonlinear system convergence tolerance"] doubleValue];
    }
    
    newtonTol = 1.0;
    newtonIter = 0;
    
    if (newtonIter == 0) _newtonLinearization = YES;
    
    if ((solution.solutionInfo)[@"nonlinear system relaxation factor"] != nil) {
        relax = [(solution.solutionInfo)[@"nonlinear system relaxation factor"] doubleValue];
    } else relax = 1.0;
    
    _transientAssembly = transient;
    found = NO;
    if ((solution.solutionInfo)[@"steady state transition time step"] != nil) {
        dt0 = [(solution.solutionInfo)[@"steady state transition time step"] doubleValue];
        found = YES;
    }
    if (found == YES && _dt > dt0) _transientAssembly = NO;
    
    saveRelax = relax;
    cumulativeTime = 0.0;
    
    firstTime = YES;
    _prevSolution = doublevec(0, _localNodes-1);
    
    FEMNumericIntegration *integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) errorfunct("FEMHeatSolution_OpenCL:fieldSolutionComputer", "Allocation error in FEMNumericIntegration!");
    
    NSString *stabilizationFlag = (solution.solutionInfo)[@"stabilization method"];
    BOOL vms = ([stabilizationFlag isEqualToString:@"vms"] == YES) ? YES : NO;
    
    // Connect to the OpenCL device
    dispatch_queue_t queue = [kernel getDispatchQueueAndInfoForDeviceType:@"cpu"];
    
    kernelDof = 26;
    
    // Set the matrix to be passed to the device
    void *matrixDiag = gcl_malloc(sizeof(cl_int)*matContainers->sizeDiag, matContainers->Diag, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    void *matrixRows = gcl_malloc(sizeof(cl_int)*matContainers->sizeRows, matContainers->Rows, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    void *matrixCols = gcl_malloc(sizeof(cl_int)*matContainers->sizeCols, matContainers->Cols, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
    void *matrixValues = gcl_malloc(sizeof(cl_double)*matContainers->sizeValues, matContainers->Values, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
    void *matrixRHS = gcl_malloc(sizeof(cl_double)*matContainers->sizeRHS, matContainers->RHS, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR);
    
    nodesInfo = doublevec(0, (kernelDof*mesh.numberOfNodes)-1);
    void *nodesInfoIn = gcl_malloc(sizeof(cl_double)*(kernelDof*mesh.numberOfNodes), NULL, CL_MEM_READ_WRITE);
    
    while (cumulativeTime < timeStep-1.0e-12 || transient == NO) {
        // The first time around this has been done by the caller...
        if (transient == YES && firstTime == NO) [kernel initializeTimeStepInSolution:solution model:model];
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
            at = cputime();
            at0 = realtime();
            
            NSLog(@"FEMHeatSolution:fieldSolutionComputer:\n");
            NSLog(@"FEMHeatSolution:fieldSolutionComputer:\n");
            NSLog(@"FEMHeatSolution:fieldSolutionComputer: -----------------------------------------------------------\n");
            NSLog(@"FEMHeatSolution:fieldSolutionComputer: TEMPERATURE ITERATION %d\n", iter);
            NSLog(@"FEMHeatSolution:fieldSolutionComputer: -----------------------------------------------------------\n");
            NSLog(@"FEMHeatSolution:fieldSolutionComputer:\n");
            NSLog(@"FEMHeatSolution:fieldSolutionComputer: Starting Assembly...\n");
            
            [kernel defaultInitializeSolution:solution model:model];
            
            body_id = -1;
            materialAtID = nil;
            bodyForceAtID = nil;
            nb = 0;
            
            for (NSMutableArray *color in mesh.colors) {

                memset( nodesInfo, 0.0, (kernelDof*mesh.numberOfNodes)*sizeof(double) );
                elementBasis = doublevec(0, ([color[0] intValue]*((mesh.maxElementNodes+1)*mesh.maxElementNodes))-1);
                elementdBasis = doublevec(0, ([color[0] intValue]*((mesh.maxElementNodes*3)*mesh.maxElementNodes))-1);
                elementPermutationStore = intvec(0, ([color[0] intValue]*mesh.maxElementDofs)-1);
                memset( elementBasis, 0.0, ([color[0] intValue]*((mesh.maxElementNodes+1)*mesh.maxElementNodes))*sizeof(double) );
                memset( elementdBasis, 0.0, ([color[0] intValue]*((mesh.maxElementNodes*3)*mesh.maxElementNodes))*sizeof(double) );
                memset( elementPermutationStore, -1, ([color[0] intValue]*mesh.maxElementDofs)*sizeof(int) );
                indx1 = 0;
                indx2 = 0;
                indx3 = 0;
                for (t=0; t<solution.numberOfActiveElements; t++) {
                    element = [kernel getActiveElement:t solution:solution model:model];
                    if (element->color.colorIndex == [color[1] intValue]) {
                        if (element->BodyID != body_id) {
                            eq_id = [kernel getEquationIDForElement:element model:model];
                            equationAtID = (model.equations)[eq_id-1];
                            convectionFlag = [listUtilities listGetString:model inArray:equationAtID.valuesList forVariable:@"convection" info:&found];
                            
                            mat_id = [kernel getMaterialIDForElement:element model:model];
                            materialAtID = (model.materials)[mat_id-1];
                            
                            compressibilityFlag = [listUtilities listGetString:model inArray:materialAtID.valuesList forVariable:@"compressibility model" info:&found];
                            if (found == NO) compressibilityModel = incompressible;
                            
                            if ([compressibilityFlag isEqualToString:@"incompressible"] == YES) {
                                compressibilityModel = incompressible;
                            } else if ([compressibilityFlag isEqualToString:@"user defined"] == YES) {
                                compressibilityModel = user_defined1;
                            } else if ([compressibilityFlag isEqualToString:@"perfect gas"] == YES || [compressibilityFlag isEqualToString:@"perfect gas equation 1"] == YES) {
                                compressibilityModel = perfect_gas1;
                            } else if ([compressibilityFlag isEqualToString:@"thermal"] == YES) {
                                compressibilityModel = thermal;
                            } else {
                                compressibilityModel = incompressible;
                            }
                        }
                        
                        n = element->Type.NumberOfNodes;
                        nn = n;
                        [kernel getNodes:solution model:model inElement:element resultNodes:_elementNodes numberOfNodes:NULL];
                        [kernel getScalarLocalField:_localTemperature sizeField:solution.mesh.maxElementDofs name:nil element:element solution:solution model:model timeStep:NULL];
                        
                        // Get element material parameters
                        memset( _heatCapacity, 0.0, n*sizeof(double) );
                        found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"heat capacity" buffer:&buffer];
                        if (found == YES) memcpy(_heatCapacity, buffer.vector, n*sizeof(double));
                        
                        memset( **_heatConductivity, 0.0, (3*3*solution.mesh.maxElementDofs)*sizeof(double) );
                        found = [listUtilities listGetRealArray:model inArray:materialAtID.valuesList forVariable:@"heat conductivity" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer];
                        if (found == YES) {
                            if (buffer.m == 1) {
                                for (i=0; i<3; i++) {
                                    for (j=0; j<n; j++) {
                                        _heatConductivity[i][i][j] = buffer.tensor[0][0][j];
                                    }
                                }
                            } else if (buffer.n == 1) {
                                for (i=0; i<min(3, buffer.m); i++) {
                                    for (j=0; j<n; j++) {
                                        _heatConductivity[i][i][j] = buffer.tensor[i][0][j];
                                    }
                                }
                            } else {
                                for (i=0; i<min(3, buffer.m); i++) {
                                    for (j=0; j<min(3, buffer.n); j++) {
                                        for (k=0; k<n; k++) {
                                            _heatConductivity[i][j][k] = buffer.tensor[i][j][k];
                                        }
                                    }
                                }
                            }
                        }
                        
                        if (compressibilityModel == perfect_gas1) {
                            // Read specific heat ratio
                            specificHeatRatio = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"specific heat ratio" info:&found minValue:NULL maxValue:NULL];
                            if (found == NO) specificHeatRatio = 5.0/3.0;
                            
                            // For an ideal gas, \gamma, c_p and R are really a constant.
                            // GasConstant is an array only since HeatCapacity formally is
                            for (i=0; i<n; i++) {
                                _gasConstant[i] = (specificHeatRatio - 1.0) * _heatCapacity[i] /  specificHeatRatio;
                            }
                            
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"pressure coefficient" buffer:&buffer];
                            if (found == YES) {
                                memcpy(_pressureCoeff, buffer.vector, n*sizeof(double));
                            } else {
                                for (i=0; i<n; i++) {
                                    _pressureCoeff[i] = 1.0;
                                }
                            }
                        } else if (compressibilityModel == thermal) {
                            memset( _referenceTemperature, 0.0, n*sizeof(double) );
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"reference temperature" buffer:&buffer];
                            if (found == YES) memcpy(_referenceTemperature, buffer.vector, n*sizeof(double));
                            
                            memset( _heatExpansionCoeff, 0.0, n*sizeof(double) );
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"heat expansion coefficient" buffer:&buffer];
                            if (found == YES) memcpy(_heatExpansionCoeff, buffer.vector, n*sizeof(double));
                            
                            memset( _density, 0.0, n*sizeof(double) );
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&buffer];
                            if (found == YES) memcpy(_density, buffer.vector, n*sizeof(double));
                            for (i=0; i<n; i++) {
                                _density[i] = _density[i] * ( 1.0 - _heatExpansionCoeff[i] * (_localTemperature[i] - _referenceTemperature[i]) );
                            }
                            
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"pressure coefficient" buffer:&buffer];
                            if (found == YES) {
                                memcpy(_pressureCoeff, buffer.vector, n*sizeof(double));
                            } else {
                                for (i=0; i<n; i++) {
                                    _pressureCoeff[i] = _localTemperature[i] *
                                    _heatExpansionCoeff[i] / (1.0 - _heatExpansionCoeff[i] * (_localTemperature[i] - _referenceTemperature[i]));
                                }
                            }
                        } else if (compressibilityModel == user_defined1) {
                            if (densitySol != nil) {
                                [kernel getScalarLocalField:_density sizeField:solution.mesh.maxElementDofs name:@"density" element:element solution:solution model:model timeStep:NULL];
                            } else {
                                memset( _density, 0.0, n*sizeof(double) );
                                found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&buffer];
                                if (found == YES) memcpy(_density, buffer.vector, n*sizeof(double));
                            }
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"pressure coefficient" buffer:&buffer];
                            if (found == YES) {
                                memcpy(_pressureCoeff, buffer.vector, n*sizeof(double));
                            } else memset( _pressureCoeff, 0.0, n*sizeof(double) );
                        } else {
                            memset( _pressureCoeff, 0.0, n*sizeof(double) );
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"pressure coefficient" buffer:&buffer];
                            if (found == YES) memcpy(_pressureCoeff, buffer.vector, n*sizeof(double));
                            
                            memset( _density, 0.0, n*sizeof(double) );
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&buffer];
                            if (found == YES) memcpy(_density, buffer.vector, n*sizeof(double));
                        }
                        
                        // Take pressure deviation p_d as the dependent variable p = p_0 + p_d.
                        // For perfect gas, read p_0
                        if (compressibilityModel != incompressible) {
                            referencePressure = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"reference pressure" info:&found minValue:NULL maxValue:NULL];
                            if (found == NO) referencePressure = 0.0;
                        }
                        
                        memset( _load, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                        memset( _pressure, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                        memset( _dPressureDt, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                        
                        // Check for convection model
                        C1 = 1.0;
                        memset( _u, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                        memset( _v, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                        memset( _w, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                        
                        memset( *_mu, 0.0, (3*solution.mesh.maxElementDofs)*sizeof(double) );
                        [kernel getVectorLocalField:_mu size1Field:3 size2Field:solution.mesh.maxElementDofs name:@"mesh velocity" element:element solution:solution model:model timeStep:NULL];
                        
                        if ([convectionFlag isEqualToString:@"constant"] == YES) {
                            
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"convection velocity 1" buffer:&buffer];
                            if (found == YES) {
                                memcpy(_u, buffer.vector, n*sizeof(double));
                            } else {
                                found = [kernel getReal:model forElement:element inArray:equationAtID.valuesList variableName:@"convection velocity 1" buffer:&buffer];
                                if (found == YES) memcpy(_u, buffer.vector, n*sizeof(double));
                            }
                            
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"convection velocity 2" buffer:&buffer];
                            if (found == YES) {
                                memcpy(_v, buffer.vector, n*sizeof(double));
                            } else {
                                [kernel getReal:model forElement:element inArray:equationAtID.valuesList variableName:@"convection velocity 2" buffer:&buffer];
                                if (found == YES) memcpy(_v, buffer.vector, n*sizeof(double));
                            }
                            
                            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"convection velocity 3" buffer:&buffer];
                            if (found == YES) {
                                memcpy(_w, buffer.vector, n*sizeof(double));
                            } else {
                                found = [kernel getReal:model forElement:element inArray:equationAtID.valuesList variableName:@"convection velocity 3" buffer:&buffer];
                                if (found == YES) memcpy(_w, buffer.vector, n*sizeof(double));
                            }
                            //                        } else if ([convectionFlag isEqualToString:@"computed"] == YES && flowSol != nil) {
                            //                            for (i=0; i<n; i++) {
                            //                                k = flowPerm[element->NodeIndexes[i]];
                            //                                if (k >= 0) {
                            //                                    _pressure[i] = flowSolution[nsdofs*(k+1)-1] + referencePressure;
                            //                                    switch (compressibilityModel) {
                            //                                        case perfect_gas1:
                            //                                            _density[i] = _pressure[i] / (_gasConstant[i] * _localTemperature[i]);
                            //                                            break;
                            //                                    }
                            //                                    if (transient == YES) {
                            //                                        _dPressureDt[i] = ( flowSolution[nsdofs*(k+1)-1] - flowSolContainers->PrevValues[nsdofs*(k+1)-1][0]) / timeStep;
                            //                                    }
                            //                                    switch (nsdofs) {
                            //                                        case 3:
                            //                                            _u[i] = flowSolution[nsdofs*k];
                            //                                            _v[i] = flowSolution[nsdofs*k+1];
                            //                                            _w[i] = 0.0;
                            //                                            break;
                            //                                        case 4:
                            //                                            _u[i] = flowSolution[nsdofs*k];
                            //                                            _v[i] = flowSolution[nsdofs*k+1];
                            //                                            _w[i] = flowSolution[nsdofs*k+2];
                            //                                            break;
                            //                                    }
                            //                                } else {
                            //                                    _u[i] = 0.0;
                            //                                    _v[i] = 0.0;
                            //                                    _w[i] = 0.0;
                            //                                }
                            //                            }
                        } else if ([convectionFlag isEqualToString:@"computed"] == YES) {
                            NSLog(@"FEMHeatSolution:fieldSolutionComputer: convection model specified but no accociated flow field present?\n");
                        } else {
                            all = YES;
                            for (i=0; i<3; i++) {
                                for (j=0; j<solution.mesh.maxElementDofs; j++) {
                                    if (_mu[i][j] != 0.0) {
                                        all = NO;
                                        break;
                                    }
                                }
                                if (all == NO) break;
                            }
                            if (all == YES) C1 = 0.0;
                        }
                        
                        // Check if modeling phase change with Eulerian approach
                        _phaseSpatial = NO;
                        for (i=0; i<n; i++) {
                            _heatCapacity[i] = _density[i] * _heatCapacity[i];
                        }
                        
                        memset( _viscosity, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                        
                        // Add body forces if any
                        bf_id = [kernel getBodyForceIDForElement:element model:model];
                        bodyForceAtID = (model.bodyForces)[bf_id-1];
                        if (bodyForceAtID != nil) {
                            // Frictional viscous heating
                            if ([listUtilities listGetLogical:model inArray:bodyForceAtID.valuesList forVariable:@"friction heat" info:&found] == YES) {
                                found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"viscosity" buffer:&buffer];
                                if (found == YES) memcpy(_viscosity, buffer.vector, n*sizeof(double));
                            }
                        }
                        
                        // Get heat source
                        found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"heat source" buffer:&buffer];
                        if (found == YES) {
                            for (i=0; i<n; i++) {
                                _load[i] = _density[i] * buffer.vector[i];
                            }
                        }
                        
                        memset( _c0, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                        
                        // Note at this point HeatCapacity = \rho * c_p or \rho * (c_p - R)
                        // and C1 = 0 (diffusion) or 1 (convection)
                        
                        // Perfusion (added as suggested by Matthias Zenker)
                        memset( _perfusionRate, 0.0, n*sizeof(double) );
                        found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"perfusion rate" buffer:&buffer];
                        if (found == YES) {
                            memcpy(_perfusionRate, buffer.vector, n*sizeof(double));
                            
                            memset( _perfusionRefTemperature, 0.0, n*sizeof(double) );
                            found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"perfusion reference temperature" buffer:&buffer];
                            if (found == YES) memcpy(_perfusionRefTemperature, buffer.vector, n*sizeof(double));
                            
                            memset( _perfusionDensity, 0.0, n*sizeof(double) );
                            found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"perfusion density" buffer:&buffer];
                            if (found == YES) memcpy(_perfusionDensity, buffer.vector, n*sizeof(double));
                            
                            memset( _perfusionHeatCapacity, 0.0, n*sizeof(double) );
                            found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"perfusion heat capacity" buffer:&buffer];
                            if (found == YES) memcpy(_perfusionHeatCapacity, buffer.vector, n*sizeof(double));
                            for (i=0; i<n; i++) {
                                _c0[i] = _perfusionHeatCapacity[i] * _perfusionRate[i] * _perfusionDensity[i];
                                _load[i] = _load[i] + _c0[i] * _perfusionRefTemperature[i];
                            }
                        }
                        
                        BOOL any = NO;
                        for (i=0; i<n; i++) {
                            if (C1 * _heatCapacity[i] != 0.0) {
                                any = YES;
                                break;
                            }
                        }
                        BOOL convection = (any == YES) ? YES : NO;
                        nBasis = n;
                        bubbles = NO;
                        if (convection == YES && !(vms == YES || stabilize == YES) && useBubbles == YES) {
                            nBasis = 2*n;
                            bubbles = YES;
                        }
                        
                        // Integration stuff
                        if (bubbles == YES) {
                            IP = GaussQuadrature(element, NULL, &element->Type.GaussPoints2);
                        } else IP = GaussQuadrature(element, NULL, NULL);
                        
                        // Copy the element wise node info to global node that each kernel can read
                        // Data are serialized as follows:
                        // [_heatCapacity, _heatConductivity[][], _pressureCoeff, _density, _localTemperature, _load, _pressure, _dPressureDt, _enthalpy, _u, _v, _w, _viscosity, heatCapacity, mux, muy, muz, _c0, phaseSpatial, next node...]
                        for (i=0; i<n; i++) {
                            j = tempContainers->Perm[element->NodeIndexes[i]];
                            nodesInfo[kernelDof*j] = _heatCapacity[i];
                            nodesInfo[kernelDof*j+1] = _heatConductivity[0][0][i];
                            nodesInfo[kernelDof*j+2] = _heatConductivity[0][1][i];
                            nodesInfo[kernelDof*j+3] = _heatConductivity[0][2][i];
                            nodesInfo[kernelDof*j+4] = _heatConductivity[1][0][i];
                            nodesInfo[kernelDof*j+5] = _heatConductivity[1][1][i];
                            nodesInfo[kernelDof*j+6] = _heatConductivity[1][2][i];
                            nodesInfo[kernelDof*j+7] = _heatConductivity[2][0][i];
                            nodesInfo[kernelDof*j+8] = _heatConductivity[2][1][i];
                            nodesInfo[kernelDof*j+9] = _heatConductivity[2][2][i];
                            nodesInfo[kernelDof*j+10] = _pressureCoeff[i];
                            nodesInfo[kernelDof*j+11] = _density[i];
                            nodesInfo[kernelDof*j+12] = _localTemperature[i];
                            nodesInfo[kernelDof*j+13] = _load[i];
                            nodesInfo[kernelDof*j+14] = _pressure[i];
                            nodesInfo[kernelDof*j+15] = _dPressureDt[i];
                            nodesInfo[kernelDof*j+16] = _enthalpy[i];
                            nodesInfo[kernelDof*j+17] = _u[i];
                            nodesInfo[kernelDof*j+18] = _v[i];
                            nodesInfo[kernelDof*j+19] = _w[i];
                            nodesInfo[kernelDof*j+20] = _viscosity[i];
                            nodesInfo[kernelDof*j+21] =  C1 * _heatCapacity[i];
                            nodesInfo[kernelDof*j+22] = _mu[0][i];
                            nodesInfo[kernelDof*j+23] = _mu[1][i];
                            nodesInfo[kernelDof*j+24] = _mu[2][i];
                            nodesInfo[kernelDof*j+25] = _c0[i];
                        }
                        
                        nnn = IP->n;
                        for (int tt=0; tt<IP->n; tt++) {
                            u = IP->u[tt];
                            v = IP->v[tt];
                            w = IP->w[tt];
                            
                            // Basis function values & derivatives at the integration point
                            stat = [integration setBasisForElement:element elementNodes:_elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:bubbles basisDegree:NULL];
                            stat = [integration setBasisFirstDerivativeForElement:element elementNodes:_elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:bubbles basisDegree:NULL];
                            stat = [integration setMetricDeterminantForElement:element elementNodes:_elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
                            detJ = integration.metricDeterminant;
                            s = detJ * IP->s[tt];
                            for (i=0; i<n; i++) {
                                elementBasis[indx1] = integration.basis[i];
                                indx1++;
                            }
                            elementBasis[indx1] = s;
                            indx1++;
                            
                            for (i=0; i<n; i++) {
                                for (j=0; j<3; j++) {
                                    elementdBasis[indx2] = integration.basisFirstDerivative[i][j];
                                    indx2++;
                                }
                            }
                        }
                        
                        memset( indexStore, -1, sizeof(indexStore) );
                        n = [kernel getElementDofsSolution:solution model:model forElement:element atIndexes:indexStore];
                        for (i=0; i<n; i++) {
                            elementPermutationStore[indx3] = tempContainers->Perm[indexStore[i]];
                            indx3++;
                        }
                        GaussQuadratureDeallocation(IP);
                    }
                }
                
                // Set data to device
                void *elementPermutationStoreIn = gcl_malloc(sizeof(cl_int)*([color[0] intValue]*mesh.maxElementDofs), elementPermutationStore, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
                void *elementBasisIn = gcl_malloc(sizeof(cl_double)*([color[0] intValue]*((mesh.maxElementNodes+1)*mesh.maxElementNodes)), elementBasis, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
                void *elementdBasisIn = gcl_malloc(sizeof(cl_double)*([color[0] intValue]*((mesh.maxElementNodes*3)*mesh.maxElementNodes)), elementdBasis, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR);
                
                // Dispatch the kernel block using the queue created earlier
                dispatch_sync(queue, ^{
                    
                    void *map = gcl_map_ptr(nodesInfoIn, CL_MEM_READ_WRITE, sizeof(cl_double)*(kernelDof*mesh.numberOfNodes));
                    memcpy(map, nodesInfo, sizeof(cl_double)*(kernelDof*mesh.numberOfNodes));
                    
                    cl_ndrange range = {
                        1,
                        {0, 0, 0},
                        {[color[0] intValue], 0, 0},
                        {0, 0, 0}
                    };
                    
                    // Call the kernel
                    heatSolutionAssembly_kernel(&range, (cl_int*)matrixDiag, (cl_int*)matrixRows, (cl_int*)matrixCols, (cl_double*)matrixValues, (cl_double*)matrixRHS, (cl_double*)nodesInfoIn, (cl_int*)elementPermutationStoreIn,
                                                (cl_double*)elementBasisIn, (cl_double*)elementdBasisIn, model.dimension, nn, n, nnn, nBasis, kernelDof, solution.variable.dofs);
                });
                
                gcl_free(elementPermutationStoreIn);
                gcl_free(elementBasisIn);
                gcl_free(elementdBasisIn);
                
                free_dvector(elementBasis, 0, ([color[0] intValue]*((mesh.maxElementNodes+1)*mesh.maxElementNodes))-1);
                free_dvector(elementdBasis, 0, ([color[0] intValue]*((mesh.maxElementNodes*3)*mesh.maxElementNodes))-1);
                free_ivector(elementPermutationStore, 0, ([color[0] intValue]*mesh.maxElementDofs)-1);
            }
            
             dispatch_sync(queue, ^{
                 // Copy data from device to host
                 gcl_memcpy(matContainers->Values, matrixValues, sizeof(cl_double)*matContainers->sizeValues);
                 gcl_memcpy(matContainers->RHS, matrixRHS, sizeof(cl_double)*matContainers->sizeRHS);
             });
            
            at = cputime() -  at;

            // Neumann & Newton boundary conditions
            for (t=0; t<solution.mesh.numberOfBoundaryElements; t++) {
                element = [kernel getBoundaryElement:solution atIndex:t];
                if ([kernel isActiveBoundaryElement:element inSolution:solution model:model] == NO) continue;
                
                n = element->Type.NumberOfNodes;
                if ([kernel getElementFamily:element] == 1) continue;
                
                bc = [kernel getBoundaryCondition:model forElement:element];
                if (bc == nil) continue;
                
                heatFluxBC = [listUtilities listGetLogical:model inArray:bc forVariable:@"heat fluc bc" info:&found];
                if (found == YES && heatFluxBC == NO) continue;
            } // Neumann & Newton BCs
            
            [kernel defaultFinishAssemblySolution:solution model:model];
            NSLog(@"FEMHeatSolution:fieldSolutionComputer: Assembly done\n");
            
            [kernel dirichletBoundaryConditions:model inSolution:solution usingOffset:NULL offDiaginalMatrix:NULL];
            
            // Solve the system and check for convergence
            st = cputime();
            
            prevNorm = norm;
            norm = [kernel findSolution:solution model:model backRorateNT:NULL];
            
            st = cputime() - st;
            totat = totat + at;
            totst = totst + st;
            NSLog(@"FEMHeatSolution:fieldSolutionComputer: iter: %d, Assembly (s): %f %f\n", iter, at, totat);
            NSLog(@"FEMHeatSolution:fieldSolutionComputer: iter: %d, Solve (s): %f %f\n", iter, st, totst);
            
            relativeChange = solution.variable.nonLinChange;
            NSLog(@"FEMHeatSolution:fieldSolutionComputer: result norm: %e\n", norm);
            NSLog(@"FEMHeatSolution:fieldSolutionComputer: relative change: %e\n", relativeChange);
            
            if (relativeChange < newtonTol || iter >= newtonIter) _newtonLinearization = YES;
            if (relativeChange < nonLinearTol) break;
            
        } // Non linear iteration loop
        
        // Compute cumulative time done by now and time remaining
        if (transient == NO) break;
        cumulativeTime = cumulativeTime + _dt;
        _dt = timeStep - cumulativeTime;
        
        if (buffer.vector != NULL) {
            free_dvector(buffer.vector, 0, buffer.m-1);
            buffer.vector = NULL;
        }
        if (buffer.tensor != NULL) {
            free_d3tensor(buffer.tensor, 0, buffer.m-1, 0, buffer.n-1, 0, buffer.p-1);
            buffer.tensor = NULL;
        }
        if (indexes != NULL) {
            free_ivector(indexes, 0, nb-1);
            indexes = NULL;
        }
    } // time interval
    
    gcl_free(matrixDiag);
    gcl_free(matrixRows);
    gcl_free(matrixCols);
    gcl_free(matrixValues);
    gcl_free(matrixRHS);
    gcl_free(nodesInfoIn);
    dispatch_release(queue);
    
    free_dvector(nodesInfo, 0, (kernelDof*mesh.numberOfNodes)-1);
    [integration deallocation:mesh];
    
    solution.dt = timeStep;
    [solution.solutionInfo setValue:@(saveRelax) forKey:@"nonlinear system relaxation factor"];
    free_dvector(_prevSolution, 0, _localNodes-1);
    
    if ([(solution.solutionInfo)[@"adaptive mesh refinement"] boolValue] == YES) {
        // TODO: implement mesh refinement
    }
}

@end

