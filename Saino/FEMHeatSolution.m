//
//  FEMHeatEquation.m
//  Saino
//
//  Created by Seddik hakime on 19/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMHeatSolution.h"
#import "FEMCore.h"
#import "FEMListUtilities.h"
#import "FEMElementUtils.h"
#import "FEMRadiation.h"
#import "FEMCoordinateSystems.h"
#import "FEMElementDescription.h"
#import "FEMDiffuseConvectiveAnisotropic.h"
#import "FEMDiffuseConvectiveGeneralAnisotropic.h"
#import "GaussIntegration.h"
#import "Utils.h"
#import "TimeProfile.h"

static int k1 = 0, n1 = 0;
static double *saveValues = NULL;
static double **stiff = NULL, **mass = NULL, **x = NULL;

enum {
  PHASE_SPATIAL_1 = 1,
  PHASE_SPATIAL_2,
  PHASE_TEMPORAL
};

@interface FEMHeatSolution ()
-(void)FEMHeatSolution_nullify;
-(void)FEMHeatSolution_findGapIndexesElement:(Element_t * __nonnull)element indexes:(int * __nonnull)indexes numberOfNodes:(int)n solution:(FEMSolution * __nonnull)solution;
-(void)FEMHeatSolution_effectiveHeatCapacityElement:(Element_t * __nonnull)element numberOfNodes:(int)n material:(FEMMaterial * __nonnull)material model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities transientSimulation:(BOOL)transient;
-(void)FEMHeatSolution_integrationOverSurfaceElement:(Element_t * __nonnull)element boundaryNumberOfNodes:(int)n radiationBoundaryOfNodes:(int)m model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution mesh:(FEMMesh * __nonnull)mesh;
-(void)FEMHeatSolution_diffuseGrayRadiationModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution core:(FEMCore * __nonnull)core element:(Element_t * __nonnull)element numberOfNodes:(int)n forceVector:(double * __nonnull)forceVector angleFraction:(double * __nonnull)angleFraction text:(double * __nonnull)text timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities;
-(void)FEMHeatSolution_addHeatGapSolution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element numberOfNodes:(int)n core:(FEMCore * __nonnull)core;
-(void)FEMHeatSolution_addHeatFluxBC:(NSArray * __nonnull)bc element:(Element_t * __nonnull)element parent:(Element_t * __nonnull)parent numberOfNodes:(int)n forceVector:(double * __nonnull)forceVector core:(FEMCore * __nonnull)core solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilites crsMatrix:(FEMMatrixCRS * __nonnull)crsMatrix bandMatrix:(FEMMatrixBand * __nonnull)bandMatrix integration:(FEMNumericIntegration * __nonnull)integration diffuseConvectiveAnisotropic:(FEMDiffuseConvectiveAnisotropic * __nonnull)diffuseConvectiveAnisotropic diffuseConvectiveGeneralAnisotropic:(FEMDiffuseConvectiveGeneralAnisotropic * __nonnull)diffuseConvectiveGeneralAnisotropic timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities;
-(void)FEMHeatSolution_addGlobalTimeSolution:(FEMSolution * __nonnull)solution;
-(BOOL)FEMHeatSolution_checkLatentHeatModel:(FEMModel * __nonnull)model;
@end

@implementation FEMHeatSolution {
    
    BOOL _allocationDone;
    BOOL _constantBulk;
    BOOL _heatGapBC;
    BOOL _newtonLinearization;
    BOOL _phaseSpatial;
    BOOL _transientAssembly;
    int _cols;
    int _doneTime;
    int _phaseChangeModel;
    int _localNodes;
    int * __nullable _indexes;
    int _rows;
    int * __nullable _saveIndexes;
    int * __nullable _tempPerm;
    double _dt;
    double _emissivity;
    double _powerScaling;
    double _prevPowerScaling;
    double _s;
    double _stefanBoltzmann;
    double _text;
    double _visibleFraction;
    double *__nullable _u, * __nullable _v, * __nullable _w;
    double *__nullable * __nullable _mu;
    double * __nullable _pressure;
    double * __nullable _dPressureDt;
    double * __nullable _pressureCoeff;
    double * __nullable _density;
    double * __nullable _work;
    double *__nullable _latentHeat;
    double * __nullable * __nullable _phaseVelocity;
    double *__nullable _electricConductivity;
    double _normal[3];
    double * __nullable _permeability;
    double * __nullable _viscosity;
    double * __nullable _c0;
    double * __nullable _heatTransferCoeff;
    double * __nullable _heatExpansionCoeff;
    double * __nullable _referenceTemperature;
    double * __nullable * __nullable _mass;
    double * __nullable _localTemperature;
    double * __nullable _heatCapacity;
    double * __nullable _enthalpy;
    double * __nullable _nodalEmissivity;
    double * __nullable _gasConstant;
    double * __nullable _aText;
    double * __nullable * __nullable * __nullable _heatConductivity;
    double * __nullable *__nullable _stiff;
    double * __nullable _load;
    double * __nullable _force;
    double * __nullable _timeForce;
    double * __nullable _perfusionRate;
    double * __nullable _perfusionDensity;
    double * __nullable _perfusionHeatCapacity;
    double * __nullable _perfusionRefTemperature;
    double * __nullable _heaterArea;
    double * __nullable _heaterDensity;
    double * __nullable _heaterSource;
    double * __nullable _heaterScaling;
    double * __nullable _heaterTarget;
    double * __nullable _xx;
    double * __nullable _yy;
    double * __nullable _forceHeater;
    double * __nullable _prevSolution;
    double * __nullable _temperature;
    double * __nullable _tSolution;
    double * __nullable _tSolution1;
    bool * __nullable _smarterHeaters;
    bool * __nullable _integralHeaters;
    Nodes_t * __nullable _elementNodes;
    NSString * __nullable _phaseModel;
    NSString * __nullable _radiationFlag;
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

-(void)FEMHeatSolution_findGapIndexesElement:(Element_t * __nonnull)element indexes:(int * __nonnull)indexes numberOfNodes:(int)n solution:(FEMSolution * __nonnull)solution {
    
    int i, j, k;
    double x0, y0, z0, x, y, z;
    BOOL any;
    Element_t *left = NULL, *parent = NULL, *right = NULL;
    Nodes_t *meshNodes = NULL;
    
    left = element->BoundaryInfo->Left;
    right = element->BoundaryInfo->Right;
    
    if (left == NULL || right == NULL) return;
    
    meshNodes = solution.mesh.getNodes;
    for (i=0; i<n; i++) {
        parent = left;
        k = element->NodeIndexes[i];
        
        any = NO;
        for (j=0; j<parent->Type.NumberOfNodes; j++) {
            if (parent->NodeIndexes[j] == k) {
                any = YES;
                break;
            }
        }
        if (any == YES) parent = right;
        
        x0 = _elementNodes->x[i];
        y0 = _elementNodes->y[i];
        z0 = _elementNodes->z[i];
        for (j=0; j<parent->Type.NumberOfNodes; j++) {
            k = parent->NodeIndexes[j];
            x = meshNodes->x[k] - x0;
            y = meshNodes->y[k] - y0;
            z = meshNodes->z[k] - z0;
            if (pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0) < AEPS) break;
        }
        indexes[i] = k;
    }
}

-(void)FEMHeatSolution_effectiveHeatCapacityElement:(Element_t * __nonnull)element numberOfNodes:(int)n material:(FEMMaterial * __nonnull)material model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities transientSimulation:(BOOL)transient {
    
    int i;
    BOOL any, found, specific;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    //---------------------------------------------------------------------
    // See if temperature gradient inside the element is large
    // enough to use c_p = sqrt( (dH/dx)^2 / (dT/dx)^2 ), otherwise
    // use c_p = dH/dT or if in time dependent simulation, use
    // c_p = (dH/dt) / (dT/dt), if requested
    //---------------------------------------------------------------------
    
    if ([_phaseModel isEqualToString:@"spatial 1"] == YES) {
        _phaseChangeModel = PHASE_SPATIAL_1;
    } else if ([_phaseModel isEqualToString:@"spatial 2"] == YES) {
        // Check if local variation of temperature is large enough to actually use the
        // Spatial 2 model. Should perhaps be scaled to element size (or actually
        // compute the gradient, but this will do for now...)
        double maxVal, minVal;
        vDSP_maxvD(_localTemperature, 1, &maxVal, n);
        vDSP_minvD(_localTemperature, 1, &minVal, n);
        _s = maxVal - minVal;
        if (_s < AEPS) {
            _phaseChangeModel = PHASE_SPATIAL_1;
        } else {
            _phaseChangeModel = PHASE_SPATIAL_2;
        }
    } else if ([_phaseModel isEqualToString:@"temporal"] == YES) {
        // Note that here HeatCapacity is misused for saving dT
        if (transient == YES && _tSolution != NULL) {
            for (i=0; i<n; i++) {
                _heatCapacity[i] = _temperature[_tempPerm[element->NodeIndexes[i]]] - _tSolution[_tempPerm[element->NodeIndexes[i]]];
            }
            any = NO;
            for (i=0; i<n; i++) {
                if (fabs(_heatCapacity[i]) < AEPS) {
                    any = YES;
                    break;
                }
            }
            if (any == YES) {
                _phaseChangeModel = PHASE_SPATIAL_1;
            } else {
                _phaseChangeModel = PHASE_TEMPORAL;
            }
        } else {
            _phaseChangeModel = PHASE_SPATIAL_1;
        }
    } else {
        _phaseChangeModel = PHASE_SPATIAL_1;
    }
    
    _phaseSpatial = (_phaseChangeModel == PHASE_SPATIAL_2) ? YES : NO;
    
    specific = [listUtilities listCheckPresentVariable:@"specific enthalpy" inArray:material.valuesList];
    
    switch (_phaseChangeModel) {
        // This phase change model is avaiable only for some type of real entries
        // that have an implemented analytical derivation rule
        case PHASE_SPATIAL_1:
            found = [listUtilities listGetReal:model inArray:material.valuesList forVariable:@"effective heat capacity" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
            if (found == YES) {
                memcpy(_heatCapacity, buffer.vector, n*sizeof(double));
            } else {
                if (specific == YES) {
                    found = [listUtilities listGetDerivativeValue:model inArray:material.valuesList forVariable:@"specific enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer];
                    if (found == YES) memcpy(_heatCapacity, buffer.vector, n*sizeof(double));
                    for (i=0; i<n; i++) {
                        _heatCapacity[i] = _density[i] * _heatCapacity[i];
                    }
                } else {
                    found = [listUtilities listGetDerivativeValue:model inArray:material.valuesList forVariable:@"enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer];
                    if (found == YES)  memcpy(_heatCapacity, buffer.vector, n*sizeof(double));
                }
            }
            break;
            
        // Note that for the 'spatial 2' mode, the evaluation of c_p is done in each integration point
        // and thus Enthalpy and PhaseSpatial flag are used instead of HeatCapacity directly
        case PHASE_SPATIAL_2:
            if (specific == YES) {
                found = [listUtilities listGetReal:model inArray:material.valuesList forVariable:@"specific enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
                if (found == YES) memcpy(_enthalpy, buffer.vector, n*sizeof(double));
                for (i=0; i<n; i++) {
                    _enthalpy[i] = _density[i] * _enthalpy[i];
                }
            } else {
                found = [listUtilities listGetReal:model inArray:material.valuesList forVariable:@"enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
                if (found == YES)  memcpy(_enthalpy, buffer.vector, n*sizeof(double));
            }
            
        case PHASE_TEMPORAL:
            memcpy(_tSolution1, _temperature, _localNodes*sizeof(double));
            
            // When retrieving the value of enthalpy on the previous timeStep, the
            // relevant entries of the Temperature solution in the global vector are
            // tampered in order to make the ListGetReal method work as wanted
            if (specific == YES) {
                found = [listUtilities listGetReal:model inArray:material.valuesList forVariable:@"specific enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
                if (found == YES) memcpy(_work, buffer.vector, n*sizeof(double));
                memcpy(_temperature, _tSolution, _localNodes*sizeof(double));
                for (i=0; i<n; i++) {
                    _work[i] = _work[i] - buffer.vector[i];
                }
            } else {
                found = [listUtilities listGetReal:model inArray:material.valuesList forVariable:@"enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
                if (found == YES) memcpy(_work, buffer.vector, n*sizeof(double));
                memcpy(_temperature, _tSolution, _localNodes*sizeof(double));
                for (i=0; i<n; i++) {
                    _work[i] = _work[i] - buffer.vector[i];
                }
            }
            memcpy(_temperature, _tSolution1, _localNodes*sizeof(double));
            for (i=0; i<n; i++) {
                _heatCapacity[i] = _work[i] / _heatCapacity[i];
            }
            if (specific == YES) {
                for (i=0; i<n; i++) {
                    _heatCapacity[i] = _density[i] * _heatCapacity[i];
                }
            }
    }
    
    if (buffer.vector != NULL) free_dvector(buffer.vector, 0, buffer.m-1);
}

-(void)FEMHeatSolution_integrationOverSurfaceElement:(Element_t * __nonnull)element boundaryNumberOfNodes:(int)n radiationBoundaryOfNodes:(int)m model:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution mesh:(FEMMesh * __nonnull)mesh {
 
    int i, p, q, t;
    double alpha, detJ, force, s, x, y, z;
    BOOL stat;
    FEMNumericIntegration *integration;
    FEMCoordinateSystems *coordinateSystem;
    GaussIntegrationPoints *IP = NULL;
    
    memset( *_stiff, 0.0, ((2*solution.mesh.maxElementDofs)*(2*solution.mesh.maxElementDofs))*sizeof(double) );
    memset( _force, 0.0, (2*solution.mesh.maxElementDofs)*sizeof(double) );
    
    // Integration stuff
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) fatal("FEMHeatSolution:FEMHeatSolution_integrationOverSurfaceElement", "Allocation error in FEMNumericIntegration.");
    IP = GaussQuadrature(element, NULL, NULL);

    coordinateSystem = [[FEMCoordinateSystems alloc] init];
    for (t=0; t<IP->n; t++) {
        // Basis function values & derivatives at the integration point
        stat = [integration setBasisForElement:element elementNodes:_elementNodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t] withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:_elementNodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t]];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
        // Coordinate system dependent info
        if (model.coordinates != cartesian) {
            x = 0.0;
            y = 0.0;
            z = 0.0;
            for (i=0; i<n; i++) {
                x = x + (_elementNodes->x[i]*integration.basis[i]);
                y = y + (_elementNodes->y[i]*integration.basis[i]);
                z = z + (_elementNodes->z[i]*integration.basis[i]);
            }
            s = s * [coordinateSystem coordinateSquareRootMetricModel:model coordX:x coordY:y coordZ:z];
        }
        
        force = 0.0;
        alpha = 0.0;
        for (i=0; i<n; i++) {
            force = force + (_load[i] * integration.basis[i]);
            alpha = alpha + (_heatTransferCoeff[i] * integration.basis[i]);
        }
        
        for (p=0; p<n; p++) {
            for (q=0; q<m; q++) {
                _stiff[p][q] = _stiff[p][q] + s * alpha * integration.basis[p] / m;
            }
        }
        for (p=0; p<n; p++) {
            _force[p] = _force[p] + s * force * integration.basis[p];
        }
    }

    [integration deallocation:mesh];
}

-(void)FEMHeatSolution_diffuseGrayRadiationModel:(FEMModel * __nonnull)model solution:(FEMSolution * __nonnull)solution core:(FEMCore * __nonnull)core element:(Element_t * __nonnull)element numberOfNodes:(int)n forceVector:(double * __nonnull)forceVector angleFraction:(double * __nonnull)angleFraction text:(double * __nonnull)text timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities {
    
    int i, j, k, k1, k2, l, m, cols, implicitFactors, perm[element->Type.NumberOfNodes], rows;
    double area, asum, s, sum;
    FEMElementUtils *elementUtils;
    FEMRadiation *radiation;
    Element_t *radiationElement = NULL;
    
    elementUtils = [[FEMElementUtils alloc] init];
    radiation = [[FEMRadiation alloc] init];
    
    asum = 0.0;
    if (_newtonLinearization == NO) {
        *text = [radiation computeRadiationLoadModel:model mesh:solution.mesh element:element temperature:_temperature reorder:_tempPerm emissivity:_emissivity angleFraction:angleFraction];
    } else { // Full Newton-Raphson solution
       
        // Go through surfaces (j) this surface (i) is getting radiation from
        area = [elementUtils elementArea:element numberOfNodes:n mesh:solution.mesh nodel:model];
        
        radiationElement = solution.mesh.getElements;
        for (j=0; j<element->BoundaryInfo->GebhardtFactors->NumberOfFactors; j++) {
            *text = [radiation computeRadiationCoeffModel:model mesh:solution.mesh element:element index:j] / area;
            asum = asum + *text;
            
            // Gebhardt factors are given element-wise at the center of the element,
            // so take average of nodal temperatures (or integrate over surface j)
            k = radiationElement[element->BoundaryInfo->GebhardtFactors->Elements[j]].Type.NumberOfNodes;
            implicitFactors = element->BoundaryInfo->GebhardtFactors->NumberOfImplicitFactors;
            if (implicitFactors == 0) implicitFactors = element->BoundaryInfo->GebhardtFactors->NumberOfFactors;
            
            if (j < implicitFactors) {
                sum = 0.0;
                for (i=0; i<k; i++) {
                    sum = sum + pow(_temperature[_tempPerm[radiationElement[element->BoundaryInfo->GebhardtFactors->Elements[j]].NodeIndexes[i]]], 4.0);
                }
                s = pow((sum/k), (1.0/4.0));
                // Linearization of G_jiT^4_j term
                for (i=0; i<n; i++) {
                    _heatTransferCoeff[i] = -4.0 * *text * pow(s, 3.0) * _stefanBoltzmann;
                    _load[i] = -3.0 * *text * pow(s, 4.0) * _stefanBoltzmann;
                }
                // Integrate the contribution of surface j over surface i and
                // add to global matrix
                [self FEMHeatSolution_integrationOverSurfaceElement:element boundaryNumberOfNodes:n radiationBoundaryOfNodes:k model:model solution:solution mesh:solution.mesh];
                if (_transientAssembly == YES) {
                    memset( *_mass, 0.0, ((2*solution.mesh.maxElementDofs)*(2*solution.mesh.maxElementDofs))*sizeof(double) );
                    cols = 2*solution.mesh.maxElementDofs;
                    rows = 2*solution.mesh.maxElementDofs;
                    for (i=0; i<element->Type.NumberOfNodes; i++) {
                        perm[i] = _tempPerm[element->NodeIndexes[i]];
                    }
                    [core addFirstOrderTimeModel:model solution:solution element:element massMatrix:_mass stiffMatrix:_stiff force:_force dt:_dt size:n dofs:1 nodeIndexes:perm rows:&rows cols:&cols timeIntegration:timeIntegration utilities:utilities];
                }
                
                for (m=0; m<n; m++) {
                    k1 = _tempPerm[element->NodeIndexes[m]];
                    for (l=0; l<k; l++) {
                        k2 = _tempPerm[radiationElement[element->BoundaryInfo->GebhardtFactors->Elements[j]].NodeIndexes[l]];
                        [core addToMatrixElementForSolution:solution atIndex:k1 andIndex:k2 value:_stiff[m][l]];
                    }
                    forceVector[k1] = forceVector[k1] + _force[m];
                }
            } else {
                sum = 0.0;
                for (i=0; i<k; i++) {
                    sum = sum + pow(_temperature[_tempPerm[radiationElement[element->BoundaryInfo->GebhardtFactors->Elements[j]].NodeIndexes[i]]], 4.0);
                }
                s = sum / k;
                memset( _heatTransferCoeff, 0.0, n*sizeof(double) );
                for (i=0; i<n; i++) {
                    _load[i] = *text * s * _stefanBoltzmann;
                }
                
                [self FEMHeatSolution_integrationOverSurfaceElement:element boundaryNumberOfNodes:n radiationBoundaryOfNodes:k model:model solution:solution mesh:solution.mesh];
                for (m=0; m<n; m++) {
                    k1 = _tempPerm[element->NodeIndexes[m]];
                    forceVector[k1] = forceVector[k1] + _force[m];
                }
            }
        }
        
        // We have already added all external temperature contributions to the matrix
        // for the Newton type iteration
        *angleFraction = asum / _emissivity;
        *text = 0.0;
    }
}

-(void)FEMHeatSolution_addHeatGapSolution:(FEMSolution * __nonnull)solution element:(Element_t * __nonnull)element numberOfNodes:(int)n core:(FEMCore * __nonnull)core {
    
    int i, j, k, l, ind[n];
    
    [self FEMHeatSolution_findGapIndexesElement:element indexes:ind numberOfNodes:n solution:solution];
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            k = _tempPerm[element->NodeIndexes[i]];
            l = _tempPerm[ind[j]];
            if (k >= 0 && l >= 0) {
                [core addToMatrixElementForSolution:solution atIndex:k andIndex:l value:-_stiff[i][j]];
            }
        }
    }
}

-(void)FEMHeatSolution_addHeatFluxBC:(NSArray * __nonnull)bc element:(Element_t * __nonnull)element parent:(Element_t * __nonnull)parent numberOfNodes:(int)n forceVector:(double * __nonnull)forceVector core:(FEMCore * __nonnull)core solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilites crsMatrix:(FEMMatrixCRS * __nonnull)crsMatrix bandMatrix:(FEMMatrixBand * __nonnull)bandMatrix integration:(FEMNumericIntegration * __nonnull)integration diffuseConvectiveAnisotropic:(FEMDiffuseConvectiveAnisotropic * __nonnull)diffuseConvectiveAnisotropic diffuseConvectiveGeneralAnisotropic:(FEMDiffuseConvectiveGeneralAnisotropic * __nonnull)diffuseConvectiveGeneralAnisotropic timeIntegration:(FEMTimeIntegration * __nonnull)timeIntegration utilities:(FEMUtilities * __nonnull)utilities {
    
    int i, j, k;
    BOOL found;
    double sum;
    FEMElementDescription *elementDescription;
    FEMMaterial *materialAtID = nil;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    [core getNodes:solution model:model inElement:element resultNodes:_elementNodes numberOfNodes:NULL mesh:nil];
    
    memset( _heatTransferCoeff, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
    memset( _load, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
    
    // BC: -k@T/@n = \epsilon\sigma(T^4 - Text^4)
    _radiationFlag = [listUtilites listGetString:model inArray:bc forVariable:@"radiation" info:&found];
    if (found == YES && [_radiationFlag isEqualToString:@"none"] == NO) {
        found = [core getReal:model forElement:element inArray:bc variableName:@"emissivity" buffer:&buffer listUtilities:listUtilites];
        memset( _nodalEmissivity, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        if (found == YES) {
            memcpy(_nodalEmissivity, buffer.vector, n*sizeof(double));
        } else {
            found = [core getParentMaterialProperty:@"emissivity" forElement:element parentElement:parent model:model listUtilities:listUtilites buffer:&buffer];
            if (found == YES) memcpy(_nodalEmissivity, buffer.vector, n*sizeof(double));
        }
        vDSP_sveD(_nodalEmissivity, 1, &sum, n);
        _emissivity = sum / n;
        
        memset( _aText, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        if ([_radiationFlag isEqualToString:@"idealized"] == YES) {
            found = [core getReal:model forElement:element inArray:bc variableName:@"radiation external temperature" buffer:&buffer listUtilities:listUtilites];
            if (found == YES) {
                memcpy(_aText, buffer.vector, n*sizeof(double));
            } else {
                found = [core getReal:model forElement:element inArray:bc variableName:@"external temperature" buffer:&buffer listUtilities:listUtilites];
                if (found == YES) memcpy(_aText, buffer.vector, n*sizeof(double));
            }
        } else {
            [self FEMHeatSolution_diffuseGrayRadiationModel:model solution:solution core:core element:element numberOfNodes:n forceVector:forceVector angleFraction:&_visibleFraction text:&_text timeIntegration:timeIntegration utilities:utilities];
            
            if ([listUtilites listGetLogical:model inArray:bc forVariable:@"radiation boundary open" info:&found] == YES) {
                found = [core getReal:model forElement:element inArray:bc variableName:@"radiation external temperature" buffer:&buffer listUtilities:listUtilites];
                if (found == YES) {
                    memcpy(_aText, buffer.vector, n*sizeof(double));
                } else {
                    found = [core getReal:model forElement:element inArray:bc variableName:@"external temperature" buffer:&buffer listUtilities:listUtilites];
                    if (found == YES) memcpy(_aText, buffer.vector, n*sizeof(double));
                }
                if (_visibleFraction >= 1.0) {
                    for (i=0; i<n; i++) {
                        _aText[i] = _text;
                    }
                } else {
                    for (i=0; i<n; i++) {
                        _aText[i] = pow( ( (1.0 - _visibleFraction) * pow(_aText[i], 4.0) + _visibleFraction * pow(_text, 4.0) ), 0.25);
                    }
                }
            } else {
                for (i=0; i<n; i++) {
                    _aText[i] = _text;
                }
            }
        }
        
        // Add our own contribution to surface temperature (and external if
        // using linear type iteration or idealized radiation)
        for (j=0; j<n; j++) {
            k = _tempPerm[element->NodeIndexes[j]];
            _text = _aText[j];
            
            if (_heatGapBC == NO && _newtonLinearization == YES) {
                _heatTransferCoeff[j] = _emissivity * 4.0 * pow(_temperature[k], 3.0) * _stefanBoltzmann;
                _load[j] = _emissivity * (3.0 * pow(_temperature[k], 4.0) + pow(_text, 4.0)) * _stefanBoltzmann;
            } else {
                _heatTransferCoeff[j] = _emissivity *
                     (pow(_temperature[k], 3.0) + pow(_temperature[k], 2.0) * _text + _temperature[k] * pow(_text, 2.0) + pow(_text, 3.0)) * _stefanBoltzmann;
                _load[j] = _heatTransferCoeff[j] * _text;
            }
        }
    } // end of radiation
    
    found = [core getReal:model forElement:element inArray:bc variableName:@"heat transfert coefficient" buffer:&buffer listUtilities:listUtilites];
    if (found == YES) {
        memset( _work, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        memset( _aText, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
        memcpy(_work, buffer.vector, n*sizeof(double));
        found = [core getReal:model forElement:element inArray:bc variableName:@"external temperature" buffer:&buffer listUtilities:listUtilites];
        if (found) memcpy(_aText, buffer.vector, n*sizeof(double));
        for (j=0; j<n; j++) {
            // BC: -k@T/@n = \alpha(T - Text)
            k = _tempPerm[element->NodeIndexes[j]];
            _load[j] = _load[j] + _work[j] * _aText[j];
            _heatTransferCoeff[j] = _heatTransferCoeff[j] + _work[j];
        }
    }
    
    // BC: -k@T/@n = (rho*L)*v.n
    // Heating related to pulling is possible only in cases where pull velocity is described
    if ([listUtilites listGetLogical:model inArray:bc forVariable:@"phase change" info:&found] == YES) {
        elementDescription = [FEMElementDescription sharedElementDescription];
        
        found = [core getReal:model forElement:element inArray:bc variableName:@"phase change 1" buffer:&buffer listUtilities:listUtilites];
        for (i=0; i<n; i++) {
            _phaseVelocity[0][i] = buffer.vector[i];
        }
        found = [core getReal:model forElement:element inArray:bc variableName:@"phase change 2" buffer:&buffer listUtilities:listUtilites];
        for (i=0; i<n; i++) {
            _phaseVelocity[1][i] = buffer.vector[i];
        }
        found = [core getReal:model forElement:element inArray:bc variableName:@"phase change 3" buffer:&buffer listUtilities:listUtilites];
        for (i=0; i<n; i++) {
            _phaseVelocity[2][i] = buffer.vector[i];
        }
        
        // Ensure that the latent heat and density come from the same side
        found = [core getParentMaterialProperty:@"latent heat" forElement:element parentElement:parent model:model listUtilities:listUtilites buffer:&buffer];
        memcpy(_latentHeat, buffer.vector, n*sizeof(double));
        if (parent == NULL) {
            NSLog(@"FEMHeatSolution:FEMHeatSolution_addHeatFluxBC: parent not associated.\n");
        } else {
            k = [(model.bodies)[parent->BodyID-1][@"material"] intValue];
            materialAtID = (model.materials)[k-1];
            found = [core getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&buffer listUtilities:listUtilites];
             if (found) memcpy(_density, buffer.vector, n*sizeof(double));
        }
        
        // This could be rather put as a new type of BC into the assembly routine and
        // then the Normal could be taken at the proper Gaussian integration points
        double u = 0.0;
        double v = 0.0;
        BOOL check = YES;
        [elementDescription normalVectorForBDElement:element boundaryNodes:_elementNodes mesh:solution.mesh paraU:&u paraV:&v check:&check normals:_normal];
        
        for (i=0; i<n; i++) {
            sum = 0.0;
            for (j=0; j<3; j++) {
                sum = sum + (_normal[j]*_phaseVelocity[j][i]);
            }
            _load[i] = _load[i] + _latentHeat[i] * _density[i] * sum;
        }
    }
    
    // BC: -k@T/@n = g
    found = [core getReal:model forElement:element inArray:bc variableName:@"heat flux" buffer:&buffer listUtilities:listUtilites];
    if (found == YES) {
        for (i=0; i<n; i++) {
            _load[i] = _load[i] + buffer.vector[i];
        }
    }
    
    // Get element matrix and rhs due to boundary conditions...
    if (model.coordinates == cartesian) {
        Dimensions_t dimensions = (Dimensions_t){.mat1=(2*solution.mesh.maxElementDofs), .mat2=(2*solution.mesh.maxElementDofs), .vec=(2*solution.mesh.maxElementDofs)};
        [diffuseConvectiveAnisotropic diffuseConvectiveBoundaryMatrix:_stiff boundaryVector:_force dimensions:dimensions loadVector:_load nodalAlpha:_heatTransferCoeff element:element numberOfNodes:n nodes:_elementNodes mesh:solution.mesh integration:integration];
    } else {
        Dimensions_t dimensions = (Dimensions_t){.mat1=(2*solution.mesh.maxElementDofs), .mat2=(2*solution.mesh.maxElementDofs), .vec=(2*solution.mesh.maxElementDofs)};
        [diffuseConvectiveGeneralAnisotropic diffuseConvectiveGeneralBoundaryMatrix:_stiff boundaryVector:_force dimensions:dimensions loadVector:_load nodalAlpha:_heatTransferCoeff element:element numberOfNodes:n nodes:_elementNodes model:model mesh:solution.mesh integration:integration];
    }
    
    // Update global matrices from local matrices
    int rows = 2*solution.mesh.maxElementDofs;
    int cols = 2*solution.mesh.maxElementDofs;
    if (_transientAssembly == YES && _constantBulk == NO) {
        memset( *_mass, 0.0, ((2*solution.mesh.maxElementDofs)*(2*solution.mesh.maxElementDofs))*sizeof(double) );
        [core defaultFirstOrderTime:model inSolution:solution forElement:element realMass:_mass realStiff:_stiff realForce:_force stiffRows:&rows stiffCols:&cols timeIntegration:timeIntegration utilities:utilities];
    }
    if (_heatGapBC == YES) {
        [self FEMHeatSolution_addHeatGapSolution:solution element:element numberOfNodes:n core:core];
    }
    [core defaultUpdateEquations:model inSolution:solution forElement:element realStiff:_stiff realForce:_force stiffRows:&rows stiffCols:&cols crsMatrix:crsMatrix bandMatrix:bandMatrix];
    
    if (buffer.vector != NULL) free_dvector(buffer.vector, 0, buffer.m-1);
}

-(void)FEMHeatSolution_addGlobalTimeSolution:(FEMSolution * __nonnull)solution {
    
    int i, j, k, n;
    double force[1];
    variableArraysContainer *varContainers = NULL;
    matrixArraysContainer *matContainers = NULL;
    
    FEMTimeIntegration *timeIntegration = [[FEMTimeIntegration alloc] init];
    
    varContainers = solution.variable.getContainers;
    matContainers = solution.matrix.getContainers;
    if (saveValues != varContainers->Values) {
        if (stiff != NULL) {
            free_dmatrix(stiff, 0, 0, 0, n1-1);
            free_dmatrix(mass, 0, 0, 0, n1-1);
            free_dmatrix(x, 0, n1-1, 0, k1-1);
        }
        n1 = 0;
        for (i=0; i<solution.matrix.numberOfRows; i++) {
            n1 = max(n1, (int)(matContainers->Rows[i+1]-matContainers->Rows[i]));
        }
        k1 = varContainers->size2PrevValues;
        stiff = doublematrix(0, 0, 0, n1-1);
        mass = doublematrix(0, 0, 0, n1-1);
        x = doublematrix(0, n1-1, 0, k1-1);
        
        saveValues = varContainers->Values;
    }
    
    int rows = 1;
    for (i=0; i<solution.matrix.numberOfRows; i++) {
        n = 0;
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            stiff[0][n] = matContainers->Values[j];
            mass[0][n] = matContainers->MassValues[j];
            for (k=0; k<varContainers->size2PrevValues; k++) {
                x[n][k] = varContainers->PrevValues[matContainers->Cols[j]][k];
            }
            n++;
        }
        force[0] = matContainers->RHS[i];
        matContainers->Force[i][0] = force[0];
        k = min(solution.doneTime, solution.order);
        [timeIntegration bdfLocalInSolution:solution numberOfNodes:n dt:_dt massMatrix:mass stiffMatrix:stiff force:force prevSolution:x order:k rows:&rows cols:&n1];
        
        n = 0;
        for (j=matContainers->Rows[i]; j<=matContainers->Rows[i+1]-1; j++) {
            matContainers->Values[j] = stiff[0][n];
            n++;
        }
        matContainers->RHS[i] = force[0];
    }
}

-(BOOL)FEMHeatSolution_checkLatentHeatModel:(FEMModel * __nonnull)model {
    
    int i, j, k, t, body_id, eq_id, n;
    FEMMesh *mesh;
    Element_t *elements = NULL;
    BOOL failure = NO, any, checkLatentHeatRelease, found, phaseChange;
    NSString *phaseModel;
    FEMListUtilities *listUtilities;
    FEMEquation *equationAtID = nil;
    FEMMaterial *materialAtID = nil;
    listBuffer phaseChangeIntervals = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    mesh = (FEMMesh *)model.mesh;
    elements = mesh.getElements;
    
    listUtilities = [[FEMListUtilities alloc] init];
    for (t=0; t<mesh.numberOfBulkElements; t++) {
        // Check if this element belongs to a body where temperature
        // has been calculated
        any = NO;
        for (i=0; i<elements[t].Type.NumberOfNodes; i++) {
            if (_tempPerm[elements[t].NodeIndexes[i]] < 0) {
                any = YES;
                break;
            }
        }
        if (any == YES) continue;
        
        body_id = elements[t].BodyID;
        eq_id = [(model.bodies)[body_id-1][@"equation"] intValue];
        if (eq_id < 1) eq_id = 1;
        if (eq_id > model.numberOfEquations) eq_id = model.numberOfEquations;
        
        equationAtID = (model.equations)[eq_id-1];
        phaseModel = [listUtilities listGetString:model inArray:equationAtID.valuesList forVariable:@"phase change model" info:&found];
        
        phaseChange = (found == YES && [phaseModel isEqualToString:@"none"] == NO) ? YES : NO;
        checkLatentHeatRelease = NO;
        if (phaseChange == YES) {
            checkLatentHeatRelease = [listUtilities listGetLogical:model inArray:equationAtID.valuesList forVariable:@"check latent heat release" info:&found];
        }
        if (!(phaseChange == YES && checkLatentHeatRelease == YES)) continue;
        
        n = elements[t].Type.NumberOfNodes;
        
        k = [(model.bodies)[body_id-1][@"material"] intValue];
        materialAtID = (model.materials)[k-1];
        found = [listUtilities listGetConstRealArray:model inArray:materialAtID.valuesList forVariable:@"phase change intervals" buffer:&phaseChangeIntervals];
        
        for (k=0; k<n; k++) {
            i = _tempPerm[elements[t].NodeIndexes[k]];
            for (j=0; j<phaseChangeIntervals.n; j++) {
                if ((_temperature[i] < phaseChangeIntervals.matrix[0][j] && _prevSolution[i] > phaseChangeIntervals.matrix[1][j]) ||
                    (_temperature[i] > phaseChangeIntervals.matrix[1][j] && _prevSolution[i] < phaseChangeIntervals.matrix[0][j])) {
                    failure = YES;
                    break;
                }
            }
            if (failure == YES) break;
        }
        if (failure == YES) break;
    }
    
    if (phaseChangeIntervals.matrix != NULL) {
        free_dmatrix(phaseChangeIntervals.matrix, 0, phaseChangeIntervals.m-1, 0, phaseChangeIntervals.n-1);
    }
    
    return failure;
}

- (id)init
{
    self = [super init];
    if (self) {
        _allocationDone = NO;
        _constantBulk = NO;
        _newtonLinearization = NO;
        _cols = 0;
        _doneTime = 0;
        _rows = 0;
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

-(void)deallocation:(FEMSolution * __nonnull)solution {

    int n = solution.mesh.maxElementDofs;
    variableArraysContainer *tempContainers = solution.variable.getContainers;
    
    if (stiff != NULL) free_dmatrix(stiff, 0, 0, 0, n1-1);
    if (mass != NULL) free_dmatrix(mass, 0, 0, 0, n1-1);
    if (x != NULL) free_dmatrix(x, 0, n1-1, 0, k1-1);
    
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

-(void)solutionComputer:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int)timeStep transientSimulation:(BOOL)transient {
    
    int i, j, k, l, n, nb=0, nd, t, bf_id, body_id, compressibilityModel=-1, eq_id, iter, mat_id, nsdofs=0, nonLinearIter, newtonIter, smartHeaterNode=-1, smartHeaterBC=-1;
    int *indexes = NULL, *flowPerm = NULL;
    double at=0.0, at0, C1, dist, dt0, jx, jy, jz, cumulativeTime, newtonTol, meltPoint=0.0, minDist, nonLinearTol, norm, powerRelax, powerSensitivity, powerTimeScale=0.0,
           prevNorm, referencePressure=0.0, relax, relativeChange, s, saveRelax, smartTol=0.0, specificHeatRatio, st, sum, totat, totst, xave=0.0, yave=0.0;
    double controlPoint[3], *forceVector, *flowSolution = NULL;
    BOOL all, bubbles, checkLatentHeatRelease=NO, found, heatFluxBC, gotIt, isRadiation, firstTime, gotMeltPoint,
         integralHeaterControl, heaterControlLocal = NO, phaseChange=NO, smartHeaterControl, smartHeaterAverage, smartTolReached=NO, stabilize = YES, transientHeaterControl, useBubbles;
    NSString *convectionField, *stabilizeFlag, *convectionFlag, *compressibilityFlag;
    NSArray *bc;
    Element_t *elements = NULL, *element = NULL, *parent = NULL;
    Nodes_t *meshNodes = NULL;
    FEMMesh *mesh;
    FEMVariable *flowSol, *densitySol;
    FEMBodyForce *bodyForceAtID = nil;
    FEMBoundaryCondition *boundaryConditionAtID = nil;
    FEMMaterial *materialAtID = nil;
    FEMEquation *equationAtID = nil;
    matrixArraysContainer *matContainers = NULL;
    variableArraysContainer *tempContainers = NULL, *flowSolContainers = NULL;
    listBuffer realWork = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer vector = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer tensor = { NULL, NULL, NULL, NULL, 0, 0, 0};
    FEMCore *core = [FEMCore sharedCore];
    FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    FEMMatrixCRS *crsMatrix = [[FEMMatrixCRS alloc] init];
    FEMMatrixBand *bandMatrix = [[FEMMatrixBand alloc] init];
    FEMDiffuseConvectiveAnisotropic *diffuseConvectiveAnisotropic = [[FEMDiffuseConvectiveAnisotropic alloc] init];
    FEMDiffuseConvectiveGeneralAnisotropic *diffuseConvectiveGeneralAnisotropic = [[FEMDiffuseConvectiveGeneralAnisotropic alloc] init];
    FEMElementUtils *elementUtils = [[FEMElementUtils alloc] init];
    FEMMaterialModels *materialModels = [[FEMMaterialModels alloc] init];
    FEMDifferentials *differentials = [[FEMDifferentials alloc] init];
    FEMCoordinateSystems *coordinatesSystems = [[FEMCoordinateSystems alloc] init];
    FEMTimeIntegration *timeIntegration;
    
    static Element_t* (*getActiveElementIMP)(id, SEL, int, FEMSolution*, FEMModel*) = nil;
    static int (*getEquationIDForElementIMP)(id, SEL, Element_t*, FEMModel*) = nil;
    static NSString* (*listGetStringIMP)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*) = nil;
    static int (*getMaterialIDForElementIMP)(id, SEL, Element_t*, FEMModel*) = nil;
    static BOOL (*listGetLogicalIMP)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*) = nil;
    static void (*getNodesIMP)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*) = nil;
    static void (*getScalarLocalFieldIMP)(id, SEL, double*, int, NSString*, Element_t*, FEMSolution*, FEMModel*, int*) = nil;
    static BOOL (*getRealIMP)(id, SEL, FEMModel*, Element_t*, NSArray*, NSString*, listBuffer*, FEMListUtilities*) = nil;
    static BOOL (*listGetRealArrayIMP)(id, SEL, FEMModel*, NSArray*, NSString*, int, int*, listBuffer*) = nil;
    static double (*listGetConstRealIMP)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*, double*, double*) = nil;
    static void (*getVectorLocalFieldIMP)(id, SEL, double**, int, int, NSString*, Element_t*, FEMSolution*, FEMModel*, int*) = nil;
    static void (*effectiveHeatCapacityElementIMP)(id, SEL, Element_t*, int, FEMMaterial*, FEMModel*, FEMListUtilities*, BOOL) = nil;
    static int (*getBodyForceIDForElementIMP)(id, SEL, Element_t*, FEMModel*) = nil;
    static void (*defaultUpdateEquationsIMP)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double*, int*, int*, FEMMatrixCRS*, FEMMatrixBand*) = nil;
    static void (*diffuseConvectiveComposeMassMatrixIMP)(id, SEL, double**, double**, double*, double*, double*, double*, double*, double***, BOOL, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, BOOL, BOOL, Element_t*, int, Nodes_t*, FEMSolution*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*) = nil;
    static void (*diffuseConvectiveGeneralComposeMassMatrixIMP)(id, SEL, double**, double**, double*, double*, double*, double*, double*, double***, BOOL, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, BOOL, Element_t*, int, Nodes_t*, FEMSolution*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMMaterialModels*, FEMDifferentials*, FEMCoordinateSystems*, FEMListUtilities*) = nil;
    static void (*defaultFirstOrderTimeIMP)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double**, double*, int*, int*, FEMTimeIntegration*, FEMUtilities*) = nil;
    static void (*updateGlobalEquationsModelIMP)(id, SEL, FEMModel*, FEMSolution*, Element_t*, double**, double*, double*, int, int, int*, int*, int*, BOOL*, FEMMatrixCRS*, FEMMatrixBand*) = nil;
    static void (*defaultUpdateMassIMP)(id, SEL, double**, Element_t*, FEMSolution*, FEMModel*) = nil;
    static void (*defaultUpdateDampIMP)(id, SEL, double**, Element_t*, FEMSolution*, FEMModel*) = nil;
    static void (*condensateStiffIMP)(id, SEL, double**, double*, int, double*) = nil;
    
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        if (!getActiveElementIMP) {
            getActiveElementIMP = (Element_t* (*)(id, SEL, int, FEMSolution*, FEMModel*))
            [core methodForSelector: @selector(getActiveElement:solution:model:)];
        }
        if (!getEquationIDForElementIMP) {
            getEquationIDForElementIMP = (int (*)(id, SEL, Element_t*, FEMModel*))
            [core methodForSelector: @selector(getEquationIDForElement:model:)];
        }
        if (!listGetStringIMP) {
            listGetStringIMP = (NSString* (*)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*))
            [listUtilities methodForSelector: @selector(listGetString:inArray:forVariable:info:)];
        }
        if (!getMaterialIDForElementIMP) {
            getMaterialIDForElementIMP = (int (*)(id, SEL, Element_t*, FEMModel*))
            [core methodForSelector: @selector(getMaterialIDForElement:model:)];
        }
        if (!listGetLogicalIMP) {
            listGetLogicalIMP = (BOOL (*)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*))
            [listUtilities methodForSelector: @selector(listGetLogical:inArray:forVariable:info:)];
        }
        if (!getNodesIMP) {
            getNodesIMP = (void (*)(id, SEL, FEMSolution*, FEMModel*, Element_t*, Nodes_t*, int*, FEMMesh*))
            [core methodForSelector: @selector(getNodes:model:inElement:resultNodes:numberOfNodes:mesh:)];
        }
        if (!getScalarLocalFieldIMP) {
            getScalarLocalFieldIMP = (void (*)(id, SEL, double*, int, NSString*, Element_t*, FEMSolution*, FEMModel*, int*))
            [core methodForSelector: @selector(getScalarLocalField:sizeField:name:element:solution:model:timeStep:)];
        }
        if (!getRealIMP) {
            getRealIMP = (BOOL (*)(id, SEL, FEMModel*, Element_t*, NSArray*, NSString*, listBuffer*, FEMListUtilities*))
            [core methodForSelector: @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:)];
        }
        if (!listGetRealArrayIMP) {
            listGetRealArrayIMP = (BOOL (*)(id, SEL, FEMModel*, NSArray*, NSString*, int, int*, listBuffer*))
            [listUtilities methodForSelector: @selector(listGetRealArray:inArray:forVariable:numberOfNodes:indexes:buffer:)];
        }
        if (!listGetConstRealIMP) {
            listGetConstRealIMP = (double (*)(id, SEL, FEMModel*, NSArray*, NSString*, BOOL*, double*, double*))
            [listUtilities methodForSelector: @selector(listGetConstReal:inArray:forVariable:info:minValue:maxValue:)];
        }
        if (!getVectorLocalFieldIMP) {
            getVectorLocalFieldIMP = (void (*)(id, SEL, double**, int, int, NSString*, Element_t*, FEMSolution*, FEMModel*, int*))
            [core methodForSelector: @selector(getVectorLocalField:size1Field:size2Field:name:element:solution:model:timeStep:)];
        }
        if (!effectiveHeatCapacityElementIMP) {
            effectiveHeatCapacityElementIMP = (void (*)(id, SEL, Element_t*, int, FEMMaterial*, FEMModel*, FEMListUtilities*, BOOL))
            [self methodForSelector: @selector(FEMHeatSolution_effectiveHeatCapacityElement:numberOfNodes:material:model:listUtilities:transientSimulation:)];
        }
        if (!getBodyForceIDForElementIMP) {
            getBodyForceIDForElementIMP = (int (*)(id, SEL, Element_t*, FEMModel*))
            [core methodForSelector: @selector(getBodyForceIDForElement:model:)];
        }
        if (!defaultUpdateEquationsIMP) {
            defaultUpdateEquationsIMP = (void (*)(id, SEL, FEMModel*, FEMSolution*, Element_t*, double**, double*, int*, int*, FEMMatrixCRS*, FEMMatrixBand*))
            [core methodForSelector: @selector(defaultUpdateEquations:inSolution:forElement:realStiff:realForce:stiffRows:stiffCols:crsMatrix:bandMatrix:)];
        }
        if (!diffuseConvectiveComposeMassMatrixIMP) {
            diffuseConvectiveComposeMassMatrixIMP = (void (*)(id, SEL, double**, double**, double*, double*, double*, double*, double*, double***, BOOL, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, BOOL, BOOL, Element_t*, int, Nodes_t*, FEMSolution*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMMaterialModels*, FEMDifferentials*, FEMListUtilities*))
            [diffuseConvectiveAnisotropic methodForSelector: @selector(diffuseConvectiveComposeMassMatrix:stiffMatrix:forceVector:loadVector:timeDerivativeTerm:zeroDegreeTerm:convectionTerm:diffusionTerm:phaseChange:nodalTemperature:enthalpy:velocityX:velocitY:velocityZ:meshVeloX:meshVeloY:meshVeloZ:nodalViscosity:nodaldensity:nodalPressure:nodalPressureDt:nodalPressureCoeff:compressible:stabilize:useBubbles:element:numberOfNodes:nodes:solution:core:mesh:model:integration:materialModels:differentials:listUtilities:)];
        }
        if (!diffuseConvectiveGeneralComposeMassMatrixIMP) {
            diffuseConvectiveGeneralComposeMassMatrixIMP = (void (*)(id, SEL, double**, double**, double*, double*, double*, double*, double*, double***, BOOL, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, BOOL, BOOL, Element_t*, int, Nodes_t*, FEMSolution*, FEMCore*, FEMMesh*, FEMModel*, FEMNumericIntegration*, FEMMaterialModels*, FEMDifferentials*, FEMCoordinateSystems*, FEMListUtilities*))
            [diffuseConvectiveGeneralAnisotropic methodForSelector: @selector(diffuseConvectiveGeneralComposeMassMatrix:stiffMatrix:forceVector:loadVector:timeDerivativeTerm:zeroDegreeTerm:convectionTerm:diffusionTerm:phaseChange:nodalTemperature:enthalpy:velocityX:velocitY:velocityZ:meshVeloX:meshVeloY:meshVeloZ:nodalViscosity:nodaldensity:nodalPressure:nodalPressureDt:nodalPressureCoeff:compressible:stabilize:element:numberOfNodes:nodes:solution:core:mesh:model:integration:materialModels:differentials:coordinatesSystems:listUtilities:)];
        }
        if (!defaultFirstOrderTimeIMP) {
            defaultFirstOrderTimeIMP = (void (*)(id, SEL, FEMModel*, FEMSolution*, Element_t *, double**, double**, double*, int*, int*, FEMTimeIntegration*, FEMUtilities*))
            [core methodForSelector: @selector(defaultFirstOrderTime:inSolution:forElement:realMass:realStiff:realForce:stiffRows:stiffCols:timeIntegration:utilities:)];
        }
        if (!updateGlobalEquationsModelIMP) {
            updateGlobalEquationsModelIMP = (void (*)(id, SEL, FEMModel*, FEMSolution*, Element_t*, double**, double*, double*, int, int, int*, int*, int*, BOOL*, FEMMatrixCRS*, FEMMatrixBand*))
            [core methodForSelector: @selector(updateGlobalEquationsModel:inSolution:element:localStiffMatrix:forceVector:localForce:size:dofs:nodeIndexes:rows:cols:rotateNT:crsMatrix:bandMatrix:)];
        }
        if (!defaultUpdateMassIMP) {
            defaultUpdateMassIMP = (void (*)(id, SEL, double**, Element_t*, FEMSolution*, FEMModel*))
            [core methodForSelector: @selector(defaultUpdateMass:element:solution:model:)];
        }
        if (!defaultUpdateDampIMP) {
            defaultUpdateDampIMP = (void (*)(id, SEL, double**, Element_t*, FEMSolution*, FEMModel*))
            [core methodForSelector: @selector(defaultUpdateDamp:element:solution:model:)];
        }
        if (!condensateStiffIMP) {
            condensateStiffIMP = (void (*)(id, SEL, double**, double*, int, double*))
            [core methodForSelector: @selector(condensateStiff:force:numberOfNodes:force1:)];
        }
    });
    
    if (transient == YES) {
        timeIntegration = [[FEMTimeIntegration alloc] init];
    }

    mesh = (FEMMesh *)model.mesh;
    elements = mesh.getElements;
    meshNodes = mesh.getNodes;
    
    isRadiation = NO;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        isRadiation = ([listUtilities listCheckPresentVariable:@"radiation" inArray:boundaryCondition.valuesList] == YES) ? YES : NO;
        if (isRadiation == YES) break;
    }
    // The view and Gebhardt factors may change. If this is necessary, this is done without this method.
    // This method is called in the start as it may affect the matrix topplogy
    // Newton lineariarization option is needed only when there is radiation
    if (isRadiation == YES) {
        //TODO: implement the radiation factors computation
    }
    
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
        
        _rows = 2*n;
        _cols = 2*n;
        
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
                _heaterArea = NULL;
                _heaterDensity = NULL;
                _heaterSource = NULL;
                _heaterScaling = NULL;
                _heaterTarget = NULL;
                _smarterHeaters = NULL;
                _integralHeaters = NULL;
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
                 fatal("FEMHeatSolution:solutionComputer", "Memory allocation error.");
            }
            memset( _smarterHeaters, false, n*sizeof(bool) );
            memset( _integralHeaters, false, n*sizeof(bool) );
        }
        
        if (smartHeaterControl == YES) {
            if (_allocationDone == YES) {
                free_dvector(_xx, 0, tempContainers->sizeValues-1);
                free_dvector(_yy, 0, tempContainers->sizeValues-1);
                free_dvector(_forceHeater, 0, tempContainers->sizeValues-1);
                _xx = NULL;
                _yy = NULL;
                _forceHeater = NULL;
            }
            
            _xx = doublevec(0, tempContainers->sizeValues-1);
            _yy = doublevec(0, tempContainers->sizeValues-1);
            _forceHeater = doublevec(0, tempContainers->sizeValues-1);
            if (_xx == NULL || _yy == NULL || _forceHeater == NULL) {
                fatal("FEMHeatSolution:solutionComputer", "Memory allocation error.");
            }
            memset( _xx, 0.0, tempContainers->sizeValues*sizeof(double) );
            memset( _yy, 0.0, tempContainers->sizeValues*sizeof(double) );
            memset( _forceHeater, 0.0, tempContainers->sizeValues*sizeof(double) );
        }
        
        _allocationDone = YES;
    }
    
    // Do some additional initialization and go for it
    _dt = timeStep;
    if (isRadiation == YES) {
        _stefanBoltzmann = [listUtilities listGetConstReal:model inArray:model.constants.valuesList forVariable:@"stefan boltzmann" info:&found minValue:NULL maxValue:NULL];
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
    
    if ((solution.solutionInfo)[@"nonlinear system maximum iterations"] != nil) {
        nonLinearIter = [(solution.solutionInfo)[@"nonlinear system maximum iterations"] intValue];
    } else nonLinearIter = 1;
    
    nonLinearTol = 0.0;
    if ((solution.solutionInfo)[@"nonlinear system convergence tolerance"] != nil) {
        nonLinearTol = [(solution.solutionInfo)[@"nonlinear system convergence tolerance"] doubleValue];
    }
    
    // Newton linearization option is only needed when there is radiation
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
    if (newtonIter == 0) _newtonLinearization = YES;
    
    if ((solution.solutionInfo)[@"nonlinear system relaxation factor"] != nil) {
        relax = [(solution.solutionInfo)[@"nonlinear system relaxation factor"] doubleValue];
    } else relax = 1.0;
    
    _transientAssembly = transient;
    found = NO;
    dt0 = 0.0;
    if ((solution.solutionInfo)[@"steady state transition time step"] != nil) {
        dt0 = [(solution.solutionInfo)[@"steady state transition time step"] doubleValue];
        found = YES;
    } else {
        if ((solution.solutionInfo)[@"smart heater time scale"] != nil) {
            dt0 = [(solution.solutionInfo)[@"smart heater time scale"] doubleValue];
            found = YES;
        }
    }
    if (found == YES && _dt > dt0) _transientAssembly = NO;
    
    transientHeaterControl = NO;
    smartHeaterAverage = NO;
    if (smartHeaterControl == YES) {
        
        // Mark the smart heaters
        memset( _smarterHeaters, false, model.numberOfBodyForces*sizeof(bool) );
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
        
        minDist = 0.0;
        smartHeaterNode = -1;
        smartHeaterNode = [listUtilities listGetInteger:model inArray:bodyForceAtID.valuesList forVariable:@"smart heater control node" info:&found minValue:NULL maxValue:NULL];
        if (found == YES) {
            smartHeaterNode = smartHeaterNode -1;
        } else if (found == NO) {
            found = [listUtilities listGetConstRealArray:model inArray:bodyForceAtID.valuesList forVariable:@"smart heater control point" buffer:&realWork];
            if (found == YES) {
                for (i=0; i<3; i++) {
                    controlPoint[i] = realWork.matrix[i][0];
                }
                minDist = DBL_MAX;
                for (l=0; l<model.numberOfNodes; l++) {
                    if (_tempPerm[l] < 0) continue;
                    
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
            NSLog(@"FEMHeatSolution:solutionComputer: found control point at distance: %f.\n", sqrt(minDist));
            NSLog(@"FEMHeatSolution:solutionComputer: control point index: %d.\n", smartHeaterNode);
            if (realWork.matrix != NULL) {
                free_dmatrix(realWork.matrix, 0, realWork.m-1, 0, realWork.n-1);
                realWork.matrix = NULL;
            }
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
                fatal("FEMHeatSolution:solutionComputer", "Smart heater boundary / Phase change is undefined.");
            }
            
            boundaryConditionAtID = (model.boundaryConditions)[smartHeaterBC];
            meltPoint = [listUtilities listGetConstReal:model inArray:boundaryConditionAtID.valuesList forVariable:@"smart heater temperature" info:&found minValue:NULL maxValue:NULL];
            if (found == NO) {
                for (FEMMaterial *material in model.materials) {
                    meltPoint = [listUtilities listGetConstReal:model inArray:material.valuesList forVariable:@"melting point" info:&found minValue:NULL maxValue:NULL];
                    if (found == YES) break;
                }
                if (found == NO) fatal("FEMHeatSolution:solutionComputer", "Smart heater temperature / melting point is undefined.");
            }
            
            // Find the node related to temperature control
            if ((solution.solutionInfo)[@"smart heater average"] != nil) {
                smartHeaterAverage = [(solution.solutionInfo)[@"smart heater average"] boolValue];
            }
            if (smartHeaterAverage == NO) {
                jx = -DBL_MAX;
                j = -1;
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
    
        if (transient == YES && _dt < powerTimeScale) {
            transientHeaterControl = YES;
            NSLog(@"FEMHeatSolution:solutionComputer: using transient heater control.\n");
        } else {
            transientHeaterControl = NO;
            NSLog(@"FEMHeatSolution:solutionComputer: using steady-state heater control.\n");
        }
    
        if (solution.doneTime != _doneTime) {
            _prevPowerScaling = _powerScaling;
            _doneTime = solution.doneTime;
        }
    }
    
    if (integralHeaterControl == YES) {
        NSLog(@"FEMHeatSolution:solutionComputer: using integral heater control.\n");
        memset( _integralHeaters, false, model.numberOfBodyForces*sizeof(bool) );
        i = 0;
        for (FEMBodyForce *bodyForce in model.bodyForces) {
            if ([listUtilities listCheckPresentVariable:@"integral heat source" inArray:bodyForce.valuesList] == YES) {
                _integralHeaters[i] = true;
            }
            i++;
        }
    }
    
    if ((solution.solutionInfo)[@"constant bulk system"] != nil) {
        _constantBulk = [(solution.solutionInfo)[@"constant bulk system"] boolValue];
    }
    
    saveRelax = relax;
    cumulativeTime = 0.0;
    
    firstTime = YES;
    _prevSolution = doublevec(0, _localNodes-1);
    
    FEMNumericIntegration *integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) fatal("FEMDiffuseConvectiveAnisotropic:diffuseConvectiveComposeyMassMatrix", "Allocation error in FEMNumericIntegration.");
    
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
            
            NSLog(@"FEMHeatSolution:solutionComputer:\n");
            NSLog(@"FEMHeatSolution:solutionComputer:\n");
            NSLog(@"FEMHeatSolution:solutionComputer: -----------------------------------------------------------\n");
            NSLog(@"FEMHeatSolution:solutionComputer: TEMPERATURE ITERATION %d.\n", iter);
            NSLog(@"FEMHeatSolution:solutionComputer: -----------------------------------------------------------\n");
            NSLog(@"FEMHeatSolution:solutionComputer:\n");
            NSLog(@"FEMHeatSolution:solutionComputer: Starting Assembly...\n");
            
            if (_constantBulk == YES && matContainers->BulkValues != NULL) {
                memcpy(matContainers->Values, matContainers->BulkValues, matContainers->sizeBulkValues*sizeof(double));
                memcpy(matContainers->RHS, matContainers->BulkRHS, matContainers->sizeBulkRHS*sizeof(double));
                goto jump;
            }
            
            [core defaultInitializeSolution:solution model:model];
            
            if (smartHeaterControl == YES || integralHeaterControl == YES) {
                if (smartHeaterControl == YES) memset( _forceHeater, 0.0, tempContainers->sizeValues*sizeof(double) );
                memset( _heaterArea, 0.0, model.numberOfBodyForces*sizeof(double) );
                memset( _heaterSource, 0.0, model.numberOfBodyForces*sizeof(double) );
                memset( _heaterScaling, 0.0, model.numberOfBodyForces*sizeof(double) );
                memset( _heaterTarget, 0.0, model.numberOfBodyForces*sizeof(double) );
                
                for (t=0; t<solution.numberOfActiveElements; t++) {
                    element = [core getActiveElement:t solution:solution model:model];
                    bf_id = [core getBodyForceIDForElement:element model:model];
                    bodyForceAtID = (model.bodyForces)[bf_id-1];
                    if (bodyForceAtID == nil) continue;
                    if (!_smarterHeaters[bf_id-1] || !_integralHeaters[bf_id-1]) continue;
                    
                    n = element->Type.NumberOfNodes;
                    
                    mat_id = [core getMaterialIDForElement:element model:model];
                    materialAtID = (model.materials)[mat_id-1];
                    
                    memset( _density, 0.0, n*sizeof(double) );
                    found = [core getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&vector listUtilities:listUtilities];
                    if (found == YES) memcpy(_density, vector.vector, n*sizeof(double));

                    
                    memset( _load, 0.0, n*sizeof(double) );
                    found = [core getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"heat source" buffer:&vector listUtilities:listUtilities];
                    if (found == YES) memcpy(_load, vector.vector, n*sizeof(double));
                    
                    _s = [elementUtils elementArea:element numberOfNodes:n mesh:solution.mesh nodel:model];
                    
                    if (model.coordinates == axis_symmetric || model.coordinates == cylindric_symmetric) _s = 2.0 * M_PI * _s;
                    
                    _heaterSource[bf_id-1] = _heaterSource[bf_id-1] + _s * cblas_ddot(n, _density, 1, _load, 1) / n;
                    _heaterArea[bf_id-1] = _heaterArea[bf_id-1] + _s;
                    vDSP_sveD(_density, 1, &sum, n);
                    _heaterDensity[bf_id-1] = _heaterDensity[bf_id-1] + _s * sum / n;
                }
                i = 0;
                for (FEMBodyForce *bodyForce in model.bodyForces) {
                    if (_integralHeaters[i] || _smarterHeaters[i]) {
                        _heaterDensity[i] = _heaterDensity[i] / _heaterArea[i];
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
            bodyForceAtID = nil;
            nb = 0;
            
            startAdvanceOutput((char *)[@"FEMHeatSolution" UTF8String], (char *)[@"Assembly:" UTF8String]);
            at = cputime();
            // Bulk elements
            for (t=0; t<solution.numberOfActiveElements; t++) {
                
                advanceOutput(t, solution.numberOfActiveElements, NULL, NULL);
                
                // Check if this element belongs to a body where temperature
                // should be calculated
                element = getActiveElementIMP(core, @selector(getActiveElement:solution:model:), t, solution, model);
                if (element->BodyID != body_id) {
                    body_id = element->BodyID;
                    
                    eq_id = getEquationIDForElementIMP(core, @selector(getEquationIDForElement:model:), element, model);
                    equationAtID = (model.equations)[eq_id-1];
                    convectionFlag = listGetStringIMP(listUtilities, @selector(listGetString:inArray:forVariable:info:), model, equationAtID.valuesList, @"convection", &found);
                    
                    mat_id = getMaterialIDForElementIMP(core, @selector(getMaterialIDForElement:model:), element, model);
                    materialAtID = (model.materials)[mat_id-1];
                    
                    compressibilityFlag = listGetStringIMP(listUtilities, @selector(listGetString:inArray:forVariable:info:), model, materialAtID.valuesList, @"compressibility model", &found);
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
                    
                    _phaseModel = listGetStringIMP(listUtilities, @selector(listGetString:inArray:forVariable:info:), model, equationAtID.valuesList, @"phase change model", &found);
                    if (found == NO) _phaseModel = listGetStringIMP(listUtilities, @selector(listGetString:inArray:forVariable:info:), model, materialAtID.valuesList, @"phase change model", &found);
                    phaseChange = (found == YES && [_phaseModel isEqualToString:@"none"] == NO) ? YES : NO;
                    if (phaseChange == YES) {
                        checkLatentHeatRelease = listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, equationAtID.valuesList, @"check latent heat release", &found);
                    }
                }
                
                n = element->Type.NumberOfNodes;
                getNodesIMP(core, @selector(getNodes:model:inElement:resultNodes:numberOfNodes:mesh:), solution, model, element, _elementNodes, NULL, nil);
                getScalarLocalFieldIMP(core, @selector(getScalarLocalField:sizeField:name:element:solution:model:timeStep:), _localTemperature, solution.mesh.maxElementDofs, nil, element, solution, model, NULL);
                
                // Get element material parameters
                memset( _heatCapacity, 0.0, n*sizeof(double) );
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"heat capacity", &vector, listUtilities);
                if (found == YES) memcpy(_heatCapacity, vector.vector, n*sizeof(double));
                
                memset( **_heatConductivity, 0.0, (3*3*solution.mesh.maxElementDofs)*sizeof(double) );
                found = listGetRealArrayIMP(listUtilities, @selector(listGetRealArray:inArray:forVariable:numberOfNodes:indexes:buffer:), model, materialAtID.valuesList, @"heat conductivity", n, element->NodeIndexes, &tensor);
                if (found == YES) {
                    if (tensor.m == 1) {
                        for (i=0; i<3; i++) {
                            for (j=0; j<n; j++) {
                                _heatConductivity[i][i][j] = tensor.tensor[0][0][j];
                            }
                        }
                    } else if (tensor.n == 1) {
                        for (i=0; i<min(3, tensor.m); i++) {
                            for (j=0; j<n; j++) {
                                _heatConductivity[i][i][j] = tensor.tensor[i][0][j];
                            }
                        }
                    } else {
                        for (i=0; i<min(3, tensor.m); i++) {
                            for (j=0; j<min(3, tensor.n); j++) {
                                for (k=0; k<n; k++) {
                                    _heatConductivity[i][j][k] = tensor.tensor[i][j][k];
                                }
                            }
                        }
                    }
                }
                
                if (compressibilityModel == perfect_gas1) {
                    // Read specific heat ratio
                    specificHeatRatio = listGetConstRealIMP(listUtilities, @selector(listGetConstReal:inArray:forVariable:info:minValue:maxValue:), model, materialAtID.valuesList, @"specific heat ratio", &found, NULL, NULL);
                    if (found == NO) specificHeatRatio = 5.0/3.0;
                    
                    // For an ideal gas, \gamma, c_p and R are really a constant.
                    // GasConstant is an array only since HeatCapacity formally is
                    for (i=0; i<n; i++) {
                        _gasConstant[i] = (specificHeatRatio - 1.0) * _heatCapacity[i] /  specificHeatRatio;
                    }
                    
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"pressure coefficient", &vector, listUtilities);
                    if (found == YES) {
                        memcpy(_pressureCoeff, vector.vector, n*sizeof(double));
                    } else {
                        for (i=0; i<n; i++) {
                            _pressureCoeff[i] = 1.0;
                        }
                    }
                } else if (compressibilityModel == thermal) {
                    memset( _referenceTemperature, 0.0, n*sizeof(double) );
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"reference temperature", &vector, listUtilities);
                    if (found == YES) memcpy(_referenceTemperature, vector.vector, n*sizeof(double));
                    
                    memset( _heatExpansionCoeff, 0.0, n*sizeof(double) );
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"heat expansion coefficient", &vector, listUtilities);
                    if (found == YES) memcpy(_heatExpansionCoeff, vector.vector, n*sizeof(double));
                    
                    memset( _density, 0.0, n*sizeof(double) );
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"density", &vector, listUtilities);
                    if (found == YES) memcpy(_density, vector.vector, n*sizeof(double));
                    for (i=0; i<n; i++) {
                        _density[i] = _density[i] * ( 1.0 - _heatExpansionCoeff[i] * (_localTemperature[i] - _referenceTemperature[i]) );
                    }
                    
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"pressure coefficient", &vector, listUtilities);
                    if (found == YES) {
                        memcpy(_pressureCoeff, vector.vector, n*sizeof(double));
                    } else {
                        for (i=0; i<n; i++) {
                            _pressureCoeff[i] = _localTemperature[i] *
                               _heatExpansionCoeff[i] / (1.0 - _heatExpansionCoeff[i] * (_localTemperature[i] - _referenceTemperature[i]));
                        }
                    }
                } else if (compressibilityModel == user_defined1) {
                    if (densitySol != nil) {
                        getScalarLocalFieldIMP(core, @selector(getScalarLocalField:sizeField:name:element:solution:model:timeStep:), _density, solution.mesh.maxElementDofs, @"density", element, solution, model, NULL);
                    } else {
                        memset( _density, 0.0, n*sizeof(double) );
                        found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"density", &vector, listUtilities);
                        if (found == YES) memcpy(_density, vector.vector, n*sizeof(double));
                    }
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"pressure coefficient", &vector, listUtilities);
                    if (found == YES) {
                        memcpy(_pressureCoeff, vector.vector, n*sizeof(double));
                    } else memset( _pressureCoeff, 0.0, n*sizeof(double) );
                } else {
                    memset( _pressureCoeff, 0.0, n*sizeof(double) );
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"pressure coefficient", &vector, listUtilities);
                    if (found == YES) memcpy(_pressureCoeff, vector.vector, n*sizeof(double));
                    
                    memset( _density, 0.0, n*sizeof(double) );
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"density", &vector, listUtilities);
                    if (found == YES) memcpy(_density, vector.vector, n*sizeof(double));
                }
                
                // Take pressure deviation p_d as the dependent variable p = p_0 + p_d.
                // For perfect gas, read p_0
                if (compressibilityModel != incompressible) {
                    referencePressure = listGetConstRealIMP(listUtilities, @selector(listGetConstReal:inArray:forVariable:info:minValue:maxValue:), model, materialAtID.valuesList, @"reference pressure", &found, NULL, NULL);
                    if (found == NO) referencePressure = 0.0;
                }
                
                heaterControlLocal = NO;
                memset( _load, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                memset( _pressure, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                memset( _dPressureDt, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                
                // Check for convection model
                C1 = 1.0;
                memset( _u, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                memset( _v, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                memset( _w, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                
                memset( *_mu, 0.0, (3*solution.mesh.maxElementDofs)*sizeof(double) );
                getVectorLocalFieldIMP(core, @selector(getVectorLocalField:size1Field:size2Field:name:element:solution:model:timeStep:), _mu, 3, solution.mesh.maxElementDofs, @"mesh velocity", element, solution, model, NULL);
                
                if ([convectionFlag isEqualToString:@"constant"] == YES) {
                    
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"convection velocity 1", &vector, listUtilities);
                    if (found == YES) {
                        memcpy(_u, vector.vector, n*sizeof(double));
                    } else {
                        found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, equationAtID.valuesList, @"convection velocity 1", &vector, listUtilities);
                        if (found == YES) memcpy(_u, vector.vector, n*sizeof(double));
                    }
                    
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"convection velocity 2", &vector, listUtilities);
                    if (found == YES) {
                        memcpy(_v, vector.vector, n*sizeof(double));
                    } else {
                        found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, equationAtID.valuesList, @"convection velocity 2", &vector, listUtilities);
                        if (found == YES) memcpy(_v, vector.vector, n*sizeof(double));
                    }
                    
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"convection velocity 3", &vector, listUtilities);
                    if (found == YES) {
                        memcpy(_w, vector.vector, n*sizeof(double));
                    } else {
                         found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, equationAtID.valuesList, @"convection velocity 3", &vector, listUtilities);
                        if (found == YES) memcpy(_w, vector.vector, n*sizeof(double));
                    }
                } else if ([convectionFlag isEqualToString:@"computed"] == YES && flowSol != nil) {
                    for (i=0; i<n; i++) {
                        k = flowPerm[element->NodeIndexes[i]];
                        if (k >= 0) {
                            _pressure[i] = flowSolution[nsdofs*(k+1)-1] + referencePressure;
                            switch (compressibilityModel) {
                                case perfect_gas1:
                                    _density[i] = _pressure[i] / (_gasConstant[i] * _localTemperature[i]);
                                    break;
                            }
                            if (transient == YES) {
                                _dPressureDt[i] = ( flowSolution[nsdofs*(k+1)-1] - flowSolContainers->PrevValues[nsdofs*(k+1)-1][0]) / timeStep;
                            }
                            switch (nsdofs) {
                                case 3:
                                    _u[i] = flowSolution[nsdofs*k];
                                    _v[i] = flowSolution[nsdofs*k+1];
                                    _w[i] = 0.0;
                                    break;
                                case 4:
                                    _u[i] = flowSolution[nsdofs*k];
                                    _v[i] = flowSolution[nsdofs*k+1];
                                    _w[i] = flowSolution[nsdofs*k+2];
                                    break;
                            }
                        } else {
                            _u[i] = 0.0;
                            _v[i] = 0.0;
                            _w[i] = 0.0;
                        }
                    }
                } else if ([convectionFlag isEqualToString:@"computed"] == YES) {
                    NSLog(@"FEMHeatSolution:solutionComputer: convection model specified but no accociated flow field present.\n");
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
                if (phaseChange == YES) {
                    effectiveHeatCapacityElementIMP(core, @selector(FEMHeatSolution_effectiveHeatCapacityElement:numberOfNodes:material:model:listUtilities:transientSimulation:), element, n, materialAtID, model, listUtilities, transient);
                } else {
                    for (i=0; i<n; i++) {
                        _heatCapacity[i] = _density[i] * _heatCapacity[i];
                    }
                }
                
                memset( _viscosity, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
                
                // Add body forces if any
                bf_id = getBodyForceIDForElementIMP(core, @selector(getBodyForceIDForElement:model:), element, model);
                bodyForceAtID = (model.bodyForces)[bf_id-1];
                if (bodyForceAtID != nil) {
                    // Frictional viscous heating
                    if (listGetLogicalIMP(listUtilities, @selector(listGetLogical:inArray:forVariable:info:), model, bodyForceAtID.valuesList, @"friction heat", &found) == YES) {
                        found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, materialAtID.valuesList, @"viscosity", &vector, listUtilities);
                        if (found == YES) memcpy(_viscosity, vector.vector, n*sizeof(double));
                    }
                }
                
                // Get heat source
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"heat source", &vector, listUtilities);
                if (found == YES) {
                    for (i=0; i<n; i++) {
                        _load[i] = _density[i] * vector.vector[i];
                    }
                }
                
                if (smartHeaterControl == YES && _newtonLinearization == YES && smartTolReached == YES) {
                    if (_smarterHeaters[bf_id-1]) {
                        heaterControlLocal = YES;
                        if (transientHeaterControl == YES) {
                            for (i=0; i<n; i++) {
                                _load[i] = _prevPowerScaling * _load[i];
                            }
                            _heaterScaling[bf_id-1] = _prevPowerScaling;
                        }
                    }
                }
                
                if (integralHeaterControl == YES) {
                    if (_integralHeaters[bf_id-1]) {
                        for (i=0; i<n; i++) {
                            _load[i] = _load[i] * _heaterScaling[bf_id-1];
                        }
                    }
                }
                
                memset( _c0, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
            
                // Note at this point HeatCapacity = \rho * c_p or \rho * (c_p - R)
                // and C1 = 0 (diffusion) or 1 (convection)
            
                // Perfusion (added as suggested by Matthias Zenker)
                memset( _perfusionRate, 0.0, n*sizeof(double) );
                found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"perfusion rate", &vector, listUtilities);
                if (found == YES) {
                    memcpy(_perfusionRate, vector.vector, n*sizeof(double));
                    
                    memset( _perfusionRefTemperature, 0.0, n*sizeof(double) );
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"perfusion reference temperature", &vector, listUtilities);
                    if (found == YES) memcpy(_perfusionRefTemperature, vector.vector, n*sizeof(double));
                    
                    memset( _perfusionDensity, 0.0, n*sizeof(double) );
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"perfusion density", &vector, listUtilities);
                    if (found == YES) memcpy(_perfusionDensity, vector.vector, n*sizeof(double));
                    
                    memset( _perfusionHeatCapacity, 0.0, n*sizeof(double) );
                    found = getRealIMP(core, @selector(getReal:forElement:inArray:variableName:buffer:listUtilities:), model, element, bodyForceAtID.valuesList, @"perfusion heat capacity", &vector, listUtilities);
                    if (found == YES) memcpy(_perfusionHeatCapacity, vector.vector, n*sizeof(double));
                    for (i=0; i<n; i++) {
                        _c0[i] = _perfusionHeatCapacity[i] * _perfusionRate[i] * _perfusionDensity[i];
                        _load[i] = _load[i] + _c0[i] * _perfusionRefTemperature[i];
                    }
                }
            
                // Get element local matrices and RHS vectors
                
                // We initialize these arrays before calling the assembly methods. Elmer does it inside the routines
                // but we do it here so that we don't need to pass the sizes of the arrays
                memset( *_stiff, 0.0, ((2*solution.mesh.maxElementDofs)*(2*solution.mesh.maxElementDofs))*sizeof(double) );
                memset( *_mass, 0.0, ((2*solution.mesh.maxElementDofs)*(2*solution.mesh.maxElementDofs))*sizeof(double) );
                memset( _force, 0.0, (2*solution.mesh.maxElementDofs)*sizeof(double) );
                if (model.coordinates == cartesian) {
                    double heatCapacity[n];
                    double mux[n], muy[n], muz[n];
                    for (i=0; i<n; i++) {
                        heatCapacity[i] = C1 * _heatCapacity[i];
                        mux[i] = _mu[0][i];
                        muy[i] = _mu[1][i];
                        muz[i] = _mu[2][i];
                    }
                    diffuseConvectiveComposeMassMatrixIMP(diffuseConvectiveAnisotropic,@selector(diffuseConvectiveComposeMassMatrix:stiffMatrix:forceVector:loadVector:timeDerivativeTerm:zeroDegreeTerm:convectionTerm:diffusionTerm:phaseChange:nodalTemperature:enthalpy:velocityX:velocitY:velocityZ:meshVeloX:meshVeloY:meshVeloZ:nodalViscosity:nodaldensity:nodalPressure:nodalPressureDt:nodalPressureCoeff:compressible:stabilize:useBubbles:element:numberOfNodes:nodes:solution:core:mesh:model:integration:materialModels:differentials:listUtilities:), _mass, _stiff, _force, _load, _heatCapacity, _c0, heatCapacity, _heatConductivity, _phaseSpatial, _localTemperature, _enthalpy, _u, _v, _w, mux, muy, muz, _viscosity, _density, _pressure, _dPressureDt, _pressureCoeff, (compressibilityModel != incompressible) ? YES : NO, stabilize, useBubbles, element, n, _elementNodes, solution, core, solution.mesh, model, integration, materialModels, differentials, listUtilities);
                } else {
                    double heatCapacity[n];
                    double mux[n], muy[n], muz[n];
                    for (i=0; i<n; i++) {
                        heatCapacity[i] = C1 * _heatCapacity[i];
                        mux[i] = _mu[0][i];
                        muy[i] = _mu[1][i];
                        muz[i] = _mu[2][i];
                    }
                    
                    diffuseConvectiveGeneralComposeMassMatrixIMP(diffuseConvectiveGeneralAnisotropic, @selector(diffuseConvectiveGeneralComposeMassMatrix:stiffMatrix:forceVector:loadVector:timeDerivativeTerm:zeroDegreeTerm:convectionTerm:diffusionTerm:phaseChange:nodalTemperature:enthalpy:velocityX:velocitY:velocityZ:meshVeloX:meshVeloY:meshVeloZ:nodalViscosity:nodaldensity:nodalPressure:nodalPressureDt:nodalPressureCoeff:compressible:stabilize:element:numberOfNodes:nodes:solution:core:mesh:model:integration:materialModels:differentials:coordinatesSystems:listUtilities:), _mass, _stiff, _force, _load, _heatCapacity, _c0, heatCapacity, _heatConductivity, _phaseSpatial, _localTemperature, _enthalpy, _u, _v, _w, mux, muy, muz, _viscosity, _density, _pressure, _dPressureDt, _pressureCoeff, (compressibilityModel != incompressible) ? YES : NO, stabilize, element, n, _elementNodes, solution, core, solution.mesh, model, integration, materialModels, differentials, coordinatesSystems, listUtilities);
                }
                
                if (heaterControlLocal == YES && transientHeaterControl == NO) {
                    if (_transientAssembly == YES && _constantBulk == NO) {
                        defaultFirstOrderTimeIMP(core, @selector(defaultFirstOrderTime:inSolution:forElement:realMass:realStiff:realForce:stiffRows:stiffCols:timeIntegration:utilities:), model, solution, element, _mass, _stiff, _force, &_rows, &_cols, timeIntegration, utilities);
                    }
                    if (indexes == NULL || n != nb) {
                        if (indexes != NULL) free_ivector(indexes, 0, nb-1);
                        indexes = intvec(0, n-1);
                        nb = n;
                    }
                    for (i=0; i<n; i++) {
                        indexes[i] = _tempPerm[element->NodeIndexes[i]];
                    }
                    updateGlobalEquationsModelIMP(core, @selector(updateGlobalEquationsModel:inSolution:element:localStiffMatrix:forceVector:localForce:size:dofs:nodeIndexes:rows:cols:rotateNT:crsMatrix:bandMatrix:), model, solution, element, _stiff, _forceHeater, _force, n, 1, indexes, &_rows, &_cols, NULL, crsMatrix, bandMatrix);
                } else {
                    bubbles = (useBubbles == YES && stabilize == NO && ([convectionFlag isEqualToString:@"computed"] == YES || [convectionFlag isEqualToString:@"constant"] == YES)) ? YES : NO;
                    
                    // If time dependent simulation, add mass matrix to stiff matrix
                    memset( _timeForce, 0.0, (2*solution.mesh.maxElementDofs)*sizeof(double) );
                    if (_transientAssembly == YES) {
                        if (_constantBulk == YES) {
                            defaultUpdateMassIMP(core, @selector(defaultUpdateMass:element:solution:model:), _mass, element, solution, model);
                        } else {
                            defaultFirstOrderTimeIMP(core, @selector(defaultFirstOrderTime:inSolution:forElement:realMass:realStiff:realForce:stiffRows:stiffCols:timeIntegration:utilities:), model, solution, element, _mass, _stiff, _force, &_rows, &_cols, timeIntegration, utilities);
                        }
                    } else if (solution.nOfEigenValues > 0) {
                        defaultUpdateDampIMP(core, @selector(defaultUpdateDamp:element:solution:model:), _mass, element, solution, model);
                    }
                    // Update global matrices and local matrices
                    if (bubbles == YES) {
                        condensateStiffIMP(core, @selector(condensateStiff:force:numberOfNodes:force1:), _stiff, _force, n, _timeForce);
                    }
                    
                    defaultUpdateEquationsIMP(core, @selector(defaultUpdateEquations:inSolution:forElement:realStiff:realForce:stiffRows:stiffCols:crsMatrix:bandMatrix:),
                                              model, solution, element, _stiff, _force, &_rows, &_cols, crsMatrix, bandMatrix);
                }
            } // Bulk elements
            
            [core defaultFinishBulkAssemblySolution:solution bulkUpdate:NULL];
            
        jump:
            at = cputime() -  at;
            
            // Neumann & Newton boundary conditions
            for (t=0; t<solution.mesh.numberOfBoundaryElements; t++) {
                 element = [core getBoundaryElement:solution atIndex:t];
                if ([core isActiveBoundaryElement:element inSolution:solution model:model] == NO) continue;
                
                n = element->Type.NumberOfNodes;
                
                // Check that the dimension of element is suitable for fluxes
                if ([core isFluxElement:element mesh:mesh] == NO) continue;
                
                bc = [core getBoundaryCondition:model forElement:element];
                if (bc == nil) continue;
                
                // This check whether there are any Dirichlet conditions on the smart heater boundary.
                // If there are, the r.h.s must be zero as there can possibly not be any effect on temperature
                if (heaterControlLocal == YES && transientHeaterControl == NO) {
                    if ([listUtilities listCheckPresentVariable:solution.variable.name inArray:bc] == YES) {
                        memset( core.indexStore, -1, core.sizeIndexStore*sizeof(int) );
                        nd = [core getElementDofsSolution:solution model:model forElement:element atIndexes:core.indexStore disableDiscontinuousGalerkin:NULL];
                        for (i=0; i<nd; i++) {
                            _forceHeater[_tempPerm[core.indexStore[i]]] = 0.0;
                        }
                    }
                }
                
                heatFluxBC = [listUtilities listGetLogical:model inArray:bc forVariable:@"heat fluc bc" info:&found];
                if (found == YES && heatFluxBC == NO) continue;
                
                _heatGapBC = [listUtilities listGetLogical:model inArray:bc forVariable:@"heat gap" info:&found];
                [self FEMHeatSolution_addHeatFluxBC:bc element:element parent:parent numberOfNodes:n forceVector:forceVector core:core solution:solution model:model listUtilities:listUtilities crsMatrix:crsMatrix bandMatrix:bandMatrix integration:integration diffuseConvectiveAnisotropic:diffuseConvectiveAnisotropic diffuseConvectiveGeneralAnisotropic:diffuseConvectiveGeneralAnisotropic timeIntegration:timeIntegration utilities:utilities];
                
                if (_heatGapBC == YES) {
                    [self FEMHeatSolution_findGapIndexesElement:element indexes:_indexes numberOfNodes:n solution:solution];
                    memcpy(_saveIndexes, element->NodeIndexes, n*sizeof(double));
                    memcpy(element->NodeIndexes, _indexes, n*sizeof(double));
                    [self FEMHeatSolution_addHeatFluxBC:bc element:element parent:parent numberOfNodes:n forceVector:forceVector core:core solution:solution model:model listUtilities:listUtilities crsMatrix:crsMatrix bandMatrix:bandMatrix integration:integration diffuseConvectiveAnisotropic:diffuseConvectiveAnisotropic diffuseConvectiveGeneralAnisotropic:diffuseConvectiveGeneralAnisotropic timeIntegration:timeIntegration utilities:utilities];
                    memcpy(element->NodeIndexes, _saveIndexes, n*sizeof(double));
                }
            } // Neumann & Newton BCs
            
            if (transient == YES && _constantBulk == YES) [self FEMHeatSolution_addGlobalTimeSolution:solution];
            
            [core defaultFinishAssemblySolution:solution model:model timeIntegration:timeIntegration utilities:utilities];
            NSLog(@"FEMHeatSolution:solutionComputer: Assembly done.\n");
            
            [core dirichletBoundaryConditions:model inSolution:solution usingOffset:NULL offDiaginalMatrix:NULL];
            
            // Solve the system and check for convergence
            st = cputime();
            
            prevNorm = norm;
            
            xave = 0.0;
            yave = 0.0;
            if (smartHeaterControl == YES && _newtonLinearization == YES && smartTolReached == YES) {
                if (transientHeaterControl == NO) {
                    
                     [solution.solutionInfo setObject:@YES forKey:@"skip compute nonlinear change"];
                    
                    if ((solution.solutionInfo)[@"nonlinear system relaxation factor"] != nil) {
                        relax = [(solution.solutionInfo)[@"nonlinear system relaxation factor"] doubleValue];
                        if (relax != 1.0) {
                            [solution.solutionInfo setValue:@(1.0) forKey:@"nonlinear system relaxation factor"];
                        }
                    } else {
                        relax = 1.0;
                    }
                    
                    [core solveSystemMatrix:solution.matrix rhs:_forceHeater result:_xx norm:&norm dofs:1 solution:solution model:model];
                    [core solveSystemMatrix:solution.matrix rhs:matContainers->RHS result:_yy norm:&norm dofs:1 solution:solution model:model];
                    
                    [solution.solutionInfo setObject:@NO forKey:@"skip compute nonlinear change"];
                } else {
                    [core solveSystemMatrix:solution.matrix rhs:matContainers->RHS result:_temperature norm:&norm dofs:1 solution:solution model:model];
                    memcpy(_yy, _temperature, tempContainers->sizeValues*sizeof(double));
                }
                
                if (smartHeaterAverage == NO) {
                    xave = _xx[_tempPerm[smartHeaterNode]];
                    yave = _yy[_tempPerm[smartHeaterNode]];
                } else {
                    j = 0;
                    double sum1 = 0.0;
                    
                    for (k=mesh.numberOfBulkElements; k<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; k++) {
                        if (elements[k].BoundaryInfo->Constraint == smartHeaterBC+1) {
                            l = elements[k].Type.NumberOfNodes;
                            j = j + l;
                            sum = 0.0;
                            sum1 = 0.0;
                            for (i=0; i<elements[k].Type.NumberOfNodes; i++) {
                                sum = sum + _xx[_tempPerm[elements[k].NodeIndexes[i]]];
                                sum1 = sum1 + _yy[_tempPerm[elements[k].NodeIndexes[i]]];
                            }
                            xave = xave + sum;
                            yave = yave + sum1;
                        }
                    }
                    xave = xave / j;
                    yave = yave / j;
                    [listUtilities addConstRealInClassList:model.simulation theVariable:@"res: smart heater temperature" withValue:&yave orUsingBlock:nil string:nil];
                }
                
                if (transientHeaterControl == NO) {
                    if (tempContainers->NonLinValues != NULL) {
                        tempContainers->NonLinValues = _temperature;
                    }
                    
                    _powerScaling = (meltPoint - yave) / xave;
                    for (i=0; i<tempContainers->sizeValues; i++) {
                        _temperature[i] = _yy[i] + _powerScaling * _xx[i];
                    }
                    
                    // The change is computed seperately for the controlled temperature field
                    [core computeChange:solution model:model isSteadyState:NO nsize:&_localNodes values:_temperature values0:NULL sizeValues0:NULL];
                    norm = solution.variable.norm;
                }
                
                if (_dt > powerTimeScale) {
                    if (relax != 0.0) [solution.solutionInfo setValue:@(relax) forKey:@"nonlinear system relaxation factor"];
                }
            } else {
                norm = [core findSolution:solution model:model backRorateNT:NULL];
            }
            
            if (smartHeaterControl == YES || integralHeaterControl == YES) {
                [listUtilities addConstRealInClassList:model.simulation theVariable:@"res: heater power scaling" withValue:&_powerScaling orUsingBlock:nil string:nil];
                NSLog(@"FEMHeatSolution:solutionComputer: heater control information.\n");
                for (i=0; i<model.numberOfBodyForces; i++) {
                    if (!(_smarterHeaters[i] || _integralHeaters[i])) continue;
                    if (_smarterHeaters[i]) _heaterScaling[i] = _powerScaling;
                    NSLog(@"FEMHeatSolution:solutionComputer: heater for body %d.\n", i+1);
                    if (_smarterHeaters[i]) NSLog(@"FEMHeatSolution:solutionComputer: heater type: smart heater.\n");
                    if (_integralHeaters[i]) NSLog(@"FEMHeatSolution:solutionComputer: heater type: integral heater.\n");
                    
                    NSLog(@"FEMHeatSolution:solutionComputer: heater volutme (m^3): %f.\n", _heaterArea[i]);
                    s = _heaterSource[i] * _heaterScaling[i];
                    NSLog(@"FEMHeatSolution:solutionComputer: heater power (W): %f.\n", s);
                    
                    NSLog(@"FEMHeatSolution:solutionComputer: heater scaling: %f.\n", _heaterScaling[i]);
                    NSLog(@"FEMHeatSolution:solutionComputer: heater power density (W/kg): %e.\n", s/(_heaterDensity[i]*_heaterArea[i]));
                    
                    if (_smarterHeaters[i]) {
                        double value = s/(_heaterDensity[i]*_heaterArea[i]);
                        [listUtilities addConstRealInClassList:model.simulation theVariable:@"res: heater power density" withValue:&value orUsingBlock:nil string:nil];
                    }
                }
            }
            
            st = cputime() - st;
            totat = totat + at;
            totst = totst + st;
            NSLog(@"FEMHeatSolution:solutionComputer: iter: %d, Assembly (s): %f %f.\n", iter, at, totat);
            NSLog(@"FEMHeatSolution:solutionComputer: iter: %d, Solve (s): %f %f.\n", iter, st, totst);
            
            // If modeling phase change (and if requested by the user), check if any
            // node has jumped over the phase interval, and if so, reduce time step
            // and or relaxation and recompute
            if (phaseChange == YES && checkLatentHeatRelease == YES) {
                if ([self FEMHeatSolution_checkLatentHeatModel:model] == YES) {
                    memcpy(_temperature, _prevSolution, _localNodes*sizeof(double));
                    norm = prevNorm;
                    
                    if (transient == YES) {
                        _dt = _dt / 2;
                        solution.dt = _dt;
                        NSLog(@"FEMHeatSolution:solutionComputer: latent heat release check: reducing time step to: %f.\n", _dt);
                    } else {
                        relax = relax / 2;
                        [solution.solutionInfo setValue:@(relax) forKey:@"nonlinear system relaxation factor"];
                        NSLog(@"FEMHeatSolution:solutionComputer: latent heat release check: reducing relaxation to: %f.\n", relax);
                    }
                    continue;
                }
                if (transient == NO) memcpy(_prevSolution, _temperature, _localNodes*sizeof(double));
            }
            
            relativeChange = solution.variable.nonLinChange;
            NSLog(@"FEMHeatSolution:solutionComputer: result norm: %e.\n", norm);
            NSLog(@"FEMHeatSolution:solutionComputer: relative change: %e.\n", relativeChange);
            
            if (relativeChange < newtonTol || iter >= newtonIter) _newtonLinearization = YES;
            if (relativeChange < nonLinearTol && (smartHeaterControl == NO || smartTolReached == YES)) break;
            
            if (smartHeaterControl == YES) {
                if (relativeChange < smartTol) {
                    smartTolReached = YES;
                    memcpy(_yy, _temperature, tempContainers->sizeValues*sizeof(double));
                }
            }
        } // Non linear iteration loop
        
        if (transientHeaterControl == YES) {
            if ((solution.solutionInfo)[@"smart heater relaxation factor"] != nil) {
                powerRelax = [(solution.solutionInfo)[@"smart heater relaxation factor"] doubleValue];
            } else powerRelax = 1.0;
            if ((solution.solutionInfo)[@"smart heater power sensitivity"] != nil) { 
                powerSensitivity = [(solution.solutionInfo)[@"smart heater power sensitivity"] doubleValue];
            } powerSensitivity = 4.0;
            _powerScaling = _powerScaling * (1.0 + powerSensitivity * powerRelax * (meltPoint/yave - 1.0));
            if ([(solution.solutionInfo)[@"smart heater transient speedup"] boolValue] == YES) {
                for (i=0; i<tempContainers->sizeValues; i++) {
                    _temperature[i] = _temperature[i] * ( 1.0 + powerRelax *(meltPoint/yave - 1.0) );
                }
            }
            memcpy(_yy, _temperature, tempContainers->sizeValues*sizeof(double));
        }
        
        // Compute cumulative time done by now and time remaining
        if (transient == NO) break;
        cumulativeTime = cumulativeTime + _dt;
        _dt = timeStep - cumulativeTime;
        
        if (vector.vector != NULL) {
            free_dvector(vector.vector, 0, vector.m-1);
            vector.vector = NULL;
        }
        if (tensor.tensor != NULL) {
            free_d3tensor(tensor.tensor, 0, tensor.m-1, 0, tensor.n-1, 0, tensor.p-1);
            tensor.tensor = NULL;
        }
        if (indexes != NULL) {
            free_ivector(indexes, 0, nb-1);
            indexes = NULL;
        }
    } // time interval
    
    solution.dt = timeStep;
    [solution.solutionInfo setValue:@(saveRelax) forKey:@"nonlinear system relaxation factor"];
    free_dvector(_prevSolution, 0, _localNodes-1);
    [integration deallocation:mesh];
    
    if ([(solution.solutionInfo)[@"adaptive mesh refinement"] boolValue] == YES) {
        // TODO: implement mesh refinement
    }
}

@end
