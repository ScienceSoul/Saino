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

static int k1 = 0, n1 = 0;
static double *saveValues = NULL;
static double **stiff = NULL, **mass = NULL, **x = NULL;

enum {
  PHASE_SPATIAL_1 = 1,
  PHASE_SPATIAL_2,
  PHASE_TEMPORAL
};

@interface FEMHeatSolution ()
-(void)FEMHeatSolution_findGapIndexesElement:(Element_t *)element indexes:(int *)indexes numberOfNodes:(int)n solution:(FEMSolution *)solution;
-(void)FEMHeatSolution_effectiveHeatCapacityElement:(Element_t *)element numberOfNodes:(int)n material:(FEMMaterial *)material model:(FEMModel *)model transientSimulation:(BOOL)transient;
-(void)FEMHeatSolution_integrationOverSurfaceElement:(Element_t*)element boundaryNumberOfNodes:(int)n radiationBoundaryOfNodes:(int)m model:(FEMModel *)model solution:(FEMSolution *)solution mesh:(FEMMesh *)mesh;
-(void)FEMHeatSolution_diffuseGrayRadiationModel:(FEMModel *)model solution:(FEMSolution *)solution kernel:(FEMKernel *)kernel element:(Element_t *)element numberOfNodes:(int)n forceVector:(double *)forceVector angleFraction:(double *)angleFraction text:(double *)text;
-(void)FEMHeatSolution_addHeatGapSolution:(FEMSolution *)solution element:(Element_t *)element numberOfNodes:(int)n kernel:(FEMKernel *)kernel;
-(void)FEMHeatSolution_addHeatFluxBC:(NSArray *)bc element:(Element_t *)element parent:(Element_t *)parent numberOfNodes:(int)n forceVector:(double *)forceVector kernel:(FEMKernel *)kernel solution:(FEMSolution *)solution model:(FEMModel *)model;
-(void)FEMHeatSolution_addGlobalTimeSolution:(FEMSolution *)solution;
-(BOOL)FEMHeatSolution_checkLatentHeatModel:(FEMModel *)model;
@end

@implementation FEMHeatSolution {
    
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

-(void)FEMHeatSolution_findGapIndexesElement:(Element_t *)element indexes:(int *)indexes numberOfNodes:(int)n solution:(FEMSolution *)solution {
    
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

-(void)FEMHeatSolution_effectiveHeatCapacityElement:(Element_t *)element numberOfNodes:(int)n material:(FEMMaterial *)material model:(FEMModel *)model transientSimulation:(BOOL)transient {
    
    int i;
    BOOL any, found, specific;
    FEMListUtilities *listUtilites;
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
        _s = max_array(_localTemperature, n) - min_array(_localTemperature, n);
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
    
    listUtilites = [[FEMListUtilities alloc] init];
    specific = [listUtilites listCheckPresentVariable:@"specific enthalpy" inArray:material.valuesList];
    
    switch (_phaseChangeModel) {
        // This phase change model is avaiable only for some type of real entries
        // that have an implemented analytical derivation rule
        case PHASE_SPATIAL_1:
            found = [listUtilites listGetReal:model inArray:material.valuesList forVariable:@"effective heat capacity" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
            if (found == YES) {
                memcpy(_heatCapacity, buffer.vector, n*sizeof(double));
            } else {
                if (specific == YES) {
                    found = [listUtilites listGetDerivativeValue:model inArray:material.valuesList forVariable:@"specific enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer];
                    if (found == YES) memcpy(_heatCapacity, buffer.vector, n*sizeof(double));
                    for (i=0; i<n; i++) {
                        _heatCapacity[i] = _density[i] * _heatCapacity[i];
                    }
                } else {
                    found = [listUtilites listGetDerivativeValue:model inArray:material.valuesList forVariable:@"enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer];
                    if (found == YES)  memcpy(_heatCapacity, buffer.vector, n*sizeof(double));
                }
            }
            break;
            
        // Note that for the 'spatial 2' mode, the evaluation of c_p is done in each integration point
        // and thus Enthalpy and PhaseSpatial flag are used instead of HeatCapacity directly
        case PHASE_SPATIAL_2:
            if (specific == YES) {
                found = [listUtilites listGetReal:model inArray:material.valuesList forVariable:@"specific enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
                if (found == YES) memcpy(_enthalpy, buffer.vector, n*sizeof(double));
                for (i=0; i<n; i++) {
                    _enthalpy[i] = _density[i] * _enthalpy[i];
                }
            } else {
                found = [listUtilites listGetReal:model inArray:material.valuesList forVariable:@"enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
                if (found == YES)  memcpy(_enthalpy, buffer.vector, n*sizeof(double));
            }
            
        case PHASE_TEMPORAL:
            memcpy(_tSolution1, _temperature, _localNodes*sizeof(double));
            
            // When retrieving the value of enthalpy on the previous timeStep, the
            // relevant entries of the Temperature solution in the global vector are
            // tampered in order to make the ListGetReal method work as wanted
            if (specific == YES) {
                found = [listUtilites listGetReal:model inArray:material.valuesList forVariable:@"specific enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
                if (found == YES) memcpy(_work, buffer.vector, n*sizeof(double));
                memcpy(_temperature, _tSolution, _localNodes*sizeof(double));
                for (i=0; i<n; i++) {
                    _work[i] = _work[i] - buffer.vector[i];
                }
            } else {
                found = [listUtilites listGetReal:model inArray:material.valuesList forVariable:@"enthalpy" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
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

-(void)FEMHeatSolution_integrationOverSurfaceElement:(Element_t*)element boundaryNumberOfNodes:(int)n radiationBoundaryOfNodes:(int)m model:(FEMModel *)model solution:(FEMSolution *)solution mesh:(FEMMesh *)mesh {
 
    int i, p, q, t;
    double alpha, detJ, force, s, x, y, z;
    BOOL stat;
    FEMNumericIntegration *integration;
    FEMCoordinateSystems *coordinateSystem;
    GaussIntegrationPoints *IP;
    
    memset( *_stiff, 0.0, ((2*solution.mesh.maxElementDofs)*(2*solution.mesh.maxElementDofs))*sizeof(double) );
    memset( _force, 0.0, (2*solution.mesh.maxElementDofs)*sizeof(double) );
    
    // Integration stuff
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) errorfunct("FEMHeatSolution_integrationOverSurfaceElement", "Allocation error in FEMNumericIntegration!");
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
    GaussQuadratureDeallocation(IP);
    [integration deallocation:mesh];
}

-(void)FEMHeatSolution_diffuseGrayRadiationModel:(FEMModel *)model solution:(FEMSolution *)solution kernel:(FEMKernel *)kernel element:(Element_t *)element numberOfNodes:(int)n forceVector:(double *)forceVector angleFraction:(double *)angleFraction text:(double *)text {
    
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
                    [kernel addFirstOrderTimeModel:model solution:solution element:element massMatrix:_mass stiffMatrix:_stiff force:_force dt:_dt size:n dofs:1 nodeIndexes:perm rows:&rows cols:&cols];
                }
                
                for (m=0; m<n; m++) {
                    k1 = _tempPerm[element->NodeIndexes[m]];
                    for (l=0; l<k; l++) {
                        k2 = _tempPerm[radiationElement[element->BoundaryInfo->GebhardtFactors->Elements[j]].NodeIndexes[l]];
                        [kernel addToMatrixElementForSolution:solution atIndex:k1 andIndex:k2 value:_stiff[m][l]];
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

-(void)FEMHeatSolution_addHeatGapSolution:(FEMSolution *)solution element:(Element_t *)element numberOfNodes:(int)n kernel:(FEMKernel *)kernel {
    
    int i, j, k, l, ind[n];
    
    [self FEMHeatSolution_findGapIndexesElement:element indexes:ind numberOfNodes:n solution:solution];
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            k = _tempPerm[element->NodeIndexes[i]];
            l = _tempPerm[ind[j]];
            if (k >= 0 && l >= 0) {
                [kernel addToMatrixElementForSolution:solution atIndex:k andIndex:l value:-_stiff[i][j]];
            }
        }
    }
}

-(void)FEMHeatSolution_addHeatFluxBC:(NSArray *)bc element:(Element_t *)element parent:(Element_t *)parent numberOfNodes:(int)n forceVector:(double *)forceVector kernel:(FEMKernel *)kernel solution:(FEMSolution *)solution model:(FEMModel *)model {
    
    int i, j, k;
    BOOL found;
    double sum;
    FEMListUtilities *listUtilities;
    FEMElementDescription *elementDescription;
    FEMMaterial *materialAtID = nil;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    listUtilities = [[FEMListUtilities alloc] init];
    
    [kernel getNodes:solution model:model inElement:element resultNodes:_elementNodes numberOfNodes:NULL];
    
    memset( _heatTransferCoeff, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
    memset( _load, 0.0, solution.mesh.maxElementDofs*sizeof(double) );
    
    // BC: -k@T/@n = \epsilon\sigma(T^4 - Text^4)
    _radiationFlag = [listUtilities listGetString:model inArray:bc forVariable:@"radiation" info:&found];
    if (found == YES && [_radiationFlag isEqualToString:@"none"] == NO) {
        found = [kernel getReal:model forElement:element inArray:bc variableName:@"emissivity" buffer:&buffer];
        if (found == YES) {
            memcpy(_nodalEmissivity, buffer.vector, n*sizeof(double));
        } else {
            found = [kernel getParentMaterialProperty:@"emissivity" forElement:element parentElement:parent model:model buffer:&buffer];
            if (found == YES) memcpy(_nodalEmissivity, buffer.vector, n*sizeof(double));
        }
        sum = 0.0;
        for (i=0; i<n; i++) {
            sum = sum + _nodalEmissivity[i];
        }
        _emissivity = sum / n;
        
        if ([_radiationFlag isEqualToString:@"idealized"] == YES) {
            found = [kernel getReal:model forElement:element inArray:bc variableName:@"radiation external temperature" buffer:&buffer];
            if (found == YES) {
                memcpy(_aText, buffer.vector, n*sizeof(double));
            } else {
                found = [kernel getReal:model forElement:element inArray:bc variableName:@"external temperature" buffer:&buffer];
                if (found == YES) memcpy(_aText, buffer.vector, n*sizeof(double));
            }
        } else {
            [self FEMHeatSolution_diffuseGrayRadiationModel:model solution:solution kernel:kernel element:element numberOfNodes:n forceVector:forceVector angleFraction:&_visibleFraction text:&_text];
            
            if ([listUtilities listGetLogical:model inArray:bc forVariable:@"radiation boundary open" info:&found] == YES) {
                found = [kernel getReal:model forElement:element inArray:bc variableName:@"radiation external temperature" buffer:&buffer];
                if (found == YES) {
                    memcpy(_aText, buffer.vector, n*sizeof(double));
                } else {
                    found = [kernel getReal:model forElement:element inArray:bc variableName:@"external temperature" buffer:&buffer];
                    if (found == YES) memcpy(_aText, buffer.vector, n*sizeof(double));
                }
                for (i=0; i<n; i++) {
                    _aText[i] = pow( ( (1.0 - _visibleFraction) * pow(_aText[i], 4.0) + _visibleFraction * pow(_text, 4.0) ), 0.25);
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
    
    found = [kernel getReal:model forElement:element inArray:bc variableName:@"heat transfert coefficient" buffer:&buffer];
    if (found == YES) {
        memcpy(_work, buffer.vector, n*sizeof(double));
        found = [kernel getReal:model forElement:element inArray:bc variableName:@"external temperature" buffer:&buffer];
        memcpy(_aText, buffer.vector, n*sizeof(double));
        for (j=0; j<n; j++) {
            // BC: -k@T/@n = \alpha(T - Text)
            k = _tempPerm[element->NodeIndexes[j]];
            _load[j] = _load[j] + _work[j] * _aText[j];
            _heatTransferCoeff[j] = _heatTransferCoeff[j] + _work[j];
        }
    }
    
    // BC: -k@T/@n = (rho*L)*v.n
    // Heating related to pulling is possible only in cases where pull velocity is described
    if ([listUtilities listGetLogical:model inArray:bc forVariable:@"phase change" info:&found] == YES) {
        elementDescription = [[FEMElementDescription alloc] init];
        
        found = [kernel getReal:model forElement:element inArray:bc variableName:@"phase change 1" buffer:&buffer];
        for (i=0; i<n; i++) {
            _phaseVelocity[0][i] = buffer.vector[i];
        }
        found = [kernel getReal:model forElement:element inArray:bc variableName:@"phase change 2" buffer:&buffer];
        for (i=0; i<n; i++) {
            _phaseVelocity[1][i] = buffer.vector[i];
        }
        found = [kernel getReal:model forElement:element inArray:bc variableName:@"phase change 3" buffer:&buffer];
        for (i=0; i<n; i++) {
            _phaseVelocity[2][i] = buffer.vector[i];
        }
        
        // Ensure that the latent heat and density come from the same side
        found = [kernel getParentMaterialProperty:@"latent heat" forElement:element parentElement:parent model:model buffer:&buffer];
        memcpy(_latentHeat, buffer.vector, n*sizeof(double));
        if (parent == NULL) {
            NSLog(@"FEMHeatSolution:FEMHeatSolution_addHeatFluxBC: parent not associated.\n");
        } else {
            k = [(model.bodies)[parent->BodyID-1][@"material"] intValue];
            materialAtID = (model.materials)[k-1];
            found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&buffer];
             memcpy(_density, buffer.vector, n*sizeof(double));
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
        [elementDescription deallocation];
    }
    
    // BC: -k@T/@n = g
    found = [kernel getReal:model forElement:element inArray:bc variableName:@"heat flux" buffer:&buffer];
    if (found == YES) {
        for (i=0; i<n; i++) {
            _load[i] = _load[i] + buffer.vector[i];
        }
    }
    
    // Get element matrix and rhs due to boundary conditions...
    if (model.coordinates == cartesian) {
        FEMDiffuseConvectiveAnisotropic *diffuseConvectiveAnisotropic = [[FEMDiffuseConvectiveAnisotropic alloc] init];
        Dimensions_t dimensions = (Dimensions_t){.mat1=(2*solution.mesh.maxElementDofs), .mat2=(2*solution.mesh.maxElementDofs), .vec=(2*solution.mesh.maxElementDofs)};
        [diffuseConvectiveAnisotropic diffuseConvectiveBoundaryMatrix:_stiff boundaryVector:_force dimensions:dimensions loadVector:_load nodalAlpha:_heatTransferCoeff element:element numberOfNodes:n nodes:_elementNodes mesh:solution.mesh];
    } else {
        FEMDiffuseConvectiveGeneralAnisotropic *diffuseConvectiveGeneralAnisotropic = [[FEMDiffuseConvectiveGeneralAnisotropic alloc] init];
        Dimensions_t dimensions = (Dimensions_t){.mat1=(2*solution.mesh.maxElementDofs), .mat2=(2*solution.mesh.maxElementDofs), .vec=(2*solution.mesh.maxElementDofs)};
        [diffuseConvectiveGeneralAnisotropic diffuseConvectiveGeneralBoundaryMatrix:_stiff boundaryVector:_force dimensions:dimensions loadVector:_load nodalAlpha:_heatTransferCoeff element:element numberOfNodes:n nodes:_elementNodes model:model mesh:solution.mesh];
    }
    
    // Update global matrices from local matrices
    int rows = 2*solution.mesh.maxElementDofs;
    int cols = 2*solution.mesh.maxElementDofs;
    if (_transientAssembly == YES && _constantBulk == NO) {
        memset( *_mass, 0.0, ((2*solution.mesh.maxElementDofs)*(2*solution.mesh.maxElementDofs))*sizeof(double) );
        [kernel defaultFirstOrderTime:model inSolution:solution forElement:element realMass:_mass realStiff:_stiff realForce:_force stiffRows:&rows stiffCols:&cols];
    }
    if (_heatGapBC == YES) {
        [self FEMHeatSolution_addHeatGapSolution:solution element:element numberOfNodes:n kernel:kernel];
    }
    [kernel defaultUpdateEquations:model inSolution:solution forElement:element realStiff:_stiff realForce:_force stiffRows:&rows stiffCols:&cols requestBulkUpdate:NULL];
    
    if (buffer.vector != NULL) free_dvector(buffer.vector, 0, buffer.m-1);
}

-(void)FEMHeatSolution_addGlobalTimeSolution:(FEMSolution *)solution {
    
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
            n1 = max(n1, (int)(matContainers->RHS[i+1]-matContainers->RHS[i]));
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
        for (j=matContainers->RHS[i]; j<=matContainers->RHS[i+1]-1; j++) {
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
        for (j=matContainers->RHS[i]; j<=matContainers->RHS[i+1]-1; j++) {
            matContainers->Values[j] = stiff[0][n];
            n++;
        }
        matContainers->RHS[i] = force[0];
    }
}

-(BOOL)FEMHeatSolution_checkLatentHeatModel:(FEMModel *)model {
    
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
        phaseChangeIntervals.matrix = NULL;
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
        _elementNodes->x = NULL;
        _elementNodes->y = NULL;
        _elementNodes->z = NULL;

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
    
    int i, j, k, l, n, nb, nd, t, bf_id, body_id, cols, compressibilityModel, eq_id, iter, mat_id, nsdofs, nonLinearIter, newtonIter, rows, smartHeaterNode, smartHeaterBC;
    int *indexes = NULL, *flowPerm;
    double at, at0, C1, dist, dt0, jx, jy, jz, cumulativeTime, newtonTol, meltPoint, minDist, nonLinearTol, norm, powerRelax, powerSensitivity, powerTimeScale,
           prevNorm, referencePressure, relax, relativeChange, s, saveRelax, smartTol, specificHeatRatio, st, sum, totat, totst, xave, yave;
    double controlPoint[3], *forceVector, *flowSolution;
    BOOL all, bubbles = YES, checkLatentHeatRelease=NO, found, heatFluxBC, gotIt, isRadiation, firstTime, gotMeltPoint,
         integralHeaterControl, heaterControlLocal, phaseChange=NO, saveBulk, smartHeaterControl, smartHeaterAverage, smartTolReached, stabilize = YES, transientHeaterControl, useBubbles;
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
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    FEMKernel *kernel = [FEMKernel sharedKernel];
    FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    
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
        
        cols = 2*n;
        rows = 2*n;
        
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
        
        smartHeaterAverage = NO;
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
            NSLog(@"FEMHeatSolution: found control point at distance: %f\n", sqrt(minDist));
            NSLog(@"FEMHeatSolution: control point index: %d\n", smartHeaterNode);
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
                jx = -DBL_MAX;
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
    saveBulk = (_constantBulk == YES || [(solution.solutionInfo)[@"save bulk system"] boolValue] == YES) ? YES : NO;
    saveBulk = (_constantBulk == YES || [(solution.solutionInfo)[@"calculate loads"] boolValue] == YES) ? YES : NO;
    
    saveRelax = relax;
    cumulativeTime = 0.0;
    
    firstTime = YES;
    _prevSolution = doublevec(0, _localNodes-1);
    
    FEMElementUtils *elementUtils = [[FEMElementUtils alloc] init];
    FEMDiffuseConvectiveAnisotropic *diffuseConvectiveAnisotropic = [[FEMDiffuseConvectiveAnisotropic alloc] init];
    FEMDiffuseConvectiveGeneralAnisotropic *diffuseConvectiveGeneralAnisotropic = [[FEMDiffuseConvectiveGeneralAnisotropic alloc] init];
    
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
            
            NSLog(@"FEMHeatSolution\n");
            NSLog(@"FEMHeatSolution\n");
            NSLog(@"FEMHeatSolution: --------------------------------------------------------\n");
            NSLog(@"FEMHeatSolution: TEMPERATURE ITERATION %d\n", iter);
            NSLog(@"FEMHeatSolution: --------------------------------------------------------\n");
            NSLog(@"FEMHeatSolution\n");
            NSLog(@"FEMHeatSolution: Starting Assembly...\n");
            
            if (_constantBulk == YES && matContainers->BulkValues != NULL) {
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
                    
                    n = element->Type.NumberOfNodes;
                    
                    mat_id = [kernel getMaterialIDForElement:element model:model];
                    materialAtID = (model.materials)[mat_id-1];
                    
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&buffer];
                    memcpy(_density, buffer.vector, n*sizeof(double));
                    
                    found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"heat source" buffer:&buffer];
                    memcpy(_load, buffer.vector, n*sizeof(double));
                    
                    _s = [elementUtils elementArea:element numberOfNodes:n mesh:solution.mesh nodel:model];
                    
                    if (model.coordinates == axis_symmetric || model.coordinates == cylindric_symmetric) _s = 2.0 * pi * _s;
                    
                    sum = 0.0;
                    for (i=0; i<n; i++) {
                        sum = sum + (_density[i] * _load[i]);
                    }
                    _heaterSource[bf_id-1] = _heaterSource[bf_id-1] + _s * sum / n;
                    _heaterArea[bf_id-1] = _heaterArea[bf_id-1] + _s;
                    sum = 0.0;
                    for (i=0; i<n; i++) {
                        sum = sum + _density[i];
                    }
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
                    
                    _phaseModel = [listUtilities listGetString:model inArray:equationAtID.valuesList forVariable:@"phase change model" info:&found];
                    if (found == NO) _phaseModel = [listUtilities listGetString:model inArray:materialAtID.valuesList forVariable:@"phase change model" info:&found];
                    phaseChange = (found == YES && [_phaseModel isEqualToString:@"none"] == NO) ? YES : NO;
                    if (phaseChange == YES) {
                        checkLatentHeatRelease = [listUtilities listGetLogical:model inArray:equationAtID.valuesList forVariable:@"check latent heat release" info:&found];
                    }
                }
                
                n = element->Type.NumberOfNodes;
                [kernel getNodes:solution model:model inElement:element resultNodes:_elementNodes numberOfNodes:NULL];
                [kernel getScalarLocalField:_localTemperature sizeField:solution.mesh.maxElementDofs name:nil element:element solution:solution model:model timeStep:NULL];
                
                // Get element material parameters
                found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"heat capacity" buffer:&buffer];
                 memcpy(_heatCapacity, buffer.vector, n*sizeof(double));
                
                found = [listUtilities listGetRealArray:model inArray:materialAtID.valuesList forVariable:@"heat conductivity" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer];
                memset( **_heatConductivity, 0.0, (3*3*solution.mesh.maxElementDofs)*sizeof(double) );
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
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"reference temperature" buffer:&buffer];
                    if (found == YES) memcpy(_referenceTemperature, buffer.vector, n*sizeof(double));
                    
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"heat expansion coefficient" buffer:&buffer];
                    if (found == YES) memcpy(_heatExpansionCoeff, buffer.vector, n*sizeof(double));
                    
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
                        found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&buffer];
                        if (found == YES) memcpy(_density, buffer.vector, n*sizeof(double));
                    }
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"pressure coefficient" buffer:&buffer];
                    if (found == YES) {
                        memcpy(_pressureCoeff, buffer.vector, n*sizeof(double));
                    } else memset( _pressureCoeff, 0.0, n*sizeof(double) );
                } else {
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"pressure coefficient" buffer:&buffer];
                    if (found == YES) {
                        memcpy(_pressureCoeff, buffer.vector, n*sizeof(double));
                    } else memset( _pressureCoeff, 0.0, n*sizeof(double) );
                    
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"density" buffer:&buffer];
                    if (found == YES) memcpy(_density, buffer.vector, n*sizeof(double));
                }
                
                // Take pressure deviation p_d as the dependent variable p = p_0 + p_d.
                // For perfect gas, read p_0
                if (compressibilityModel != incompressible) {
                    referencePressure = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"reference pressure" info:&found minValue:NULL maxValue:NULL];
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
                [kernel getVectorLocalField:_mu size1Field:3 size2Field:solution.mesh.maxElementDofs name:@"mesh velocity" element:element solution:solution model:model timeStep:NULL];
                
                if ([convectionFlag isEqualToString:@"constant"] == YES) {
                    
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"convection velocity 1" buffer:&buffer];
                    if (found == YES) {
                        memcpy(_u, buffer.vector, n*sizeof(double));
                    } else {
                        found = [kernel getReal:model forElement:element inArray:equationAtID.valuesList variableName:@"convection velocity 1" buffer:&buffer];
                        memcpy(_u, buffer.vector, n*sizeof(double));
                    }
                    
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"convection velocity 2" buffer:&buffer];
                    if (found == YES) {
                        memcpy(_v, buffer.vector, n*sizeof(double));
                    } else {
                        [kernel getReal:model forElement:element inArray:equationAtID.valuesList variableName:@"convection velocity 2" buffer:&buffer];
                        memcpy(_v, buffer.vector, n*sizeof(double));
                    }
                    
                    found = [kernel getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"convection velocity 3" buffer:&buffer];
                    if (found == YES) {
                        memcpy(_w, buffer.vector, n*sizeof(double));
                    } else {
                         found = [kernel getReal:model forElement:element inArray:equationAtID.valuesList variableName:@"convection velocity 3" buffer:&buffer];
                        memcpy(_w, buffer.vector, n*sizeof(double));
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
                    NSLog(@"FEMHeatSolution: convection model specified but no accociated flow field present?\n");
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
                    [self FEMHeatSolution_effectiveHeatCapacityElement:element numberOfNodes:n material:materialAtID model:model transientSimulation:transient];
                } else {
                    for (i=0; i<n; i++) {
                        _heatCapacity[i] = _density[i] * _heatCapacity[i];
                    }
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
                found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"perfusion rate" buffer:&buffer];
                if (found == YES) {
                    memcpy(_perfusionRate, buffer.vector, n*sizeof(double));
                    
                    found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"perfusion reference temperature" buffer:&buffer];
                    if (found == YES) memcpy(_perfusionRefTemperature, buffer.vector, n*sizeof(double));
                    found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"perfusion density" buffer:&buffer];
                    if (found == YES) memcpy(_perfusionDensity, buffer.vector, n*sizeof(double));
                    found = [kernel getReal:model forElement:element inArray:bodyForceAtID.valuesList variableName:@"perfusion heat capacity" buffer:&buffer];
                    if (found == YES) memcpy(_perfusionHeatCapacity, buffer.vector, n*sizeof(double));
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
                   [diffuseConvectiveAnisotropic diffuseConvectiveComposeMassMatrix:_mass
                                                                        stiffMatrix:_stiff
                                                                        forceVector:_force
                                                                         loadVector:_load
                                                                 timeDerivativeTerm:_heatCapacity
                                                                     zeroDegreeTerm:_c0
                                                                     convectionTerm:heatCapacity
                                                                      diffusionTerm:_heatConductivity
                                                                        phaseChange:_phaseSpatial
                                                                   nodalTemperature:_localTemperature
                                                                           enthalpy:_enthalpy
                                                                          velocityX:_u
                                                                           velocitY:_v
                                                                          velocityZ:_w
                                                                          meshVeloX:mux
                                                                          meshVeloY:muy
                                                                          meshVeloZ:muz
                                                                     nodalViscosity:_viscosity
                                                                       nodaldensity:_density
                                                                      nodalPressure:_pressure
                                                                    nodalPressureDt:_dPressureDt
                                                                 nodalPressureCoeff:_pressureCoeff
                                                                       compressible:(compressibilityModel != incompressible) ? YES : NO
                                                                          stabilize:stabilize
                                                                         useBubbles:useBubbles
                                                                            element:element
                                                                      numberOfNodes:n
                                                                              nodes:_elementNodes
                                                                           solution:solution
                                                                               mesh:solution.mesh
                                                                              model:model];
                } else {
                    double heatCapacity[n];
                    double mux[n], muy[n], muz[n];
                    for (i=0; i<n; i++) {
                        heatCapacity[i] = C1 * _heatCapacity[i];
                        mux[i] = _mu[0][i];
                        muy[i] = _mu[1][i];
                        muz[i] = _mu[2][i];
                    }
                    [diffuseConvectiveGeneralAnisotropic diffuseConvectiveGeneralComposeMassMatrix:_mass
                                                                                       stiffMatrix:_stiff
                                                                                       forceVector:_force
                                                                                        loadVector:_load
                                                                                timeDerivativeTerm:_heatCapacity
                                                                                    zeroDegreeTerm:_c0
                                                                                    convectionTerm:heatCapacity
                                                                                     diffusionTerm:_heatConductivity
                                                                                       phaseChange:_phaseSpatial
                                                                                  nodalTemperature:_localTemperature
                                                                                          enthalpy:_enthalpy
                                                                                         velocityX:_u
                                                                                          velocitY:_v
                                                                                         velocityZ:_w
                                                                                         meshVeloX:mux
                                                                                         meshVeloY:muy
                                                                                         meshVeloZ:muz
                                                                                    nodalViscosity:_viscosity
                                                                                      nodaldensity:_density
                                                                                     nodalPressure:_pressure
                                                                                   nodalPressureDt:_dPressureDt
                                                                                nodalPressureCoeff:_pressureCoeff
                                                                                      compressible:(compressibilityModel != incompressible) ? YES : NO
                                                                                         stabilize:stabilize
                                                                                           element:element
                                                                                     numberOfNodes:n
                                                                                             nodes:_elementNodes
                                                                                          solution:solution
                                                                                              mesh:solution.mesh
                                                                                             model:model];
                }
                
                if (heaterControlLocal == YES && transientHeaterControl == NO) {
                    if (_transientAssembly == YES && _constantBulk == NO) {
                        [kernel defaultFirstOrderTime:model inSolution:solution forElement:element realMass:_mass realStiff:_stiff realForce:_force stiffRows:&rows stiffCols:&cols];
                    }
                    if (indexes == NULL || n != nb) {
                        if (indexes != NULL) free_ivector(indexes, 0, nb-1);
                        indexes = intvec(0, n-1);
                        nb = n;
                    }
                    for (i=0; i<n; i++) {
                        indexes[i] = _tempPerm[element->NodeIndexes[i]];
                    }
                    [kernel updateGlobalEquationsModel:model inSolution:solution element:element localStiffMatrix:_stiff forceVector:_forceHeater localForce:_force size:n dofs:1 nodeIndexes:indexes rows:&rows cols:&cols rotateNT:NULL];
                } else {
                    bubbles = (useBubbles == YES && stabilize == NO && ([convectionFlag isEqualToString:@"computed"] == YES || [convectionFlag isEqualToString:@"constant"] == YES)) ? YES : NO;
                    
                    // If time dependent simulation, add mass matrix to stiff matrix
                    memset( _timeForce, 0.0, (2*solution.mesh.maxElementDofs)*sizeof(double) );
                    if (_transientAssembly == YES) {
                        if (_constantBulk == YES) {
                            [kernel defaultUpdateMass:_mass element:element solution:solution model:model];
                        } else {
                            [kernel defaultFirstOrderTime:model inSolution:solution forElement:element realMass:_mass realStiff:_stiff realForce:_force stiffRows:&rows stiffCols:&cols];
                        }
                    } else if (solution.nOfEigenValues > 0) {
                        [kernel defaultUpdateDamp:_mass element:element solution:solution model:model];
                    }
                    // Update global matrices and local matrices
                    if (bubbles == YES) {
                        [kernel condensateStiff:_stiff force:_force numberOfNodes:n force1:_timeForce];
                    }
                    [kernel defaultUpdateEquations:model inSolution:solution forElement:element realStiff:_stiff realForce:_force stiffRows:&rows stiffCols:&cols requestBulkUpdate:&saveBulk];
                }
            } // Bulk elements
            
        jump:
            // Neumann & Newton boundary conditions
            for (t=0; t<solution.mesh.numberOfBoundaryElements; t++) {
                 element = [kernel getBoundaryElement:solution atIndex:t];
                if ([kernel isActiveBoundaryElement:element inSolution:solution model:model] == NO) continue;
                
                n = element->Type.NumberOfNodes;
                if ([kernel getElementFamily:element] == 1) continue;
                
                bc = [kernel getBoundaryCondition:model forElement:&elements[t]];
                if (bc == nil) continue;
                
                // This check whether there are any Dirichlet conditions on the smart heater boundary.
                // If there are, the r.h.s must be zero as there can possibly not be any effect on temperature
                if (heaterControlLocal == YES && transientHeaterControl == NO) {
                    if ([listUtilities listCheckPresentVariable:solution.variable.name inArray:bc] == YES) {
                        memset( kernel.indexStore, -1, kernel.sizeIndexStore*sizeof(int) );
                        nd = [kernel getElementDofsSolution:solution model:model forElement:element atIndexes:kernel.indexStore];
                        for (i=0; i<nd; i++) {
                            _forceHeater[_tempPerm[kernel.indexStore[i]]] = 0.0;
                        }
                    }
                }
                
                heatFluxBC = [listUtilities listGetLogical:model inArray:bc forVariable:@"heat fluc bc" info:&found];
                if (found == YES && heatFluxBC == NO) continue;
                
                _heatGapBC = [listUtilities listGetLogical:model inArray:bc forVariable:@"heat gap" info:&found];
                [self FEMHeatSolution_addHeatFluxBC:bc element:element parent:parent numberOfNodes:n forceVector:forceVector kernel:kernel solution:solution model:model];
                
                if (_heatGapBC == YES) {
                    [self FEMHeatSolution_findGapIndexesElement:element indexes:_indexes numberOfNodes:n solution:solution];
                    memcpy(_saveIndexes, element->NodeIndexes, n*sizeof(double));
                    memcpy(element->NodeIndexes, _indexes, n*sizeof(double));
                    [self FEMHeatSolution_addHeatFluxBC:bc element:element parent:parent numberOfNodes:n forceVector:forceVector kernel:kernel solution:solution model:model];
                    memcpy(element->NodeIndexes, _saveIndexes, n*sizeof(double));
                }
            } // Neumann & Newton BCs
            
            if (transient == YES && _constantBulk == YES) [self FEMHeatSolution_addGlobalTimeSolution:solution];
            
            [kernel defaultFinishAssemblySolution:solution model:model];
            NSLog(@"FEMHeatSolution: Assembly done.\n");
            
            [kernel dirichletBoundaryConditions:model inSolution:solution usingOffset:NULL];
            
            // Solve the system and check for convergence
            at = cputime() -  at;
            st = cputime();
            
            prevNorm = norm;
            
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
                    
                    [kernel solveSystemMatrix:solution.matrix rhs:_forceHeater result:_xx norm:&norm dofs:1 solution:solution model:model];
                    [kernel solveSystemMatrix:solution.matrix rhs:matContainers->RHS result:_yy norm:&norm dofs:1 solution:solution model:model];
                    
                    [solution.solutionInfo setObject:@NO forKey:@"skip compute nonlinear change"];
                } else {
                    [kernel solveSystemMatrix:solution.matrix rhs:matContainers->RHS result:_temperature norm:&norm dofs:1 solution:solution model:model];
                }
                
                if (smartHeaterAverage == NO) {
                    xave = _xx[_tempPerm[smartHeaterNode]];
                    yave = _yy[_tempPerm[smartHeaterNode]];
                } else {
                    xave = 0.0;
                    yave = 0.0;
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
                    [listUtilities addConstRealInClassList:model.simulation theVariable:@"res: smart heater temperature" withValue:yave string:nil];
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
                    [kernel computeChange:solution model:model isSteadyState:NO nsize:&_localNodes values:_temperature values0:NULL];
                    norm = solution.variable.norm;
                }
                
                if (_dt > powerTimeScale) {
                    if (relax != 0.0) [solution.solutionInfo setValue:@(relax) forKey:@"nonlinear system relaxation factor"];
                }
            } else {
                norm = [kernel findSolution:solution model:model backRorateNT:NULL];
            }
            
            if (smartHeaterControl == YES || integralHeaterControl == YES) {
                [listUtilities addConstRealInClassList:model.simulation theVariable:@"res: heater power scaling" withValue:_powerScaling string:nil];
                NSLog(@"FEMHeatSolution: Heater control information\n");
                i = 0;
                for (FEMBodyForce *bodyForce in model.bodyForces) {
                    if (!(_smarterHeaters[i] || _integralHeaters[i])) continue;
                    if (_smarterHeaters[i]) _heaterScaling[i] = _powerScaling;
                    NSLog(@"FEMHeatSolution: heater for body %d\n", i+1);
                    if (_smarterHeaters[i]) NSLog(@"FEMHeatSolution: heater type: smart heater\n");
                    if (_integralHeaters[i]) NSLog(@"FEMHeatSolution: hater type: integral heater\n");
                    
                    NSLog(@"FEMHeatSolution: heater volutme (m^3): %f\n", _heaterArea[i]);
                    s = _heaterSource[i] * _heaterScaling[i];
                    NSLog(@"FEMHeatSolution: heater power (W): %f\n", s);
                    
                    NSLog(@"FEMHeatSolution: heater scaling: %f\n", _heaterScaling[i]);
                    NSLog(@"FEMHeatSolution: heater power density (W/kg): %e\n", s/(_heaterDensity[i]*_heaterArea[i]));
                    
                    if (_smarterHeaters[i]) [listUtilities addConstRealInClassList:model.simulation theVariable:@"res: heater power density" withValue:s/(_heaterDensity[i]*_heaterArea[i]) string:nil];
                    i++;
                }
            }
            
            st = cputime() - st;
            totat = totat + at;
            totst = totst + st;
            NSLog(@"FEMHeatSolution: iter: %d, Assembly (s): %f %f\n", iter, at, totat);
            NSLog(@"FEMHeatSolution: iter: %d, Solve (s): %f %f\n", iter, st, totst);
            
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
                        NSLog(@"FEMHeatSolution: latent heat release check: reducing time step to: %f\n", _dt);
                    } else {
                        relax = relax / 2;
                        [solution.solutionInfo setValue:@(relax) forKey:@"nonlinear system relaxation factor"];
                        NSLog(@"FEMHeatSolution: latent heat release check: reducing relaxation to: %f\n", relax);
                    }
                    continue;
                }
                if (transient == NO) memcpy(_prevSolution, _temperature, _localNodes*sizeof(double));
            }
            
            relativeChange = solution.variable.nonLinChange;
            NSLog(@"FEMHeatSolution: result norm: %f\n", norm);
            NSLog(@"FEMHeatSolution: relative change: %f\n", relativeChange);
            
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
        if (transient == NO) {
            cumulativeTime = cumulativeTime + _dt;
            _dt = timeStep - cumulativeTime;
        }
        
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
    
    solution.dt = timeStep;
    [solution.solutionInfo setValue:@(saveRelax) forKey:@"nonlinear system relaxation factor"];
    free_dvector(_prevSolution, 0, _localNodes-1);
    
    if ([(solution.solutionInfo)[@"adaptive mesh refinement"] boolValue] == YES) {
        // TODO: implement mesh refinement
    }
}

@end
