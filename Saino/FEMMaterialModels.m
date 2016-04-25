//
//  FEMMaterialModels.m
//  Saino
//
//  Created by Seddik hakime on 13/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMaterialModels.h"
#import "FEMCore.h"
#import "FEMListUtilities.h"
#import "FEMCoordinateSystems.h"
#import "FEMLinearAlgebra.h"
#include "FEMElementDescription.h"
#include "Utils.h"

@implementation FEMMaterialModels

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
    }
    
    return self;
}

/******************************************************************************
 
    Return second invariant.
    Note: Actually square of the second invariant of the velocity is returned
 
******************************************************************************/
-(double)secondInvariantVelo:(double[3])velo dVelodx:(double[][3])dVelodx crtMatrix:(double[][3])crtMatrix symbols:(double[][3][3])symbols model:(FEMModel* __nonnull)model {
    
    int i, j, k, l;
    double s, t;
    double secInv;
    
    secInv = 0.0;
    if (model.coordinates == cartesian) {
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) {
                s = dVelodx[i][j] + dVelodx[j][i];
                secInv = secInv + s * s;
            }
        }
    } else if (model.coordinates == axis_symmetric) {
        secInv = pow((2.0*dVelodx[0][0]), 2.0) + pow((2.0*dVelodx[1][1]), 2.0) + 2.0*pow((dVelodx[0][1] + dVelodx[1][0]), 2.0)
                 + pow((2.0*velo[0]*symbols[0][2][2]), 2.0);
    } else {
        FEMLinearAlgebra *linearAlgebra = [[FEMLinearAlgebra alloc] init];
        double **covMetric = doublematrix(0, 2, 0, 2);
        memcpy(*covMetric, *crtMatrix, (3*3)*sizeof(double));
        [linearAlgebra invertMatrix:covMetric ofSize:3];
        
        for (i=0; i<3; i++) {
            for (j=0; j<3; j++) {
                s = 0.0;
                t = 0.0;
                for (k=0; k<3; k++) {
                    s = s + covMetric[i][k] * dVelodx[k][j] + covMetric[j][k] * dVelodx[k][i];
                    t = t + crtMatrix[j][k] *dVelodx[i][k] + crtMatrix[i][k] * dVelodx[j][k];
                    for (l=0; l<3; l++) {
                        s = s - covMetric[i][k] * symbols[l][j][k] * velo[l];
                        s = s - covMetric[j][k] * symbols[l][i][k] * velo[l];
                        
                        t = t - crtMatrix[j][k] * symbols[l][k][i] * velo[l];
                        t = t - crtMatrix[i][k] * symbols[l][k][j] * velo[l];
                    }
                }
                secInv = secInv + s * t;
            }
        }
        free_dmatrix(covMetric, 0, 2, 0, 2);
    }
    
    return secInv;
}

/*************************************************************************************************
 
    Returns effective viscosity for Navier-Stokes equation.
    The viscosity model may be either some non-newtonian material law or from turbulence models, 
    but not from both at the same time.
 
*************************************************************************************************/
-(double)effectiveViscosity:(double)viscosity density:(double)density velocityX:(double * __nonnull)ux velocitY:(double * __nonnull)uy velocityZ:(double * __nonnull)uz element:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int)n numberOfPoints:(int)nd integrationU:(double)u integrationV:(double)v integrationW:(double)w muder:(double * __nullable)muder mesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration {
    
    int i, j, k;
    double a1, a2, arrheniusFactor, c1, c2, c3, c4, cmu, ct, dVelodx[3][3], dSymb[3][3][3][3], ehf, h, ke_k, ke_e, ke_z, metric[3][3], mu, q1, q2,
           sqrtElementMetric, sqrtMetric, r, s, ss, symb[3][3][3], temp, timeScale, tLimit, velo[3], x, y, z;
    NSString *viscosityFlag;
    BOOL found, setArrheniusFactor, stat;
    FEMMaterial *materialAtID = nil;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    mu = viscosity;
    if (muder != NULL) *muder = 0;
    
    k = [(model.bodies)[element->BodyID-1][@"material"] intValue];
    materialAtID = (model.materials)[k-1];
    
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    viscosityFlag = [listUtilities listGetString:model inArray:materialAtID.valuesList forVariable:@"viscosity model" info:&found];
    if (found == NO) return mu;
    
    // Basis function values & derivatives at the calculation point
    stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
    stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
    stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
    sqrtElementMetric = integration.metricDeterminant;
    // Coordinate system dependent information
    x = cblas_ddot(n, nodes->x, 1, integration.basis, 1);
    y = cblas_ddot(n, nodes->y, 1, integration.basis, 1);
    z = cblas_ddot(n, nodes->z, 1, integration.basis, 1);
    FEMCoordinateSystems *coordinateSystems = [[FEMCoordinateSystems alloc] init];
    [coordinateSystems coordinateSystemInfoModel:model metric:metric sqrtMetric:&sqrtMetric symbols:symb dSymbols:dSymb coordX:x coordY:y coordZ:z];
    
    memset( *dVelodx, 0.0, (3*3)*sizeof(double) );
    for (j=0; j<3; j++) {
        for (i=0; i<nd; i++) {
            dVelodx[0][j] = dVelodx[0][j] + ux[i]*integration.basisFirstDerivative[i][j];
            dVelodx[1][j] = dVelodx[1][j] + uy[i]*integration.basisFirstDerivative[i][j];
            dVelodx[2][j] = dVelodx[2][j] + uz[i]*integration.basisFirstDerivative[i][j];
        }
    }
    
    memset( velo, 0.0, sizeof(velo) );
    velo[0] = cblas_ddot(nd, integration.basis, 1, ux, 1);
    velo[1] = cblas_ddot(nd, integration.basis, 1, uy, 1);
    velo[2] = cblas_ddot(nd, integration.basis, 1, uz, 1);
    
    // This is the square of shear rate which results to 1/2 in exponent
    // Also the derivative is taken with respect to the square
    ss = [self secondInvariantVelo:velo dVelodx:dVelodx crtMatrix:metric symbols:symb model:model] / 2.0;
    
    if ([viscosityFlag isEqualToString:@"glen"] == YES) {
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"glen exponent" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL]; // This is the real exponent, n, not 1/n
        if (found == YES) {
            c2 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        } else {
            c2 = 0.0;
            for (i=0; i<n; i++) {
                c2 = c2 + integration.basis[i]*3.0;
            }
        }
        
        s = ss / 4.0; // The second invariant is not taken from the strain rate tensor, but rather 2*strain rate tensor (that's why we divide by 4 = 2^2)
        
        setArrheniusFactor = [listUtilities listGetLogical:model inArray:materialAtID.valuesList forVariable:@"set arrhenius factor" info:&found];
        if (found == NO || setArrheniusFactor == NO) {
            found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"constant temperature" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
            temp = 0.0;
            if (found == NO) { // Need to find a temperature field
                NSString *temperatureName = [listUtilities listGetString:model inArray:materialAtID.valuesList forVariable:@"temperature field variable" info:&found];
                if (found == NO) temperatureName = @"temperature";
                FEMUtilities *utilities = [[FEMUtilities alloc] init];
                FEMVariable *tempSol = [utilities getVariableFrom:model.variables model:model name:temperatureName onlySearch:NULL maskName:nil info:&found];
                if (found == YES) {
                    variableArraysContainer *tempContainers = tempSol.getContainers;
                    for (i=0; i<n; i++) {
                        temp = temp + integration.basis[i] * tempContainers->Values[tempContainers->Perm[element->NodeIndexes[i]]];
                    }
                } else {
                    fprintf(stderr, "FEMMaterialModels:effectiveViscosity: can't find variable %s to inquire temperature field for Glen.\n", [temperatureName UTF8String]);
                    fatal("FEMMaterialModels:effectiveViscosity", "Saino will abort the simulation now...");
                }
            } else {
                temp = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
            }
            r = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"gas constant" info:&found minValue:NULL maxValue:NULL];
            if (found == NO) r = 8.314;
            
            // Hard coded limit for the time being
            tLimit = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"limit temperature" info:&found minValue:NULL maxValue:NULL];
            if (found == NO) {
                tLimit = -10.0;
                fprintf(stdout, "FEMMaterialModels:effectiveViscosity: limit temperature not found. Setting -10.0.\n");
            }
            
            a1 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"rate factor 1" info:&found minValue:NULL maxValue:NULL];
            if (found == NO) {
                a1 = 3.985e-13;
                fprintf(stdout, "FEMMaterialModels:effectiveViscosity: rate factor 1 not found. Setting 3.985e-13.\n");
            }
            
            a2 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"rate factor 2" info:&found minValue:NULL maxValue:NULL];
            if (found == NO) {
                a2 = 1.916e03;
                fprintf(stdout, "FEMMaterialModels:effectiveViscosity: rate factor 2 not found. Setting 1.916e03.\n");
            }
            
            q1 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"activation energy 1" info:&found minValue:NULL maxValue:NULL];
            if (found == NO) {
                a1 = 60.0e03;
                fprintf(stdout, "FEMMaterialModels:effectiveViscosity: activation energy 1 not found. Setting 60.0e03.\n");
            }
            
            q2 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"activation energy 2" info:&found minValue:NULL maxValue:NULL];
            if (found == NO) {
                a2 = 139.0e03;
                fprintf(stdout, "FEMMaterialModels:effectiveViscosity: activation energy 2 not found. Setting 139.0e03.\n");
            }

            if (temp <= tLimit) {
                arrheniusFactor = a1 * exp( -q1/(r * (273.15 + temp)) );
            } else if (tLimit < temp && temp <= 0.0) {
                arrheniusFactor = a2 * exp( -q2/(r * (273.15 + temp)) );
            } else {
                arrheniusFactor = a2 * exp( -q2/(r * 273.15) );
                fprintf(stdout, "FEMMaterialModels:effectiveViscosity: positive temperature detected in Glen - limiting to zero.\n");
            }
            
        } else {
            arrheniusFactor = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"arrhenius factor" info:&found minValue:NULL maxValue:NULL];
            if (found == NO) {
                fatal("FEMMaterialModels:effectiveViscosity", "'set arrhenius factor' parameter is true but not value for 'arrhenius factor'.");
            }
        }
        
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"glen enhancement factor" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        if (found == YES) {
            ehf = cblas_ddot(n, integration.basis, 1,  buffer.vector, 1);
        } else {
            ehf = 0.0;
            for (i=0; i<n; i++) {
                ehf = ehf + integration.basis[i] * 1.0;
            }
        }
        
        if (muder != NULL) *muder = 0.5 * pow(ehf*arrheniusFactor, -1.0/c2) * ( (1.0/c2)-1.0 )/2.0 * pow(s, ((1.0/c2)-1.0)/2.0 - 1.0)/4.0;
        
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"critical shear rate" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        if (found == YES) {
            c3 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
            if (s < pow(c3, 2.0)) {
                s = pow(c3, 2.0);
                if (muder != NULL) *muder = 0.0;
            }
        }
        // Compute the effective viscosity
        mu = 0.5 * pow(ehf * arrheniusFactor, -1.0/c2) * pow(s, ((1.0/c2)-1.0)/2.0 );
        
    } else if ([viscosityFlag isEqualToString:@"power law"] == YES) {
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity exponent" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c2 = 0.0;
        if (found == YES) {
            c2 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        
        s = ss;
        if (muder != NULL) *muder = viscosity * (c2-1.0)/2.0 * pow(s, (c2-1.0)/2.0 - 1.0);
        
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"critical shear rate" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        if (found == YES) {
            c3 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
            if (s < pow(c3, 2.0)) {
                s = pow(c3, 2.0);
                if (muder != NULL) *muder = 0.0;
            }
        }
        mu = viscosity * pow(s, (c2-1.0)/2.0);
        
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"nominal shear rate" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        if (found == YES) {
            c4 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
            mu = mu / pow(c4, c2-1.0);
            if (muder != NULL) *muder = *muder / pow(c4, c2-1.0);
        }
        
    } else if ([viscosityFlag isEqualToString:@"power law too"] == YES) {
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity exponent" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c2 = 0.0;
        if (found == YES) {
            c2 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        mu = pow(viscosity, -1.0/c2) * pow(ss, -(c2-1.0)/(2.0*c2)) / 2.0;
        if (muder != NULL) *muder = pow(viscosity, -1.0/c2) * ( -(c2-1.0)/(2.0*c2) ) * pow(ss, -(c2-1.0)/(2.0*c2)-1.0) / 2.0;
        
    } else if ([viscosityFlag isEqualToString:@"carreau"] == YES) {
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity difference" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c1 = 0.0;
        if (found == YES) {
            c1 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity exponent" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c2 = 0.0;
        if (found == YES) {
            c2 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity transition" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c3 = 0.0;
        if (found == YES) {
            c3 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        c4 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"yasuda exponent" info:&found minValue:NULL maxValue:NULL];
        if (found == YES) {
            s = sqrt(ss);
            mu = viscosity + c1 * pow(1.0 + pow(c3, c4) * pow(ss, c4/2.0), (c2-1.0)/c4);
            if (muder != NULL) *muder = c1 * pow(1.0 + pow(c3, c4) * pow(ss, c4/2.0), (c2-1.0)/c4-1.0) * (c2-1.0)/2.0 * pow(c3, c4)*pow(ss, c4/2.0-1.0);
        } else {
            mu = viscosity + c1 * pow(1.0 + c3*c3*ss, (c2-1.0)/2.0);
            if (muder != NULL) *muder = c1 * (c2-1.0)/2.0 * pow(c3, 2.0) * pow(1.0 + pow(c3, 2.0) * ss, (c2-1.0)/2.0-1.0);
        }
    } else if ([viscosityFlag isEqualToString:@"cross"] == YES) {
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity difference" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c1 = 0.0;
        if (found == YES) {
            c1 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity exponent" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c2 = 0.0;
        if (found == YES) {
            c2 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity transition" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c3 = 0.0;
        if (found == YES) {
            c3 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        mu = viscosity + c1 / (1.0 + c3 * pow(ss, c2/2.0));
        if (muder != NULL) *muder = -c1 * c3 * pow(ss, c2/2.0) * c2 / (2.0*pow(1.0+c3*pow(ss, c2/2.0), 2.0) * ss);
    } else if ([viscosityFlag isEqualToString:@"powell eyring"] == YES) {
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity difference" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c1 = 0.0;
        if (found == YES) {
            c1 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        c2 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"viscosity transition" info:&found minValue:NULL maxValue:NULL];
        s = sqrt(ss);
        if (c2*s < 1.0e-5) {
            mu = viscosity + c1;
        } else {
            mu = viscosity + c1 * log(c2 * s + sqrt(c2*c2*ss+1.0)) / (c2*s);
            if (muder != NULL) *muder = c1 * (c2/(2.0*s)+pow(c2, 2.0)/(2.0*sqrt(pow(c2, 2.0)*ss+1.0)))/((c2*s+sqrt(c2*ss+1.0))*c2*s) -
                                        c1 * log(c2*s+sqrt(pow(c2, 2.0)*ss+1.0))/(c2*pow(s, 3.0))/2.0;
        }
    } else if ([viscosityFlag isEqualToString:@"smagorinsky"] == YES) {
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"smagorinsky constant" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c2 = 0.0;
        if (found == YES) {
            c2 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
        h = [elementDescription elementDiameter:element nodes:nodes];
        viscosity = viscosity + density * c2 * pow(h, 2.0) * sqrt(2.0*ss) / 2.0;
        if (muder != NULL) *muder = density * c2 * pow(h, 2.0) * sqrt(2.0)/(4.0*sqrt(ss));
    } else if ([viscosityFlag isEqualToString:@"ke"] == YES || [viscosityFlag isEqualToString:@"k-epsilon"] == YES) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *var = [utilities getVariableFrom:model.variables model:model name:@"kinetic energy" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The kinetic energy variable is not defined.");
        variableArraysContainer *varContainers = var.getContainers;
        ke_k = 0.0;
        for (i=0; i<n; i++) {
            ke_k = ke_k + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        var = nil;
        varContainers = NULL;
        var = [utilities getVariableFrom:model.variables model:model name:@"kinetic dissipation" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The kinetic dissipation variable is not defined.");
        varContainers = var.getContainers;
        ke_e = 0.0;
        for (i=0; i<n; i++) {
            ke_e = ke_e + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        if (![[listUtilities listGetString:model inArray:materialAtID.valuesList forVariable:@"ke model" info:&found] isEqualToString:@"v2-f"]) {
            found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"ke cmu" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
            if (found == YES) {
                cmu = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
            } else {
                cmu = 0.09;
            }
            mu = viscosity + cmu  * density * pow(ke_k, 2.0) / ke_e;
        } else {
            var = nil;
            varContainers = NULL;
            var = [utilities getVariableFrom:model.variables model:model name:@"v2" onlySearch:NULL maskName:nil info:&found];
            if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The v2 variable is not defined.");
            varContainers = var.getContainers;
            ke_z = 0.0;
            for (i=0; i<n; i++) {
                ke_z = ke_z + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
            }
            
            found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"v2-f ct" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
            ct = 0.0;
            if (found == YES) {
                ct = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
            }
            timeScale = max(ke_k/ke_e, ct*sqrt(viscosity/density/ke_e));
            
            found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"ke cmu" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
            cmu = 0.0;
            if (found == YES) {
                cmu = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
            }
            mu = viscosity + cmu*density*ke_z*timeScale;
        }
    } else if ([viscosityFlag isEqualToString:@"rng k-epsilon"] == YES) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *var = [utilities getVariableFrom:model.variables model:model name:@"effective viscosity" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The effective viscosity variable is not defined.");
        variableArraysContainer *varContainers = var.getContainers;
        for (i=0; i<n; i++) {
            mu = mu + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
    } else if ([viscosityFlag isEqualToString:@"spalart-allmaras"] == YES) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *var = [utilities getVariableFrom:model.variables model:model name:@"turbulent viscosity" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The turbulent viscosity variable is not defined.");
        variableArraysContainer *varContainers = var.getContainers;
        for (i=0; i<n; i++) {
            mu = mu + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        c1 = mu / (viscosity/density);
        c1 = pow(c1, 3.0) / (pow(c1, 3.0)+pow(7.1, 3.0));
        mu = viscosity + mu*density*c1;
    } else if ([viscosityFlag isEqualToString:@"k-omega"] == YES) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *var = [utilities getVariableFrom:model.variables model:model name:@"kinetic energy" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The kinetic energy variable is not defined.");
        variableArraysContainer *varContainers = var.getContainers;
        ke_k = 0.0;
        for (i=0; i<n; i++) {
            ke_k = ke_k + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        var = nil;
        varContainers = NULL;
        var = [utilities getVariableFrom:model.variables model:model name:@"kinetic dissipation" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The kinetic dissipation variable is not defined.");
        varContainers = var.getContainers;
        ke_e = 0.0;
        for (i=0; i<n; i++) {
            ke_e = ke_e + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        mu = viscosity + density * ke_k/ke_e;
    } else if ([viscosityFlag isEqualToString:@"sst k-omega"] == YES) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *var = [utilities getVariableFrom:model.variables model:model name:@"kinetic energy" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The kinetic energy variable is not defined.");
        variableArraysContainer *varContainers = var.getContainers;
        ke_k = 0.0;
        for (i=0; i<n; i++) {
            ke_k = ke_k + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        var = nil;
        varContainers = NULL;
        var = [utilities getVariableFrom:model.variables model:model name:@"kinetic dissipation" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The kinetic dissipation variable is not defined.");
        varContainers = var.getContainers;
        ke_e = 0.0;
        for (i=0; i<n; i++) {
            ke_e = ke_e + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        var = nil;
        varContainers = NULL;
        var = [utilities getVariableFrom:model.variables model:model name:@"wall distance" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "The wall distance variable is not defined.");
        varContainers = var.getContainers;
        double dist = 0.0;
        for (i=0; i<n; i++) {
            dist = dist + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        double f2 = tanh(pow(max(2.0+sqrt(ke_k)/(0.09*ke_e*dist), 500.0*viscosity/(density*ke_e*pow(dist, 2.0))), 2.0));
        double f3 = 1.0;
        mu = viscosity+0.31*density*ke_k/max(0.31*ke_e, sqrt(ss)*f2*f3);
    } else if ([viscosityFlag isEqualToString:@"levelset"] == YES) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *var = [utilities getVariableFrom:model.variables model:model name:@"surface" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "Variable surface needed in levelset viscosity model but not found.");
        variableArraysContainer *varContainers = var.getContainers;
        
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"viscosity difference" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        c1 = 0.0;
        if (found == YES) {
            c1 = cblas_ddot(n, integration.basis, 1, buffer.vector, 1);
        }
        c2 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"levelset bandwidth" info:&found minValue:NULL maxValue:NULL];
        double temp = 0.0;
        for (i=0; i<n; i++) {
            temp = temp + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        temp = temp / c2;
        if (temp < -1.0) {
            c3 = 0.0;
        } else if (temp > 1.0) {
            c3 = 1.0;
        } else {
            c3 = 0.75 * (temp - pow(temp, 3)/3.0) + 0.5;
        }
        mu = viscosity + c3 * c1;
    } else if ([viscosityFlag isEqualToString:@"user function"] == YES) {
        // TODO: implement this later if we need it
    } else {
        fprintf(stdout, "FEMMaterialModels:effectiveViscosity: unknown material model.\n");
    }
    
    if (buffer.vector != NULL) {
        free_dvector(buffer.vector, 0, buffer.m-1);
        buffer.vector = NULL;
    }
    
    // Add a generic temperature coefficient at the integration point for backward
    // compatibility, this is activated by an existing keyword
    c1 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"viscosity temp exp" info:&found minValue:NULL maxValue:NULL];
    if (found == YES) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *var = [utilities getVariableFrom:model.variables model:model name:@"temperature" onlySearch:NULL maskName:nil info:&found];
        if (var == nil) fatal("FEMMaterialModels:effectiveViscosity", "Variable temperature needed for thermal model but not found.");
        variableArraysContainer *varContainers = var.getContainers;
        c2 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"viscosity temp ref" info:&found minValue:NULL maxValue:NULL];
        c3 = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"viscosity temp offset" info:&found minValue:NULL maxValue:NULL];
        double temp = 0.0;
        for (i=0; i<n; i++) {
            temp = temp + integration.basis[i]*varContainers->Values[varContainers->Perm[element->NodeIndexes[i]]];
        }
        double tempCoeff = exp(c1*(1.0/(temp+c3)-1.0/c2));
        mu = tempCoeff * mu;
    }
    
    return mu;
}

/***************************************************************************
 
    Returns effective heat conductivity mainly related to turbulent models.
 
***************************************************************************/
-(double)effectiveConductivity:(double)conductivity density:(double)density element:(Element_t * __nonnull)element temperature:(double * __nullable)temperature velocityX:(double * __nonnull)ux velocitY:(double * __nonnull)uy velocityZ:(double * __nonnull)uz nodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int)n numberOfPoints:(int)nd integrationU:(double)u integrationV:(double)v integrationW:(double)w conductivityFlag:(NSString * __nonnull)conductivityFlag core:(FEMCore * __nonnull)core mesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    int i;
    double c_p, mu, pCond, pr_t, tmu;
    BOOL found, stat;
    FEMMaterial *materialAtID = nil;
    listBuffer c1n = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    pCond = conductivity;
    
    if ([conductivityFlag isEqualToString:@"ke"] == YES || [conductivityFlag isEqualToString:@"k-epsilon"] == YES || [conductivityFlag isEqualToString:@"turbulent"] == YES) {
        
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
     
        found = [core getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"heat capacity" buffer:&c1n listUtilities:listUtilities];
        c_p = 0.0;
        for (i=0; i<n; i++) {
            c_p = c_p + integration.basis[i]*c1n.vector[i];
        }
        
        found = [core getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"viscosity" buffer:&c1n listUtilities:listUtilities];
        mu = 0.0;
        for (i=0; i<n; i++) {
            mu = mu + integration.basis[i]*c1n.vector[i];
        }
        
        found = [core getReal:model forElement:element inArray:materialAtID.valuesList variableName:@"turbulent prandtl number" buffer:&c1n listUtilities:listUtilities];
        if (found == YES) {
            pr_t = 0.0;
            for (i=0; i<n; i++) {
                pr_t = pr_t + integration.basis[i]*c1n.vector[i];
            }
        } else pr_t = 0.85;
        
        tmu = [self effectiveViscosity:mu density:density velocityX:ux velocitY:uy velocityZ:uz element:element nodes:nodes numberOfNodes:n numberOfPoints:nd integrationU:u integrationV:v integrationW:w muder:NULL mesh:mesh model:model integration:integration] - mu;
        pCond = conductivity + c_p * tmu/pr_t;
        
        if (c1n.vector != NULL) {
            free_dvector(c1n.vector, 0, c1n.m-1);
            c1n.vector = NULL;
        }
    } else if ([conductivityFlag isEqualToString:@"user function"] == YES) {
        // TODO: implement this later if we need it
        // In that case change the method argument temperature from __nullable to __nonnull
    } else {
        fprintf(stdout, "FEMMaterialModels:effectiveConductivity: unknown material model.\n");
    }
    
    return pCond;
}

@end
