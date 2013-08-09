//
//  FEMDiffuseConvectiveGeneralAnisotropic.m
//  Saino
//
//  Created by Seddik hakime on 04/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMDiffuseConvectiveGeneralAnisotropic.h"
#import "FEMListUtilities.h"
#import "FEMNumericIntegration.h"
#import "FEMMaterial.h"
#import "FEMBodyForce.h"
#import "FEMMaterialModels.h"
#import "GaussIntegration.h"
#import "FEMDifferentials.h"
#import "FEMCoordinateSystems.h"
#import "Utils.h"

@implementation FEMDiffuseConvectiveGeneralAnisotropic

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
    }
    
    return self;
}

/********************************************************************************************************************************************
    Return element local matrices and RHS vector for diffusion-convection equation (general euclidian coordinates)
 
    Arguments:
 
        double **massMatrix             -> output: time derivative coefficient matrix
        double **stiffMatrix            -> output: rest of the equation coefficients
        double *forceVector             -> output: RHS vector
        double *loadVector              -> load vector
        double NodalCT,NodalC0,NodalC1  -> coefficient of the time derivative term, 0 degree term, and the convection term respectively
        double ***nodalC2               -> nodal values of the diffusion term coefficient tensor
        BOOL phaseChange                -> if model phase change
        double *modalTemperature        -> nodalTemperature from previous iteration
        double *enthalpy                -> enthalpy from previous iteration, needed if we model phase change
        double *ux, *uy, *uz            -> nodal values of velocity components from previous iteration used only if coefficient of
                                           the convection term (C1) is nonzero
        double *nodalViscosity          -> nodal viscosity
        double *nodaldensity            -> nodal density
        double *nodalPressure           -> nodal pressure
        double *nodalPressureDt         -> nodal pressure time derivative
        double *nodalPressureCoeff      -> nodal pressure coefficients
        BOOL compressible               -> if compressible flow
        BOOL stabilize                  -> should stabilzation be used ? Used only if coefficient of the convection term (C1) is nonzero
        Element_t *element              -> structure describing the element (dimension,nof nodes, interpolation degree, etc...)
        int numberOfNodes               -> number of element nodes
        Nodes_t *nodes                  -> element node coordinates
********************************************************************************************************************************************/
-(void)diffuseConvectiveGeneralComposeMassMatrix:(double **)massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double *)loadVector timeDerivativeTerm:(double *)nodalCT zeroDegreeTerm:(double *)nodalC0 convectionTerm:(double *)nodalC1 diffusionTerm:(double ***)nodalC2 phaseChange:(BOOL)phaseChange nodalTemperature:(double *)nodalTemperature enthalpy:(double *)enthalpy velocityX:(double *)ux velocitY:(double *)uy velocityZ:(double *)uz meshVeloX:(double *)mux meshVeloY:(double *)muy meshVeloZ:(double *)muz nodalViscosity:(double *)nodalviscosity nodaldensity:(double *)nodalDensity nodalPressure:(double *)nodalPressure nodalPressureDt:(double *)nodalPressureDt nodalPressureCoeff:(double *)nodalPressureCoeff compressible:(BOOL)compressible stabilize:(BOOL)stabilize element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution mesh:(FEMMesh *)mesh model:(FEMModel *)model {
    
    int i, j, k, l, p, q, t, dim, body_id, nBasis;
    double a, c0, c1, ct, dEnth, dTemp, force, hk, mk, dc2dx[3][3][3], density, detJ, divVelo, dNodalBasisdx[n][n][3], dsymb[3][3][3][3],
           dVelodx[3][3], load, m, metric[3][3], pe, pressure, s, su[n], symb[3][3][3], sw[n], tau, sqrtMetric, sum, u, v, velo[3], viscosity,
           vnorm, x, y, z, w;
    BOOL any, bubbles, convection, convectiveAndStabilize, cylindrincSymmetry, found, frictionHeat, stat;
    FEMKernel *kernel;
    FEMListUtilities *listUtilities;
    FEMNumericIntegration *integration;
    FEMBodyForce *bodyForceAtID = nil;
    GaussIntegrationPoints *IP;
    
    cylindrincSymmetry = (model.coordinates == cylindric_symmetric || model.coordinates == axis_symmetric) ? YES : NO;
    
    if (cylindrincSymmetry) {
        dim = 3;
    } else {
        dim = model.dimension;
    }
    load = 0.0;
    
    any = NO;
    for (i=0; i<n; i++) {
        if (nodalC1[i] != 0.0) {
            any = YES;
            break;
        }
    }
    convection = (any == YES) ? YES : NO;
    nBasis = n;
    bubbles = NO;
    if (convection == YES && stabilize == NO) {
        nBasis = 2*n;
        bubbles = YES;
    }
    
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) errorfunct("FEMDiffuseConvectiveGeneralAnisotropic:diffuseConvectiveGeneralComposeyMassMatrix", "Allocation error in FEMNumericIntegration!");
    // Integration stuff
    if (bubbles == YES) {
        IP = GaussQuadrature(element, NULL, &element->Type.GaussPoints2);
    } else IP = GaussQuadrature(element, NULL, NULL);

    // Stabilization parameters: hk, mk (Franca et al.)
    // If there is no convection term, we don't need stabilization
    convectiveAndStabilize = NO;
    if (stabilize == YES && any == YES) {
        convectiveAndStabilize = YES;
        hk = element->hK;
        mk = element->StabilizationMK;
        memset( **dNodalBasisdx, 0.0, (n*n*3)*sizeof(double) );
        for (p=0; p<n; p++) {
            u = element->Type.NodeU[p];
            v = element->Type.NodeV[p];
            w = element->Type.NodeW[p];
            stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
            for (i=0; i<n; i++) {
                for (j=0; j<3; j++) {
                    dNodalBasisdx[i][p][j] = integration.basisFirstDerivative[i][j];
                }
            }
        }
    }
    
    frictionHeat = NO;
    body_id = [kernel getBodyForceIDForElement:element model:model];
    bodyForceAtID = (model.bodyForces)[body_id-1];
    if (bodyForceAtID != nil) frictionHeat = [listUtilities listGetLogical:model inArray:bodyForceAtID.valuesList forVariable:@"friction heat" info:&found];

    // Now we start integrating
    FEMMaterialModels *materialModels = [[FEMMaterialModels alloc] init];
    FEMDifferentials *differentials = [[FEMDifferentials alloc] init];
    FEMCoordinateSystems *coordinatesSystems = [[FEMCoordinateSystems alloc] init];
    double **c2 = doublematrix(0, 2, 0, 2);
    for (t=0; t<IP->n; t++) {
    
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];

        // Basis function values & derivatives at the integration point
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:bubbles basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        
        // Coordinate system dependent info
        if (model.coordinates != cartesian) {
            x = 0.0;
            y = 0.0;
            z = 0.0;
            for (i=0; i<n; i++) {
                x = x + nodes->x[i]*integration.basis[i];
                y = y + nodes->y[i]*integration.basis[i];
                z = z + nodes->z[i]*integration.basis[i];
            }
        }
        [coordinatesSystems coordinateSystemInfoModel:model metric:metric sqrtMetric:&sqrtMetric symbols:symb dSymbols:dsymb coordX:x coordY:y coordZ:z];
        s = sqrtMetric * detJ * IP->s[t];
        
        // Coefficient of the convection and time derivative terms at the integration point
        c0 = 0.0;
        for (i=0; i<n; i++) {
            c0 = c0 + nodalC0[i]*integration.basis[i];
        }
        c1 = 0.0;
        for (i=0; i<n; i++) {
            c1 = c1 + nodalC1[i]*integration.basis[i];
        }
        ct = 0.0;
        for (i=0; i<n; i++) {
            ct = ct + nodalCT[i]*integration.basis[i];
        }
        
        // Compute effective heat capacity, if modeling phase change, at the integration point
        // Note: this is for heat equation only, not generally for diff.conv.equ.
        if (phaseChange == YES) {
            dEnth = 0.0;
            dTemp = 0.0;
            for (i=0; i<3; i++) {
                sum = 0.0;
                for (j=0; j<n; j++) {
                    sum = sum + enthalpy[j]*integration.basisFirstDerivative[j][i];
                }
                dEnth = dEnth + pow(sum, 2.0);
                sum = 0.0;
                for (j=0; j<n; j++) {
                    sum = sum + nodalTemperature[j]*integration.basisFirstDerivative[j][i];
                }
                dTemp = dTemp + pow(sum, 2.0);
            }
            ct = sqrt(dEnth/dTemp);
        }
        
        // Coefficient of the diffusion term & its derivatives at the integration point
        density = 0.0;
        for (i=0; i<n; i++) {
            density = density + nodalDensity[i]*integration.basis[i];
        }
        for (i=0; i<dim; i++) {
            for (j=0; j<dim; j++) {
                sum = 0.0;
                for (k=0; k<n; k++) {
                    sum = sum + nodalC2[i][j][k]*integration.basis[k];
                }
                c2[i][j] = sqrt(metric[i][i]) * sqrt(metric[j][j]) * sum;
            }
        }
        for (i=0; i<dim; i++) {
            c2[i][i] = [materialModels effectiveConductivity:c2[i][i] density:density element:element temperature:nodalTemperature velocityX:ux velocitY:uy velocityZ:uz nodes:nodes numberOfNodes:n numberOfPoints:n integrationU:u integrationV:v integrationW:w mesh:mesh model:model];
        }
        
        // If there's no convection term we don't need the velocities and also no need for stabilzation
        convection = NO;
        if (c1 != 0.0) {
            convection = YES;
            if (phaseChange == YES) c1 = ct;
            // Velocity from previous iteration at the integration point
            memset( velo, 0.0, sizeof(int) );
            for (i=0; i<n; i++) {
                velo[0] = velo[0] + (ux[i]-mux[i])*integration.basis[i];
                velo[1] = velo[1] + (uy[i]-muy[i])*integration.basis[i];
                if (dim > 2 && model.coordinates != axis_symmetric) velo[2] = velo[2] + (uz[i]-muz[i])*integration.basis[i];
            }
            
            if (compressible == YES) {
                pressure = 0.0;
                for (i=0; i<n; i++) {
                    pressure = pressure + nodalPressure[i]*integration.basis[i];
                }
                
                memset( *dVelodx, 0.0, (3*3)*sizeof(double) );
                for (i=0; i<3; i++) {
                    for (j=0; j<n; j++) {
                        dVelodx[0][i] = dVelodx[0][i] + ux[j] * integration.basisFirstDerivative[j][i];
                        dVelodx[1][i] = dVelodx[1][i] + uy[j] * integration.basisFirstDerivative[j][i];
                        if (dim > 2 && model.coordinates != axis_symmetric) dVelodx[2][i] = dVelodx[2][i] + uz[j] * integration.basisFirstDerivative[j][i];
                    }
                }
                
                divVelo = 0.0;
                for (i=0; i<dim; i++) {
                    divVelo = divVelo + dVelodx[i][i];
                }
                if (model.coordinates >= cylindric && model.coordinates <= axis_symmetric) { // Cylindric coordinates
                    divVelo = divVelo + velo[0]/x;
                } else { // General coordinate system
                    for (i=0; i<dim; i++) {
                        for (j=0; j<dim; j++) {
                            divVelo = divVelo + velo[j]*symb[i][j][i];
                        }
                    }
                }
            }
            
            // Stabilization parameters...
            if (stabilize == YES) {
                vnorm = 0.0;
                for (i=0; i<dim; i++) {
                    vnorm = vnorm + velo[i] * velo[i]/metric[i][i];
                }
                vnorm = sqrt(vnorm);
                
                pe = min(1.0, mk*hk*c1*vnorm/(2.0*fabs(c2[0][0])));
                tau = 0.0;
                if (vnorm != 0.0) {
                    tau = hk * pe / (2.0 * c1 * vnorm);
                }
                
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        for (k=0; k<dim; k++) {
                            sum = 0.0;
                            for (l=0; l<n; l++) {
                                sum = sum + nodalC2[i][j][l]*integration.basisFirstDerivative[l][k];
                            }
                            dc2dx[i][j][k] = sqrt(metric[i][i]) * sqrt(metric[j][j]) * sum;
                        }
                    }
                }
                
                // Compute residual & stabilization weight vectors
                for (p=0; p<n; p++) {
                    su[p] = c0 * integration.basis[p];
                    for (i=0; i<dim; i++) {
                        su[p] = su[p] + c1 * integration.basisFirstDerivative[p][i]*velo[i];
                        if (element->Type.BasisFunctionDegree <= 1) continue;
                        for (j=0; j<dim; j++) {
                            su[p] = su[p] - dc2dx[i][j][j]*integration.basisFirstDerivative[p][i];
                            sum = 0.0;
                            for (l=0; l<n; l++) {
                                sum = sum + dNodalBasisdx[p][l][i]*integration.basisFirstDerivative[l][j];
                            }
                            su[p] = su[p] - c2[i][j]*sum;
                            for (k=0; k<dim; k++) {
                                su[p] = su[p] + c2[i][j] * symb[i][j][k] * integration.basisFirstDerivative[p][k];
                                su[p] = su[p] - c2[i][k] * symb[k][j][j] + integration.basisFirstDerivative[p][i];
                                su[p] = su[p] - c2[k][j] * symb[k][j][i] * integration.basisFirstDerivative[p][i];
                            }
                        }
                    }
                    sw[p] = c0 * integration.basis[p];
                    for (i=0; i<dim; i++) {
                        sw[p] = sw[p] + c1 * integration.basisFirstDerivative[p][i]*velo[i];
                        if (element->Type.BasisFunctionDegree <= 1) continue;
                        for (j=0; j<dim; j++) {
                            sw[p] = sw[p] - dc2dx[i][j][j]*integration.basisFirstDerivative[p][i];
                            sum = 0.0;
                            for (l=0; l<n; l++) {
                                sum = sum + dNodalBasisdx[p][l][i]*integration.basisFirstDerivative[l][j];
                            }
                            sw[p] = sw[p] - c2[i][j]*sum;
                            for (k=0; k<dim; k++) {
                                sw[p] = sw[p] + c2[i][j] * symb[i][j][k] * integration.basisFirstDerivative[p][k];
                                sw[p] = sw[p] - c2[i][k] * symb[k][j][j] * integration.basisFirstDerivative[p][i];
                                sw[p] = sw[p] - c2[k][j] * symb[k][j][i] * integration.basisFirstDerivative[p][i];
                                
                            }
                        }
                    }
                }
            }
        }
        
        // Loop over basis functions of both unknowns and weights
        for (p=0; p<nBasis; p++) {
            for (q=0; q<nBasis; q++) {
                // The diffusive-convective equation without stabilization
                m = ct * integration.basis[q] * integration.basis[p];
                a = c0 * integration.basis[q] * integration.basis[p];
                
                // The diffusion term
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        a = a + c2[i][j] * integration.basisFirstDerivative[q][i] * integration.basisFirstDerivative[p][j];
                    }
                }
                
                if (convection == YES) {
                    // The convection term
                    for (i=0; i<dim; i++) {
                        a = a + c1 * velo[i] * integration.basisFirstDerivative[q][i] * integration.basis[p];
                    }
                    // Next we add the stabilization...
                    if (stabilize == YES) {
                        a = a + tau * su[q] * sw[p];
                        m = m + tau * ct * integration.basis[q] * sw[p];
                    }
                }
                stiffMatrix[p][q] = stiffMatrix[p][q] + s * a;
                massMatrix[p][q] = massMatrix[p][q] + s * m;
            }
        }
        
        // Force at the integration point
        sum = 0.0;
        for (i=0; i<n; i++) {
            sum = sum + loadVector[i]*integration.basis[i];
        }
        force = sum + [differentials jouleHeatElement:element nodes:nodes numberOfNodes:n integrationU:u integrationV:v integrationW:w mesh:mesh model:model];
        
        if (convection == YES) {
            double pcoeff = 0.0;
            for (i=0; i<n; i++) {
                pcoeff = pcoeff + nodalPressureCoeff[i]+integration.basis[i];
            }
            if (pcoeff != 0.0) {
                sum = 0.0;
                for (i=0; i<n; i++) {
                    sum = sum + nodalPressureDt[i]*integration.basis[i];
                }
                force = force + pcoeff * sum;
                for (i=0; i<dim; i++) {
                    sum = 0.0;
                    for (j=0; j<n; j++) {
                        sum = sum + nodalPressure[j]*integration.basisFirstDerivative[j][i];
                    }
                    force = force + pcoeff*velo[i]*sum;
                }
            }
            
            if (frictionHeat == YES) {
                viscosity = 0.0;
                for (i=0; i<n; i++) {
                    viscosity = viscosity + nodalviscosity[i]*integration.basis[i];
                }
                viscosity = [materialModels effectiveViscosity:viscosity density:density velocityX:ux velocitY:uy velocityZ:uz element:element nodes:nodes numberOfNodes:n numberOfPoints:n integrationU:u integrationV:v integrationW:w muder:NULL mesh:mesh model:model];
                if (viscosity > 0.0) {
                    if (compressible == NO) {
                        memset( *dVelodx, 0.0, (3*3)*sizeof(double) );
                        for (i=0; i<3; i++) {
                            for (j=0; j<n; j++) {
                                dVelodx[0][i] = dVelodx[0][i] + ux[j]*integration.basisFirstDerivative[j][i];
                                dVelodx[1][i] = dVelodx[1][i] + uy[j]*integration.basisFirstDerivative[j][i];
                                if (dim > 2 && model.coordinates != axis_symmetric)
                                                   dVelodx[2][i] = dVelodx[2][i] + uz[j]*integration.basisFirstDerivative[j][i];
                            }
                        }
                    }
                    force = force + 0.5*viscosity*[materialModels secondInvariantVelo:velo dVelodx:dVelodx crtMatrix:NULL symbols:NULL model:model];
                }
            }
        }
        
        // The righthand side...
        for (p=0; p<nBasis; p++) {
            load = integration.basis[p];
            if (convectiveAndStabilize == YES) {
                load = load + tau * sw[p];
            }
            forceVector[p] = forceVector[p] + s * load * force;
        }
    }
    free_dmatrix(c2, 0, 2, 0, 2);
    
    GaussQuadratureDeallocation(IP);
    [integration deallocation:mesh];    
}

/*********************************************************************************************************************************
    Return element local matrices and RHS vector for boundary conditions of diffusion convection equation.
 
    Arguments:
        - double **boundaryMatrix   ->  coefficient matrix if equations
        - double *boundaryVector    ->  RHS vector
        - double *loadVector        ->  coefficient of the force term
        - double *nodalAlpha        ->  coefficient for temperature dependent term
        - Element_t *element        ->  Structure describing the element (dimension, nb of nodes, interpolation degree, etc...)
        - int n                     ->  number of element nodes
        - Nodes_t *nodes            ->  Element node coordinates
*********************************************************************************************************************************/
-(void)diffuseConvectiveGeneralBoundaryMatrix:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector dimensions:(Dimensions_t)dimensions loadVector:(double *)loadVector nodalAlpha:(double *)nodalAlpha element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes model:(FEMModel *)model mesh:(FEMMesh *)mesh {
    
    int i, p, q, t;
    double alpha, detJ, force, s, x, y, z;
    BOOL stat;
    FEMCoordinateSystems *coordinateSystem;
    FEMNumericIntegration *integration;
    GaussIntegrationPoints *IP;
    
    memset( *boundaryMatrix, 0.0, (dimensions.mat1*dimensions.mat2)*sizeof(double) );
    memset( boundaryVector, 0.0, dimensions.vec*sizeof(double) );
    
    // Integration stuff
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) errorfunct("FEMDiffuseConvectiveGeneralAnisotropic:diffuseConvectiveGeneralBoundaryMatrix", "Allocation error in FEMNumericIntegration!");
    IP = GaussQuadrature(element, NULL, NULL);

    coordinateSystem = [[FEMCoordinateSystems alloc] init];
    for (t=0; t<IP->n; t++) {
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t] withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t]];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
        // Coordinate system dependent info
        if (model.coordinates != cartesian) {
            x = 0.0;
            y = 0.0;
            z = 0.0;
            for (i=0; i<n; i++) {
                x = x + (nodes->x[i]*integration.basis[i]);
                y = y + (nodes->y[i]*integration.basis[i]);
                z = z + (nodes->z[i]*integration.basis[i]);
            }
            s = s * [coordinateSystem coordinateSquareRootMetricModel:model coordX:x coordY:y coordZ:z];
        }
        
        force = 0.0;
        alpha = 0.0;
        for (i=0; i<n; i++) {
            force = force + (loadVector[i] * integration.basis[i]);
            alpha = alpha + (nodalAlpha[i] * integration.basis[i]);
        }

        for (p=0; p<n; p++) {
            for (q=0; q<n; q++) {
                boundaryMatrix[p][q] = boundaryMatrix[p][q] + s * alpha * integration.basis[q] * integration.basis[p];
            }
        }
        for (q=0; q<n; q++) {
            boundaryVector[q] = boundaryVector[q] + s * integration.basis[q] * force;
        }
    }
    GaussQuadratureDeallocation(IP);
    [integration deallocation:mesh];
}

@end
