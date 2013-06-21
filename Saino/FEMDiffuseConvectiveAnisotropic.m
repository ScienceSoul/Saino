//
//  FEMDiffuseConvectiveAnisotropic.m
//  Saino
//
//  Created by Seddik hakime on 04/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#include <Accelerate/Accelerate.h>

#import "FEMDiffuseConvectiveAnisotropic.h"
#import "FEMListUtilities.h"
#import "FEMNumericIntegration.h"
#import "FEMMaterial.h"
#import "FEMBodyForce.h"
#import "FEMMaterialModels.h"
#import "GaussIntegration.h"
#import "FEMDifferentials.h"
#import "Utils.h"

@implementation FEMDiffuseConvectiveAnisotropic

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
    }
    
    return self;
}

/********************************************************************************************************************************************
    Return element local matrices and RHS vector for diffusion-convection equation (cartesian coordinates)
 
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
        BOOL useBubbles                 -> if use bubbles
        Element_t *element              -> structure describing the element (dimension,nof nodes, interpolation degree, etc...)
        int numberOfNodes               -> number of element nodes
        Nodes_t *nodes                  -> element node coordinates
********************************************************************************************************************************************/
-(void)diffuseConvectiveComposeMassMatrix:(double **)massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double *)loadVector timeDerivativeTerm:(double *)nodalCT zeroDegreeTerm:(double *)nodalC0 convectionTerm:(double *)nodalC1 diffusionTerm:(double ***)nodalC2 phaseChange:(BOOL)phaseChange nodalTemperature:(double *)nodalTemperature enthalpy:(double *)enthalpy velocityX:(double *)ux velocitY:(double *)uy velocityZ:(double *)uz meshVeloX:(double *)mux meshVeloY:(double *)muy meshVeloZ:(double *)muz nodalViscosity:(double *)nodalviscosity nodaldensity:(double *)nodalDensity nodalPressure:(double *)nodalPressure nodalPressureDt:(double *)nodalPressureDt nodalPressureCoeff:(double *)nodalPressureCoeff compressible:(BOOL)compressible stabilize:(BOOL)stabilize useBubbles:(BOOL)useBubbles element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution mesh:(FEMMesh *)mesh model:(FEMModel *)model {
    
    int i, j, k, l, p, q, t, dim, body_id, mat_id, nBasis, order, tStep;
    double a, c0, c1, ct, dEnth, dTemp, expc, force, hk, mk, dc2dx[3][3][3], detJ, divVelo, dNodalBasisdx[n][n][3], dt=0.0, gmat[3][3], grad[3][3],
           gradNodal[n][3][3], gradP[n], grav[3], gvec[3], lc[3][n], lc1, load, m, mu, **nodalPVelo, **nodalVelo, pe, pressure, reft, rho, rm[n], s,
           su[n], sw[n], tau, tau_m, temperature, sum, sum2, u, v, velo[3], vnorm, vrm[3], y[3], w;
    NSString *stabilizationFlag;
    BOOL any, bubbles, convection, convectiveAndStabilize, found, frictionHeat, stat, transient, vms;
    FEMKernel *kernel;
    FEMListUtilities *listUtilities;
    FEMNumericIntegration *integration;
    FEMMaterial *materialAtID = nil;
    FEMBodyForce *bodyForceAtID = nil;
    GaussIntegrationPoints *IP;
    listBuffer gwrk = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    kernel = [FEMKernel sharedKernel];
        
    stabilizationFlag = (solution.solutionInfo)[@"stabilization method"];
    vms = ([stabilizationFlag isEqualToString:@"vms"] == YES) ? YES : NO;
    
    listUtilities = [[FEMListUtilities alloc] init];
    transient = ([[listUtilities listGetString:model inArray:model.simulation.valuesList forVariable:@"simulation type" info:&found] isEqualToString:@"transient"] == YES) ? YES : NO;
    
    dim = model.dimension;
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
    if (convection == YES && !(vms == YES || stabilize == YES) && useBubbles == YES) {
        nBasis = 2*n;
        bubbles = YES;
    }
    
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) errorfunct("diffuseConvectiveComposeyMassMatrix", "Allocation error in FEMNumericIntegration!");
    // Integration stuff
    if (bubbles == YES) {
        IP = GaussQuadrature(element, NULL, &element->Type.GaussPoints2);
    } else IP = GaussQuadrature(element, NULL, NULL);
    
    // Stabilization parameters: hk, mk (Franca et al.)
    // If there is no convection term, we don't need stabilization
    hk = element->hK;
    mk = element->StabilizationMK;
    
    convectiveAndStabilize = NO;
    if (vms == YES) {
        nodalVelo = doublematrix(0, 3, 0, n-1);
        nodalPVelo = doublematrix(0, 3, 0, n-1);
        for (i=0; i<n; i++) {
            nodalVelo[0][i] = ux[i];
            nodalVelo[1][i] = uy[i];
            nodalVelo[2][i] = uz[i];
            nodalVelo[dim][i] = nodalPressure[i];
        }
        memset( **dNodalBasisdx, 0.0, (n*n*3)*sizeof(double) );
        memset( **gradNodal, 0.0, (n*3*3)*sizeof(double) );
        double c1[dim][dim];
        double yy[dim];
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
            memset( *c1, 0.0, (dim*dim)*sizeof(double) );
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, n, 1.0, *nodalVelo, n, *integration.basisFirstDerivative, 3, 0.0, (double *)c1, dim);
            for (i=0; i<dim; i++) {
                for (j=0; j<dim; j++) {
                    gradNodal[p][i][j] = c1[i][j];
                }
            }
            cblas_dgemv(CblasRowMajor, CblasNoTrans, n, dim, 1.0, *integration.basisFirstDerivative, 3, *nodalVelo+(dim*n), 1.0, 0.0, yy, 1.0);
            for (i=0; i<dim; i++) {
                gradNodal[p][dim][i] = yy[i];
            }
        }
        memset( *nodalPVelo, 0.0, (4*n)*sizeof(double) );
        if (transient) {
            dt = solution.dt;
            order = min(solution.doneTime, solution.order);
            tStep = -1;
            [kernel getVectorLocalField:nodalPVelo size1Field:4 size2Field:n name:@"flow solution" element:element solution:solution model:model timeStep:&tStep];
            if (order < 2) {
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        nodalPVelo[i][j] = nodalVelo[i][j] - nodalPVelo[i][j]/dt;
                    }
                }
            } else {
                double **work = doublematrix(0, 2, 0, n-1);
                tStep = -2;
                [kernel getVectorLocalField:work size1Field:3 size2Field:n name:@"flow solution" element:element solution:solution model:model timeStep:&tStep];
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        nodalPVelo[i][j] = (1.5*nodalVelo[i][j]) - 2.0*nodalPVelo[i][j] + 0.5*work[i][j]/dt;
                    }
                }
            }
        }
        
        mat_id = [kernel getMaterialIDForElement:element model:model];
        materialAtID = (model.materials)[mat_id-1];

        expc = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"heat expansion coefficient" info:&found minValue:NULL maxValue:NULL];
        reft = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"reference temperature" info:&found minValue:NULL maxValue:NULL];
        found = [listUtilities listGetConstRealArray:model inArray:model.constants.valuesList forVariable:@"grav" buffer:&gwrk];
        if (found == YES) {
            for (i=0; i<3; i++) {
                grav[i] = gwrk.matrix[i][0]*gwrk.matrix[3][0];
            }
        } else {
            memset( grav, 0.0, sizeof(double) );
            grav[1] = -9.81;
        }
        lc1 = 2.0/mk;
        for (i=0; i<n; i++) {
            lc[0][i] = element->Type.NodeU[i];
            lc[1][i] = element->Type.NodeV[i];
            lc[2][i] = element->Type.NodeW[i];
        }
        double maxl, minl;
        for (i=0; i<element->Type.dimension; i++) {
            minl = HUGE_VAL;
            for (j=0; j<n; j++) {
                if (lc[i][j]<minl) {
                    minl = lc[i][j];
                }
            }
            maxl = -HUGE_VAL;
            for (j=0; j<n; j++) {
                if (lc[i][j]>maxl) {
                    maxl = lc[i][j];
                }
            }
            for (j=0; j<n; j++) {
                lc[i][j] = 2.0* (lc[i][j]-minl)/(maxl-minl) - 1.0;
            }
        }
        
        if (gwrk.matrix != NULL) {
            free_dmatrix(gwrk.matrix, 0, gwrk.m-1, 0, gwrk.n-1);
            gwrk.matrix = NULL;
        }
        free_dmatrix(nodalVelo, 0, 3, 0, n-1);
        free_dmatrix(nodalPVelo, 0, 3, 0, n-1);
    } else if (stabilize == YES && convection == YES) {
        convectiveAndStabilize = YES;
        hk = element->hK;
        mk = element->StabilizationMK;
        memset( **dNodalBasisdx, 0.0, (n*n*3)*sizeof(double) );
        for (p=0; p<n; p++) {
            u = element->Type.NodeU[p];
            v = element->Type.NodeV[p];
            w = element->Type.NodeW[p];
        }
        stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        for (i=0; i<n; i++) {
            for (j=0; j<3; j++) {
                dNodalBasisdx[i][p][j] = integration.basisFirstDerivative[i][j];
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
    double **c2 = doublematrix(0, 2, 0, 2);
    for (t=0; t<IP->n; t++) {
        
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];
        
        // Basis function values & derivatives at the integration point
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
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
        rho = 0.0;
        for (i=0; i<n; i++) {
            rho = rho + nodalDensity[i]*integration.basis[i];
        }
        
        for (i=0; i<dim; i++) {
            for (j=0; j<dim; j++) {
                sum = 0.0;
                for (k=0; k<n; k++) {
                    sum = sum + nodalC2[i][j][k]*integration.basis[k];
                }
                c2[i][j] = sum;
            }
        }
        for (i=0; i<dim; i++) {
            c2[i][i] = [materialModels effectiveConductivity:c2[i][i] density:rho element:element temperature:nodalTemperature velocityX:ux velocitY:uy velocityZ:uz nodes:nodes numberOfNodes:n numberOfPoints:n integrationU:u integrationV:v integrationW:w mesh:mesh model:model];
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
                if (dim > 2) velo[2] = velo[2] + (uz[i]-muz[i])*integration.basis[i];
            }
            
            if (compressible == YES) {
                memset( *grad, 0.0, (3*3)*sizeof(double) );
                for (i=0; i<3; i++) {
                    for (j=0; j<n; j++) {
                        grad[0][i] = grad[0][i] + ux[j] * integration.basisFirstDerivative[j][i];
                        grad[1][i] = grad[1][i] + uy[j] * integration.basisFirstDerivative[j][i];
                        if (dim > 2) grad[2][i] = grad[2][i] + uz[j] * integration.basisFirstDerivative[j][i];
                    }
                }
                
                pressure = 0.0;
                for (i=0; i<n; i++) {
                    pressure = pressure + nodalPressure[i]*integration.basis[i];
                }
                divVelo = 0.0;
                for (i=0; i<dim; i++) {
                    divVelo = divVelo + grad[i][i];
                }
            }
            
            if (vms == YES) {
                mu = [listUtilities listGetConstReal:model inArray:materialAtID.valuesList forVariable:@"viscosity" info:&found minValue:NULL maxValue:NULL];
                mu = [materialModels effectiveViscosity:mu density:rho velocityX:ux velocitY:uy velocityZ:uz element:element nodes:nodes numberOfNodes:n numberOfPoints:n integrationU:u integrationV:v integrationW:w muder:NULL mesh:mesh model:model];
                
                memset( *grad, 0.0, (3*3)*sizeof(double) );
                for (i=0; i<3; i++) {
                    for (j=0; j<n; j++) {
                        grad[0][i] = grad[0][i] + ux[j] * integration.basisFirstDerivative[j][i];
                        grad[1][i] = grad[1][i] + uy[j] * integration.basisFirstDerivative[j][i];
                        if (dim > 2) grad[2][i] = grad[2][i] + uz[j] * integration.basisFirstDerivative[j][i];
                    }
                }
                sum = 0.0;
                for (i=0; i<dim; i++) {
                    sum = sum + pow(velo[i], 2.0);
                }
                vnorm = sqrt(sum);
                temperature = 0.0;
                for (i=0; i<n; i++) {
                    temperature = temperature + integration.basis[i]*nodalTemperature[i];
                }
                memset( gradP, 0.0, sizeof(double) );
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        gradP[i] = gradP[i] + nodalPressure[j]*integration.basisFirstDerivative[j][i];
                    }
                }
                
                memset( *gmat, 0.0, (3*3)*sizeof(double) );
                memset( gvec, 0.0, sizeof(double) );
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        sum = 0.0;
                        for (l=0; l<n; l++) {
                            sum = sum + lc[j][l]*integration.basisFirstDerivative[l][i];
                        }
                        gvec[i] = gvec[i] + sum;
                        for (k=0; k<dim; k++) {
                            sum = 0.0;
                            sum2 = 0.0;
                            for (l=0; l<n; l++) {
                                sum = sum + lc[k][l]*integration.basisFirstDerivative[l][i];
                                sum2 = sum2 + lc[k][l]*integration.basisFirstDerivative[l][j];
                            }
                            gmat[i][j] = gmat[i][j] + sum * sum2;
                        }
                    }
                }
                
                sum = 0.0;
                memset( y, 0.0, sizeof(double) );
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, (double *)gmat, 3, velo, 1, 0.0, y, 1);
                for (i=0; i<3; i++) {
                    sum = sum + velo[i]*y[i];
                }
                sum2 = 0.0;
                for (i=0; i<3; i++) {
                    for (j=0; j<3; j++) {
                        sum2 = sum2 + gmat[i][j]*gmat[i][j];
                    }
                }
                if (transient == YES) {
                    tau_m = 1.0 / sqrt( sum + pow(lc1, 2.0) * pow((mu/rho), 2.0) * sum2/dim + 4.0/pow(dt, 2.0));
                } else {
                    tau_m = 1.0 / sqrt( sum + pow(lc1, 2.0) * pow((mu/rho), 2.0) * sum2/dim);
                }
                
                pe = min(1.0, mk*hk*c1*vnorm/(2.0*fabs(c2[0][0])));
                if (vnorm != 0.0) {
                    tau = hk * pe / (2.0 * c1 * vnorm);
                }
                
                memset( rm, 0.0, sizeof(double) );
                for (p=0; p<n; p++) {
                    rm[p] = c0 * integration.basis[p];
                    for (i=0; i<dim; i++) {
                        rm[p] = rm[p] + c1 * velo[i] * integration.basisFirstDerivative[p][i];
                        for (j=0; j<dim; j++) {
                            sum = 0.0;
                            for (l=0; l<n; l++) {
                                sum = sum + dNodalBasisdx[p][l][i]*integration.basisFirstDerivative[l][j];
                            }
                            rm[p] = rm[p] - c2[i][j]*sum;
                        }
                    }
                }
                
                memset( vrm, 0.0, sizeof(double) );
                for (i=0; i<dim; i++) {
                    for (l=0; l<n; l++) {
                        vrm[i] = vrm[i] + nodalPVelo[i][l]*integration.basis[l];
                    }
                    for (j=0; j<dim; j++) {
                        vrm[i] = vrm[i] + velo[j] * grad[i][j];
                        sum = 0.0;
                        for (l=0; l<n; l++) {
                            sum = sum + gradNodal[l][i][j]*integration.basisFirstDerivative[l][j];
                        }
                        vrm[i] = vrm[i] - (mu/rho)*sum;
                    }
                    vrm[i] = vrm[i] + gradP[i];
                    vrm[i] = vrm[i] + grav[i] * expc * (temperature - reft);
                }
            } else if (stabilize == YES) {
                // Stabilization parameter tau
                sum = 0.0;
                for (i=0; i<dim; i++) {
                    sum = sum + pow(velo[i], 2.0);
                }
                vnorm = sqrt(sum);
                
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
                            dc2dx[i][j][k] = sum;
                        }
                    }
                }
                
                // Compute residual & stabilization vectors
                for (p=0; p<n; p++) {
                    su[p] = c0 * integration.basis[p];
                    for (i=0; i<dim; i++) {
                        su[p] = su[p] + c1 * integration.basisFirstDerivative[p][i]*velo[i];
                        for (j=0; j<dim; j++) {
                            su[p] = su[p] - dc2dx[i][j][j]*integration.basisFirstDerivative[p][i];
                            sum = 0.0;
                            for (l=0; l<n; l++) {
                                sum = sum + dNodalBasisdx[p][l][i]*integration.basisFirstDerivative[l][j];
                            }
                            su[p] = su[p] - c2[i][j]*sum;
                        }
                    }
                    sw[p] = c0 * integration.basis[p];
                    for (i=0; i<dim; i++) {
                        sw[p] = sw[p] + c1 * integration.basisFirstDerivative[p][i]*velo[i];
                        for (j=0; j<dim; j++) {
                            sw[p] = sw[p] - dc2dx[i][j][j]*integration.basisFirstDerivative[p][i];
                            sum = 0.0;
                            for (l=0; l<n; l++) {
                                sum = sum + dNodalBasisdx[p][l][i]*integration.basisFirstDerivative[l][j];
                            }
                            sw[p] = sw[p] - c2[i][j]*sum;
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
                    if (vms == YES) {
                        for (i=0; i<dim; i++) {
                            a = a - c1 * tau_m * vrm[i] * integration.basisFirstDerivative[q][i] * integration.basis[p];
                            
                            a = a + c1 * velo[i] * tau * rm[q] * integration.basisFirstDerivative[p][i];
                            m = m + c1 * velo[i] * tau * ct * integration.basis[q] * integration.basisFirstDerivative[p][i];
                            
                            a = a - c1 * tau_m * vrm[i] * tau * rm[q] * integration.basisFirstDerivative[p][i];
                            m = m - c1 * tau_m * vrm[i] * tau * ct * integration.basis[q] * integration.basisFirstDerivative[p][i];
                        }
                    } else if (stabilize == YES) {
                        a = a + tau * su[q] * sw[p];
                        m = m + tau * ct * integration.basis[q] * sw[p];
                    }
                }
                stiffMatrix[p][q] = stiffMatrix[p][q] + s * a;
                massMatrix[p][q] = massMatrix[p][q] + s * m;
            }
        }

        // The righthand side....
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
                mu = 0.0;
                for (i=0; i<n; i++) {
                    mu = mu + nodalviscosity[i]*integration.basis[i];
                }
                mu = [materialModels effectiveViscosity:mu density:rho velocityX:ux velocitY:uy velocityZ:uz element:element nodes:nodes numberOfNodes:n numberOfPoints:n integrationU:u integrationV:v integrationW:w muder:NULL mesh:mesh model:model];
                if (mu > 0.0) {
                    if (compressible == NO) {
                        memset( *grad, 0.0, (3*3)*sizeof(double) );
                        for (i=0; i<3; i++) {
                            for (j=0; j<n; j++) {
                                grad[0][i] = grad[0][i] + ux[j]*integration.basisFirstDerivative[j][i];
                                grad[1][i] = grad[1][i] + uy[j]*integration.basisFirstDerivative[j][i];
                                if (dim > 2) grad[2][i] = grad[2][i] + uz[j]*integration.basisFirstDerivative[j][i];
                            }
                        }
                    }
                    force = force + 0.5*mu*[materialModels secondInvariantVelo:velo dVelodx:grad crtMatrix:NULL symbols:NULL model:model];
                }
            }
        }
        
        for (p=0; p<nBasis; p++) {
            load = force * integration.basis[p];
            if (vms == YES) {
                for (i=0; i<dim; i++) {
                    load = load + c1 * velo[i] * tau * force * integration.basisFirstDerivative[p][i];
                    load = load - c1 * tau_m * vrm[i] * tau * force * integration.basisFirstDerivative[p][i];
                }
            } else if (convectiveAndStabilize == YES) {
                load = load + tau * force * sw[p];
            }
            forceVector[p] = forceVector[p] + s * load;
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
-(void)diffuseConvectiveBoundaryMatrix:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector dimensions:(Dimensions_t)dimensions loadVector:(double *)loadVector nodalAlpha:(double *)nodalAlpha element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh; {
    
    int i, p, q, t;
    double alpha, detJ, force, s;
    BOOL stat;
    FEMNumericIntegration *integration;
    GaussIntegrationPoints *IP;
    
    memset( *boundaryMatrix, 0.0, (dimensions.mat1*dimensions.mat2)*sizeof(double) );
    memset( boundaryVector, 0.0, dimensions.vec*sizeof(double) );
    
    // Integration stuff
    integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) errorfunct("diffuseConvectiveBoundaryMatrix", "Allocation error in FEMNumericIntegration!");
    IP = GaussQuadrature(element, NULL, NULL);

    for (t=0; t<IP->n; t++) {
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t] withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:IP->u[t] secondEvaluationPoint:IP->v[t] thirdEvaluationPoint:IP->w[t]];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
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
