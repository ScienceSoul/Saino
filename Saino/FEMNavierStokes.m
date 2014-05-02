//
//  FEMNavierStokes.m
//  Saino
//
//  Created by Seddik hakime on 01/05/2014.
//  Copyright (c) 2014 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMNavierStokes.h"

#import "GaussIntegration.h"
#import "Utils.h"

@implementation FEMNavierStokes {
    
    int _linearCode[8];
}

- (id)init
{
    self = [super init];
    if (self) {
        memset(_linearCode, 0, sizeof(_linearCode) );
        _linearCode[2] = 303;
        _linearCode[3] = 404;
        _linearCode[4] = 504;
        _linearCode[5] = 605;
        _linearCode[6] = 706;
        _linearCode[7] = 808;
    }
    
    return self;
}

/*********************************************************************************************************************
    Return element local matrices and RHS vector for Navier-Stokes-Equations in cartesian coordinates.
 
    Arguments:
        
        double **massMatrix                        -> output: time derivative coefficient matrix
        double **stiffMatrix                       -> output: rest of the equation coefficients
        double *forceVector                        -> output: RHS vector
        double **loadVector                         -> loadVector vector
        double *nodalViscosity                     -> nodal values for viscosity (i.e. if turbulence model or
                                                      power-law viscosity is used, the values vary in space)
        double *nodalDensity                       -> nodal values of density
        double *velocityX, *velocityY, *velocityZ  -> nodal values of velocity components from previous iteration
        double *nodalPressure                      -> nodal values of total pressure from previous iteration
        BOOL stabilize                             -> should stabilization be used ?
        BOOL pseudoCompressible                    -> should artificial compressibility be added ?
        double *nodalCompressibility               -> artificial compressibility for the nodes
        BOOL magneticForce                         -> should Lorentz force for magneto-hydrodynamics be included
        BOOL rotating                              -> is the coordinate system rotating
        double *omega                              -> if previous is YES, components of angular velocity
        BOOL newtonLinearization                   -> Picard or Newton linearization of the convetion term ?
        Element_t *element                         -> structure describing the element (dimension,nof nodes,
                                                      interpolation degree, etc...)
        int n                                      -> number of element nodes
        Nodes_t *nodes                             -> element node coordinates
*********************************************************************************************************************/
-(void)navierStokesComposeMassMatrix:(double **) massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double **)loadVector nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)velocityX velocityY:(double *)velocityY velocityZ:(double *)velocityZ meshVelocityX:(double *)meshVelocityX meshVelocityY:(double *)meshVelocityY meshVelocityZ:(double *)meshVelocityZ nodalPressure:(double *)nodalPressure nodalTemperature:(double *)nodalTemperature convect:(BOOL)convect stabilizeFlag:(NSString *)stabilizeFlag compressibilityModel:(int)compressibilityModel pseudoCompressible:(BOOL)pseudoCompressible nodalCompressibility:(double *)nodalCompressibility nodalGasConstant:(double *)nodalGasConstant porous:(BOOL)porous nodalDrag:(double **)nodalDrag potentialForce:(BOOL)potentialForce potentialField:(double *)potentialField potentialCoefficient:(BOOL)potentialCoefficient magneticForce:(BOOL)magneticForce rotating:(BOOL)rotating omega:(double *)omega divDiscretization:(BOOL)divDiscretization gradDriscretization:(BOOL)gradDriscretization newtonLinearization:(BOOL)newtonLinearization transient:(BOOL)transient element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution core:(FEMCore *)core mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration material:(FEMMaterial *)material elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems differentials:(FEMDifferentials *)differentials listUtilities:(FEMListUtilities *)listUtilities utilities:(FEMUtilities *)utilities {
    
    int c, i, j, k, p, t, dim, linearBasis, nBasis, order, tStep;
    double c1, compress, detJ, drhodp, dt=0.0, gasC, hk, hScale, mk, pressure, rho, s, sum, temperature, u, v, viscConstantCondition, w;
    double dNodalBasisdx[n][n][3], dPressuredx[3], dRhodx[3], dTemperaturedx[3], force[4], grad[3][3], gradNodal[n][4][3], lc[3][n], lrf[3], jacM[8*n][8*n],
           pBasis[n], pdBasisdx[n][3], velo[3];
    double *basis = NULL, **basisFirstDerivative = NULL;
    BOOL bubbles, compressible, found, isotropic, laplaceDiscretization, pBubbles, p2p1, stat, stabilize, viscNonnewtonian, vms;
    ElementType_t *linearType = NULL, *saveType = NULL;
    GaussIntegrationPoints *IP = NULL;
    listBuffer tensor = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer drhodp_n = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    dim = model.dimension;
    
    basis = integration.basis;
    basisFirstDerivative = integration.basisFirstDerivative;
    
    if ((solution.solutionInfo)[@"h scale"] != nil) {
        hScale = [(solution.solutionInfo)[@"h scale"] doubleValue];
    } else {
        hScale = 1.0;
    }
    
    c = dim + 1;
    if (newtonLinearization == YES) memset( *jacM, 0.0, ((8*n)*(8*n))*sizeof(double) );
    
    // Integration stuff
    nBasis = n;
    bubbles = NO;
    pBubbles = NO;
    p2p1 = NO;
    vms = ([stabilizeFlag isEqualToString:@"vms"]) ? YES : NO;
    
    isotropic = ([[listUtilities listGetString:model inArray:material.valuesList forVariable:@"viscosity model" info:&viscNonnewtonian] isEqualToString:@"anisotropic"] == NO) ? YES : NO;
    if (isotropic == NO) {
        found = [listUtilities listGetRealArray:model inArray:material.valuesList forVariable:@"viscosity" numberOfNodes:n indexes:element->NodeIndexes buffer:&tensor];
    }
    
    if (viscNonnewtonian == YES) {
        if ((solution.solutionInfo)[@"newtonian viscosity condition"] != nil) {
            viscConstantCondition = [(solution.solutionInfo)[@"newtonian viscosity condition"] doubleValue];
            if (viscConstantCondition > 0.0) viscNonnewtonian = NO;
        }
    }
    
    laplaceDiscretization = [(solution.solutionInfo)[@"laplace discretization"] boolValue];
    
    if (vms == YES && transient == YES) {
        dt = solution.dt;
        order = min(solution.doneTime, solution.order);
    }
    
    compressible = (compressibilityModel != incompressible) ? YES : NO;
    found = [core getReal:model forElement:element inArray:material.valuesList variableName:@"drho/dp" buffer:&drhodp_n listUtilities:listUtilities];
    
    stabilize = ([stabilizeFlag isEqualToString:@"stabilized"] == YES) ? YES : NO;
    if (!(vms == YES || stabilize == YES) || compressible == YES) {
        pBubbles = ([core isPElement:element] == YES && element->BDOFs > 0) ? YES : NO;
        if (pBubbles == YES) {
            nBasis = n + element->BDOFs;
        } else if ([stabilizeFlag isEqualToString:@"bubbles"] == YES || element->Type.BasisFunctionDegree <= 1) {
            nBasis = 2 * n;
            bubbles = YES;
        } else {
            p2p1 = YES;
        }
        stabilize = NO;
    }
    
    linearBasis = 0;
    if (p2p1 == YES) {
        saveType = &element->Type;
        j = [core getElementFamily:element];
        linearType = [elementDescription getElementType:_linearCode[j-1] inMesh:mesh stabilization:NULL];
        linearBasis = linearType->NumberOfNodes;
    }
    
    if (bubbles) {
        IP = GaussQuadrature(element, NULL, &element->Type.GaussPoints2);
    } else IP = GaussQuadrature(element, NULL, NULL);
    
    // Stabilization parameters: hk, mk (Franca et al.)
    hk = element->hK * hScale;
    mk = element->StabilizationMK;
    
    if (stabilize == YES) {
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
    } else if (vms == YES) {
        double **nodalVelo = doublematrix(0, 3, 0, n-1);
        double **nodalPVelo = doublematrix(0, 3, 0, n-1);
        for (i=0; i<n; i++) {
            nodalVelo[0][i] = velocityX[i];
            nodalVelo[1][i] = velocityY[i];
            nodalVelo[2][i] = velocityZ[i];
            nodalVelo[dim][i] = nodalPressure[i];
        }
        
        memset( **dNodalBasisdx, 0.0, (n*n*3)*sizeof(double) );
        memset( **gradNodal, 0.0, (n*4*3)*sizeof(double) );
        double c1[dim][dim];
        double yy[dim];
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
        
        memset( *nodalPVelo, 0.0, (4*n)*sizeof(double) );
        if (transient == YES) {
            tStep = -1;
            [core getVectorLocalField:nodalPVelo size1Field:4 size2Field:n name:nil element:element solution:solution model:model timeStep:&tStep];
            
            if (order < 2) {
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        nodalPVelo[i][j] = nodalVelo[i][j] - nodalPVelo[i][j]/dt;
                    }
                }
            } else {
                double **work = doublematrix(0, 2, 0, n-1);
                tStep = -2;
                [core getVectorLocalField:work size1Field:3 size2Field:n name:nil element:element solution:solution model:model timeStep:&tStep];
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        nodalPVelo[i][j] = (1.5*nodalVelo[i][j]) - 2.0*nodalPVelo[i][j] + 0.5*work[i][j]/dt;
                    }
                }
            }
        }
        free_dmatrix(nodalVelo, 0, 3, 0, n-1);
        free_dmatrix(nodalPVelo, 0, 3, 0, n-1);
    }
    
    c1 = 2.0 / mk;
    for (i=0; i<n; i++) {
        lc[0][i] = element->Type.NodeU[i];
        lc[1][i] = element->Type.NodeV[i];
        lc[2][i] = element->Type.NodeW[i];
    }
    
    double maxl = -HUGE_VAL;
    double minl = HUGE_VAL;
    for (i=0; i<element->Type.dimension; i++) {
        for (j=0; j<n; j++) {
            if (lc[i][j]>maxl) {
                maxl = lc[i][j];
            }
        }
        for (j=0; j<n; j++) {
            if (lc[i][j]<minl) {
                minl = lc[i][j];
            }
        }
        for (j=0; j<n; j++) {
            lc[i][j] = 2.0 * (lc[i][j]-minl) / (maxl-minl) - 1.0;
        }
    }
    
    viscNonnewtonian = NO;
    
    // Now we start integrating
    for (t=0; t<IP->n; t++) {
        
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];

        // Basis function values and derivatives at the integration point
        if (p2p1 == YES) {
            element->Type = *linearType;
            stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:bubbles basisDegree:NULL];
            stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:bubbles basisDegree:NULL];
            for (i=0; i<n; i++) {
                pBasis[i] = integration.basis[i];
                for (j=0; j<3; j++) {
                    pdBasisdx[i][j] = integration.basisFirstDerivative[i][j];
                }
            }
            element->Type = *saveType;
        }
        
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:bubbles basisDegree:NULL];
        stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:bubbles basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
        // Density at the integration point
        rho = 0.0;
        for (i=0; i<n; i++) {
            rho = rho + nodalDensity[i]*basis[i];
        }
        switch (compressibilityModel) {
            case perfect_gas1:
                if (p2p1 == YES) {
                    k = linearBasis;
                    pressure = 0.0;
                    for (i=0; i<k; i++) {
                        pressure = pressure + nodalPressure[i]*pBasis[i];
                    }
                    temperature = 0.0;
                    for (i=0; i<k; i++) {
                        temperature = temperature + nodalTemperature[i]*pBasis[i];
                    }
                    memset(dPressuredx, 0.0, sizeof(dPressuredx) );
                    memset(dTemperaturedx, 0.0, sizeof(dTemperaturedx) );
                    for (i=0; i<dim; i++) {
                        for (j=0; j<k; j++) {
                            dPressuredx[i] = dPressuredx[i] + nodalPressure[j]*pdBasisdx[j][i];
                            dTemperaturedx[i] = dTemperaturedx[i] + nodalTemperature[j]*pdBasisdx[j][i];
                        }
                    }
                } else {
                    pressure = 0.0;
                    for (i=0; i<n; i++) {
                        pressure = pressure + nodalPressure[i]*basis[i];
                    }
                    temperature = 0.0;
                    for (i=0; i<n; i++) {
                        temperature = temperature + nodalTemperature[i]*basis[i];
                    }
                    for (i=0; i<dim; i++) {
                        for (j=0; j<n; j++) {
                            dPressuredx[i] = dPressuredx[i] + nodalPressure[j]*basisFirstDerivative[j][i];
                            dTemperaturedx[i] = dTemperaturedx[i] + nodalTemperature[j]*basisFirstDerivative[j][i];
                        }
                    }
                }
                gasC = 0.0;
                for (i=0; i<n; i++) {
                    gasC = gasC + nodalGasConstant[i]*basis[i];
                }
                rho = pressure / (gasC * temperature);
                break;
                
            case user_defined1:
            case thermal:
                memset(dRhodx, 0.0, sizeof(dRhodx) );
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        dRhodx[i] = dRhodx[i] + nodalDensity[j]*basisFirstDerivative[j][i];
                    }
                }
                break;
                
            case user_defined2:
                memset(dPressuredx, 0.0, sizeof(dPressuredx) );
                for (i=0; i<n; i++) {
                    for (j=0; j<n; j++) {
                        dPressuredx[i] = dPressuredx[i] + nodalPressure[j]*basisFirstDerivative[j][i];
                    }
                }
                drhodp = 0.0;
                for (i=0; i<n; i++) {
                    drhodp = drhodp + drhodp_n.vector[i]*basis[i];
                }
                break;
        }
        
        if (pseudoCompressible == YES) {
            pressure = 0.0;
            for (i=0; i<n; i++) {
                pressure = pressure + nodalPressure[i]*basis[i];
            }
            sum = 0.0;
            for (i=0; i<n; i++) {
                sum = sum + nodalCompressibility[i]*basis[i];
            }
            compress = rho * sum;
        }
        
        // Velocity from previous iteration (relative to mesh veclocity) at the
        // integration point
        memset(velo, 0.0, sizeof(velo) );
        for (i=0; i<n; i++) {
            velo[0] = velo[0] + (velocityX[i]-meshVelocityX[i])*basis[i];
            velo[1] = velo[1] + (velocityY[i]-meshVelocityY[i])*basis[i];
            if (dim > 2) velo[2] = velo[2] + (velocityZ[i]-meshVelocityZ[i])*basis[i];
        }
        
        memset( *grad, 0.0, (3*3)*sizeof(double) );
        for (i=0; i<3; i++) {
            for (j=0; j<n; j++) {
                grad[0][i] = grad[0][i] + velocityX[j]*basisFirstDerivative[j][i];
                grad[1][i] = grad[0][i] + velocityY[j]*basisFirstDerivative[j][i];
                if (dim > 2)  grad[2][i] = grad[2][i] + velocityZ[j]*basisFirstDerivative[j][i];
            }
        }
        
        //Force at integration point
        memset(force, 0.0, sizeof(force) );
        for (i=0; i<c; i++) {
            for (j=0; j<n; j++) {
                force[i] = force[i] + loadVector[i][j]*basis[j];
            }
        }
        
        if (magneticForce == YES) {
            memset(lrf, 0.0, sizeof(lrf) );
            [differentials lorentzForceElement:element nodes:nodes numberOfNodes:n integrationU:u integrationV:v integrationW:w lorentzForce:lrf mesh:mesh model:model integration:integration coordinateSystems:coordinateSystems listUtilities:listUtilities utilities:utilities];
            for (i=0; i<dim; i++) {
                force[i] = force[i] + lrf[i] / rho;
            }
        }

    }
}

@end
