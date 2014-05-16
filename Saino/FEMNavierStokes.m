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
-(void)navierStokesComposeMassMatrix:(double **) massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double **)loadVector nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz meshVelocityX:(double *)mux meshVelocityY:(double *)muy meshVelocityZ:(double *)muz nodalPressure:(double *)nodalPressure nodalTemperature:(double *)nodalTemperature convect:(BOOL)convect stabilizeFlag:(NSString *)stabilizeFlag compressibilityModel:(int)compressibilityModel pseudoCompressible:(BOOL)pseudoCompressible nodalCompressibility:(double *)nodalCompressibility nodalGasConstant:(double *)nodalGasConstant porous:(BOOL)porous nodalDrag:(double **)nodalDrag potentialForce:(BOOL)potentialForce potentialField:(double *)potentialField potentialCoefficient:(double *)potentialCoefficient magneticForce:(BOOL)magneticForce rotating:(BOOL)rotating omega:(double *)omega divDiscretization:(BOOL)divDiscretization gradDriscretization:(BOOL)gradDriscretization newtonLinearization:(BOOL)newtonLinearization transient:(BOOL)transient element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution core:(FEMCore *)core mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration material:(FEMMaterial *)material elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems materialModels:(FEMMaterialModels*)materialModels differentials:(FEMDifferentials *)differentials listUtilities:(FEMListUtilities *)listUtilities utilities:(FEMUtilities *)utilities {
    
    int c, i, j, k, l, p, q, t, dim, linearBasis, nBasis, order, tStep;
    double baseP, c1, compress, delta, detJ, drhodp, dt=0.0, gasC, hk, hScale, lambda=1.0, massCoeff, mk, mu, muder, muder0, pressure, re, rho, s, sum, sum1,
           tau, tau_c, tau_m, temperature, u, v, viscConstantCondition, vnorm, w;
    double b[6][3], coord[3], dNodalBasisdx[n][n][3], dPressuredx[3], drag[3], dRhodx[3], dTemperaturedx[3], dmudx[3], force[4], g[3][6], gmat[3][3],
           gvec[3], grad[3][3], gradNodal[n][4][3], gradP[n], gradT[3][3], lc[3][n], lrf[3], jacM[8*n][8*n], pBasis[n], pdBasisdx[n][3], prm[4], pVelo[3],
           rm[n][4][4], strain[3][3], su[n][4][4], sw[n][4][4], vc[6][6], velo[3], uVelo[3];
    double **nodalVelo = NULL, **nodalPVelo = NULL;
    double *basis = NULL, **basisFirstDerivative = NULL;
    BOOL bubbles, compressible, found, isotropic, laplaceDiscretization, pBubbles, p2p1, stat, stabilize, viscNewtonLin, viscNonnewtonian, vms;
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
    
    c = dim;
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
        nodalVelo = doublematrix(0, 3, 0, n-1);
        nodalPVelo = doublematrix(0, 3, 0, n-1);
        for (i=0; i<n; i++) {
            nodalVelo[0][i] = ux[i];
            nodalVelo[1][i] = uy[i];
            nodalVelo[2][i] = uz[i];
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
            velo[0] = velo[0] + (ux[i]-mux[i])*basis[i];
            velo[1] = velo[1] + (uy[i]-muy[i])*basis[i];
            if (dim > 2) velo[2] = velo[2] + (uz[i]-muz[i])*basis[i];
        }
        
        memset( *grad, 0.0, (3*3)*sizeof(double) );
        for (i=0; i<3; i++) {
            for (j=0; j<n; j++) {
                grad[0][i] = grad[0][i] + ux[j]*basisFirstDerivative[j][i];
                grad[1][i] = grad[0][i] + uy[j]*basisFirstDerivative[j][i];
                if (dim > 2)  grad[2][i] = grad[2][i] + uz[j]*basisFirstDerivative[j][i];
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
        
        if (rotating == YES) {
            memset(coord, 0.0, sizeof(coord) );
            for (i=0; i<n; i++) {
                coord[0] = coord[0] + basis[i]*nodes->x[i];
                coord[1] = coord[1] = basis[i]*nodes->y[i];
                coord[2] = coord[2] = basis[i]*nodes->z[i];
            }
            
            // Langranges formula is used to simplify the triple product
            // omega x (omega x coord) = omega(omega.coord) - coord(omega.omega)
            
            // This is will be multiplied by density later on
            sum = 0.0;
            for (i=0; i<dim; i++) {
                sum = sum + omega[i]*coord[i];
            }
            for (i=0; i<dim; i++) {
                force[i] = force[i] - omega[i] * sum;
            }
            sum = 0.0;
            for (i=0; i<dim; i++) {
                sum = sum + omega[i]*omega[i];
            }
            for (i=0; i<dim; i++) {
                force[i] = force[i] + coord[i] * sum;
            }
        }
        
        // Additional forces due to gradient forces (electrokinetic flow) and viscous flag in porous media
        if (potentialForce == YES) {
            sum = 0.0;
            for (i=0; i<n; i++) {
                sum = sum + potentialCoefficient[i]*basis[i];
            }
            for (i=0; i<dim; i++) {
                sum1 = 0.0;
                for (j=0; j<n; j++) {
                    sum1 = sum1 + potentialField[j]*basisFirstDerivative[j][i];
                }
                force[i] = force[i] -  sum * sum1;
            }
        }
        
        if (porous == YES) {
            memset(drag, 0.0, sizeof(drag) );
            for (i=0; i<dim; i++) {
                for (j=0; j<n; j++) {
                    drag[i] = drag[i] + nodalDrag[i][j]*basis[j];
                }
            }
        }
        
        if (convect == YES && newtonLinearization == YES) {
            memset(uVelo, 0.0, sizeof(uVelo) );
            for (i=0; i<n; i++) {
                uVelo[0] = uVelo[0] + basis[i]*ux[i];
                uVelo[1] = uVelo[1] + basis[i]*uy[i];
                if (dim > 2) uVelo[2] = uVelo[2] + basis[i]*uz[i];
            }
            
            for (i=0; i<dim; i++) {
                for (j=0; j<dim; j++) {
                    force[i] = force[i] + grad[i][j] + uVelo[j];
                }
            }
        }
        
        // Effective viscosity and derivatives at integration point
        if (isotropic == YES) {
            mu = 0.0;
            for (i=0; i<n; i++) {
                mu = mu + nodalViscosity[i]*basis[i];
            }
            
            if (viscNonnewtonian == YES) {
                mu = [materialModels effectiveViscosity:mu density:rho velocityX:ux velocitY:uy velocityZ:uz element:element nodes:nodes numberOfNodes:n numberOfPoints:n integrationU:u integrationV:v integrationW:w muder:&muder0 mesh:mesh model:model integration:integration];
                
                viscNewtonLin = (newtonLinearization == YES && muder0 != 0.0) ? YES : NO;
                if (viscNewtonLin == YES) {
                    // Transpose
                    for (i=0; i<3; i++) {
                        for (j=0; j<3; j++) {
                            gradT[j][i] = grad[i][j];
                        }
                    }
                    for (i=0; i<3; i++) {
                        for (j=0; j<3; j++) {
                            strain[i][j] = (grad[i][j] + gradT[i][j]) / 2.0;
                        }
                    }
                }
            }
        } else {
            memset(*vc, 0.0, (6*6)*sizeof(double) );
            for (i=0; i<6; i++) {
                for (j=0; j<6; j++) {
                    for (k=0; k<n; k++) {
                        vc[i][j] = vc[i][j] + tensor.tensor[i][j][k]*basis[k];
                    }
                }
            }
        }
        
        if (stabilize == YES) {
            memset(dmudx, 0.0, sizeof(dmudx) );
            for (i=0; i<3; i++) {
                for (j=0; j<n; j++) {
                    dmudx[i] = dmudx[i] + nodalViscosity[j]*basisFirstDerivative[j][i];
                }
            }
            // Stabilization parameters Tau and Delta
            if (convect == YES) {
                sum = 0.0;
                for (i=0; i<dim; i++) {
                    sum = sum + pow(velo[i], 2.0);
                }
                vnorm = max(sqrt(sum), 1.0e-12);
                re = min(1.0, rho * mk * hk *vnorm / (4.0 * mu));
                
                tau = hk * re / (2.0 * rho * vnorm);
                delta = rho * lambda * re * hk * vnorm;
            } else {
                delta = 0.0;
                tau = mk * pow(hk, 2.0) / (8.0 * mu);
            }
            // su will contain residual of NS-equations (except for the time derivative and force terms).
            // sw will contain the wright function values
            memset( **su, 0.0, (n*4*4)*sizeof(double) );
            memset( **sw, 0.0, (n*4*4)*sizeof(double) );
            for (p=0; p<n; p++) {
                for (i=0; i<dim; i++) {
                    su[p][i][c] = su[p][i][c] + basisFirstDerivative[p][i];
                    if (porous == YES) {
                        su[p][i][i] = su[p][i][i] + mu *drag[i] * basis[p];
                    }
                    
                    if (convect == YES) {
                        for (j=0; j<dim; j++) {
                            su[p][i][i] = su[p][i][i] + rho * basisFirstDerivative[p][j] * velo[j];
                        }
                    }
                    
                    for (j=0; j<dim; j++) {
                        su[p][i][i] = su[p][i][i] - dmudx[j] * basisFirstDerivative[p][j];
                        su[p][i][j] = su[p][i][j] - dmudx[j] * basisFirstDerivative[p][i];
                        sum = 0.0;
                        for (k=0; k<n; k++) {
                            sum = sum + dNodalBasisdx[p][k][j] * basisFirstDerivative[k][j];
                        }
                        su[p][i][i] = su[p][i][i] - mu * sum;
                        sum = 0.0;
                        for (k=0; k<n; k++) {
                            sum = sum + dNodalBasisdx[p][k][i] * basisFirstDerivative[k][j];
                        }
                        su[p][i][j] = su[p][i][j] - mu * sum;
                    }
                    
                    if (convect == YES && newtonLinearization == YES) {
                        for (j=0; j<dim; j++) {
                            su[p][i][j] = su[p][i][j] + rho * grad[i][j] * basis[p];
                        }
                    }
                    
                    if (convect == YES) {
                        sw[p][c][i] = sw[p][c][i] + rho * basisFirstDerivative[p][i];
                        for (j=0; j<dim; j++) {
                            sw[p][i][i] = sw[p][i][i] + rho * basisFirstDerivative[p][j] * velo[j];
                        }
                    } else {
                        sw[p][c][i] = sw[p][c][i] + basisFirstDerivative[p][i];
                    }
                    
                    for (j=0; j<dim; j++) {
                        sw[p][i][i] = sw[p][i][i] - dmudx[j] * basisFirstDerivative[p][j];
                        sw[p][j][i] = sw[p][j][i] - dmudx[j] * basisFirstDerivative[p][i];
                        sum = 0.0;
                        for (k=0; k<n; k++) {
                            sum = sum + dNodalBasisdx[p][k][j] * basisFirstDerivative[k][j];
                        }
                        sw[p][i][i] = sw[p][i][i] - mu * sum;
                        sum = 0.0;
                        for (k=0; k<n; k++) {
                            sum = sum + dNodalBasisdx[p][k][i] * basisFirstDerivative[k][j];
                        }
                        sw[p][j][i] = sw[p][j][i] - mu * sum;
                    }
                }
            }
        } else if (vms == YES) {
            mu = mu / rho;
            rho = 1.0;
            for (i=0; i<dim; i++) {
                sum = 0.0;
                for (j=0; j<n; j++) {
                    sum = sum + nodalPVelo[i][j]*basis[j];
                }
                pVelo[i] = rho * sum;
                memset(gradP, 0.0, sizeof(gradP) );
                for (j=0; j<n; j++) {
                    gradP[i] = gradP[i] + nodalVelo[dim][j]*basisFirstDerivative[j][i];
                }
            }
            
            memset(*gmat, 0.0, (3*3)*sizeof(double) );
            memset(gvec, 0.0, sizeof(gvec) );
            for (i=0; i<dim; i++) {
                for (j=0; j<dim; j++) {
                    sum = 0.0;
                    for (k=0; k<n; k++) {
                        sum = sum + lc[j][k] * basisFirstDerivative[k][i];
                    }
                    gvec[i] = gvec[i] + sum;
                    
                    for (k=0; k<dim; k++) {
                        sum = 0.0;
                        for (l=0; l<n; l++) {
                            sum = sum + lc[k][l] * basisFirstDerivative[l][i];
                        }
                        sum1 = 0.0;
                        for (l=0; l<n; l++) {
                            sum1 = sum1 + lc[k][l] * basisFirstDerivative[l][j];
                        }
                        gmat[i][j] = gmat[i][j] + sum * sum1;
                    }
                }
            }
            
            double yy[3];
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, (double *)gmat, 3, velo, 1.0, 0.0, yy, 1.0);
            sum = 0.0;
            for (i=0; i<3; i++) {
                sum = sum + velo[i]*yy[i];
            }
            sum1 = 0.0;
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    sum1 = sum1 + gmat[i][j] * gmat[i][j];
                }
            }
            if (transient == YES) {
                tau_m = 1.0 / sqrt( sum + pow(c1, 2.0) * pow((mu/rho), 2.0) * sum1 / dim + 4.0 / pow(dt, 2.0) );
            } else {
                tau_m = 1.0 / sqrt( sum + pow(c1, 2.0) * pow((mu/rho), 2.0) * sum1 / dim );
            }
            sum = 0.0;
            for (i=0; i<3; i++) {
                sum = sum + gvec[i]*gvec[i];
            }
            tau_c = 1.0 / (tau_m * sum);
            
            memset( **rm, 0.0, (n*4*4)*sizeof(double) );
            for (p=0; p<n; p++) {
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        rm[p][i][i] = rm[p][i][i] + rho * velo[j] * basisFirstDerivative[p][j];
                        sum = 0.0;
                        for (k=0; k<n; k++) {
                            sum = sum + dNodalBasisdx[p][k][j] * basisFirstDerivative[k][j];
                        }
                        rm[p][i][i] = rm[p][i][i] - mu * sum;
                    }
                    rm[p][i][dim] = rm[p][i][dim] + basisFirstDerivative[p][i];
                    rm[p][dim][i] = rm[p][dim][i] + rho * basisFirstDerivative[p][i];
                }
            }
            
            memset(prm, 0.0, sizeof(prm) );
            for (i=0; i<dim; i++) {
                for (j=0; j<dim; j++) {
                    prm[i] = prm[i] + rho * velo[j] * grad[i][j];
                    sum = 0.0;
                    for (k=0; k<n; k++) {
                        sum = sum + gradNodal[k][i][j] * basisFirstDerivative[k][j];
                    }
                    prm[i] = prm[i] - mu * sum;
                }
                prm[i] = prm[i] + gradP[i];
                prm[dim] = prm[dim] + rho * grad[i][i];
            }
        }
        
        // Loop over basis functions (of both unknowns and weights)
        for (p=0; p<nBasis; p++) {
            if (isotropic == NO) {
                memset( *g, 0.0, (3*6)*sizeof(double) );
                g[0][0] = basisFirstDerivative[p][0];
                g[1][1] = basisFirstDerivative[p][1];
                g[2][2] = basisFirstDerivative[p][2];
                g[0][3] = basisFirstDerivative[p][1];
                g[1][3] = basisFirstDerivative[p][0];
                g[1][4] = basisFirstDerivative[p][2];
                g[2][4] = basisFirstDerivative[p][1];
                g[0][5] = basisFirstDerivative[p][2];
                g[2][5] = basisFirstDerivative[p][0];
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 6, 6, 1.0, (double *)g, 6, (double *)vc, 6, 0.0, (double *)g, 6);
            }
            
            for (q=0; q<nBasis; q++) {
                i = (c+1) * p;
                j = (c+1) * q;
                double m[c+1][c+1];
                double a[c+1][c+1];
                double jac[c+1][c+1];
                int kk = 0;
                int ll = 0;
                for (k=i; k<=i+c; k++) {
                    for (l=j; l<=j+c; l++) {
                        m[kk][ll] = massMatrix[k][l];
                        a[kk][ll] = stiffMatrix[k][l];
                        if (viscNewtonLin) jac[kk][ll] = jacM[k][l];
                        ll++;
                    }
                    kk++;
                }
                
                // First plain Navier-Stokes
                if (p2p1 == YES && p < linearBasis) {
                    baseP = pBasis[p];
                } else {
                    baseP = basis[p];
                }
                
                // Mass matrix
                // Momentum equations
                for (i=0; i<dim; i++) {
                    m[i][i] = m[i][i] + s * rho * basis[q] * basis[p];
                }
                
                // Mass for the continuity equation (in terms of pressure)
                if (compressibilityModel == perfect_gas1) {
                    m[c][c] = m[c][c] + s * (rho / pressure) * basis[q] * baseP;
                } else if (compressibilityModel == user_defined2) {
                    m[c][c] = m[c][c] + s * drhodp * basis[q] * baseP;
                }
                
                // Stiffness matrix
                // Rotating coordinates
                if (rotating == YES) {
                    massCoeff = 2.0 * s * rho *basis[q] * basis[p];
                    a[0][1] = a[0][1] - massCoeff * omega[2];
                    a[1][0] = a[1][0] + massCoeff * omega[2];
                    if (dim == 3) {
                        a[0][2] = a[0][2] + massCoeff * omega[1];
                        a[1][2] = a[1][2] - massCoeff * omega[0];
                        a[2][1] = a[2][1] + massCoeff * omega[0];
                        a[2][0] = a[2][0] - massCoeff * omega[1];
                    }
                }
                
                // Possible Porous media effects
                if (porous == YES) {
                    for (i=0; i<dim; i++) {
                        a[i][i] = a[i][i] + s * mu * drag[i] * basis[q] * basis[p];
                    }
                }
                
                // Diffusive terms
                // Convection terms, Picard linearization
                if (isotropic == NO) {
                    memset( *b, 0.0, (6*3)*sizeof(double) );
                    double c[3][3];
                    memset( *c, 0.0, (3*3)*sizeof(double) );
                    b[0][0] = basisFirstDerivative[q][0];
                    b[1][1] = basisFirstDerivative[q][1];
                    b[2][2] = basisFirstDerivative[q][2];
                    b[3][0] = basisFirstDerivative[q][1];
                    b[3][1] = basisFirstDerivative[q][0];
                    b[4][1] = basisFirstDerivative[q][2];
                    b[4][2] = basisFirstDerivative[q][1];
                    b[5][0] = basisFirstDerivative[q][2];
                    b[5][2] = basisFirstDerivative[q][0];
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 6, 1.0, (double *)g, 6, (double *)b, 3, 0.0, (double *)c, 3);
                    for (i=0; i<3; i++) {
                        for (j=0; j<3; j++) {
                            a[i][j] = a[i][j] + s * c[i][j];
                        }
                    }
                }
                
                if (viscNewtonLin) {
                    for (i=0; i<dim; i++) {
                        sum = 0.0;
                        for (k=0; k<3; k++) {
                            sum = sum + strain[i][k] * basisFirstDerivative[q][k];
                        }
                        muder = muder0 * 4.0 * sum;
                        for (j=0; j<dim; j++) {
                            sum = 0.0;
                            for (k=0; k<3; k++) {
                                sum = sum + strain[j][k] * basisFirstDerivative[p][k];
                            }
                            jac[j][i] = jac[j][i] + s * 2.0 * muder * sum;
                        }
                    }
                }
                
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        if (isotropic == YES) {
                            a[i][i] = a[i][i] + s * mu * basisFirstDerivative[q][j] * basisFirstDerivative[p][j];
                            if (divDiscretization == YES) {
                                a[i][j] = a[i][j] + s * mu * basisFirstDerivative[q][j] * basisFirstDerivative[p][i];
                            } else if (laplaceDiscretization == NO) {
                                a[i][j] = a[i][j] + s * mu * basisFirstDerivative[q][i] * basisFirstDerivative[p][j];
                            }
                            
                            // For compressible flow add grad((2/3) \ mu div(u))
                            if (compressible == YES) {
                                a[i][j] = a[i][j] - s * (2.0/3.0) * mu * basisFirstDerivative[q][j] * basisFirstDerivative[p][i];
                            }
                        }
                        
                        if (convect == YES) {
                            a[i][i] = a[i][i] + s * rho * basisFirstDerivative[q][j] * velo[j] * basis[p];
                        }
                    }
                    // Pressure terms
                    if (gradDriscretization == YES) {
                        a[i][c] = a[i][c] + s * basisFirstDerivative[q][i] * basis[p];
                    } else {
                        a[i][c] = a[i][c] - s * basis[q] * basisFirstDerivative[p][i];
                    }
                    
                    // Continuity equation
                    if (gradDriscretization == YES) {
                        a[c][i] = a[c][i] - s * rho * basis[q] * basisFirstDerivative[p][i];
                    } else {
                        if (compressible == YES || convect == YES) {
                            a[c][i] = a[c][i] + s * rho * basisFirstDerivative[q][i] * baseP;
                        } else {
                            a[c][i] = a[c][i] + s * basisFirstDerivative[q][i] * baseP;
                        }
                        
                        switch (compressibilityModel) {
                            case perfect_gas1:
                                a[c][i] = a[c][i] + s * (rho / pressure) * basis[q] * dPressuredx[i] * baseP / 2.0;
                                a[c][c] = a[c][c] + s * (rho / pressure) * velo[i] * basisFirstDerivative[q][i] * baseP / 2.0;
                                a[c][c] = a[c][c] - s * (rho / (temperature * pressure) ) * dTemperaturedx[i] * basis[q] * baseP;
                                break;
                            case user_defined1:
                            case thermal:
                                a[c][i] = a[c][i] + s * dRhodx[i] * basis[q] * baseP;
                                break;
                            case user_defined2:
                                a[c][c] = a[c][c] + s * drhodp * basisFirstDerivative[q][i] * velo[i] * baseP / 2.0;
                                a[c][i] = a[c][i] + s * drhodp * dPressuredx[i] * basis[q] * baseP / 2.0;
                                break;
                        }
                    }
                }
                
                // Artificial compressibility, affects only the continuity equation
                if (pseudoCompressible == YES) {
                    a[c][c] = a[c][c] + s * compress * basis[q] * baseP;
                }
                
                // Convection, Newton linearization
                if (convect == YES && newtonLinearization == YES) {
                    for (i=0; i<dim; i++) {
                        for (j=0; j<dim; j++) {
                            a[i][j] = a[i][j] + s * rho + grad[i][j] * basis[q] * basis[p];
                        }
                    }
                }
                
                // Add stabilization
                if (stabilize == YES) {
                    for (i=0; i<dim; i++) {
                        for (j=0; j<c+1; j++) {
                            m[j][i] = m[j][i] + s * tau * rho * basis[q] * sw[p][j][i];
                        }
                        for (j=0; j<dim; j++) {
                            a[j][i] = a[j][i] + s * delta * basisFirstDerivative[q][i] * basisFirstDerivative[p][j];
                        }
                    }
                    double aa[c+1][dim];
                    double bb[dim][c+1];
                    double cc[c+1][c+1];
                    for (i=0; i<c+1; i++) {
                        for (j=0; j<dim; j++) {
                            aa[i][j] = sw[p][i][j];
                        }
                    }
                    for (i=0; i<dim; i++) {
                        for (j=0; j<c1; j++) {
                            bb[i][j] = sw[q][i][j];
                        }
                    }
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, c+1, c+1, dim, 1.0, (double *)aa, dim, (double *)bb, c+1, 0.0, (double *)cc, c+1);
                    for (i=0; i<c+1; i++) {
                        for (j=0; j<c+1; j++) {
                            a[i][j] = a[i][j] + s * tau * cc[i][j];
                        }
                    }
                } else if (vms == YES) {
                    
                }
            }
        }
    }
    
    if (nodalVelo != NULL) free_dmatrix(nodalVelo, 0, 3, 0, n-1);
    if (nodalPVelo != NULL)free_dmatrix(nodalPVelo, 0, 3, 0, n-1);
}

@end
