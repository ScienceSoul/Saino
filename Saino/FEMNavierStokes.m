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
#include "Walls.h"

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
        double **loadVector                        -> loadVector vector
        double *nodalViscosity                     -> nodal values for viscosity (i.e. if turbulence model or
                                                      power-law viscosity is used, the values vary in space)
        double *nodalDensity                       -> nodal values of density
        double *velocityX, *velocityY, *velocityZ  -> nodal values of velocity components from previous iteration
        double *nodalPressure                      -> nodal values of total pressure from previous iteration
        NSString *stabilizeFlag                    -> should stabilization be used?
        BOOL pseudoCompressible                    -> should artificial compressibility be added?
        double *nodalCompressibility               -> artificial compressibility for the nodes
        BOOL magneticForce                         -> should Lorentz force for magneto-hydrodynamics be included
        BOOL rotating                              -> is the coordinate system rotating
        double *omega                              -> if previous is YES, components of angular velocity
        BOOL newtonLinearization                   -> Picard or Newton linearization of the convetion term?
        Element_t *element                         -> structure describing the element (dimension,nof nodes,
                                                      interpolation degree, etc...)
        int n                                      -> number of element nodes
        Nodes_t *nodes                             -> element node coordinates
*********************************************************************************************************************/
-(void)navierStokesComposeMassMatrix:(double **) massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double **)loadVector nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz meshVelocityX:(double *)mux meshVelocityY:(double *)muy meshVelocityZ:(double *)muz nodalPressure:(double *)nodalPressure nodalTemperature:(double *)nodalTemperature isConvect:(BOOL)convect stabilizeFlag:(NSString *)stabilizeFlag compressibilityModel:(int)compressibilityModel isPseudoCompressible:(BOOL)pseudoCompressible nodalCompressibility:(double *)nodalCompressibility nodalGasConstant:(double *)nodalGasConstant isPorous:(BOOL)porous nodalDrag:(double **)nodalDrag isPotentialForce:(BOOL)potentialForce potentialField:(double *)potentialField potentialCoefficient:(double *)potentialCoefficient isMagneticForce:(BOOL)magneticForce isRotating:(BOOL)rotating omega:(double *)omega isDivDiscretization:(BOOL)divDiscretization isGradPDriscretization:(BOOL)gradPDriscretization isNewtonLinearization:(BOOL)newtonLinearization isTransient:(BOOL)transient element:(Element_t *)element numberOfNodes:(int)n rows:(int)rows cols:(int)cols nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution core:(FEMCore *)core mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration material:(FEMMaterial *)material elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems materialModels:(FEMMaterialModels *)materialModels differentials:(FEMDifferentials *)differentials listUtilities:(FEMListUtilities *)listUtilities utilities:(FEMUtilities *)utilities {
    
    int c, i, j, k, l, p, q, t, dim, linearBasis, nBasis, order=0, tStep;
    double baseP, c1=0.0, compress=0.0, delta=0.0, detJ, drhodp=0.0, dt=0.0, gasC, hk, hScale, lambda=1.0, massCoeff, mk, mu=0.0, muder, muder0=0.0, pressure=0.0, re, rho, s, sum, sum1,
           tau=0.0, tau_c=0.0, tau_m=0.0, temperature=0.0, u, v, viscConstantCondition, vnorm, w;
    double b[6][3], coord[3], diffux[n], diffuy[n], diffuz[n], dNodalBasisdx[n][n][3], dPressuredx[3], drag[3], dRhodx[3], dTemperaturedx[3], dmudx[3], force[4], g[3][6], gmat[3][3],
           gvec[3], grad[3][3], gradNodal[n][4][3], gradP[n], gradT[3][3], lc[3][n], lrf[3], jacM[8*n][8*n], pBasis[n], pdBasisdx[n][3], prm[4], pVelo[3],
           rm[n][4][4], strain[3][3], su[n][4][4], sw[n][4][4], vc[6][6], velo[3], uVelo[3];
    double **nodalVelo = NULL, **nodalPVelo = NULL;
    double *basis = NULL, **basisFirstDerivative = NULL;
    BOOL bubbles, compressible, found, isotropic, laplaceDiscretization, pBubbles, p2p1, stat, stabilize, viscNewtonLin=NO, viscNonnewtonian, vms;
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
    vms = ([stabilizeFlag isEqualToString:@"vms"] == YES) ? YES : NO;
    
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
    
    if ((solution.solutionInfo)[@"laplace discretization"] != nil) {
        laplaceDiscretization = [(solution.solutionInfo)[@"laplace discretization"] boolValue];
    } else laplaceDiscretization = NO;
    
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
                    dNodalBasisdx[i][p][j] = basisFirstDerivative[i][j];
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
        double c[dim][dim];
        double yy[dim];
        for (p=0; p<n; p++) {
            u = element->Type.NodeU[p];
            v = element->Type.NodeV[p];
            w = element->Type.NodeW[p];
            stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
            for (i=0; i<n; i++) {
                for (j=0; j<3; j++) {
                    dNodalBasisdx[i][p][j] = basisFirstDerivative[i][j];
                }
            }
            memset( *c, 0.0, (dim*dim)*sizeof(double) );
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, n, 1.0, *nodalVelo, n, *basisFirstDerivative, 3, 0.0, (double *)c, dim);
            for (i=0; i<dim; i++) {
                for (j=0; j<dim; j++) {
                    gradNodal[p][i][j] = c[i][j];
                }
            }
            memset(yy, 0.0, sizeof(yy) );
            cblas_dgemv(CblasRowMajor, CblasNoTrans, n, dim, 1.0, *basisFirstDerivative, 3, *nodalVelo+(dim*n), 1.0, 0.0, yy, 1.0);
            for (i=0; i<dim; i++) {
                gradNodal[p][dim][i] = yy[i];
            }
        }
        
        memset( *nodalPVelo, 0.0, (4*n)*sizeof(double) );
        if (transient == YES) {
            tStep = -1;
            [core getVectorLocalField:nodalPVelo size1Field:4 size2Field:n name:nil element:element solution:solution model:model timeStep:&tStep];
            
            if (order < 2) {
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        nodalPVelo[i][j] = nodalVelo[i][j] - nodalPVelo[i][j] / dt;
                    }
                }
            } else {
                double **work = doublematrix(0, 2, 0, n-1);
                tStep = -2;
                [core getVectorLocalField:work size1Field:3 size2Field:n name:nil element:element solution:solution model:model timeStep:&tStep];
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        nodalPVelo[i][j] = (1.5 * nodalVelo[i][j]) - 2.0 * nodalPVelo[i][j] + 0.5 * work[i][j] / dt;
                    }
                }
                free_dmatrix(work, 0, 2, 0, n-1);
            }
        }
        
        c1 = 2.0 / mk;
        for (i=0; i<n; i++) {
            lc[0][i] = element->Type.NodeU[i];
            lc[1][i] = element->Type.NodeV[i];
            lc[2][i] = element->Type.NodeW[i];
        }
        
        double maxv = -HUGE_VAL;
        double minv = HUGE_VAL;
        for (i=0; i<element->Type.dimension; i++) {
            for (j=0; j<n; j++) {
                if (lc[i][j] > maxv) {
                    maxv = lc[i][j];
                }
            }
            for (j=0; j<n; j++) {
                if (lc[i][j] < minv) {
                    minv = lc[i][j];
                }
            }
            for (j=0; j<n; j++) {
                lc[i][j] = 2.0 * (lc[i][j] - minv) / (maxv - minv) - 1.0;
            }
        }
    }
    // Arrays of pointers to matrices elements
    double ***m, ***a, ***jac;
    m = (double ***)malloc( (c+1) * sizeof ( double ** ));
    for (i=0; i<(c+1); i++) {
        m[i] = (double **)malloc( (c+1) * sizeof ( double * ));
    }
    a =  (double ***)malloc( (c+1) * sizeof ( double ** ));
    for (i=0; i<(c+1); i++) {
        a[i] = (double **)malloc( (c+1) * sizeof ( double * ));
    }
    jac =  (double ***)malloc( (c+1) * sizeof ( double ** ));
    for (i=0; i<(c+1); i++) {
        jac[i] = (double **)malloc( (c+1) * sizeof ( double * ));
    }
    // Array of pointers to vector elements
    double **load;
    load = (double **)malloc( (c+1) * sizeof ( double * ));

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
            memcpy(pBasis, basis, n*sizeof(double));
            memcpy(*pdBasisdx, *basisFirstDerivative, (n*3)*sizeof(double));
            element->Type = *saveType;
        }
        
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:bubbles basisDegree:NULL];
        stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:bubbles basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
        // Density at the integration point
        rho = cblas_ddot(n, nodalDensity, 1, basis, 1);
        switch (compressibilityModel) {
            case perfect_gas1:
                memset(dPressuredx, 0.0, sizeof(dPressuredx) );
                memset(dTemperaturedx, 0.0, sizeof(dTemperaturedx) );
                if (p2p1 == YES) {
                    k = linearBasis;
                    pressure = cblas_ddot(k, nodalPressure, 1, pBasis, 1);
                    temperature = cblas_ddot(k, nodalTemperature, 1, pBasis, 1);
                    for (i=0; i<dim; i++) {
                        for (j=0; j<k; j++) {
                            dPressuredx[i] = dPressuredx[i] + nodalPressure[j] * pdBasisdx[j][i];
                            dTemperaturedx[i] = dTemperaturedx[i] + nodalTemperature[j] * pdBasisdx[j][i];
                        }
                    }
                } else {
                    pressure = cblas_ddot(n, nodalPressure, 1, basis, 1);
                    temperature = cblas_ddot(n, nodalTemperature, 1, basis, 1);
                    for (i=0; i<dim; i++) {
                        for (j=0; j<n; j++) {
                            dPressuredx[i] = dPressuredx[i] + nodalPressure[j] * basisFirstDerivative[j][i];
                            dTemperaturedx[i] = dTemperaturedx[i] + nodalTemperature[j] * basisFirstDerivative[j][i];
                        }
                    }
                }
                gasC = cblas_ddot(n, nodalGasConstant, 1, basis, 1);
                rho = pressure / (gasC * temperature);
                break;
                
            case user_defined1:
            case thermal:
                memset(dRhodx, 0.0, sizeof(dRhodx) );
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        dRhodx[i] = dRhodx[i] + nodalDensity[j] * basisFirstDerivative[j][i];
                    }
                }
                break;
                
            case user_defined2:
                memset(dPressuredx, 0.0, sizeof(dPressuredx) );
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        dPressuredx[i] = dPressuredx[i] + nodalPressure[j] * basisFirstDerivative[j][i];
                    }
                }
                drhodp = cblas_ddot(n, drhodp_n.vector, 1, basis, 1);
                break;
        }
        
        if (pseudoCompressible == YES) {
            pressure = cblas_ddot(n, nodalPressure, 1, basis, 1);
            compress = rho * cblas_ddot(n, nodalCompressibility, 1, basis, 1);
        }
        
        // Velocity from previous iteration (relative to mesh veclocity) at the
        // integration point
        memset(velo, 0.0, sizeof(velo) );
        memset(diffux, 0.0, sizeof(diffux));
        memset(diffuy, 0.0, sizeof(diffuy));
        memset(diffuz, 0.0, sizeof(diffuz));
        vDSP_vsubD(mux, 1, ux, 1, diffux, 1, n);
        vDSP_vsubD(muy, 1, uy, 1, diffuy, 1, n);
        if (dim > 2) vDSP_vsubD(muz, 1, uz, 1, diffuz, 1, n);
        velo[0] = cblas_ddot(n, diffux, 1, basis, 1);
        velo[1] = cblas_ddot(n, diffuy, 1, basis, 1);
        if (dim > 2) velo[2] = cblas_ddot(n, diffuz, 1, basis, 1);
        
        memset( *grad, 0.0, (3*3)*sizeof(double) );
        for (i=0; i<3; i++) {
            for (j=0; j<n; j++) {
                grad[0][i] = grad[0][i] + ux[j] * basisFirstDerivative[j][i];
                grad[1][i] = grad[1][i] + uy[j] * basisFirstDerivative[j][i];
                if (dim > 2)  grad[2][i] = grad[2][i] + uz[j] * basisFirstDerivative[j][i];
            }
        }
        
        //Force at integration point
        memset(force, 0.0, sizeof(force) );
        for (i=0; i<c+1; i++) {
            for (j=0; j<n; j++) {
                force[i] = force[i] + loadVector[i][j] * basis[j];
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
            coord[0] = cblas_ddot(n, nodes->x, 1, basis, 1);
            coord[1] = cblas_ddot(n, nodes->y, 1, basis, 1);
            coord[2] = cblas_ddot(n, nodes->z, 1, basis, 1);
            
            // Langranges formula is used to simplify the triple product
            // omega x (omega x coord) = omega(omega.coord) - coord(omega.omega)
            
            // This is will be multiplied by density later on
            for (i=0; i<dim; i++) {
                force[i] = force[i] - omega[i] * cblas_ddot(dim, omega, 1, coord, 1);
            }
            for (i=0; i<dim; i++) {
                force[i] = force[i] + coord[i] * cblas_ddot(dim, omega, 1, omega, 1);
            }
        }
        
        // Additional forces due to gradient forces (electrokinetic flow) and viscous flag in porous media
        if (potentialForce == YES) {
            for (i=0; i<dim; i++) {
                sum = 0.0;
                for (j=0; j<n; j++) {
                    sum = sum + potentialField[j] * basisFirstDerivative[j][i];
                }
                force[i] = force[i] - cblas_ddot(n, potentialCoefficient, 1, basis, 1) * sum;
            }
        }
        
        if (porous == YES) {
            memset(drag, 0.0, sizeof(drag) );
            for (i=0; i<dim; i++) {
                for (j=0; j<n; j++) {
                    drag[i] = drag[i] + nodalDrag[i][j] * basis[j];
                }
            }
        }
        
        if (convect == YES && newtonLinearization == YES) {
            memset(uVelo, 0.0, sizeof(uVelo) );
            uVelo[0] = cblas_ddot(n, ux, 1, basis, 1);
            uVelo[1] = cblas_ddot(n, uy, 1, basis, 1);
            if (dim > 2) uVelo[2] = cblas_ddot(n, uz, 1, basis, 1);
            
            for (i=0; i<dim; i++) {
                for (j=0; j<dim; j++) {
                    force[i] = force[i] + grad[i][j] * uVelo[j];
                }
            }
        }
        
        // Effective viscosity and derivatives at integration point
        if (isotropic == YES) {
            mu = cblas_ddot(n, nodalViscosity, 1, basis, 1);
            
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
                        vc[i][j] = vc[i][j] + tensor.tensor[i][j][k] * basis[k];
                    }
                }
            }
        }
        
        if (stabilize == YES) {
            memset(dmudx, 0.0, sizeof(dmudx) );
            for (i=0; i<3; i++) {
                for (j=0; j<n; j++) {
                    dmudx[i] = dmudx[i] + nodalViscosity[j] * basisFirstDerivative[j][i];
                }
            }
            // Stabilization parameters Tau and Delta
            if (convect == YES) {
                vDSP_svesqD(velo, 1, &sum, dim);
                vnorm = max(sqrt(sum), 1.0e-12);
                re = min(1.0, rho * mk * hk * vnorm / (4.0 * mu));
                
                tau = hk * re / (2.0 * rho * vnorm);
                delta = rho * lambda * re * hk * vnorm;
            } else {
                delta = 0.0;
                tau = mk * pow(hk, 2.0) / (8.0 * mu);
            }
            // su will contain residual of NS-equations (except for the time derivative and force terms).
            // sw will contain the weight function values
            memset( **su, 0.0, (n*4*4)*sizeof(double) );
            memset( **sw, 0.0, (n*4*4)*sizeof(double) );
            for (p=0; p<n; p++) {
                for (i=0; i<dim; i++) {
                    su[p][i][c] = su[p][i][c] + basisFirstDerivative[p][i];
                    if (porous == YES) {
                        su[p][i][i] = su[p][i][i] + mu * drag[i] * basis[p];
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
                    sum = sum + nodalPVelo[i][j] * basis[j];
                }
                pVelo[i] = rho * sum;
                memset(gradP, 0.0, sizeof(gradP) );
                for (j=0; j<n; j++) {
                    gradP[i] = gradP[i] + nodalVelo[dim][j] * basisFirstDerivative[j][i];
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
            memset(yy, 0.0, sizeof(yy) );
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, (double *)gmat, 3, velo, 1.0, 0.0, yy, 1.0);
            sum = 0.0;
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    sum = sum + gmat[i][j] * gmat[i][j];
                }
            }
            if (transient == YES) {
                tau_m = 1.0 / sqrt( cblas_ddot(3, velo, 1, yy, 1) + pow(c1, 2.0) * pow((mu/rho), 2.0) * sum / dim + 4.0 / pow(dt, 2.0) );
            } else {
                tau_m = 1.0 / sqrt( cblas_ddot(3, velo, 1, yy, 1) + pow(c1, 2.0) * pow((mu/rho), 2.0) * sum / dim );
            }
            tau_c = 1.0 / (tau_m * cblas_ddot(3, gvec, 1, gvec, 1));
            
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
                int ii = (c+1) * p;
                int jj = (c+1) * q;
                int kk = 0;
                for (k=ii; k<=ii+c; k++) {
                    int ll = 0;
                    for (l=jj; l<=jj+c; l++) {
                        m[kk][ll] = &massMatrix[k][l];
                        a[kk][ll] = &stiffMatrix[k][l];
                        if (viscNewtonLin == YES) jac[kk][ll] = &jacM[k][l];
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
                    *(m[i][i]) = *(m[i][i]) + s * rho * basis[q] * basis[p];
                }
                
                // Mass for the continuity equation (in terms of pressure)
                if (compressibilityModel == perfect_gas1) {
                    *(m[c][c]) = *(m[c][c]) + s * (rho / pressure) * basis[q] * baseP;
                } else if (compressibilityModel == user_defined2) {
                    *(m[c][c]) = *(m[c][c]) + s * drhodp * basis[q] * baseP;
                }
                
                // Stiffness matrix
                // Rotating coordinates
                if (rotating == YES) {
                    massCoeff = 2.0 * s * rho *basis[q] * basis[p];
                    *(a[0][1]) = *(a[0][1]) - massCoeff * omega[2];
                    *(a[1][0]) = *(a[1][0]) + massCoeff * omega[2];
                    if (dim == 3) {
                        *(a[0][2]) = *(a[0][2]) + massCoeff * omega[1];
                        *(a[1][2]) = *(a[1][2]) - massCoeff * omega[0];
                        *(a[2][1]) = *(a[2][1]) + massCoeff * omega[0];
                        *(a[2][0]) = *(a[2][0]) - massCoeff * omega[1];
                    }
                }
                
                // Possible Porous media effects
                if (porous == YES) {
                    for (i=0; i<dim; i++) {
                        *(a[i][i]) = *(a[i][i]) + s * mu * drag[i] * basis[q] * basis[p];
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
                            *(a[i][j]) = *(a[i][j]) + s * c[i][j];
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
                            *(jac[j][i]) = *(jac[j][i]) + s * 2.0 * muder * sum;
                        }
                    }
                }
                
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        if (isotropic == YES) {
                            *(a[i][i]) = *(a[i][i]) + s * mu * basisFirstDerivative[q][j] * basisFirstDerivative[p][j];
                            if (divDiscretization == YES) {
                                *(a[i][j]) = *(a[i][j]) + s * mu * basisFirstDerivative[q][j] * basisFirstDerivative[p][i];
                            } else if (laplaceDiscretization == NO) {
                                *(a[i][j]) = *(a[i][j]) + s * mu * basisFirstDerivative[q][i] * basisFirstDerivative[p][j];
                            }
                            
                            // For compressible flow add grad((2/3) \mu div(u))
                            if (compressible == YES) {
                                *(a[i][j]) = *(a[i][j]) - s * (2.0 / 3.0) * mu * basisFirstDerivative[q][j] * basisFirstDerivative[p][i];
                            }
                        }
                        
                        if (convect == YES) {
                            *(a[i][i]) = *(a[i][i]) + s * rho * basisFirstDerivative[q][j] * velo[j] * basis[p];
                        }
                    }
                    // Pressure terms
                    if (gradPDriscretization == YES) {
                        *(a[i][c]) = *(a[i][c]) + s * basisFirstDerivative[q][i] * basis[p];
                    } else {
                        *(a[i][c]) = *(a[i][c]) - s * basis[q] * basisFirstDerivative[p][i];
                    }
                    
                    // Continuity equation
                    if (gradPDriscretization == YES) {
                        *(a[c][i]) = *(a[c][i]) - s * rho * basis[q] * basisFirstDerivative[p][i];
                    } else {
                        if (compressible == YES || convect == YES) {
                            *(a[c][i]) = *(a[c][i]) + s * rho * basisFirstDerivative[q][i] * baseP;
                        } else {
                            *(a[c][i]) = *(a[c][i]) + s * basisFirstDerivative[q][i] * baseP;
                        }
                        
                        switch (compressibilityModel) {
                            case perfect_gas1:
                                *(a[c][i]) = *(a[c][i]) + s * (rho / pressure) * basis[q] * dPressuredx[i] * baseP / 2.0;
                                *(a[c][c]) = *(a[c][c]) + s * (rho / pressure) * velo[i] * basisFirstDerivative[q][i] * baseP / 2.0;
                                *(a[c][c]) = *(a[c][c]) - s * (rho / (temperature * pressure) ) * dTemperaturedx[i] * basis[q] * baseP;
                                break;
                            case user_defined1:
                            case thermal:
                                *(a[c][i]) = *(a[c][i]) + s * dRhodx[i] * basis[q] * baseP;
                                break;
                            case user_defined2:
                                *(a[c][c]) = *(a[c][c]) + s * drhodp * basisFirstDerivative[q][i] * velo[i] * baseP / 2.0;
                                *(a[c][i]) = *(a[c][i]) + s * drhodp * dPressuredx[i] * basis[q] * baseP / 2.0;
                                break;
                        }
                    }
                }
                
                // Artificial compressibility, affects only the continuity equation
                if (pseudoCompressible == YES) {
                    *(a[c][c]) = *(a[c][c]) + s * compress * basis[q] * baseP;
                }
                
                // Convection, Newton linearization
                if (convect == YES && newtonLinearization == YES) {
                    for (i=0; i<dim; i++) {
                        for (j=0; j<dim; j++) {
                            *(a[i][j]) = *(a[i][j]) + s * rho * grad[i][j] * basis[q] * basis[p];
                        }
                    }
                }
                
                // Add stabilization
                if (stabilize == YES) {
                    for (i=0; i<dim; i++) {
                        for (j=0; j<c+1; j++) {
                            *(m[j][i]) = *(m[j][i]) + s * tau * rho * basis[q] * sw[p][j][i];
                        }
                        for (j=0; j<dim; j++) {
                            *(a[j][i]) = *(a[j][i]) + s * delta * basisFirstDerivative[q][i] * basisFirstDerivative[p][j];
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
                        for (j=0; j<c+1; j++) {
                            bb[i][j] = su[q][i][j];
                        }
                    }
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, c+1, c+1, dim, 1.0, (double *)aa, dim, (double *)bb, c+1, 0.0, (double *)cc, c+1);
                    for (i=0; i<c+1; i++) {
                        for (j=0; j<c+1; j++) {
                            *(a[i][j]) = *(a[i][j]) + s * tau * cc[i][j];
                        }
                    }
                } else if (vms == YES) {
                    for (i=0; i<dim; i++) {
                        // (rho*u', grad(q))
                        *(m[dim][i]) = *(m[dim][i]) + s * rho * tau_m * rho * basis[q] * basisFirstDerivative[p][i];
                        for (k=0; k<dim+1; k++) {
                            *(a[dim][k]) = *(a[dim][k]) + s * rho * tau_m * rm[q][i][k] * basisFirstDerivative[p][i];
                        }
                        
                        for (j=0; j<dim; j++) {
                            // -(rho*u'*grad(u), w)
                            *(m[i][j]) = *(m[i][j]) - s * rho *tau_m * rho * basis[q] * grad[i][j] * basis[p];
                            for (k=0; k<dim+1; k++) {
                                *(a[i][k]) = *(a[i][k]) - s * rho * tau_m * rm[q][j][k] * grad[i][j] * basis[p] / 2.0;
                            }
                            *(a[i][i]) = *(a[i][i]) - s * rho * tau_m * prm[j] * basisFirstDerivative[q][j] * basis[p] / 2.0;
                            *(a[i][i]) = *(a[i][i]) + s * rho * tau_m * rho * force[j] * basisFirstDerivative[q][j] * basis[p];
                            
                            // (rho*u', u.grad(w))
                            *(m[i][i]) = *(m[i][i]) + s * rho * tau_m * rho * basis[q] * velo[j] * basisFirstDerivative[p][j];
                            for (k=0; k<dim+1; k++) {
                                *(a[i][k]) = *(a[i][k]) + s * rho * tau_m * rm[q][i][k] * velo[j] * basisFirstDerivative[p][j] / 2.0;
                            }
                            *(a[i][j]) = *(a[i][j]) + s * rho * tau_m * prm[i] * basis[q] * basisFirstDerivative[p][j] / 2.0;
                            *(a[i][j]) = *(a[i][j]) - s * rho * tau_m * rho * force[i] * basis[q] * basisFirstDerivative[p][j];
                        }
                        
                        // (rho*div(u'), div(w))
                        for (j=0; j<dim; j++) {
                            *(a[i][j]) = *(a[i][j]) + s * rho * tau_c * rm[q][dim][j] * basisFirstDerivative[p][i];
                        }
                        
                        // -(rho*u'*u', grad(w))
                        for (j=0; j<dim; j++) {
                            for (k=0; k<dim+1; k++) {
                                *(a[i][k]) = *(a[i][k]) - s * rho * pow(tau_m, 2.0) * rm[q][i][k] * prm[j] * basisFirstDerivative[p][j] / 2.0;
                                *(a[i][k]) = *(a[i][k]) - s * rho * pow(tau_m, 2.0) * prm[i] * rm[q][j][k] * basisFirstDerivative[p][j] / 2.0;
                                
                                *(a[i][k]) = *(a[i][k]) + s * rho * pow(tau_m, 2.0) * rm[q][i][k] * rho * force[j] * basisFirstDerivative[p][j];
                                *(a[i][k]) = *(a[i][k]) + s * rho * pow(tau_m, 2.0) * rho * force[i] * rm[q][j][k] * basisFirstDerivative[p][j];
                            }
                            
                            *(m[i][i]) = *(m[i][i]) - s * rho * pow(tau_m, 2.0) * rho * basis[q] * prm[j] * basisFirstDerivative[p][j];
                            *(m[i][j]) = *(m[i][j]) - s * rho * pow(tau_m, 2.0) * prm[i] * rho * basis[q] * basisFirstDerivative[p][j];
                            
                            *(m[i][i]) = *(m[i][i]) + s * rho * pow(tau_m, 2.0) * rho * basis[q] * rho * force[j] * basisFirstDerivative[p][j];
                            *(m[i][j]) = *(m[i][j]) + s * rho * pow(tau_m, 2.0) * rho * basis[q] * rho * force[i] * basisFirstDerivative[p][j];
                            
                            *(m[i][i]) = *(m[i][i]) - s * rho * pow(tau_m, 2.0) * rho * basis[q] * pVelo[j] * basisFirstDerivative[p][j];
                            *(m[i][j]) = *(m[i][j]) - s * rho * pow(tau_m, 2.0) * pVelo[i] * rho * basis[q] * basisFirstDerivative[p][j];
                        }
                    }
                }
            }
        }
        
        // The right hand side...
        for (p=0; p<nBasis; p++) {
            int ii = (c+1) * p;
            int kk = 0.0;
            for (k=ii; k<=ii+c; k++) {
                load[kk] = &forceVector[k];
                kk++;
            }
            
            for (i=0; i<c+1; i++) {
                *(load[i]) = *(load[i]) + s * rho * force[i] * basis[p];
            }
            if (compressibilityModel == perfect_gas1) *(load[c]) = *(load[c]) / temperature;
            
            if (pseudoCompressible == YES) {
                *(load[c]) = *(load[c]) + s * pressure * basis[p] * compress;
            }
            
            if (stabilize == YES) {
                for (i=0; i<dim; i++) {
                    for (j=0; j<c+1; j++) {
                        *(load[j]) = *(load[j]) + s * tau * rho * force[i] * sw[p][j][i];
                    }
                }
            } else if (vms == YES) {
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        *(load[i]) = *(load[i]) - s * rho * pow(tau_m, 2.0) * rho * force[i] * rho * force[j] * basisFirstDerivative[p][j];
                    }
                    *(load[dim]) = *(load[dim]) + s * rho * tau_m * rho * force[i] * basisFirstDerivative[p][i];
                }
            }
        }
    }
    free(m);
    free(a);
    free(jac);
    free(load);
    
    if (viscNewtonLin == YES) {
        double sol[8*n];
        memset(sol, 0.0, sizeof(sol) );
        int kk = 0;
        for (i=0; i<(c+1)*n; i+=(c+1)) {
            sol[i] = ux[kk];
            kk++;
        }
        kk = 0;
        for (i=1; i<(c+1)*n; i+=(c+1)) {
            sol[i] = uy[kk];
            kk++;
        }
        if (dim > 2) {
            kk = 0;
            for (i=2; i<(c+1)*n; i+=(c+1)) {
                sol[i] = uz[kk];
                kk++;
            }
        }
        p = (c+1) * nBasis;
        for (i=0; i<p; i++) {
            for (j=0; j<p; j++) {
                stiffMatrix[i][j] = stiffMatrix[i][j] + jacM[i][j];
            }
        }
        double yy[p];
        memset(yy, 0.0, sizeof(yy) );
        cblas_dgemv(CblasRowMajor, CblasNoTrans, p, p, 1.0, (double *)jacM, 8*n, sol, 1, 0.0, yy, 1);
        for (i=0; i<p; i++) {
            forceVector[i] = forceVector[i] + yy[i];
        }
    }
    
    if (p2p1 == YES) {
        int size;
        int **edgeMap = NULL;
        j = [core getElementFamily:element];
        edgeMap = [core getEdgeMap:j mapSize:&size];
        for (i=j+1; i<=j+size; i++) {
            p = edgeMap[(i-j)-1][0];
            q = edgeMap[(i-j)-1][1];
            for (k=0; k<cols; k++) {
                stiffMatrix[((c+1)*i)-1][k] = 0.0;
                massMatrix[((c+1)*i)-1][k] = 0.0;
            }
            forceVector[((c+1)*i)-1] = 0.0;
            stiffMatrix[((c+1)*i)-1][((c+1)*i)-1] = 1.0;
            stiffMatrix[((c+1)*i)-1][((c+1)*(p+1))-1] = -0.5;
            stiffMatrix[((c+1)*i)-1][((c+1)*(q+1))-1] = -0.5;
        }
    }
    
    if (pBubbles == YES) {
        for (i=n+1; i<=nBasis; i++) {
            for (k=0; k<cols; k++) {
                stiffMatrix[((c+1)*i)-1][k] = 0.0;
                massMatrix[((c+1)*i)-1][k] = 0.0;
            }
            for (k=0; k<rows; k++) {
                stiffMatrix[k][((c+1)*i)-1] = 0.0;
                massMatrix[k][((c+1)*i)-1] = 0.0;
            }
            forceVector[((c+1)*i)-1] = 0.0;
            stiffMatrix[((c+1)*i)-1][((c+1)*i)-1] = 1.0;
        }
    }
    
    if (nodalVelo != NULL) free_dmatrix(nodalVelo, 0, 3, 0, n-1);
    if (nodalPVelo != NULL) free_dmatrix(nodalPVelo, 0, 3, 0, n-1);
    
    if (tensor.tensor != NULL) {
        free_d3tensor(tensor.tensor, 0, tensor.m-1, 0, tensor.n-1, 0, tensor.p-1);
    }
    if (drhodp_n.vector != NULL) {
        free_dvector(drhodp_n.vector, 0, drhodp_n.m-1);
    }
}

/*****************************************************************************************************************************
    Return element local matrices and RHS vector for Navier-Stokes equations boundary conditions in cartesian coordinates.
 
     Arguments:
 
        double **boundaryMatrix     -> output: time derivative coefficient matrix
        double *boundaryVector      -> output: RHS vector
        double **loadVector         -> nodal values force in coordinate directives
        double *nodalAlpha          -> nodal values for force in normal direction
        double *nodalBeta           -> nodal values of something which will be taken derivative in
                                       tangential direction and added to force
        Element_t *element          -> structure describing the element (dimension,nof nodes,
                                       interpolation degree, etc...)
        int n                       -> number of boundary element nodes
        Nodes_t *nodes              -> element node coordinates
*****************************************************************************************************************************/
-(void)navierStokesBoundary:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector loadVector:(double **)loadVector nodalAlpha:(double *)nodalAlpha nodalBeta:(double *)nodalBeta nodalExtPressure:(double *)nodalExtPressure nodalSlipCoefficient:(double **)nodalSlipCoefficient isNormalTangential:(BOOL)normalTangential element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh  model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription elementUtils:(FEMElementUtils *)elementUtils {
    
    int c, i, j, k, l, p, q, t, dim;
    double alpha, detJ, massFlux, slipCoeff, s, u, v, w;
    double force[3], normals[3], tangents[3], tangents2[3], tangentForce[3], vect[3];
    double *basis = NULL, **basisFirstDerivative = NULL;
    BOOL stat;
    GaussIntegrationPoints *IP = NULL;
    
    dim = model.dimension;
    c = dim + 1;
    
    basis = integration.basis;
    basisFirstDerivative = integration.basisFirstDerivative;
    
    IP = GaussQuadrature(element, NULL, NULL);
    
    // Start integrating
    for (t=0; t<IP->n; t++) {
        
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];
        
        // Basis function values and derivatives at the integration point
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
        // Add to load: tangential derivatives of something
        memset(tangentForce, 0.0, sizeof(tangentForce) );
        for (i=0; i<dim; i++) {
            for (j=0; j<n; j++) {
                tangentForce[i] = tangentForce[i] + nodalBeta[j] * basisFirstDerivative[j][i];
            }
        }
        
        // Add to load:  given force in coordinate directions
        memset(force, 0.0, sizeof(force) );
        for (i=0; i<dim; i++) {
            for (j=0; j<n; j++) {
                force[i] = force[i] + loadVector[i][j] * basis[j];
            }
        }
        
        // Add to load: given force in normal direction
        BOOL check = YES;
        [elementDescription normalVectorForBDElement:element boundaryNodes:nodes mesh:mesh paraU:&u paraV:&v check:&check normals:normals];
        alpha = cblas_ddot(n, nodalExtPressure, 1, basis, 1);
        if (normalTangential == YES) {
            force[0] = force[0] + alpha;
        } else {
            for (i=0; i<dim; i++) {
                force[i] = force[i] + alpha * normals[i];
            }
        }
        
        alpha = cblas_ddot(n, nodalAlpha, 1, basis, 1);
        
        massFlux = 0.0;
        for (i=0; i<n; i++) {
            massFlux = massFlux + loadVector[3][i] * basis[i];
        }
        
        switch (element->Type.dimension) {
            case 1:
                tangents[0] = normals[1];
                tangents[1] = -normals[0];
                tangents[2] = 0.0;
                memset(tangents2, 0.0, sizeof(tangents2) );
                break;
            case 2:
                [elementUtils tangentDirectionsForNormal:normals tangent1:tangents tangent2:tangents2];
                break;
        }
        
        BOOL any = NO;
        for (i=0; i<3; i++) {
            for (j=0; j<n; j++) {
                if (nodalSlipCoefficient[i][j] != 0.0) {
                    any = YES;
                    break;
                }
            }
        }
        if (any == YES) {
            for (p=0; p<n; p++) {
                for (q=0; q<n; q++) {
                    for (i=0; i<dim; i++) {
                        slipCoeff = 0.0;
                        for (j=0; j<n; j++) {
                            slipCoeff = slipCoeff + nodalSlipCoefficient[i][j] * basis[j];
                        }
                        
                        if (normalTangential == YES) {
                            switch (i) {
                                case 0:
                                    memcpy(vect, normals, sizeof(normals));
                                    break;
                                case 1:
                                    memcpy(vect, tangents, sizeof(tangents));
                                    break;
                                case 2:
                                    memcpy(vect, tangents2, sizeof(tangents2));
                                    break;
                            }
                            
                            for (j=0; j<dim; j++) {
                                for (k=0; k<dim; k++) {
                                    boundaryMatrix[p*c+j][q*c+k] = boundaryMatrix[p*c+j][q*c+k] + s * slipCoeff * basis[q] * basis[p] * vect[j] * vect[k];
                                }
                            }
                        } else {
                            boundaryMatrix[p*c+i][q*c+i] = boundaryMatrix[p*c+i][q*c+i] + s * slipCoeff * basis[q] * basis[p];
                        }
                    }
                    
                }
            }
        }
        
        for (q=0; q<n; q++) {
            for (i=0; i<dim; i++) {
                k = q * c + i;
                if (normalTangential == YES) {
                    switch (i) {
                        case 0:
                            memcpy(vect, normals, sizeof(normals));
                            break;
                        case 1:
                            memcpy(vect, tangents, sizeof(tangents));
                            break;
                        case 2:
                            memcpy(vect, tangents2, sizeof(tangents2));
                            break;
                    }
                    for (j=0; j<dim; j++) {
                        l = q * c + j;
                        boundaryVector[l] = boundaryVector[l] + s * basis[q] * force[i] * vect[j];
                    }
                } else {
                    boundaryVector[k] = boundaryVector[k] + s * basis[q] * force[i];
                }
                boundaryVector[k] = boundaryVector[k] - s * alpha * basisFirstDerivative[q][i];
                boundaryVector[k] = boundaryVector[k] + s * tangentForce[i] * basis[q];
            }
            boundaryVector[((q+1)*c)-1] = boundaryVector[((q+1)*c)-1] + s * massFlux * basis[q];
        }
    }
}

/*****************************************************************************************************************************
 
        Arguments:
 
            double **boundaryMatrix     -> output: time derivative coefficient matrix
            double *boundaryVector      -> output: RHS vector
            double *layerThickness      -> boundary layer thickness
            double *surfaceRoughness    -> measure of surface roughness
            double *nodalViscosity      -> nodal values of viscosity
            double *nodalDensity        -> nodal values of density
            double *ux, *uy, *uz         -> nodal values of velocity from previous iteration
            Element_t *element          -> structure describing the element (dimension,nof nodes,
                                           interpolation degree, etc...)
            int n                       -> number of boundary element nodes
            Nodes_t *nodes              -> element node coordinates
 *****************************************************************************************************************************/
-(void)vmsWallsBoundary:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector layerThickness:(double *)layerThickness surfaceRoughness:(double *)surfaceRoughness nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes solution:(FEMSolution *)solution mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription elementUtils:(FEMElementUtils *)elementUtils {
    
    int c, i, j, k, k1, k2, p, q, t, dim;
    int *ind;
    double detJ, dkerr, dist, dfx, frictionVelocity, h, mu, rho, roughness, s, sum, u, v, vabs, x, y, w, z;
    double normals[3], tangents[3], tangents2[3], tangentialVelocity[2], velo[3];
    double *basis = NULL, **basisFirstDerivative = NULL;
    BOOL stat;
    Nodes_t *meshNodes = NULL;
    GaussIntegrationPoints *IP = NULL;
    
    if ([solution.solutionInfo[@"stabilization method"] isEqualToString:@"vms"] == NO) return;
    
    dim = model.dimension;
    c = dim + 1;
    
    basis = integration.basis;
    basisFirstDerivative = integration.basisFirstDerivative;
    
    ind = element->BoundaryInfo->Left->NodeIndexes;
    k = 0;
    x = 0.0; y = 0.0; z = 0.0;
    for (i=0; i<element->BoundaryInfo->Left->sizeNodeIndexes; i++) {
        for (j=0; j<n; j++) {
            if (ind[i] == element->NodeIndexes[j]) break;
        }
        if (j >= n) {
            meshNodes = mesh.getNodes;
            x = x + meshNodes->x[ind[i]];
            y = y + meshNodes->y[ind[i]];
            z = z + meshNodes->z[ind[i]];
            k++;
        }
    }
    
    if (k > 0) {
        vDSP_sveD(nodes->x, 1, &sum, n);
        x = x / k - sum / n;
        
        vDSP_sveD(nodes->y, 1, &sum, n);
        y = y / k - sum / n;
        
        vDSP_sveD(nodes->z, 1, &sum, n);
        z = z / k - sum / n;
    }
    h = sqrt(pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0));
    
    IP = GaussQuadrature(element, NULL, NULL);
    
    // Start integrating
    for (t=0; t<IP->n; t++) {
        
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];
        
        // Basis function values and derivatives at the integration point
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
        // Density and viscosity at integration point
        rho = cblas_ddot(n, nodalDensity, 1, basis, 1);
        mu = cblas_ddot(n, nodalViscosity, 1, basis, 1);
        
        // Velocity from previous iteration at the integration point
        memset(velo, 0.0, sizeof(velo) );
        velo[0] = cblas_ddot(n, ux, 1, basis, 1);
        velo[1] = cblas_ddot(n, uy, 1, basis, 1);
        if (dim > 2) velo[2] = cblas_ddot(n, uz, 1, basis, 1);
        
        // Normal and tangent directions
        BOOL check = NO;
        [elementDescription normalVectorForBDElement:element boundaryNodes:nodes mesh:mesh paraU:&u paraV:&v check:&check normals:normals];
        
        if (dim <= 2) {
            tangents[0] = normals[1];
            tangents[1] = -normals[0];
            tangents[2] = 0.0;
        } else {
             [elementUtils tangentDirectionsForNormal:normals tangent1:tangents tangent2:tangents2];
        }
        memset(tangentialVelocity, 0.0, sizeof(tangentialVelocity) );
        tangentialVelocity[0] = cblas_ddot(dim, velo, 1, tangents, 1);
        if (dim == 3) {
            tangentialVelocity[1] = cblas_ddot(dim, velo, 1, tangents2, 1);
        }
        
        dist = cblas_ddot(n, layerThickness, 1, basis, 1);
        if (dist == 0.0) dist = h;
        
        roughness = cblas_ddot(n, surfaceRoughness, 1, basis, 1);
        
        // Solve friction velocity and its derivative with respect to the
        // tangential velocity
        frictionVelocity = 0.0;
        dkerr = 0.0;
        vDSP_svesqD(velo, 1, &sum, 3);
        vabs = max(sqrt(sum), 1.0e-08);
        solve_ufric(rho, mu, dist, roughness, vabs, &frictionVelocity, &dfx);
        dkerr = rho * pow(frictionVelocity, 2.0);
        
        
        for (p=0; p<n; p++) {
            for (q=0; q<n; q++) {
                for (i=0; i<dim; i++) {
                    k1 = c * p + i;
                    k2 = c * q + i;
                    sum = 0.0;
                    for (j=0; j<3; j++) {
                        sum = sum + basisFirstDerivative[q][j] * normals[j];
                    }
                    boundaryMatrix[k1][k2] = boundaryMatrix[k1][k2] - s * 2.0 * mu * sum *basis[p];
                    sum = 0.0;
                    for (j=0; j<3; j++) {
                        sum = sum + basisFirstDerivative[p][j] * normals[j];
                    }
                    boundaryMatrix[k1][k2] = boundaryMatrix[k1][k2] - s * 2.0 * mu * sum * basis[q];
                }
            }
        }
        
        for (p=0; p<n; p++) {
            for (q=0; q<n; q++) {
                for (i=0; i<dim; i++) {
                    k1 = c * p + i;
                    k2 = c * q + i;
                    boundaryMatrix[k1][k2] = boundaryMatrix[k1][k2] + s * dkerr * basis[q] * basis[p];
                }
            }
        }
    }
}

/*****************************************************************************************************************************
    Return the local matrices and RHS contribution from the wall law
 
    Arguments:
 
        double **boundaryMatrix     -> output: time derivative coefficient matrix
        double *boundaryVector      -> output: RHS vector
        double *layerThickness      -> boundary layer thickness
        double *surfaceRoughness    -> measure of surface roughness
        double *nodalViscosity      -> nodal values of viscosity
        double *nodalDensity        -> nodal values of density
        double *ux, *uy, *uz         -> nodal values of velocity from previous iteration
        Element_t *element          -> structure describing the element (dimension,nof nodes,
                                       interpolation degree, etc...)
        int n                       -> number of boundary element nodes
        Nodes_t *nodes              -> element node coordinates
 *****************************************************************************************************************************/
-(void)navierStokesWallLawBoundary:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector layerThickness:(double *)layerThickness surfaceRoughness:(double *)surfaceRoughness nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription elementUtils:(FEMElementUtils *)elementUtils {
    
    int c, i, j, k1, k2, p, q, t, dim;
    double detJ, dist, mu, rho, roughness, s, sum, u, v, vabs, w;
    double dfx[2], dkerr[2], frictionVelocity[2], normals[3], tangents[3], tangents2[3], tangentialVelocity[2], velo[3];
    double *basis = NULL;
    BOOL stat;
    GaussIntegrationPoints *IP = NULL;
    
    dim = model.coordinates;
    c = dim + 1;
    
    basis = integration.basis;
    
    IP = GaussQuadrature(element, NULL, NULL);
    
    // Start integrating
    for (t=0; t<IP->n; t++) {
        
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];
        
        // Basis function values and derivatives at the integration point
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        s = detJ * IP->s[t];
        
        // Density and viscosity at integration point
        rho = cblas_ddot(n, nodalDensity, 1, basis, 1);
        mu = cblas_ddot(n, nodalViscosity, 1, basis, 1);
        
        // Velocity from previous iteration at the integration point
        memset(velo, 0.0, sizeof(velo) );
        velo[0] = cblas_ddot(n, ux, 1, basis, 1);
        velo[1] = cblas_ddot(n, uy, 1, basis, 1);
        if (dim > 2) velo[2] = cblas_ddot(n, uz, 1, basis, 1);
        
        vDSP_svesqD(velo, 1, &sum, 3);
        vabs = max(1.0e-09, sqrt(sum));
        
        // Normal and tangent directions
        BOOL check = NO;
        [elementDescription normalVectorForBDElement:element boundaryNodes:nodes mesh:mesh paraU:&u paraV:&v check:&check normals:normals];
        
        if (dim <= 2) {
            tangents[0] = normals[1];
            tangents[1] = -normals[0];
            tangents[2] = 0.0;
        } else {
            [elementUtils tangentDirectionsForNormal:normals tangent1:tangents tangent2:tangents2];
        }
        
        memset(tangentialVelocity, 0.0, sizeof(tangentialVelocity) );
        tangentialVelocity[0] = cblas_ddot(dim, velo, 1, tangents, 1);
        if (tangentialVelocity[0] < 0) {
            vDSP_vnegD(tangents, 1, tangents, 1, 3);
            tangentialVelocity[0] = -tangentialVelocity[0];
        }
        
        if (dim == 3) {
            tangentialVelocity[1] = cblas_ddot(dim, velo, 1, tangents2, 1);
            if (tangentialVelocity[1] < 0) {
                vDSP_vnegD(tangents2, 1, tangents2, 1, 3);
                tangentialVelocity[1] = -tangentialVelocity[1];
            }
        }
        
        dist = cblas_ddot(n, layerThickness, 1, basis, 1);
        roughness = cblas_ddot(n, surfaceRoughness, 1, basis, 1);
        
        // Solve friction velocity and its derivative with respect to
        // the tangential velocity
        memset(frictionVelocity, 0.0, sizeof(frictionVelocity) );
        memset(dkerr, 0.0, sizeof(dkerr) );
        if (tangentialVelocity[0] > 1.0e-09) {
            solve_ufric(rho, mu, dist, roughness, tangentialVelocity[0], &frictionVelocity[0], &dfx[0]);
            dkerr[0] = 2.0 * rho * frictionVelocity[0] / dfx[0];
        }
        
        if (dim == 3) {
            if (tangentialVelocity[1] > 1.0e-09) {
                solve_ufric(rho, mu, dist, roughness, tangentialVelocity[1], &frictionVelocity[1], &dfx[1]);
                dkerr[1] = 2.0 * rho * frictionVelocity[1] / dfx[1];
            }
        }
        
        for (p=0; p<n; p++) {
            for (q=0; q<n; q++) {
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        k1 = p * c + i;
                        k2 = q * c + j;
                        boundaryMatrix[k1][k2] = boundaryMatrix[k1][k2] + s * dkerr[0] * tangents[i] * tangents[j] * basis[q] * basis[p];
                        if (dim == 3) {
                            boundaryMatrix[k1][k2] = boundaryMatrix[k1][k2] + s * dkerr[1] * tangents2[i] * tangents2[j] * basis[q] * basis[p];
                        }
                    }
                }
            }
        }
        
        for (q=0; q<n; q++) {
            for (i=0; i<dim; i++) {
                k1 = q * c + i;
                boundaryVector[k1] = boundaryVector[k1] + s * (dkerr[0] * tangentialVelocity[0] - rho * pow(frictionVelocity[0], 2.0)) * tangents[i] * basis[q];
                if (dim == 3) {
                    boundaryVector[k1] = boundaryVector[k1] + s * (dkerr[1] * tangentialVelocity[1] - rho * pow(frictionVelocity[1], 2.0))
                    * tangents2[i] * basis[q];
                }
            }
        }
    }
}

@end
