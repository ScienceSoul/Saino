//
//  FEMNavierStokesGeneral.m
//  Saino
//
//  Created by Seddik hakime on 28/05/2014.
//  Copyright (c) 2014 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMNavierStokesGeneral.h"

@implementation FEMNavierStokesGeneral

/*********************************************************************************************************************
    Return element local matrices and RHS vector for Navier-Stokes-Equations in cartesian coordinates.
 
    Arguments:
 
        double **massMatrix                        -> output: time derivative coefficient matrix
        double **stiffMatrix                       -> output: rest of the equation coefficients
        double *forceVector                        -> output: RHS vector
        double **loadVector                         > loadVector vector
        double *nodalViscosity                     -> nodal values for viscosity (i.e. if turbulence model or
                                                      power-law viscosity is used, the values vary in space)
        double *nodalDensity                       -> nodal values of density
        double *velocityX, *velocityY, *velocityZ  -> nodal values of velocity components from previous iteration
        BOOL stabilize                             -> should stabilization be used?
        BOOL newtonLinearization                   -> Picard or Newton linearization of the convetion term ?
        Element_t *element                         -> structure describing the element (dimension,nof nodes,
                                                      interpolation degree, etc...)
        int n                                      -> number of element nodes
        Nodes_t *nodes                             -> element node coordinates
*********************************************************************************************************************/
-(void)navierStokesGeneralComposeMassMatrix:(double **) massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double **)loadVector nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz meshVelocityX:(double *)mux meshVelocityY:(double *)muy meshVelocityZ:(double *)muz isStabilize:(BOOL)stabilize isNewtonLinearization:(BOOL)newtonLinearization element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes core:(FEMCore *)core mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems materialModels:(FEMMaterialModels *)materialModels differentials:(FEMDifferentials *)differentials listUtilities:(FEMListUtilities *)listUtilities utilities:(FEMUtilities *)utilities {
    
    int c, i, j, k, l, m, p, q, t, coordinates, dim;
    int imap[3] = {0, 1, 3};
    double delta=0.0, density, detJ, hk=0.0, lambda=1.0, mk=0.0, re, s, sqrtMetric, sum, tau=0.0, u, v, viscosity, vNorm, x, y, w, z;
    double a[4][4], dSymb[3][3][3][3], dVelodx[3][3], dViscositydx[3], force[4], load[4], lrf[3], mass[4][4], metric[3][3], symb[3][3][3], velo[3],
           uVelo[3], su[n][4][4], sw[n][4][4];
    double *basis = NULL, **basisFirstDerivative = NULL, ***basisSecondDerivative = NULL;
    BOOL cylindricSymmetry, stat;
    GaussIntegrationPoints *IP = NULL;
    
    coordinates = coordinateSystems.coordinates;
    cylindricSymmetry = (coordinates == cylindric_symmetric || coordinates == axis_symmetric) ? YES : NO;
    
    if (cylindricSymmetry) {
        dim = 3;
    } else {
        dim = model.dimension;
    }
    c = dim;
    
    basis = integration.basis;
    basisFirstDerivative = integration.basisFirstDerivative;
    basisSecondDerivative = integration.basisSecondDerivative;

    IP = GaussQuadrature(element, NULL, NULL);
    
    if (stabilize) {
        // Stabilization parameters: hk, mk (Franca et al.)
        hk = element->hK;
        mk = element->StabilizationMK;
    }
    
    // Start integrating
    for (t=0; t<IP->n; t++) {
        
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];
        
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        
        // Coordinate system dependent info
        x = 0.0;
        y = 0.0;
        z = 0.0;
        if (coordinates != cartesian) {
            x = cblas_ddot(n, nodes->x, 1, basis, 1);
            y = cblas_ddot(n, nodes->y, 1, basis, 1);
            z = cblas_ddot(n, nodes->z, 1, basis, 1);
        }
        
        [coordinateSystems coordinateSystemInfoModel:model metric:metric sqrtMetric:&sqrtMetric symbols:symb dSymbols:dSymb coordX:x coordY:y coordZ:z];
        s = sqrtMetric * detJ * IP->s[t];
        
        // Density at the integration point
        density = cblas_ddot(n, nodalDensity, 1, basis, 1);
        
        // Velocity from previous iteration at the integration point
        memset(velo, 0.0, sizeof(velo) );
        double diffux[n];
        double diffuy[n];
        double diffuz[n];
        vDSP_vsubD(mux, 1, ux, 1, diffux, 1, n);
        vDSP_vsubD(muy, 1, uy, 1, diffuy, 1, n);
        if (dim > 2 && coordinates != axis_symmetric) vDSP_vsubD(muz, 1, uz, 1, diffuz, 1, n);
        velo[0] = cblas_ddot(n, diffux, 1, basis, 1);
        velo[1] = cblas_ddot(n, diffuy, 1, basis, 1);
        if (dim > 2 && coordinates != axis_symmetric) velo[2] = cblas_ddot(n, diffuz, 1, basis, 1);
        
        if (newtonLinearization == YES) {
            memset(uVelo, 0.0, sizeof(uVelo) );
            uVelo[0] = cblas_ddot(n, ux, 1, basis, 1);
            uVelo[1] = cblas_ddot(n, uy, 1, basis, 1);
            if (dim > 2 && coordinates != axis_symmetric) uVelo[2] = cblas_ddot(n, uz, 1, basis, 1);
            
            memset( *dVelodx, 0.0, (3*3)*sizeof(double) );
            for (i=0; i<3; i++) {
                for (j=0; j<n; j++) {
                    dVelodx[0][i] = dVelodx[0][i] + ux[j] * basisFirstDerivative[j][0];
                    dVelodx[1][i] = dVelodx[1][i] + uy[j] * basisFirstDerivative[j][1];
                    if (dim > 2 && coordinates != axis_symmetric) dVelodx[2][i] = dVelodx[2][i] + uz[j] * basisFirstDerivative[j][2];
                }
            }
        }
        
        // Force at the integration point
        memset(lrf, 0.0, sizeof(lrf) );
        [differentials lorentzForceElement:element nodes:nodes numberOfNodes:n integrationU:u integrationV:v integrationW:w lorentzForce:lrf mesh:mesh model:model integration:integration coordinateSystems:coordinateSystems listUtilities:listUtilities utilities:utilities];
        
        memset(force, 0.0, sizeof(force) );
        for (i=0; i<n; i++) {
            force[0] = force[0] + loadVector[0][i] * basis[i];
            force[1] = force[1] + loadVector[1][i] * basis[i];
            force[2] = force[2] + loadVector[2][i] * basis[i];
        }
        if (dim > 2 && coordinates != axis_symmetric) {
            for (i=0; i<n; i++) {
                force[3] = force[3] + loadVector[3][i] * basis[i];
            }
            for (i=0; i<3; i++) {
                force[i] = force[i] + lrf[i] / density;
            }
        } else {
            for (i=0; i<2; i++) {
                force[i] = force[i] + lrf[i] / density;
            }
        }
        
        // Effective viscosity and derivatives at integration point
        viscosity = cblas_ddot(n, nodalViscosity, 1, basis, 1);
        viscosity = [materialModels effectiveViscosity:viscosity density:density velocityX:ux velocitY:uy velocityZ:uz element:element nodes:nodes numberOfNodes:n numberOfPoints:n integrationU:u integrationV:v integrationW:w muder:NULL mesh:mesh model:model integration:integration];
        
        // Stabilization parameters tau and delta
        if (stabilize == YES) {
            memset(dViscositydx, 0.0, sizeof(dViscositydx) );
            for (i=0; i<3; i++) {
                for (j=0; j<n; j++) {
                    dViscositydx[i] = dViscositydx[i] + nodalViscosity[j] * basisFirstDerivative[j][i];
                }
            }
            
            vDSP_svesqD(velo, 1, &sum, dim);
            vNorm = max(sqrt(sum), 1.0e-12);
            re = min(1.0, density * mk * hk * vNorm / (4.0 * viscosity));
            
            tau = 0.0;
            if (vNorm != 0.0) {
                tau = hk * re / (2.0 * density * vNorm);
            }
            
            delta = density * lambda * vNorm * re * hk;
            
            // su will contain residual of NS-equations (except for the time derivative and force terms).
            // sw will contain the weight function values
            memset( **su, 0.0, (n*4*4)*sizeof(double) );
            memset( **sw, 0.0, (n*4*4)*sizeof(double) );
            for (p=0; p<n; p++) {
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        su[p][i][c] = su[p][i][c] + metric[i][j] * basisFirstDerivative[p][j];
                        su[p][i][i] = su[p][i][i] + density * basisFirstDerivative[p][j] * velo[j];
                        if (newtonLinearization == YES) {
                            su[p][i][j] = su[p][i][j] + density * basisFirstDerivative[i][j] * basis[p];
                        }
                        for (k=0; k<dim; k++) {
                            su[p][i][i] = su[p][i][i] - viscosity * metric[j][k] * basisSecondDerivative[p][j][k];
                            su[p][i][i] = su[p][i][i] - dViscositydx[j] * metric[j][k] * basisFirstDerivative[p][k];
                            su[p][i][j] = su[p][i][j] - viscosity * metric[i][k] * basisSecondDerivative[p][j][k];
                            su[p][i][j] = su[p][i][j] - dViscositydx[j] * metric[i][k] * basisFirstDerivative[p][k];
                            if (coordinates != cartesian) {
                                su[p][i][k] = su[p][i][k] + density * symb[k][j][i] * basis[p] * velo[j];
                                if (newtonLinearization == YES) {
                                    su[p][i][j] = su[p][i][j] + density * symb[k][j][i] * velo[k] * basis[p];
                                }
                                for (l=0; l<dim; l++) {
                                    su[p][i][i] = su[p][i][i] + viscosity * metric[j][k] * symb[j][k][l] * basisFirstDerivative[p][l];
                                    su[p][i][l] = su[p][i][l] - viscosity * metric[j][k] * symb[l][j][i] * basisFirstDerivative[p][k];
                                    su[p][i][l] = su[p][i][l] - viscosity * metric[j][k] * symb[l][k][i] * basisFirstDerivative[p][j];
                                    su[p][i][l] = su[p][i][l] - viscosity * metric[j][k] * dSymb[l][j][i][k] * basis[p];
                                    su[p][i][j] = su[p][i][j] + viscosity * metric[i][k] * symb[j][k][l] * basisFirstDerivative[p][l];
                                    su[p][i][l] = su[p][i][l] - viscosity * metric[i][k] * symb[l][j][j] * basisFirstDerivative[p][k];
                                    su[p][i][l] = su[p][i][l] - viscosity * metric[i][k] * symb[l][k][j] * basisFirstDerivative[p][j];
                                    su[p][i][l] = su[p][i][l] - viscosity * metric[i][k] * dSymb[l][j][j][k] * basis[p];
                                    for (m=0; m<dim; m++) {
                                        su[p][i][l] = su[p][i][l] - viscosity * metric[j][k] * symb[m][k][i] * symb[l][j][m] * basis[p];
                                        su[p][i][l] = su[p][i][l] + viscosity * metric[j][k] * symb[j][k][m] * symb[l][m][i] * basis[p];
                                        su[p][i][l] = su[p][i][l] - viscosity * metric[i][k] * symb[m][k][j] * symb[l][j][m] * basis[p];
                                        su[p][i][l] = su[p][i][l] + viscosity * metric[i][k] * symb[j][k][m] * symb[l][m][j] * basis[p];
                                    }
                                }
                            }
                        }
                    }
                    
                    sw[p][i][c] = sw[p][i][c] + density * basisFirstDerivative[p][i];
                    for (j=0; j<dim; j++) {
                        sw[p][i][i] = sw[p][i][i] + density * basisFirstDerivative[p][j] * velo[j];
                        for (k=0; k<dim; k++) {
                            sw[p][i][i] = sw[p][i][i] + viscosity * metric[j][k] * basisSecondDerivative[p][j][k];
                            sw[p][i][i] = sw[p][i][i] + dViscositydx[j] * metric[j][k] * basisFirstDerivative[p][j];
                            sw[p][i][j] = sw[p][i][j] + viscosity * metric[j][k] * basisSecondDerivative[p][i][k];
                            sw[p][i][j] = sw[p][i][j] + dViscositydx[j] * metric[j][k] * basisFirstDerivative[p][i];
                            if (coordinates != cartesian) {
                                sw[p][i][k] = sw[p][i][k] - density * symb[i][j][k] * basis[p] * velo[j];
                                for (l=0; l<dim; l++) {
                                    sw[p][i][i] = sw[p][i][i] - viscosity * metric[j][k] * symb[j][k][l] * basisFirstDerivative[p][l];
                                    sw[p][i][j] = sw[p][i][j] - viscosity * metric[j][k] * symb[i][k][l] * basisFirstDerivative[p][l];
                                    sw[p][i][l] = sw[p][i][l] - viscosity * metric[j][k] * symb[i][j][l] * basisFirstDerivative[p][k];
                                    sw[p][i][l] = sw[p][i][l] - viscosity * metric[j][k] * symb[i][j][l] * basisFirstDerivative[p][k];
                                    sw[p][i][l] = sw[p][i][l] - viscosity * metric[j][k] * symb[i][k][l] * basisFirstDerivative[p][j];
                                    sw[p][i][l] = sw[p][i][l] - viscosity * metric[j][k] * symb[j][k][l] * basisFirstDerivative[p][i];
                                    sw[p][i][l] = sw[p][i][l] - viscosity * metric[j][k] * dSymb[i][j][l][k] * basis[p];
                                    sw[p][i][l] = sw[p][i][l] - viscosity * metric[j][k] * dSymb[i][j][l][k] * basis[p];
                                    for (m=0; m<dim; m++) {
                                        sw[p][i][l] = sw[p][i][l] + viscosity * metric[j][k] * symb[i][k][m] * symb[m][j][l] * basis[p];
                                        sw[p][i][l] = sw[p][i][l] + viscosity * metric[j][k] * symb[j][k][m] * symb[m][i][l] * basis[p];
                                        sw[p][i][l] = sw[p][i][l] + viscosity * metric[j][k] * symb[i][k][m] * symb[m][j][l] * basis[p];
                                        sw[p][i][l] = sw[p][i][l] + viscosity * metric[j][k] * symb[j][k][m] * symb[m][i][l] * basis[p];
                                    }
                                }
                            }
                        }
                    }
                    if (coordinates == axis_symmetric) {
                        su[p][i][2] = 0.0;
                        sw[p][i][2] = 0.0;
                    }
                }
            }
        }
        
        // Loop over basis functions (of both unknowns and weights)
        for (p=0; p<n; p++) {
            for (q=0; q<n; q++) {
                
                // First plain Navier-Stokes
                
                // Mass metrix
                memset( *mass, 0.0, (4*4)*sizeof(double) );
                for (i=0; i<dim; i++) {
                    mass[i][i] = density * basis[q] * basis[p];
                }

                // Stiffness matrix
                memset( *a, 0.0, (4*4)*sizeof(double) );
                
                // Diffusive terms
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        for (k=0; k<dim; k++) {
                            a[i][i] = a[i][i] + viscosity * metric[j][k] * basisFirstDerivative[q][k] * basisFirstDerivative[p][j];
                            a[i][j] = a[i][j] + viscosity * metric[i][k] * basisFirstDerivative[q][k] * basisFirstDerivative[p][j];
                            if (coordinates != cartesian) {
                                for (l=0; l<dim; l++) {
                                    a[i][l] = a[i][l] + viscosity * metric[j][k] * symb[l][k][i] * basis[q] * basisFirstDerivative[p][j];
                                    a[i][l] = a[i][l] + viscosity * metric[i][k] * symb[l][k][j] * basis[q] * basisFirstDerivative[p][j];
                                    a[l][i] = a[l][i] - viscosity * metric[j][k] * basisFirstDerivative[q][k] * symb[i][j][l] * basis[p];
                                    a[l][j] = a[l][j] - viscosity * metric[i][k] * basisFirstDerivative[q][k] * symb[i][j][l] * basis[p];
                                    for (m=0; m<dim; m++) {
                                        a[l][m] = a[l][m] - viscosity * metric[j][k] * symb[m][k][i] * basis[q] * symb[i][j][l] * basis[p];
                                        a[l][m] = a[l][m] - viscosity * metric[i][k] * symb[m][k][j] * basis[q] * symb[i][j][l] * basis[p];
                                    }
                                }
                            }
                        }
                    }
                }
                
                // Convective terms, Picard linearization
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        a[i][i] = a[i][i] + density * basisFirstDerivative[q][j] * velo[j] * basis[p];
                        if (coordinates != cartesian) {
                            for (k=0; k<dim; k++) {
                                a[i][k] = a[i][k] + density * symb[k][j][i] * basis[q] * velo[j] * basis[p];
                            }
                        }
                    }
                }
                
                // Newton linearization
                if (newtonLinearization == YES) {
                    for (i=0; i<dim; i++) {
                        for (j=0; j<dim; j++) {
                            a[i][j] = a[i][j] + density * basisFirstDerivative[i][j] * basis[q] * basis[p];
                            if (coordinates != cartesian) {
                                for (k=0; k<dim; k++) {
                                    a[i][j] = a[i][j] + density * symb[k][j][i] * velo[k] * basis[q] * basis[p];
                                }
                            }
                        }
                    }
                }
                
                // Pressure terms
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        a[i][c] = a[i][c] - metric[i][j] * basis[q] * basisFirstDerivative[p][j];
                        if (coordinates != cartesian) {
                            for (k=0; k<dim; k++) {
                                a[k][c] = a[k][c] + metric[i][j] * basis[q] * symb[i][j][k] * basis[p];
                            }
                        }
                    }
                }
                
                // Continuity equation
                for (i=0; i<dim; i++) {
                    a[c][i] = a[c][i] + density * basisFirstDerivative[q][i] * basis[p];
                    if (coordinates != cartesian) {
                        for (j=0; j<dim; j++) {
                            a[c][j] = a[c][j] + density * symb[j][i][i] * basis[q] * basis[p];
                        }
                    }
                }
                
                // Stabilization
                if (stabilize == YES) {
                    for (i=0; i<dim; i++) {
                        for (j=0; j<c+1; j++) {
                            mass[j][i] = mass[j][i] + tau * density * basis[q] * sw[p][i][j];
                            for (k=0; k<c+1; k++) {
                                a[j][k] = a[j][k] + tau * su[q][i][k] * sw[p][i][j];
                            }
                        }
                        for (j=0; j<dim; j++) {
                            for (k=0; k<dim; k++) {
                                a[j][i] = a[j][i] + delta * basisFirstDerivative[q][i] * metric[j][k] * basisFirstDerivative[p][k];
                                if (coordinates != cartesian) {
                                    for (l=0; l<dim; l++) {
                                        a[l][i] = a[l][i] - delta * basisFirstDerivative[q][i] * metric[j][k] * symb[j][k][l] * basis[p];
                                        a[j][l] = a[j][l] + delta * symb[l][i][i] * basis[q] * metric[j][k] * basisFirstDerivative[p][k];
                                        for (m=0; m<dim; m++) {
                                            a[m][l] = a[m][l] - delta * symb[l][i][i] * basis[q] * metric[j][k] * symb[j][k][m] * basis[p];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
                // Add nodal matrix to element matrix
                if (coordinates == axis_symmetric) {
                    for (i=0; i<3; i++) {
                        for (j=0; j<3; j++) {
                            stiffMatrix[3*p+i][3*q+j] = stiffMatrix[3*p+i][3*q+j] + s * a[imap[i]][imap[j]];
                            massMatrix[3*p+i][3*q+j] = massMatrix[3*p+i][3*q+j] + s * mass[imap[i]][imap[j]];
                        }
                    }
                } else {
                    for (i=0; i<c+1; i++) {
                        for (j=0; j<c+1; j++) {
                            stiffMatrix[(c+1)*p+i][(c+1)*q+j] = stiffMatrix[(c+1)*p+i][(c+1)*q+j] + s * a[i][j];
                            massMatrix[(c+1)*p+i][(c+1)*q+j] = massMatrix[(c+1)*p+i][(c+1)*q+j] + s * mass[i][j];
                        }
                    }
                }
            }
        }
        
        // The righthand side
        if (newtonLinearization == YES) {
            for (i=0; i<dim; i++) {
                for (j=0; j<dim; j++) {
                    force[i] = force[i] + dVelodx[i][j] * uVelo[j];
                    if (coordinates != cartesian) {
                        for (k=0; k<dim; k++) {
                            force[i] = force[i] + symb[k][j][i] * uVelo[k] * uVelo[j];
                        }
                    }
                }
            }
        }
        
        for (p=0; p<n; p++) {
            memset(load, 0.0, sizeof(load) );
            for (i=0; i<c+1; i++) {
                load[i] = load[i] + density * force[i] * basis[p];
            }
            if (stabilize == YES) {
                for (i=0; i<c+1; i++) {
                    for (j=0; j<c+1; j++) {
                        load[j] = load[j] + tau * density * force[i] * sw[p][i][j];
                    }
                }
            }
            
            if (coordinates == axis_symmetric) {
                for (i=0; i<3; i++) {
                    forceVector[3*p+i] = forceVector[3*p+i] + s * load[imap[i]];
                }
            } else {
                for (i=0; i<c+1; i++) {
                    forceVector[(c+1)*p+i] = forceVector[(c+1)*p+i] + s * load[i];
                }
            }
        }
    }
}

/*****************************************************************************************************************************
    Return element local matrices and RHS vector for Navier-Stokes equations boundary conditions (no velocity dependent velocity
    BCs (Newton BCs) at the moement, so boundaryMatrix will contain only zeros at exit).
 
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
-(void)navierStokesGeneralBoundary:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector loadVector:(double **)loadVector nodalAlpha:(double *)nodalAlpha nodalBeta:(double *)nodalBeta nodalExtPressure:(double *)nodalExtPressure nodalSlipCoefficient:(double **)nodalSlipCoefficient element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh  model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems {
    
    int c, i, j, p, q, t, coordinates, dim;
    double alpha, detJ, s, slipCoeff, sqrtMetric, sum, u, v, x, y, w, z;
    double dSymb[3][3][3][3], force[3], metric[3][3], normals[3], symb[3][3][3];
    double *basis = NULL, **basisFirstDerivative = NULL;
    BOOL stat;
    GaussIntegrationPoints *IP = NULL;
    
    basis = integration.basis;
    basisFirstDerivative = integration.basisFirstDerivative;
    
    coordinates = coordinateSystems.coordinates;
    if (coordinates == cylindric_symmetric) {
        dim = 3;
    } else {
        dim = model.dimension;
    }
    c = dim + 1;
    
    IP = GaussQuadrature(element, NULL, NULL);
    
    for (t=0; t<IP->n; t++) {
        
        u = IP->u[t];
        v = IP->v[t];
        w = IP->w[t];
        
        // Basis function values and derivatives at the integration point
        stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
        stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
        detJ = integration.metricDeterminant;
        
        // Coordinate system dependent info
        x = 0.0;
        y = 0.0;
        z = 0.0;
        if (coordinates != cartesian) {
            x = cblas_ddot(n, nodes->x, 1, basis, 1);
            y = cblas_ddot(n, nodes->y, 1, basis, 1);
            z = cblas_ddot(n, nodes->z, 1, basis, 1);
        }
        
        [coordinateSystems coordinateSystemInfoModel:model metric:metric sqrtMetric:&sqrtMetric symbols:symb dSymbols:dSymb coordX:x coordY:y coordZ:z];
        s = sqrtMetric * detJ * IP->s[t];
        
        // Add to load: tangential derivative of something
        memset(force, 0.0, sizeof(force) );
        for (i=0; i<3; i++) {
            for (j=0; j<n; j++) {
                force[i] = force[i] + nodalBeta[j] * basisFirstDerivative[j][i];
            }
        }
        
        // Add to load: given force in normal direction
        BOOL check = YES;
        [elementDescription normalVectorForBDElement:element boundaryNodes:nodes mesh:mesh paraU:&u paraV:&v check:&check normals:normals];
        alpha = cblas_ddot(n, nodalExtPressure, 1, basis, 1);
        for (i=0; i<dim; i++) {
            for (j=0; j<dim; j++) {
                force[i] = force[i] + alpha * metric[i][j] * normals[j];
            }
        }
        
        alpha = cblas_ddot(n, nodalAlpha, 1, basis, 1);
        
        // Add to load: given force in coordinate directions
        for (i=0; i<3; i++) {
            sum = 0.0;
            for (j=0; j<n; j++) {
                sum = sum + loadVector[i][j];
            }
            force[i] = force[i] + sum;
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
                        boundaryMatrix[p*c+i][q*c+i] = boundaryMatrix[p*c+i][q*c+i] + s * slipCoeff * basis[q] * basis[p];
                    }
                }
            }
        }
        
        for (q=0; q<n; q++) {
            for (i=0; i<dim; i++) {
                boundaryVector[q*c+i] = boundaryVector[q*c+i] + s * basis[q] * force[i];
            }
        }
    }
}

@end
