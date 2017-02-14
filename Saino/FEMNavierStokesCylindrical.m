//===----------------------------------------------------------------------===//
//  NavierStokesCylindrical.m
//  Saino
//
//  Created by Seddik hakime on 26/05/2014.
//  Copyright (c) 2014 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

#import "FEMNavierStokesCylindrical.h"

@implementation FEMNavierStokesCylindrical {
    
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
        double **loadVector                         > loadVector vector
        double *nodalViscosity                     -> nodal values for viscosity (i.e. if turbulence model or
                                                      power-law viscosity is used, the values vary in space)
        double *nodalDensity                       -> nodal values of density
        double *velocityX, *velocityY, *velocityZ  -> nodal values of velocity components from previous iteration
        double *nodalPressure                      -> nodal values of total pressure from previous iteration
        NSString *stabilizeFlag                    -> should stabilization be used?
        BOOL compressible                          -> should compressible terms be added?
        BOOL pseudoCompressible                    -> should artificial compressibility be added ?
        double *nodalCompressibility               -> artificial compressibility for the nodes
        BOOL magneticForce                         -> should Lorentz force for magneto-hydrodynamics be included
        BOOL newtonLinearization                   -> Picard or Newton linearization of the convetion term ?
        Element_t *element                         -> structure describing the element (dimension,nof nodes,
                                                      interpolation degree, etc...)
        int n                                      -> number of element nodes
        Nodes_t *nodes                             -> element node coordinates
 *********************************************************************************************************************/
-(void)navierStokesCylindricalComposeMassMatrix:(double **) massMatrix stiffMatrix:(double **)stiffMatrix forceVector:(double *)forceVector loadVector:(double **)loadVector nodalViscosity:(double *)nodalViscosity nodalDensity:(double *)nodalDensity velocityX:(double *)ux velocityY:(double *)uy velocityZ:(double *)uz meshVelocityX:(double *)mux meshVelocityY:(double *)muy meshVelocityZ:(double *)muz nodalPressure:(double *)nodalPressure nodalTemperature:(double *)nodalTemperature isConvect:(BOOL)convect stabilizeFlag:(NSString *)stabilizeFlag isCompressible:(BOOL)compressible isPseudoCompressible:(BOOL)pseudoCompressible nodalCompressibility:(double *)nodalCompressibility nodalGasConstant:(double *)nodalGasConstant isPorous:(BOOL)porous nodalDrag:(double **)nodalDrag isPotentialForce:(BOOL)potentialForce potentialField:(double *)potentialField potentialCoefficient:(double *)potentialCoefficient isMagneticForce:(BOOL)magneticForce isDivDiscretization:(BOOL)divDiscretization isGradPDriscretization:(BOOL)gradPDriscretization isNewtonLinearization:(BOOL)newtonLinearization element:(Element_t *)element numberOfNodes:(int)n rows:(int)rows cols:(int)cols nodes:(Nodes_t *)nodes core:(FEMCore *)core mesh:(FEMMesh *)mesh model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription coordinateSystems:(FEMCoordinateSystems *)coordinateSystems materialModels:(FEMMaterialModels *)materialModels differentials:(FEMDifferentials *)differentials listUtilities:(FEMListUtilities *)listUtilities utilities:(FEMUtilities *)utilities {
    
    int c, i, j, k, l, m, p, q, t, coordinates, dim, linearBasis, nBasis;
    int imap[3] = {0, 1, 3};
    double baseP, compress=0.0, delta=0.0, density, detJ, gasConstant, hk=0.0, lambda=1.0, mk=0.0, pressure=0.0, re, s, sqrtMetric, sum, tau=0.0, temperature=0.0, x, u, v,
           vNorm, viscosity, y, w, z;
    double a[4][4], dDensitydx[3], dNodalBasisdx[n][n][3], dPressuredx[3], drag[3], dSymb[3][3][3][3], dTemperaturedx[3], dVelodx[3][3], dViscositydx[3],
           force[4], gMeric[3][3], load[4], lrf[3], mass[4][4], metric[3], pBasis[n], pdBasisdx[n][3], su[n][4][4], sw[n][4][4], symb[3][3][3], uVelo[3],
           velo[3];
    double *basis = NULL, **basisFirstDerivative = NULL;
    BOOL bubbles, cylindricSymmetry, pBubbles, p2p1, stabilize, stat;
    ElementType_t *linearType = NULL, *saveType = NULL;
    GaussIntegrationPoints *IP = NULL;
    
    basis = integration.basis;
    basisFirstDerivative = integration.basisFirstDerivative;
    
    coordinates = coordinateSystems.coordinates;
    cylindricSymmetry = (coordinates == cylindric_symmetric || coordinates == axis_symmetric) ? YES : NO;
    
    if (cylindricSymmetry == YES) {
        dim = 3;
    } else {
        dim = model.dimension;
    }
    c = dim;
    
    nBasis = n;
    bubbles = NO;
    pBubbles = NO;
    p2p1 = NO;
    stabilize = ([stabilizeFlag isEqualToString:@"stabilized"] == YES) ? YES : NO;
    if (stabilize == NO || compressible == YES) {
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
    
    if (stabilize == YES) {
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
                    dNodalBasisdx[i][p][j] = basisFirstDerivative[i][j];
                }
            }
        }
    }
    
    // Start integrating
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
        
        // Coordinate system dependent info
        x = 0.0;
        y = 0.0;
        z = 0.0;
        if (coordinates != cartesian) {
            x = cblas_ddot(n, nodes->x, 1, basis, 1);;
            y = cblas_ddot(n, nodes->y, 1, basis, 1);
            z = cblas_ddot(n, nodes->z, 1, basis, 1);
        }
        [coordinateSystems coordinateSystemInfoModel:model metric:gMeric sqrtMetric:&sqrtMetric symbols:symb dSymbols:dSymb coordX:x coordY:y coordZ:z];
        for (i=0; i<3; i++) {
            metric[i] = gMeric[i][i];
        }
        s = sqrtMetric * detJ * IP->s[t];
        
        // Density at the integration point
        density = cblas_ddot(n, nodalDensity, 1, basis, 1);

        if (compressible == YES) {
            memset(dDensitydx, 0.0, sizeof(dDensitydx) );
            memset(dPressuredx, 0.0, sizeof(dPressuredx) );
            memset(dTemperaturedx, 0.0, sizeof(dTemperaturedx) );
            if (p2p1 == YES) {
                k = linearBasis;
                temperature = cblas_ddot(k, nodalTemperature, 1, pBasis, 1);
                pressure = cblas_ddot(k, nodalPressure, 1, pBasis, 1);
                for (i=0; i<dim; i++) {
                    for (j=0; j<k; j++) {
                        dDensitydx[i] = dDensitydx[i] + nodalDensity[j] * pdBasisdx[j][i];
                        dPressuredx[i] = dPressuredx[i] + nodalPressure[j] * pdBasisdx[j][i];
                        dTemperaturedx[i] = dTemperaturedx[i] + nodalTemperature[j] * pdBasisdx[j][i];
                    }
                }
            } else {
                temperature = cblas_ddot(n, nodalTemperature, 1, basis, 1);
                pressure = cblas_ddot(n, nodalPressure, 1, basis, 1);
                for (i=0; i<dim; i++) {
                    for (j=0; j<n; j++) {
                        dDensitydx[i] = dDensitydx[i] + nodalDensity[j] * basisFirstDerivative[j][i];
                        dPressuredx[i] = dPressuredx[i] + nodalPressure[j] * basisFirstDerivative[j][i];
                        dTemperaturedx[i] = dTemperaturedx[i] + nodalTemperature[j] * basisFirstDerivative[j][i];
                    }
                }
            }
            gasConstant = cblas_ddot(n, nodalGasConstant, 1, basis, 1);
            density = pressure / (temperature * gasConstant);
        }
        
        if (pseudoCompressible == YES) {
            pressure = cblas_ddot(n, nodalPressure, 1, basis, 1);
            compress = density * cblas_ddot(n, nodalCompressibility, 1, basis, 1);
        }
        
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
            memset( *dVelodx, 0.0, (3*3)*sizeof(double) );
            for (i=0; i<3; i++) {
                for (j=0; j<n; j++) {
                    dVelodx[0][i] = dVelodx[0][i] + ux[j] * basisFirstDerivative[j][i];
                    dVelodx[1][i] = dVelodx[1][i] + uy[j] * basisFirstDerivative[j][i];
                    if (dim > 2 && coordinates != axis_symmetric) dVelodx[2][i] = dVelodx[2][i] + uz[j] * basisFirstDerivative[j][i];
                }
            }
        }
        
        // Force at the integration point
        if (magneticForce == YES) {
            memset(lrf, 0.0, sizeof(lrf) );
            [differentials lorentzForceElement:element nodes:nodes numberOfNodes:n integrationU:u integrationV:v integrationW:w lorentzForce:lrf mesh:mesh model:model integration:integration coordinateSystems:coordinateSystems listUtilities:listUtilities utilities:utilities];
        }
        
        memset(force, 0.0, sizeof(force) );
        for (i=0; i<n; i++) {
            force[0] = force[0] + loadVector[0][i] * basis[i];
            force[1] = force[1] + loadVector[1][i] * basis[i];
            force[2] = force[2] + loadVector[2][i] * basis[i];
        }
        if (dim > 2 && coordinates != axis_symmetric) {
            for (i=0; i<3; i++) {
                force[i] = force[i] + lrf[i] / density;
            }
            for (i=0; i<n; i++) {
                force[3] = force[3] + loadVector[3][i] * basis[i];
            }
        } else {
            for (i=0; i<2; i++) {
                force[i] = force[i] + lrf[i] / density;
            }
        }
        
        // Additional forces due to gradient forces (electrokinetic flow) and viscous drag in porous media
        if (potentialForce == YES) {
            sum = 0.0;
            for (i=0; i<n; i++) {
                sum = sum + potentialField[i] * basisFirstDerivative[i][0];
            }
            force[0] = force[0] - cblas_ddot(n, potentialCoefficient, 1, basis, 1) * sum;
            sum = 0.0;
            for (i=0; i<n; i++) {
                sum = sum + potentialField[i] * basisFirstDerivative[i][1];
            }
            force[1] = force[1] - cblas_ddot(n, potentialCoefficient, 1, basis, 1) * sum / x;
            if (dim > 2 && coordinates != axis_symmetric) {
                sum = 0.0;
                for (i=0; i<n; i++) {
                    sum = sum + potentialField[i] * basisFirstDerivative[i][2];
                }
                force[2] = force[2] - cblas_ddot(n, potentialCoefficient, 1, basis, 1) * sum;
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
            
            if (convect == YES) {
                vNorm = 0.0;
                for (i=0; i<dim; i++) {
                    vNorm = vNorm + velo[i] * velo[i] / metric[i];
                }
                vNorm = max(sqrt(vNorm), 1.0e-12);
                re = min(1.0, density * mk * hk * vNorm / (4.0 * viscosity));
                tau = 0.0;
                if (vNorm != 0.0) {
                    tau = hk * re / (2.0 * density * vNorm);
                }
                delta = density * lambda * re * hk * vNorm;
            } else {
                delta = 0.0;
                tau = mk * pow(hk, 2.0) / (8.0 * viscosity);
            }
            
            // su will contain residual of NS-equations (except for the time derivative and force terms).
            // sw will contain the weight function values
            memset( **su, 0.0, (n*4*4)*sizeof(double) );
            memset( **sw, 0.0, (n*4*4)*sizeof(double) );
            for (p=0; p<n; p++) {
                for (i=0; i<dim; i++) {
                    su[p][i][c] = su[p][i][c] + metric[i] * basisFirstDerivative[p][i];
                    if (porous == YES) {
                        su[p][i][i] = su[p][i][i] + viscosity * drag[i] * basis[p];
                    }
                    
                    if (convect == YES) {
                        for (j=0; j<dim; j++) {
                            su[p][i][i] = su[p][i][i] + density * basisFirstDerivative[p][j] * velo[j];
                            if (newtonLinearization == YES) {
                                su[p][i][j] = su[p][i][j] + density * dVelodx[i][j] * basis[p];
                            }
                            for (k=0; k<dim; k+=2) {
                                su[p][i][k] = su[p][i][k] + density * symb[k][j][i] * basis[p] * velo[j];
                                if (newtonLinearization) {
                                    su[p][i][j] = su[p][i][j] + density * symb[k][j][i] * velo[k] * basis[p];
                                }
                            }
                        }
                    }
                    
                    // Diffusion
                    for (j=0; j<dim; j++) {
                        su[p][i][i] = su[p][i][i] - dViscositydx[j] * metric[j] * basisFirstDerivative[p][j];
                        su[p][i][j] = su[p][i][j] - dViscositydx[j] * metric[i] * basisFirstDerivative[p][i];
                        sum = 0.0;
                        for (k=0; k<n; k++) {
                            sum = sum + dNodalBasisdx[p][k][j] * basisFirstDerivative[k][j];
                        }
                        su[p][i][i] = su[p][i][i] - viscosity * metric[j] * sum;
                        sum = 0.0;
                        for (k=0; k<n; k++) {
                            sum = sum + dNodalBasisdx[p][k][j] * basisFirstDerivative[k][i];
                        }
                        su[p][i][j] = su[p][i][j] - viscosity * metric[i] * sum;
                    }
                    
                    for (j=0; j<dim; j+=2) {
                        for (l=0; l<dim; l+=2) {
                            su[p][i][i] = su[p][i][i] + viscosity * metric[j] * symb[j][j][l] * basisFirstDerivative[p][l];
                            su[p][i][l] = su[p][i][l] - viscosity * metric[j] * symb[l][j][i] * basisFirstDerivative[p][j];
                            su[p][i][l] = su[p][i][l] - viscosity * metric[j] * symb[l][j][i] * basisFirstDerivative[p][j];
                            su[p][i][l] = su[p][i][l] - viscosity * metric[j] * dSymb[l][j][i][j] * basis[p];
                            su[p][i][j] = su[p][i][j] + viscosity * metric[i] * symb[j][i][l] * basisFirstDerivative[p][l];
                            su[p][i][l] = su[p][i][l] - viscosity * metric[i] * symb[l][j][j] * basisFirstDerivative[p][i];
                            su[p][i][l] = su[p][i][l] - viscosity * metric[i] * symb[l][i][j] * basisFirstDerivative[p][j];
                            su[p][i][l] = su[p][i][l] - viscosity * metric[i] * dSymb[l][j][j][i] * basis[p];
                            for (m=0; m<dim; m+=2) {
                                su[p][i][l] = su[p][i][l] - viscosity * metric[j] * symb[j][m][i] * symb[l][j][m] * basis[p];
                                su[p][i][l] = su[p][i][l] + viscosity * metric[j] * symb[j][j][m] * symb[l][m][i] * basis[p];
                                su[p][i][l] = su[p][i][l] - viscosity * metric[i] * symb[m][i][j] * symb[l][j][m] * basis[p];
                                su[p][i][l] = su[p][i][l] + viscosity * metric[i] * symb[j][i][m] * symb[l][m][j] * basis[p];
                            }
                        }
                        // Then -mu,_j (g^{jk} U^i_,k + g^{ik} U^j_,k)
                        for (l=0; l<dim; l+=2) {
                            su[p][i][l] = su[p][i][l] - dViscositydx[j] * ( metric[j] * basis[p] * symb[j][l][i] + metric[i] * basis[p] * symb[i][l][j] );
                        }
                    }
                    
                    if (convect == YES) {
                        sw[p][i][c] = sw[p][i][c] + density * basisFirstDerivative[p][i];
                        for (j=0; j<dim; j++) {
                            sw[p][i][i] = sw[p][i][i] + density * basisFirstDerivative[p][j] * velo[j];
                            for (k=0; k<dim; k+=2) {
                                sw[p][i][k] = sw[p][i][k] - density * symb[i][j][k] * basis[p] * velo[j];
                            }
                        }
                    } else {
                        sw[p][i][c] = sw[p][i][c] + basisFirstDerivative[p][i];
                    }
                    
                    // Diffusion
                    for (j=0; j<dim; j++) {
                        sw[p][i][i] = sw[p][i][i] + dViscositydx[j] * metric[j] * basisFirstDerivative[p][j];
                        sw[p][i][j] = sw[p][i][j] + dViscositydx[j] * metric[j] * basisFirstDerivative[p][i];
                        sum = 0.0;
                        for (k=0; k<n; k++) {
                            sum = sum + dNodalBasisdx[p][k][j] * basisFirstDerivative[k][j];
                        }
                        sw[p][i][i] = sw[p][i][i] - viscosity * metric[j] * sum;
                        sum = 0.0;
                        for (k=0; k<n; k++) {
                            sum = sum + dNodalBasisdx[p][k][j] * basisFirstDerivative[k][i];
                        }
                        sw[p][i][j] = sw[p][i][j] - viscosity * metric[j] * sum;
                    }
                    
                    for (j=0; j<dim; j+=2) {
                        for (l=0; l<dim; l+=2) {
                            sw[p][i][i] = sw[p][i][i] - viscosity * metric[j] * symb[j][j][l] * basisFirstDerivative[p][l];
                            sw[p][i][l] = sw[p][i][l] - viscosity * metric[j] * symb[i][j][l] * basisFirstDerivative[p][j];
                            sw[p][i][l] = sw[p][i][l] - viscosity * metric[j] * symb[i][j][l] * basisFirstDerivative[p][j];
                            sw[p][i][l] = sw[p][i][l] - viscosity * metric[j] * dSymb[i][j][l][j] * basis[p];
                            sw[p][i][j] = sw[p][i][j] - viscosity * metric[j] * symb[i][j][l] * basisFirstDerivative[p][l];
                            sw[p][i][l] = sw[p][i][l] - viscosity * metric[j] * symb[i][j][l] * basisFirstDerivative[p][j];
                            sw[p][i][l] = sw[p][i][l] - viscosity * metric[j] * symb[j][j][l] * basisFirstDerivative[p][i];
                            sw[p][i][l] = sw[p][i][l] - viscosity * metric[j] * dSymb[i][j][l][j] * basis[p];
                            for (m=0; m<dim; m+=2) {
                                sw[p][i][l] = sw[p][i][l] + viscosity * metric[j] * symb[i][j][m] * symb[m][j][l] * basis[p];
                                sw[p][i][l] = sw[p][i][l] + viscosity * metric[j] * symb[j][j][m] * symb[m][i][l] * basis[p];
                                sw[p][i][l] = sw[p][i][l] + viscosity * metric[j] * symb[j][j][m] * symb[m][i][l] * basis[p];
                                sw[p][i][l] = sw[p][i][l] + viscosity * metric[j] * symb[i][j][m] * symb[m][j][l] * basis[p];
                            }
                        }
                        // Then -mu,_j g^{jk} (w_i,_k + w_k,_i)
                        for (l=0; l<dim; l+=2) {
                            sw[p][i][l] = sw[p][i][l] + dViscositydx[j] * ( metric[j] * basis[p] * symb[i][j][l] + metric[i] * basis[p] * symb[j][i][l] );
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
        for (p=0; p<nBasis; p++) {
            for (q=0; q<nBasis; q++) {
                
                // First plain Navier-Stokes
                
                if (p2p1 == YES && p < linearBasis) {
                    baseP = pBasis[p];
                } else {
                    baseP = basis[p];
                }
                
                // Mass metrix
                memset( *mass, 0.0, (4*4)*sizeof(double) );
                for (i=0; i<dim; i++) {
                    mass[i][i] = density * basis[q] * basis[p];
                }
                
                // Continuity equation
                if (compressible == YES) {
                    mass[c][c] = (density / pressure) * basis[q] * baseP;
                }
                
                // Stiffness matrix
                memset( *a, 0.0, (4*4)*sizeof(double) );
                
                // Possible Porous media effects
                if (porous == YES) {
                    for (i=0; i<dim; i++) {
                        a[i][i] = a[i][i] + viscosity * drag[i] * basis[q] * basis[p];
                    }
                }
                
                // Diffusive terms
                for (i=0; i<dim; i++) {
                    for (j=0; j<dim; j++) {
                        a[i][i] = a[i][i] + viscosity * metric[j] * basisFirstDerivative[q][j] * basisFirstDerivative[p][j];
                        a[i][j] = a[i][j] + viscosity * metric[i] * basisFirstDerivative[q][i] * basisFirstDerivative[p][j];
                        if (compressible == YES) {
                            // For compressible flows, add (2/3) \mu \nabla \cdot u
                            // Partial derivative terms only here
                            a[i][j] = a[i][j] - (2.0 / 3.0) * viscosity * metric[i] * basisFirstDerivative[q][j] * basisFirstDerivative[p][i];
                        }
                    }
                }
                
                if (compressible == YES) {
                    // For compressible flows, add (2/3) \mu \nabla \cdot u
                    // Terms involving Christoffel symbols
                    for (i=0; i<dim; i++) {
                        a[i][0] = a[i][0] - (2.0 / 3.0) * viscosity * metric[i] * basis[q] * symb[2][0][2] * basisFirstDerivative[p][i];
                    }
                    for (j=0; j<dim; j++) {
                        a[0][j] = a[0][j] + (2.0 / 3.0) * viscosity * basisFirstDerivative[q][j] * metric[2] * symb[2][2][0] * basis[p];
                    }
                    a[0][0] = a[0][0] + (2.0 / 3.0) * viscosity * symb[2][0][2] * basis[q] * metric[2] * symb[2][2][0] * basis[p];
                }
                
                if (cylindricSymmetry == NO) {
                    a[0][2] = a[0][2] - 2.0 * viscosity * metric[2] * symb[2][2][0] * basis[p] * basisFirstDerivative[q][2];
                    a[2][0] = a[2][0] + 2.0 * viscosity * metric[2] * symb[0][2][2] * basis[q] * basisFirstDerivative[p][2];
                    a[2][0] = a[2][0] - 2.0 * viscosity * metric[2] * symb[0][2][2] * basis[p] * basisFirstDerivative[q][2];
                }
                a[0][0] = a[0][0] + 2.0 * viscosity * metric[2] * basis[q] * basis[p];
                a[2][2] = a[2][2] - 2.0 * viscosity * symb[2][0][2] * basisFirstDerivative[q][0] * basis[p];
                
                // Convection terms, Picard linearization
                if (convect == YES) {
                    for (i=0; i<dim; i++) {
                        for (j=0; j<dim; j++) {
                            a[i][i] = a[i][i] + density * basisFirstDerivative[q][j] * velo[j] * basis[p];
                        }
                    }
                    a[0][2] = a[0][2] + density * symb[2][2][0] * basis[q] * velo[2] * basis[p];
                    a[2][0] = a[2][0] + density * symb[0][2][2] * basis[q] * velo[2] * basis[p];
                    a[2][2] = a[2][2] + density * symb[2][0][2] * basis[q] * velo[0] * basis[p];
                    
                    // Convection terns, Newton linearization
                    if (newtonLinearization == YES) {
                        for (i=0; i<dim; i++) {
                            for (j=0; j<dim; j++) {
                                a[i][j] = a[i][j] + density * dVelodx[i][j] * basis[q] * basis[p];
                            }
                        }
                        a[0][2] = a[0][2] + density * symb[2][2][0] * basis[q] * velo[2] * basis[p];
                        a[2][0] = a[2][0] + density * symb[2][0][2] * basis[q] * velo[2] * basis[p];
                        a[2][2] = a[2][2] + density * symb[0][2][2] * basis[q] * velo[0] * basis[p];
                    }
                }
                
                // Pressure terms
                if (gradPDriscretization == YES) {
                    for (i=0; i<dim; i++) {
                        a[i][c] = a[i][c] + metric[i] * basisFirstDerivative[q][i] * basis[p];
                    }
                } else {
                    for (i=0; i<dim; i++) {
                        a[i][c] = a[i][c] - metric[i] * basis[q] * basisFirstDerivative[p][i];
                    }
                    a[0][c] = a[0][c] + metric[2] * basis[q] * symb[2][2][0] * basis[p];
                }
                
                // Continuity equation
                for (i=0; i<dim; i++) {
                    if (compressible == YES) {
                        a[c][i] = a[c][i] + (density / pressure) * basis[q] * dPressuredx[i] * baseP / 2.0;
                        a[c][c] = a[c][c] + (density / pressure) * velo[i] * basisFirstDerivative[q][i] * baseP / 2.0;
                        a[c][c] = a[c][c] - (density / (temperature * pressure)) * velo[i] * dTemperaturedx[i] * basis[q] * baseP;
                        a[c][i] = a[c][i] + density * basisFirstDerivative[q][i] * baseP;
                    } else {
                        if (convect == YES) {
                            if (gradPDriscretization == YES) {
                                a[c][i] = a[c][i] - density * basis[q] * basisFirstDerivative[p][i];
                            } else {
                                a[c][i] = a[c][i] + density * basisFirstDerivative[q][i] * baseP;
                            }
                        } else {
                            if (gradPDriscretization == YES) {
                                a[c][i] = a[c][i] - basis[q] * basisFirstDerivative[p][i];
                            } else {
                                a[c][i] = a[c][i] + basisFirstDerivative[q][i] * baseP;
                            }
                        }
                    }
                }
                
                if (gradPDriscretization == NO) {
                    if (compressible == YES || convect == YES) {
                        a[c][0] = a[c][0] + density * symb[0][2][2] * basis[q] * baseP;
                    } else {
                        a[c][0] = a[c][0] + symb[0][2][2] * basis[q] * baseP;
                    }
                }
                
                // Artificial compressibility, affects only the continuity equation
                if (pseudoCompressible == YES) {
                    a[c][c] = a[c][c] + compress * basis[q] * basis[p];
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
                            a[j][i] = a[j][i] + delta * basisFirstDerivative[q][i] * metric[j] * basisFirstDerivative[p][j];
                            for (l=0; l<dim; l+=2) {
                                a[l][i] = a[l][i] - delta * basisFirstDerivative[q][i] * metric[j] * symb[j][j][l] * basis[p];
                                a[j][l] = a[j][l] + delta * symb[l][i][i] * basis[q] * metric[j] * basisFirstDerivative[p][j];
                                for (m=0; m<dim; m+=2) {
                                    a[m][l] = a[m][l] - delta * symb[l][i][i] * basis[q] * metric[j] * symb[j][j][m] * basis[p];
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
        if (convect == YES && newtonLinearization == YES) {
            memset(uVelo, 0.0, sizeof(uVelo) );
            uVelo[0] = cblas_ddot(n, ux, 1, basis, 1);
            uVelo[1] = cblas_ddot(n, uy, 1, basis, 1);
            if (dim > 2) uVelo[2] = cblas_ddot(n, uz, 1, basis, 1);
            
            for (i=0; i<dim; i++) {
                for (j=0; j<dim; j++) {
                    force[i] = force[i] + dVelodx[i][j] * uVelo[j];
                }
            }
            
            if (coordinates != axis_symmetric) {
                force[0] = force[0] + symb[2][2][0] * uVelo[2] * uVelo[2];
                force[2] = force[2] + symb[2][0][2] * uVelo[0] * uVelo[2];
                force[2] = force[2] + symb[0][2][2] * uVelo[2] * uVelo[0];
            }
        }
        
        for (p=0; p<nBasis; p++) {
            memset(load, 0.0, sizeof(load) );
            for (i=0; i<c+1; i++) {
                load[i] = load[i] + density * force[i] * basis[p];
            }
            if (compressible == YES) load[c] = load[c] /  temperature;
            
            if (pseudoCompressible == YES) {
                load[c] = load[c] + pressure * basis[p] * compress;
            }
            if (stabilize == YES) {
                for (i=0; i<dim; i++) {
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
    
    k = c+1;
    if (coordinates == axis_symmetric) k = 3;
    if (p2p1 == YES) {
        int size;
        int **edgeMap = NULL;
        j = [core getElementFamily:element];
        edgeMap = [core getEdgeMap:j mapSize:&size];
        for (i=j+1; i<=j+size; i++) {
            p = edgeMap[(i-j)-1][0];
            q = edgeMap[(i-j)-1][1];
            for (l=0; l<cols; l++) {
                stiffMatrix[(k*i)-1][l] = 0.0;
                massMatrix[(k*i)-1][l] = 0.0;
            }
            forceVector[(k*i)-1] = 0.0;
            stiffMatrix[(k*i)-1][(k*i)-1] = 1.0;
            stiffMatrix[(k*i)-1][(k*(p+1))-1] = -0.5;
            stiffMatrix[(k*i)-1][(k*(q+1))-1] = -0.5;
        }
    }
    
    if (pBubbles == YES) {
        for (i=n+1; i<=nBasis; i++) {
            for (l=0; l<cols; l++) {
                stiffMatrix[(k*i)-1][l] = 0.0;
                massMatrix[(k*i)-1][l] = 0.0;
            }
            for (l=0; l<rows; l++) {
                stiffMatrix[l][(k*i)-1] = 0.0;
                massMatrix[l][(k*i)-1] = 0.0;
            }
            forceVector[(k*i)-1] = 0.0;
            stiffMatrix[(k*i)-1][(k*i)-1] = 1.0;
        }
    }
}

/*****************************************************************************************************************************
    Return element local matrices and RHS vector for Navier-Stokes equations boundary conditions in cylindrical coordinates
    (no velocity dependent velocity BCs (Newton BCs) at the moement, so boundaryMatrix will contain only zeros at exit).
 
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
-(void)navierStokesCylindricalBoundary:(double **)boundaryMatrix boundaryVector:(double *)boundaryVector loadVector:(double **)loadVector nodalAlpha:(double *)nodalAlpha nodalBeta:(double *)nodalBeta nodalExtPressure:(double *)nodalExtPressure nodalSlipCoefficient:(double **)nodalSlipCoefficient isNormalTangential:(BOOL)normalTangential element:(Element_t *)element numberOfNodes:(int)n nodes:(Nodes_t *)nodes mesh:(FEMMesh *)mesh  model:(FEMModel *)model integration:(FEMNumericIntegration *)integration elementDescription:(FEMElementDescription *)elementDescription elementUtils:(FEMElementUtils *)elementUtils coordinateSystems:(FEMCoordinateSystems *)coordinateSystems {
    
    int c, i, j, k, l, p, q, t, coordinates, dim;
    double alpha, detJ, s, slipCoeff, sqrtMetric, sum, u, v, x, y, w, z;
    double dSymb[3][3][3][3], force[4], metric[3][3], normals[3], tangents[3], tangents2[3], tangentForce[3], symb[3][3][3], vect[3];
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
        memset(tangentForce, 0.0, sizeof(tangentForce) );
        for (i=0; i<dim; i++) {
            for (j=0; j<n; j++) {
                tangentForce[i] = tangentForce[i] + nodalBeta[j] * basisFirstDerivative[j][i];
            }
        }
        
        // Add to load: given force in coordinate directions
        memset(force, 0.0, sizeof(force) );
        for (i=0; i<c; i++) {
            sum = 0.0;
            for (j=0; j<n; j++) {
                sum = sum + loadVector[i][j] * basis[j];
            }
            force[i] = force[i] + sum;
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
            k = q * c;
            boundaryVector[k] = boundaryVector[k] - s * alpha * basis[q] * symb[2][0][2];
        }
    }
}

@end
