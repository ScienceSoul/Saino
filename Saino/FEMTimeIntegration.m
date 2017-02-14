//===----------------------------------------------------------------------===//
//  FEMTimeIntegration.m
//  Saino
//
//  Created by Hakime Seddik on 14/03/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
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

#import "FEMTimeIntegration.h"
#import "Utils.h"

@implementation FEMTimeIntegration

/*************************************************
    Fractional step integration scheme.
*************************************************/
-(void)fractionalStepInSolution:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * __nonnull * __nonnull )massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution:(double * __nonnull)prevSolution rows:(int * __nonnull)rows {
    
    int i, j, nb;
    double s, fsStep, fsTheta, fsdTheta, fsAlpha, fsBeta, massCoeff=0.0, forceCoeff=0.0;

    if (rows != NULL) {
        nb = *rows;
    } else {
        fprintf(stdout, "fractionalStepInSolution: dimension(s) missing.");
        return;
    }
    
    fsStep = 0.0;
    if (solution.solutionInfo[@"fsstep"] != nil) {
        fsStep = [solution.solutionInfo[@"fsstep"] doubleValue];
    }
    fsTheta = 0.0;
    if (solution.solutionInfo[@"fstheta"] != nil) {
        fsTheta = [solution.solutionInfo[@"fstheta"] doubleValue];
    }
    fsdTheta = 0.0;
    if (solution.solutionInfo[@"fsdtheta"] != nil) {
        fsdTheta = [solution.solutionInfo[@"fsdtheta"] doubleValue];
    }
    fsAlpha = 0.0;
    if (solution.solutionInfo[@"fsalpha"] != nil) {
        fsAlpha = [solution.solutionInfo[@"fsalpha"] doubleValue];
    }
    fsBeta = 0.0;
    if (solution.solutionInfo[@"fsbeta"] != nil) {
        fsBeta = [solution.solutionInfo[@"fsbeta"] doubleValue];
    }
    
    switch ((int)fsStep) {
        case 1:
            massCoeff = fsAlpha * fsTheta;
            forceCoeff = fsBeta * fsTheta;
            break;
        case 2:
            massCoeff = fsBeta * fsdTheta;
            forceCoeff = fsAlpha * fsdTheta;
            break;
        case 3:
            massCoeff = fsAlpha * fsTheta;
            forceCoeff = fsBeta * fsTheta;
            break;
    }
    
    for (i=0; i<nb; i++) {
        s = 0.0;
        for (j=0; j<n; j++) {
            s = s + (1.00/dt) * massMatrix[i][j] * prevSolution[j] - forceCoeff * stiffMatrix[i][j] * prevSolution[j];
        }
        force[i] = force[i] + s;
        
        for (j=0; j<nb; j++) {
            stiffMatrix[i][j] = massCoeff * stiffMatrix[i][j] + (1.0/dt)*massMatrix[i][j];
        }
    }
    
}

/**************************************************************************
    Constant time step BDF integration scheme to elementwise matrix entry.
**************************************************************************/
-(void)bdfLocalInSolution:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * __nonnull * __nonnull)massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution:(double * __nonnull * __nonnull)prevSolution order:(int)order rows:(int * __nonnull)rows cols:(int * __nonnull)cols {
    
    int i, j, nb1, nb2;
    double s;
    
    if (rows != NULL && cols != NULL) {
        nb1 = *rows;
        nb2 = *cols;
    } else {
        fprintf(stdout, "bdfLocalInSolution: dimension(s) missing.");
        return;
    }
    
    switch (order) {
        case 1:
            for (i=0; i<nb1; i++) {
                s = 0.0;
                for (j=0; j<n; j++) {
                    s = s + (1.0/dt)*massMatrix[i][j]*prevSolution[j][0];
                }
                force[i] = force[i] + s;
                
                for (j=0; j<nb2; j++) {
                    stiffMatrix[i][j] = (1.0/dt)*massMatrix[i][j] + stiffMatrix[i][j];
                }
            }
            break;
        case 2:
            for (i=0; i<nb1; i++) {
                s = 0.0;
                for (j=0; j<n; j++) {
                    s = s + (1.0/dt)*massMatrix[i][j] * (2.0*prevSolution[j][0] - 0.5*prevSolution[j][1]);
                }
                force[i] = force[i] + s;
                
                for (j=0; j<nb2; j++) {
                    stiffMatrix[i][j] = (1.5/dt)*massMatrix[i][j] + stiffMatrix[i][j];
                }
            }
            break;
        case 3:
            for (i=0; i<nb1; i++) {
                s = 0.0;
                for (j=0; j<n; j++) {
                    s = s + (1.0/dt)*massMatrix[i][j] * (3.0*prevSolution[j][0] - (3.0/2.0)*prevSolution[j][1] + (1.0/3.0)*prevSolution[j][2]);
                }
                force[i] = force[i] + s;
                
                for (j=0; j<nb2; j++) {
                    stiffMatrix[i][j] = (11.0/(6.0*dt))*massMatrix[i][j] + stiffMatrix[i][j];
                }
            }
            break;
        case 4:
            for (i=0; i<nb1; i++) {
                s = 0.0;
                for (j=0; j<n; j++) {
                    s = s + (1.0/dt)*massMatrix[i][j] * (4.0*prevSolution[j][0] - 3.0*prevSolution[j][1] + (4.0/3.0)*prevSolution[j][2] - (1.0/4.0)*prevSolution[j][3]);
                }
                force[i] = force[i] + s;
                
                for (j=0; j<nb2; j++) {
                    stiffMatrix[i][j] = (25.0/(12.0*dt))*massMatrix[i][j] + stiffMatrix[i][j];
                }
            }
            break;
        case 5:
            for (i=0; i<nb1; i++) {
                s = 0.0;
                for (j=0; j<n; j++) {
                    s = s + (1.0/dt)*massMatrix[i][j] * (5.0*prevSolution[j][0] - 5.0*prevSolution[j][1] + (10.0/3.0)*prevSolution[j][2] - (5.0/4.0)*prevSolution[j][3] + (1.0/5.0)*prevSolution[j][4]);
                }
                force[i] = force[i] + s;
                
                for (j=0; j<nb2; j++) {
                    stiffMatrix[i][j] = (137.0/(60.0*dt))*massMatrix[i][j] + stiffMatrix[i][j];
                }
            }
            break;
        default:
            fatal("FEMTimeIntegration:bdfLocalInSolution", "Invalid order BDF.");
            break;
    }
}

/***************************************************************************
    Variable time step BDF integration scheme to elementwise matrix entry.
***************************************************************************/
-(void)vbdfLocalInSolution:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dts:(double * __nonnull)dts massMatrix:(double * __nonnull * __nonnull)massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution:(double * __nonnull * __nonnull)prevSolution order:(int)order rows:(int * __nonnull)rows cols:(int * __nonnull)cols {
    
    int i, j, k, nb1, nb2;
    double s, a[4];
    
    if (rows != NULL && cols != NULL) {
        nb1 = *rows;
        nb2 = *cols;
    } else {
        fprintf(stdout, "vbdfLocalInSolution: dimension(s) missing.");
        return;
    }

    memset( a, 0.0, sizeof(a) );
    
    a[0] = 1.0 / dts[0];
    a[1] = -1.0 / dts[0];
    if (order >= 2) {
        a[0] = a[0] + 1.0 / (dts[0]+dts[1]);
        a[1] = a[1] - (1.0 + dts[0]/dts[1]) / (dts[0]+dts[1]);
        a[2] = (dts[0]/dts[1]) / (dts[0]+dts[1]);
    }
    if (order >= 3) {
        a[0] = a[0] + 1.0 / (dts[0]+dts[1]+dts[2]);
        a[1] = a[1] - ( 1.0 + dts[0]/dts[1] * (1.0+(dts[0]+dts[1])/(dts[1]+dts[2])) ) / (dts[0]+dts[1]+dts[2]);
        a[2] = a[2] + ( dts[0]/dts[1]*(1.0+(dts[0]+dts[1])/(dts[1]+dts[2]))+dts[0]/dts[2]*(dts[0]+dts[1])/(dts[1]+dts[2]) ) / (dts[0]+dts[1]+dts[2]);
        a[3] = -(dts[0]/dts[2]) * (dts[0]+dts[1]) / (dts[1]+dts[2]) / (dts[0]+dts[1]+dts[2]);
    }
    if (order > 3) {
        fprintf(stdout, "FEMTimeIntegration:vbdfLocalInSolution: variable time step BDF implemented only up to order 3.\n");
    }
    
    for (i=0; i<nb1; i++) {
        s = 0.0;
        for (k=0; k<min(order, 3); k++) {
            for (j=0; j<n; j++) {
                s = s - a[k+1]*massMatrix[i][j] * prevSolution[j][k];
            }
        }
        force[i] = force[i] + s;
        
        for (j=0; j<nb2; j++) {
            stiffMatrix[i][j] = a[0]*massMatrix[i][j] + stiffMatrix[i][j];
        }
    }
}

/****************************************************************************************************************
    Apply newmark beta scheme to local matrix equation. This is used also for Implicit Euler with beta=1.0 and
    Crank-Nicolson with beta=0.5.
****************************************************************************************************************/
-(void)newMarkBetaInSolution:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * __nonnull * __nonnull)massMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution:(double * __nonnull)prevSolution beta:(double)beta rows:(int * __nonnull)rows cols:(int * __nonnull)cols {
    
    int i, j, nb1, nb2;
    double s;
  
    if (rows != NULL && cols != NULL) {
        nb1 = *rows;
        nb2 = *cols;
    } else {
        fprintf(stdout, "newMarkBetaInSolution: dimension(s) missing.");
        return;
    }

    for (i=0; i<nb1; i++) {
        s = 0.0;
        for (j=0; j<n; j++) {
            s = s + (1.0/dt) * massMatrix[i][j] * prevSolution[j] - (1.0-beta) * stiffMatrix[i][j] * prevSolution[j];
        }
        for (j=0; j<nb2; j++) {
            stiffMatrix[i][j] = beta * stiffMatrix[i][j] + (1.0/dt)*massMatrix[i][j];
        }
        force[i] = force[i] + s;
    }
}

/*********************************************************************************
    Second order Bossok time integration scheme to the elememtwise matrix entry.
*********************************************************************************/
-(void)bossakSecondOrder:(FEMSolution * __nonnull)solution numberOfNodes:(int)n dt:(double)dt massMatrix:(double * __nonnull * __nonnull)massMatrix dampMatrix:(double * __nonnull * __nonnull)dampMatrix stiffMatrix:(double * __nonnull * __nonnull)stiffMatrix force:(double * __nonnull)force prevSolution1:(double * __nonnull)x prevSolution2:(double * __nonnull)v prevSolution3:(double * __nonnull)a alpha:(double)alpha rows:(int * __nonnull)rows cols:(int * __nonnull)cols {
    
    int i, j, nb1, nb2;
    double beta, gamma, s;
    
    if (rows != NULL && cols != NULL) {
        nb1 = *rows;
        nb2 = *cols;
    } else {
        fprintf(stdout, "bossakSecondOrder: dimension(s) missing.");
        return;
    }
    
    nb1 = min(n, nb1);
    nb2 = min(n, nb2);
    
    gamma = 0.5 - alpha;
    beta = pow((1.0 - alpha), 2.0) / 4.0;
    for (i=0; i<nb1; i++) {
        s = 0.0;
        for (j=0; j<nb2; j++) {
            s = s + ( (1.0 - alpha) / (beta * pow(dt, 2.0)) ) * massMatrix[i][j] * x[j];
            s = s + ( (1.0 - alpha) / (beta*dt) ) * massMatrix[i][j] * v[j];
            s = s - ( (1.0 - alpha) * (1.0 - 1.0 / (2.0*beta)) + alpha ) * massMatrix[i][j] * a[j];
            
            s = s + ( gamma / (beta*dt) ) * dampMatrix[i][j] * x[j];
            s = s + ( gamma/beta - 1.0) * dampMatrix[i][j] * v[j];
            s = s - ( (1.0 - gamma) + gamma * (1.0 - 1.0 / (2.0*beta)) ) * dt * dampMatrix[i][j] * a[j];
            
            stiffMatrix[i][j] = stiffMatrix[i][j] + ( (1.0 - alpha) / (beta * pow(dt, 2.0)) ) * massMatrix[i][j]
            + (gamma / (beta*dt)) * dampMatrix[i][j];
        }
        force[i] = force[i] + s;
    }
}

@end
