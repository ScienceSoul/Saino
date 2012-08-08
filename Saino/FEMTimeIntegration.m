//
//  FEMTimeIntegration.m
//  Saino
//
//  Created by Hakime Seddik on 14/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMTimeIntegration.h"

@implementation FEMTimeIntegration

-(void)fractionStep:(FEMSolution *)solution: (int)n: (double)dt: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double *)prevSolution: (int *)rows {
    
    int i, j, nb;
    double s, fsStep, fsTheta, fsdTheta, fsAlpha, fsBeta, massCoeff, forceCoeff;

    if (rows != NULL) {
        nb = *rows;
    } else {
        nb = solution.mesh.maxElementNodes;
    }
    
    if ((solution.solutionInfo)[@"fsstep"] != nil) {
        fsStep = [(solution.solutionInfo)[@"fsstep"] doubleValue];
    }
    if ((solution.solutionInfo)[@"fstheta"] != nil) {
        fsTheta = [(solution.solutionInfo)[@"fstheta"] doubleValue];
    }
    if ((solution.solutionInfo)[@"fsdtheta"] != nil) {
        fsdTheta = [(solution.solutionInfo)[@"fsdtheta"] doubleValue];
    }
    if ((solution.solutionInfo)[@"fsalpha"] != nil) {
        fsAlpha = [(solution.solutionInfo)[@"fsalpha"] doubleValue];
    }
    if ((solution.solutionInfo)[@"fsbeta"] != nil) {
        fsBeta = [(solution.solutionInfo)[@"fsbeta"] doubleValue];
    }
    
    switch ((int)fsStep) {
        case 1:
            massCoeff = fsAlpha * fsTheta;
            forceCoeff = fsBeta + fsTheta;
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

-(void)bdfLocal:(FEMSolution *)solution: (int)n: (double)dt: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double **)prevSolution: (int)order: (int *)rows: (int *)cols {
    
    int i, j, nb1, nb2;
    double s;
    
    if (rows != NULL) {
        nb1 = *rows;
    } else {
        nb1 = solution.mesh.maxElementNodes;
    }
    
    if (cols != NULL) {
        nb2 = *cols;
    } else {
        nb2 = solution.mesh.maxElementNodes;
    }
    
    switch (order) {
        case 1:
            for (i=0; i<nb1; i++) {
                s = 0.0;
                for (j=0; j<n; i++) {
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
                    s = s + (1.0/dt)*massMatrix[i][j] + (2.0*prevSolution[j][0] - 0.5*prevSolution[j][1]); 
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
            errorfunct("bdfLocal", "Invalid order BDF.");
            break;
    }
}

-(void)vbdfLocal:(FEMSolution *)solution: (int)n: (double *)dts: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double **)prevSolution: (int)order: (int *)rows: (int *)cols {
/***********************************************************************************************************************************************
    Variable time step BDF
***********************************************************************************************************************************************/
    
    int i, j, k, nb1, nb2;
    double s, *a;
    
    if (rows != NULL) {
        nb1 = *rows;
    } else {
        nb1 = solution.mesh.maxElementNodes;
    }
    
    if (cols != NULL) {
        nb2 = *cols;
    } else {
        nb2 = solution.mesh.maxElementNodes;
    }
    
    a = doublevec(0, 3);
    memset( a, 0.0, (4*sizeof(a)) );
    
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
        warnfunct("vbdfLocal", "Variable time step BDF implemented only up to order 3.");
    }
    
    for (i=0; i<nb1; i++) {
        s = 0.0;
        for (k=0; k<min(order, 3); k++) {
            for (j=0; j<n; j++) {
                s = s - a[k+1]*massMatrix[i][k] * prevSolution[j][k];
            }
        }
        force[i] = force[i] + s;
        
        for (j=0; j<nb2; j++) {
            stiffMatrix[i][j] = a[0]*massMatrix[i][j] + stiffMatrix[i][j];
        }
    }
    
    free_dvector(a, 0, 3);
}

-(void)newMarkBeta:(FEMSolution *)solution: (int)n: (double)dt: (double **)massMatrix: (double **)stiffMatrix: (double *)force: (double *)prevSolution: (double)beta: (int *)rows {
    
    int i, j, nb;
    double s;
    
    if (rows != NULL) {
        nb = *rows;
    } else {
        nb = solution.mesh.maxElementNodes;
    }
    
    for (i=0; i<nb; i++) {
        s = 0.0;
        for (j=0; j<n; j++) {
            s = s + (1.0/dt) * massMatrix[i][j] * prevSolution[j] - (1.0-beta) * stiffMatrix[i][j] * prevSolution[j];
        }
        for (j=0; j<nb; j++) {
            stiffMatrix[i][j] = beta * stiffMatrix[i][j] + (1.0/dt)*massMatrix[i][j];
        }
        force[i] = force[i] + s;
    }
}

@end
