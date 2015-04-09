//
//  Walls.h
//  Saino
//
//  Created by Seddik hakime on 23/05/2014.
//  Copyright (c) 2014 Institute of Low Temperature Science. All rights reserved.
//

#include <stdio.h>
#include <math.h>

#include "Utils.h"


/*************************************************************************************************************
    Give difference between the tangential velocity given by Reichardts's wall law and the tangential 
    velocity of the previous iteration
 
    Inputs:
        ufric     -> friction velocity
        ut        -> tangential velocity
        density   -> density
        viscos    -> viscosity
        dist      -> distance from the wall
 *************************************************************************************************************/
inline double wall_law(double ufric, double ut, double densit, double viscos, double dist, double rough) {
    
    double yplus;
    double kappa = 0.41;
    
    yplus = densit * ufric * dist / viscos;
    
    // Reichardts law
    return (ufric / kappa) * log(1.0 + 0.4 * yplus) + ufric * 7.8 *
    ( 1.0 - exp(-yplus / 11.0) - (yplus / 11.0) * exp(-0.33 * yplus) ) - ut;
}

/*************************************************************************************************************
    Calculate derivative of the wall law
 
    Inputs:
        ufric     -> friction velocity
        ut        -> tangential velocity
        density   -> density
        viscos    -> viscosity
        dist      -> distance from the wall
 *************************************************************************************************************/
inline double d_wall_law(double ufric, double ut, double densit, double viscos, double dist, double rough) {
    
    double yplus;
    double kappa = 0.41;
    
    yplus = densit * ufric * dist / viscos;
    
    // Reichardts law
    return log(1.0 + 0.4 * yplus) /  kappa + (0.4 / kappa) * yplus / (1.0 + 0.4 * yplus)
           + 7.8 * ( 1.0 - exp(-yplus / 11/0) - (yplus / 11.0) * exp(-0.33 * yplus) )
           + 7.8 * (yplus / 11.0) * ( exp(-yplus / 11.0) - exp(-0.33 * yplus) + 0.33 * yplus * exp(-0.33 * yplus) );
}

/*************************************************************************************************************
    Solve the friction velocity of the previous iteration based on the wall law
 
    Inputs:
        density -> density
        viscos  -> viscosity
        dist    -> distance from the wall
        ut      -> tangential velocity of the previous iteration
 
    Outputs:
        ufric   -> friction velocity
        dfx     -> derivative of the wall law
*************************************************************************************************************/
inline void solve_ufric(double densit, double viscos, double dist, double rough, double ut, double *ufric, double *dfx) {
    
    int iter, maxiter = 100;
    double fx, tauw, yplus;
    double tol = 1.0e-14;
    
    // Default value
    tauw = ut / dist;
    *ufric = sqrt(tauw / densit);
    
    for (iter=1; iter<maxiter; iter++) {
        fx = wall_law(*ufric, ut, densit, viscos, dist, rough);
        *dfx = d_wall_law(*ufric, ut, densit, viscos, dist, rough);
        
        // Newton steps
        if (*dfx == 0.0) errorfunct("solve_ufric", "dfx = 0");
        *ufric = *ufric - fx / *dfx;
        yplus = densit * *ufric * dist / viscos;
        if (fabs(fx) <= tol) break;
    }
    
    if (fabs(fx) > 1.0e-09) {
        printf("fx: %f\n", fx);
        warnfunct("solve_ufric", "problems in value of fx");
    }
}


/*************************************************************************************************************
    Calculate the boundary values of turbulent kinetic energy and its dissipation bases on the wall law
 
    Inputs:
        ut      -> tangential velocity of the previous iteration
        dist    -> distance from the wall
        viscos  -> viscosity
        density -> density
 
    Outputs:
        tx      -> turbulent kinetic energy
        teps    -> turbulent kinetic energy dissipation
        tomg    -> ...
*************************************************************************************************************/
inline void kewall(double *tk, double *teps, double *tomg, double ut, double dist, double rough, double viscos, double densit) {
    
    double alpha, dfx, omegaPlus, tomgl, tomgt, ufric, utlocal, yplus;
    double cmyy   = 0.09;
    double karman = 0.41;
    double small  = 1.0e-10;
    
    utlocal = max(ut, small);
    solve_ufric(densit, viscos, dist, rough, utlocal, &ufric, &dfx);
    
    yplus = densit * ufric * dist / viscos;
    alpha = min(1.0, yplus / 10.0);
    
    *tk = pow(ufric, 2.0) / sqrt(cmyy) * alpha;
    *teps = pow(ufric, 3.0) / (karman * dist) * min( 1.0, alpha + 0.2 * karman * (1.0 - pow(alpha, 2.0)) / sqrt(cmyy) );
    
    omegaPlus = 6.0 / (0.072 * pow(yplus, 2.0));
    tomgl = densit * pow(ufric, 2.0) * omegaPlus / viscos;
    tomgt = ufric / (sqrt(cmyy) * karman * dist);
    
    if (yplus < 4) {
        *tomg = tomgl;
    } else if (yplus < 32) {
        *tomg = sqrt(pow(tomgl, 2.0) + pow(tomgt, 2.0));
    } else {
        *tomg = tomgt;
    }
}