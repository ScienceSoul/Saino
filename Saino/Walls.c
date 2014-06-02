//
//  Walls.c
//  Saino
//
//  Created by Seddik hakime on 23/05/2014.
//  Copyright (c) 2014 Institute of Low Temperature Science. All rights reserved.
//

#include "Walls.h"

extern inline double wall_law(double ufric, double ut, double densit, double viscos, double dist, double rough);
extern inline double d_wall_law(double ufric, double ut, double densit, double viscos, double dist, double rough);
extern inline void solve_ufric(double densit, double viscos, double dist, double rough, double ut, double *ufric, double *dfx);
extern inline void kewall(double *tk, double *teps, double *tomg, double ut, double dist, double rough, double viscos, double densit);