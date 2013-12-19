//
//  GaussIntegration.h
//  Saino
//
//  Created by Seddik hakime on 16/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Numerics.h"
#include "Constructors.h"
#include "memory.h"

#define MAXN 13
#define MAX_INTEGRATION_POINTS MAXN*MAXN*MAXN

void DerivPoly(int n, double *Q, double *P);
double EvalPoly(int n, double *P, double x);
void RefineRoots(int n, double *P, double *Q, double *Points);
void GaussQuadraturePoints1D(int n);
void GaussQuadratureInit(void);
void GaussQuadratureDeallocation(void);
GaussIntegrationPoints* GaussQuadrature0D(int np);
GaussIntegrationPoints* GaussQuadrature1D(int np);
GaussIntegrationPoints* GaussQuadratureTriangle(int np);
GaussIntegrationPoints* GaussQuadratureQuad(int np);
GaussIntegrationPoints* GaussQuadratureTetra(int np);
GaussIntegrationPoints* GaussQuadraturePyramid(int np);
GaussIntegrationPoints* GaussQuadratureWedge(int np);
GaussIntegrationPoints* GaussQuadratureBrick(int np);
GaussIntegrationPoints* GaussQuadrature(Element_t *element, int *np, int *relOrder);

