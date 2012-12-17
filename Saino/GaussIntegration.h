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

#define MAXN 14
#define MAX_INTEGRATION_POINTS MAXN*MAXN*MAXN

void DerivPoly(int n, double *Q, double *P);
double EvalPoly(int n, double *P, double x);
void RefineRoots(int n, double *P, double *Q, double *Points);
void GaussQuadraturePoints1D(int n);
void GaussQuadratureInit(GaussIntegrationPoints *pt);
void GaussQuadrature0D(int n, GaussIntegrationPoints *pt);
void GaussQuadrature1D(int n, GaussIntegrationPoints *pt);
void GaussQuadratureTriangle(int n, GaussIntegrationPoints *pt);
void GaussQuadratureQuad(int np, GaussIntegrationPoints *pt);
void GaussQuadratureTetra(int n, GaussIntegrationPoints *pt);
void GaussQuadraturePyramid(int np, GaussIntegrationPoints *pt);
void GaussQuadratureWedge(int np, GaussIntegrationPoints *pt);
void GaussQuadratureBrick(int np, GaussIntegrationPoints *pt);
GaussIntegrationPoints* GaussQuadrature(Element_t *elm, ...);

