//
//  NodalBasisFunctions.h
//  Saino
//
//  Created by Seddik hakime on 30/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Constructors.h"

void NodalBasisFunctions1D (double *y, Element_t *element, double u);
void NodalBasisFunctions2D(double *y, Element_t *element, double u, double v);
double InterpolateInElement3D(Element_t *element, double *x, double u, double v, double w);
void NodalBasisFunctions3D(double *y, Element_t *element, double u, double v, double w);
void NodalBasisFunctions(int n, double *Basis, Element_t *element, double u, double v,double w);

void NodalFirstDerivatives1D(double **y, Element_t *element, double u);
void NodalFirstDerivatives2D(double **y, Element_t *element, double u, double v);
double FirstDerivativeInU3D(Element_t *element, double *x, double u, double v, double w);
double FirstDerivativeInV3D(Element_t *element, double *x, double u, double v, double w);
double FirstDerivativeInW3D(Element_t *element, double *x, double u, double v, double w);
void NodalFirstDerivatives3D(double **y, Element_t *element, double u, double v, double w);
void NodalFirstDerivatives(int n, double **dLBasisdx, Element_t *element, double u, double v, double w);
double SecondDerivatives1D(Element_t* element, double *nodes, double u);
void SecondDerivatives2D(double **ddx, Element_t* element, double *nodes, double u, double v);
void SecondDerivatives3D(double **ddx, Element_t* element, double *nodes, double u, double v, double w);
