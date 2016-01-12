//
//  NodalBasisFunction.c
//  Saino
//
//  Created by Seddik hakime on 30/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#include "NodalBasisFunctions.h"

extern inline void NodalBasisFunctions1D (double * __nonnull y, Element_t * __nonnull element, double u);
extern inline void NodalBasisFunctions2D(double * __nonnull y, Element_t * __nonnull element, double u, double v);
extern inline double InterpolateInElement3D(Element_t * __nonnull element, double * __nonnull x, double u, double v, double w);
extern inline void NodalBasisFunctions3D(double * __nonnull y, Element_t * __nonnull element, double u, double v, double w);
extern inline void NodalBasisFunctions(int n, double * __nonnull Basis, Element_t * __nonnull element, double u, double v,double w);

extern inline void NodalFirstDerivatives1D(double * __nonnull y, Element_t * __nonnull element, double u);
extern inline void NodalFirstDerivatives2D(double * __nonnull y, Element_t * __nonnull element, double u, double v);
extern inline double FirstDerivativeInU3D(Element_t * __nonnull element, double * __nonnull x, double u, double v, double w);
extern inline double FirstDerivativeInV3D(Element_t * __nonnull element, double * __nonnull x, double u, double v, double w);
extern inline double FirstDerivativeInW3D(Element_t * __nonnull element, double * __nonnull x, double u, double v, double w);
extern inline void NodalFirstDerivatives3D(double * __nonnull y, Element_t * __nonnull element, double u, double v, double w);
extern inline void NodalFirstDerivatives(int n, double * __nonnull dLBasisdx, Element_t * __nonnull element, double u, double v, double w);

extern inline double SecondDerivatives1D(Element_t * __nonnull element, double * __nonnull nodes, double u);
extern inline void SecondDerivatives2D(double * __nonnull ddx, Element_t * __nonnull element, double * __nonnull nodes, double u, double v);
extern inline void SecondDerivatives3D(double * __nonnull ddx, Element_t * __nonnull element, double * __nonnull nodes, double u, double v, double w);