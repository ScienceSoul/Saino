//===----------------------------------------------------------------------===//
//  NodalBasisFunction.c
//  Saino
//
//  Created by Seddik hakime on 30/06/11.
//  Copyright (c) 2011 ScienceSoul. All rights reserved.
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

#include "NodalBasisFunctions.h"

extern inline void NodalBasisFunctions1D (double * _Nonnull y, Element_t * _Nonnull element, double u);
extern inline void NodalBasisFunctions2D(double * _Nonnull y, Element_t * _Nonnull element, double u, double v);
extern inline double InterpolateInElement3D(Element_t * _Nonnull element, double * _Nonnull x, double u, double v, double w);
extern inline void NodalBasisFunctions3D(double * _Nonnull y, Element_t * _Nonnull element, double u, double v, double w);
extern inline void NodalBasisFunctions(int n, double * _Nonnull Basis, Element_t * _Nonnull element, double u, double v,double w);

extern inline void NodalFirstDerivatives1D(double * _Nonnull y, Element_t * _Nonnull element, double u);
extern inline void NodalFirstDerivatives2D(double * _Nonnull y, Element_t * _Nonnull element, double u, double v);
extern inline double FirstDerivativeInU3D(Element_t * _Nonnull element, double * _Nonnull x, double u, double v, double w);
extern inline double FirstDerivativeInV3D(Element_t * _Nonnull element, double * _Nonnull x, double u, double v, double w);
extern inline double FirstDerivativeInW3D(Element_t * _Nonnull element, double * _Nonnull x, double u, double v, double w);
extern inline void NodalFirstDerivatives3D(double * _Nonnull y, Element_t * _Nonnull element, double u, double v, double w);
extern inline void NodalFirstDerivatives(int n, double * _Nonnull dLBasisdx, Element_t * _Nonnull element, double u, double v, double w);

extern inline double SecondDerivatives1D(Element_t * _Nonnull element, double * _Nonnull nodes, double u);
extern inline void SecondDerivatives2D(double * _Nonnull ddx, Element_t * _Nonnull element, double * _Nonnull nodes, double u, double v);
extern inline void SecondDerivatives3D(double * _Nonnull ddx, Element_t * _Nonnull element, double * _Nonnull nodes, double u, double v, double w);
