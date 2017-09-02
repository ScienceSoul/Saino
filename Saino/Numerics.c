//===----------------------------------------------------------------------===//
//  Numerics.c
//  Saino
//
//  Created by Seddik hakime on 26/12/2013.
//  Copyright (c) 2013 ScienceSoul. All rights reserved.
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

#include "Numerics.h"

extern inline void invertMatrix3x3(double * _Nonnull GI, double * _Nonnull covariantMetricTensor, double detG);
extern inline double det3x3(double * _Nonnull * _Nonnull a);
extern inline void derivatives(double * _Nonnull dx, Element_t * _Nonnull element, int nDOFs, Nodes_t * _Nonnull nodes, double * _Nonnull dLBasisdx);
extern inline void covariantMetric(double * _Nonnull covariantMetricTensor, double * _Nonnull dx, Element_t * _Nonnull element, int nDOFs, Nodes_t * _Nonnull nodes, int meshDimension, double * _Nonnull dLBasisdx);
extern inline double detJ(double * _Nonnull covariantMetricTensor, double * _Nonnull dx, Element_t * _Nonnull element, Nodes_t * _Nonnull nodes, int meshDimension, double u, double v, double w, double *  _Nonnull dLBasisdx);
extern inline bool contravariantMetric(double * _Nonnull elementMetric, double * _Nonnull dx, Element_t * _Nonnull element, Nodes_t * _Nonnull nodes, int meshDimension, double u, double v, double w, double * _Nonnull dLBasisdx);
extern inline bool localtoGlobalMap(double * _Nonnull ltoGMap, Element_t * _Nonnull element, Nodes_t * _Nonnull nodes, int meshDimension, double u, double v, double w, double * _Nonnull dLBasisdx);
extern inline void globalSecondDerivatives(double * _Nonnull elementMetric, Element_t * _Nonnull element, Nodes_t * _Nonnull nodes, int meshDimension, double u, double v, double w, double * _Nonnull f, double * _Nonnull dLBasisdx, double * _Nonnull values);

