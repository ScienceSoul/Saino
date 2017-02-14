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

extern inline void invertMatrix3x3(double * __nonnull GI, double * __nonnull covariantMetricTensor, double detG);
extern inline double det3x3(double * __nonnull * __nonnull a);
extern inline void derivatives(double * __nonnull dx, Element_t * __nonnull element, int nDOFs, Nodes_t * __nonnull nodes, double * __nonnull dLBasisdx);
extern inline void covariantMetric(double * __nonnull covariantMetricTensor, double * __nonnull dx, Element_t * __nonnull element, int nDOFs, Nodes_t * __nonnull nodes, int meshDimension, double * __nonnull dLBasisdx);
extern inline double detJ(double * __nonnull covariantMetricTensor, double * __nonnull dx, Element_t * __nonnull element, Nodes_t * __nonnull nodes, int meshDimension, double u, double v, double w, double *  __nonnull dLBasisdx);
extern inline bool contravariantMetric(double * __nonnull elementMetric, double * __nonnull dx, Element_t * __nonnull element, Nodes_t * __nonnull nodes, int meshDimension, double u, double v, double w, double * __nonnull dLBasisdx);
extern inline bool localtoGlobalMap(double * __nonnull ltoGMap, Element_t * __nonnull element, Nodes_t * __nonnull nodes, int meshDimension, double u, double v, double w, double * __nonnull dLBasisdx);
extern inline void globalSecondDerivatives(double * __nonnull elementMetric, Element_t * __nonnull element, Nodes_t * __nonnull nodes, int meshDimension, double u, double v, double w, double * __nonnull f, double * __nonnull dLBasisdx, double * __nonnull values);

