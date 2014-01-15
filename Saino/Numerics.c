//
//  Numerics.c
//  Saino
//
//  Created by Seddik hakime on 26/12/2013.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#include "Numerics.h"

extern inline void invertMatrix3x3(double *GI, double *covariantMetricTensor, double detG);
extern inline void derivatives(double *dx, Element_t *element, int nDOFs, Nodes_t *nodes, double *dLBasisdx);
extern inline void covariantMetric(double *covariantMetricTensor, double *dx, Element_t *element, int nDOFs, Nodes_t *nodes, int meshDimension, double *dLBasisdx);
extern inline double detJ(double *covariantMetricTensor, double *dx, Element_t *element, Nodes_t *nodes, int meshDimension, double u, double v, double w, double *dLBasisdx);
extern inline bool contravariantMetric(double *elementMetric, double *dx, Element_t *element, Nodes_t *nodes, int meshDimension, double u, double v, double w, double *dLBasisdx);
extern inline bool localtoGlobalMap(double *ltoGMap, Element_t *element, Nodes_t* nodes, int meshDimension, double u, double v, double w, double *dLBasisdx);
extern inline void globalSecondDerivatives(double *elementMetric, Element_t *element, Nodes_t *nodes, int meshDimension, double u, double v, double w, double *f, double *dLBasisdx, double *values);

