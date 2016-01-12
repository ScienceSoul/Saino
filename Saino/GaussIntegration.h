//
//  GaussIntegration.h
//  Saino
//
//  Created by Seddik hakime on 16/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#include <Accelerate/Accelerate.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>

#include "memory.h"
#include "Numerics.h"
#include "Constructors.h"

#define MAXN 13
#define MAX_INTEGRATION_POINTS MAXN*MAXN*MAXN

inline void DerivPoly(int n, double * __nonnull Q, double * __nonnull P) {
    
    int i;
    
    for(i=0;i<n;i++) {
        Q[i] = P[i] * (n-(i+1)+1);
    }
}

inline double EvalPoly(int n, double * __nonnull P, double x) {
    int i;
    double s;
    
    s = 0.0;
    for(i=0;i<n+1;i++) {
        s = s * x + P[i];
    }
    
    return s;
}

inline void RefineRoots(int n, double * __nonnull P, double * __nonnull Q, double * __nonnull Points) {
    
    int i, j;
    double x, s;
    const int MaxIter = 100;
    
    for(i=0;i<n;i++) {
        x = Points[i];
        for(j=0;j<MaxIter;j++) {
            s = EvalPoly(n,P,x) / EvalPoly(n-1, Q, x);
            x = x - s;
            if( fabs(s) <= fabs(x)*DBL_EPSILON ) break;
        }
        if( fabs(EvalPoly(n,P,x)) < fabs(EvalPoly(n,P,Points[i])) ) {
            if( fabs(x-Points[i]) < 1.0e-8 ) Points[i] = x;
        }
    }
}

void GaussQuadraturePoints1D(int n);
void GaussQuadratureInit(void);
void GaussQuadratureDeallocation(void);
GaussIntegrationPoints * __nonnull GaussQuadrature0D(int np);
GaussIntegrationPoints * __nonnull GaussQuadrature1D(int np);
GaussIntegrationPoints * __nonnull GaussQuadratureTriangle(int np);
GaussIntegrationPoints * __nonnull GaussQuadratureQuad(int np);
GaussIntegrationPoints * __nonnull GaussQuadratureTetra(int np);
GaussIntegrationPoints * __nonnull GaussQuadraturePyramid(int np);
GaussIntegrationPoints * __nonnull GaussQuadratureWedge(int np);
GaussIntegrationPoints * __nonnull GaussQuadratureBrick(int np);

inline GaussIntegrationPoints * __nonnull GaussQuadrature(Element_t * __nonnull element, int * __nullable np, int * __nullable relOrder) {
    
    int n=0, eldim, p1d;
    GaussIntegrationPoints *IntegPoint = NULL;
    bool pElement;
    
    pElement = (element->Pdefs != NULL) ? true : false;
    
    if (np != NULL) {
        n = *np;
    } else if (relOrder != NULL) {
        if (pElement == true) {
            n = element->Pdefs->GaussPoints;
            if (*relOrder == 0) {
                //No ops
            } else {
                eldim = element->Type.dimension;
                p1d = (int)round(pow((double)n, (1.0/(double)eldim))) + *relOrder;
                if (p1d < 1) {
                    fatal("GaussIntegrationPoints", "Number of integration points must remain positive.");
                }
                n = (int)pow((double)p1d, (double)eldim);
            }
        } else {
            if (*relOrder == 0) {
                n = element->Type.GaussPoints;
            } else if (*relOrder == 1) {
                n = element->Type.GaussPoints2;
            } else if (*relOrder == -1) {
                n = element->Type.GaussPoints0;
            } else {
                warning("GaussQuadrature", "RelOrder can only be {-1, 0, 1}.");
            }
        }
    } else {
        if (pElement == true) {
            n = element->Pdefs->GaussPoints;
        } else {
            n = element->Type.GaussPoints;
        }
    }
    
    switch (element->Type.ElementCode / 100) {
            
        case 1:
            IntegPoint = GaussQuadrature0D(n);
            break;
            
        case 2:
            IntegPoint = GaussQuadrature1D(n);
            break;
            
        case 3:
            // TODO: add support for p element
            IntegPoint = GaussQuadratureTriangle(n);
            break;
            
        case 4:
            IntegPoint = GaussQuadratureQuad(n);
            break;
            
        case 5:
            // TODO: add support for p element
            IntegPoint = GaussQuadratureTetra(n);
            break;
            
        case 6:
            // TODO: add support for p element
            IntegPoint = GaussQuadraturePyramid(n);
            break;
            
        case 7:
            // TODO: add support for p element
            IntegPoint = GaussQuadratureWedge(n);
            break;
            
        case 8:
            IntegPoint = GaussQuadratureBrick(n);
            break;
            
    }
    
    return IntegPoint;
}

