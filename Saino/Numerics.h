//
//  Numerics.h
//  Saino
//
//  Created by Seddik hakime on 16/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Constructors.h"
#import "NodalBasisFunctions.h"
#include "Utils.h"

typedef struct {
    
    int n;
    double *u, *v, *w, *s;
    
} GaussIntegrationPoints;

inline void invertMatrix3x3(double *GI, double *covariantMetricTensor, double detG) {
    
    double s;
    
    s = 1.0 / detG;
    
    GI[0] = s * ( covariantMetricTensor[4]*covariantMetricTensor[8] - covariantMetricTensor[7]*covariantMetricTensor[5] );
    GI[3] = -s * ( covariantMetricTensor[3]*covariantMetricTensor[8] - covariantMetricTensor[6]*covariantMetricTensor[5] );
    GI[6] = s * ( covariantMetricTensor[3]*covariantMetricTensor[7] - covariantMetricTensor[6]*covariantMetricTensor[4] );
    
    GI[1] = -s * ( covariantMetricTensor[1]*covariantMetricTensor[8] - covariantMetricTensor[7]*covariantMetricTensor[2] );
    GI[4] = s * ( covariantMetricTensor[0]*covariantMetricTensor[8] - covariantMetricTensor[6]*covariantMetricTensor[2] );
    GI[7] = -s * ( covariantMetricTensor[0]*covariantMetricTensor[7] - covariantMetricTensor[6]*covariantMetricTensor[1] );
    
    GI[2] = s * ( covariantMetricTensor[1]*covariantMetricTensor[5] - covariantMetricTensor[4]*covariantMetricTensor[2] );
    GI[5] = -s * ( covariantMetricTensor[0]*covariantMetricTensor[5] - covariantMetricTensor[3]*covariantMetricTensor[2] );
    GI[8] = s * ( covariantMetricTensor[0]*covariantMetricTensor[4] - covariantMetricTensor[3]*covariantMetricTensor[1] );
}

/***********************************************************************************
    Partial derivatives of global coordinates with respect to local coordinates
 
    Arguments:
        NDOFs              -> Number of DOFs
 
        Nodes_t *nodes     -> element nodal cooridnates
 
        double **dLBasisdx -> Derivatives of element basis function with respect
                              to local coordinates
 
        int el             -> Element number
 
***********************************************************************************/
inline void derivatives(double *dx, Element_t *element, int nDOFs, Nodes_t *nodes, double *dLBasisdx) {
    
    int i, j, n, dim;
    double accum1, accum2, accum3;
    
    dim = element->Type.dimension;
    n = min(element->Type.NumberOfNodes, nDOFs);
    
    for (i=0; i<dim; i++) {
        accum1 = 0.0;
        accum2 = 0.0;
        accum3 = 0.0;
        for (j=0; j<n; j++) {
            accum1 = accum1 + (nodes->x[j] * dLBasisdx[3*j+i]);
            accum2 = accum2 + (nodes->y[j] * dLBasisdx[3*j+i]);
            accum3 = accum3 + (nodes->z[j] * dLBasisdx[3*j+i]);
        }
        dx[i] = accum1;
        dx[3+i] = accum2;
        dx[6+i] = accum3;
    }
}

/******************************************************************************************************
    Compute the covariant metric tensor of the element coordinate system
 
    Arguments:
        NDOFs              -> Number of DOFs
 
        Nodes_t *nodes     -> element nodal cooridnates
 
        double **dLBasisdx -> Derivatives of element basis function with respect to
                              local coordinates
 
******************************************************************************************************/
inline void covariantMetric(double *covariantMetricTensor, double *dx, Element_t *element, int nDOFs, Nodes_t *nodes, int meshDimension, double *dLBasisdx) {
    int i, j, k, dim, cdim;
    double s;
    
    dim = element->Type.dimension;
    cdim = meshDimension;
    
    derivatives(dx, element, nDOFs, nodes, dLBasisdx);
    
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<cdim; k++) {
                s = s + ( dx[k*3+i] * dx[k*3+j]);
            }
            covariantMetricTensor[3*i+j] = s;
        }
    }
}

/**********************************************************************************************
    Compute determinant of covariant metric tensor det(J^TJ) and determinant
    of covariant metric tensor det(J^TJ).
 
    Arguments:
        Element_t *element  -> element structure
 
        Nodes_t *nodes      -> element nodal cooridnates
 
        double u, v, w      -> Points at which evaluate the value
 
        int el              -> Element number
 
    Function return value:
        If function failure, element is degenerate
 
**********************************************************************************************/
inline double detJ(double *covariantMetricTensor, double *dx, Element_t *element, Nodes_t *nodes, int meshDimension, double u, double v, double w, double *dLBasisdx) {
    
    int n, dim, cdim;
    double detG;
    bool success = true;
    
    n = element->Type.NumberOfNodes;
    dim = element->Type.dimension;
    cdim = meshDimension;
    
    covariantMetric(covariantMetricTensor, dx, element, n, nodes, meshDimension, dLBasisdx);
    
    detG = 0.0;
    switch (dim) {
            
        case 1:  // Line elements
            detG = covariantMetricTensor[0];
            if (detG <= DBL_MIN) success = false;
            break;
            
        case 2: // Surface elements
            detG = ( covariantMetricTensor[0]*covariantMetricTensor[4] - covariantMetricTensor[1]*covariantMetricTensor[3] );
            if (detG <= DBL_MIN) success = false;
            break;
            
        case 3: // Volume elements
            detG = covariantMetricTensor[0] * ( covariantMetricTensor[4]*covariantMetricTensor[7] - covariantMetricTensor[5]*covariantMetricTensor[7] )
            + covariantMetricTensor[1] * ( covariantMetricTensor[4]*covariantMetricTensor[6] - covariantMetricTensor[3]*covariantMetricTensor[8] )
            + covariantMetricTensor[2] * ( covariantMetricTensor[3]*covariantMetricTensor[7] - covariantMetricTensor[4]*covariantMetricTensor[6] );
            if (detG <= DBL_MIN) success = false;
            break;
            
        default:
            errorfunct("FEMNumericIntegration:detJForElement", "Dimension not supported!!");
            break;
    }
    
    if (success == false) {
        printf("detJForElement: Degenerate %d D element: %d.\n", dim, element->ElementIndex);
        printf("detJForElement: detJ: %f\n", detG);
        for (int i=0; i<cdim; i++) {
            printf("detJForElement: Dir: %d, Coord: %f %f %f\n", cdim, nodes->x[i], nodes->y[i], nodes->z[i]);
        }
        if (cdim < dim) {
            printf("detJForElement: element dimension larger than mesh dimension: %d vs %d\n", dim, cdim);
        }
        errorfunct("detJForElement", "Programm terminating now...");
    }
    return detG;
}

/********************************************************************************************************
 
    Compute contravariant metric tensor (=J^TJ)^-1 of element coordinate system
 
*********************************************************************************************************/
inline bool contravariantMetric(double *elementMetric, double *dx, Element_t *element, Nodes_t *nodes, int meshDimension, double u, double v, double w, double *dLBasisdx) {
    
    int i, j, dim, cdim;
    double detG, GI[9], covariantMetricTensor[9];
    
    dim = element->Type.dimension;
    cdim = meshDimension;
    
    memset(covariantMetricTensor, 0.0, sizeof(covariantMetricTensor) );
    detG = detJ(covariantMetricTensor, dx, element, nodes, meshDimension, u, v, w, dLBasisdx);
    
    // Convert the metric to contravariant base
    switch (dim) {
        case 1:
            elementMetric[0] = 1.0 / detG;
            break;
        case 2:
            elementMetric[0] = covariantMetricTensor[4] / detG;
            elementMetric[1] = -covariantMetricTensor[1] / detG;
            elementMetric[3] = -covariantMetricTensor[3] / detG;
            elementMetric[4] = covariantMetricTensor[0] / detG;
            break;
        case 3:
            invertMatrix3x3(GI, covariantMetricTensor, detG);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    elementMetric[3*i+j] = GI[3*i+j];
                }
            }
            break;
            
        default:
            errorfunct("FEMNumericIntegration:setMetricForElement", "Dimension not supported!!");
            break;
    }
    
    return true;
}

inline bool localtoGlobalMap(double *ltoGMap, Element_t *element, Nodes_t* nodes, int meshDimension, double u, double v, double w, double *dLBasisdx) {
    
    int i, j, k, dim, cdim;
    double dx[9], elementMetric[9], s;
    bool success;
    
    memset(dx, 0.0, sizeof(dx));
    memset(elementMetric, 0.0, sizeof(elementMetric));
    
    dim = element->Type.dimension;
    cdim = meshDimension;
    
    success = contravariantMetric(elementMetric, dx, element, nodes, meshDimension, u, v, w, dLBasisdx);
    
    for (i=0; i<cdim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                s = s + dx[3*i+k] * elementMetric[k*3+j];
            }
            ltoGMap[3*i+j] = s;
        }
    }
    
    return true;
}

/******************************************************************************************************
    Arguments:
        Element_t* element -> structure describing the element
        Nodes_t* nodes     -> nodal coordinates
        FEMMesh *mesh      -> mesh
        double u,v,w       -> point at which to evaluate
        double *f          -> nodal values of the quantity
 
    Output: 3x3 matrix (values) of partial derivatives
 
******************************************************************************************************/
inline void globalSecondDerivatives(double *elementMetric, Element_t *element, Nodes_t *nodes, int meshDimension, double u, double v, double w, double *f, double *dLBasisdx, double *values) {
    
    int i, j, k, l, n, dim, cdim;
    double C1[3][3][3], C2[3][3][3], ddx[3][3][3];
    double df[3];
    double cddf[3][3], ddf[3][3], dxx[3][3], bf[4], bff[9];
    double accum1, accum2, accum3, accum4;
    double s;
    
    // Actually not quite correct...
    if (element->Type.BasisFunctionDegree <= 1 ) return;
    
    n = element->Type.NumberOfNodes;
    dim = element->Type.dimension;
    cdim = meshDimension;
    
    // Partial derivatives of the basis functions are given, just
    // sum for the first partial derivatives.
    memset( df, 0.0, sizeof(df) );
    memset( *dxx, 0.0, (3*3)*sizeof(double) );
    
    switch (cdim) {
        case 1:
            for (i=0; i<dim; i++) {
                accum1 = 0.0;
                accum2 = 0.0;
                for (j=0; j<n; j++) {
                    accum1 = accum1 + nodes->x[j] * dLBasisdx[3*j+i];
                    accum2 = accum2 + f[j] * dLBasisdx[3*j+i];
                }
                dxx[0][i] = accum1;
                df[i] = accum2;
            }
            break;
            
        case 2:
            for (i=0; i<dim; i++) {
                accum1 = 0.0;
                accum2 = 0.0;
                accum3 = 0.0;
                for (j=0; j<n; j++) {
                    accum1 = accum1 + nodes->x[j] * dLBasisdx[3*j+i];
                    accum2 = accum2 + nodes->y[j] * dLBasisdx[3*j+i];
                    accum3 = accum3 + f[j] * dLBasisdx[3*j+i];
                }
                dxx[0][i] = accum1;
                dxx[1][i] = accum2;
                df[i] = accum3;
            }
            break;
            
        case 3:
            for (i=0; i<dim; i++) {
                accum1 = 0.0;
                accum2 = 0.0;
                accum3 = 0.0;
                accum4 = 0.0;
                for (j=0; j<n; j++) {
                    accum1 = accum1 + nodes->x[j] * dLBasisdx[3*j+i];
                    accum2 = accum2 + nodes->y[j] * dLBasisdx[3*j+i];
                    accum3 = accum3 + nodes->z[j] * dLBasisdx[3*j+i];
                    accum4 = accum4 + f[j] * dLBasisdx[3*j+i];
                }
                dxx[0][i] = accum1;
                dxx[1][i] = accum2;
                dxx[2][i] = accum3;
                df[i] = accum4;
            }
            break;
            
        default:
            errorfunct("FEMNumericIntegration:globalSecondDerivativesForElement", "Coordinate dimension not supported!!");
            break;
    }
    
    // Get second partial derivatives with respect to local coordinates
    switch (dim) {
        case 1:
            // Line elements
            ddx[0][0][0] = SecondDerivatives1D(element, nodes->x, u);
            ddx[1][0][0] = SecondDerivatives1D(element, nodes->y, u);
            ddx[2][0][0] = SecondDerivatives1D(element, nodes->z, u);
            break;
            
        case 2:
            // Surface elements
            memset( bf, 0.0, sizeof(bf) );
            SecondDerivatives2D(bf, element, nodes->x, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddx[0][i][j] = bf[2*i+j];
                }
            }
            
            memset( bf, 0.0, sizeof(bf) );
            SecondDerivatives2D(bf, element, nodes->y, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddx[1][i][j] = bf[2*i+j];
                }
            }
            
            memset( bf, 0.0, sizeof(bf) );
            SecondDerivatives2D(bf, element, nodes->z, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddx[2][i][j] = bf[2*i+j];
                }
            }
            break;
            
        case 3:
            // Volume elements
            memset( bff, 0.0, sizeof(bff) );
            SecondDerivatives3D(bff, element, nodes->x, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddx[0][i][j] = bff[3*i+j];
                }
            }
            
            memset( bff, 0.0, sizeof(bff) );
            SecondDerivatives3D(bff, element, nodes->y, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddx[1][i][j] = bff[3*i+j];
                }
            }
            
            memset( bff, 0.0, sizeof(bff) );
            SecondDerivatives3D(bff, element, nodes->z, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddx[2][i][j] = bff[3*i+j];
                }
            }
            break;
            
        default:
            errorfunct("FEMNumericIntegration:globalSecondDerivativesForElement", "Element dimension not supported");
            break;
    }
    
    // Christoffel symbols of the second kind of the element coordinate system
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            for (k=0; k<dim; k++) {
                s = 0.0;
                for (l=0; l<cdim; l++) {
                    s = s + ddx[l][i][j]*dxx[l][k];
                }
                C2[i][j][k] = s;
            }
        }
    }
    
    // Christoffel symbols of the first kind
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            for (k=0; k<dim; k++) {
                s = 0.0;
                for (l=0; l<dim; l++) {
                    s = s + elementMetric[3*k+l]*C2[i][j][l];
                }
                C1[i][j][k] = s;
            }
        }
    }
    
    // First add ordinary partials (change of the quantity with coordinates)...
    switch (dim) {
        case 1:
            ddf[0][0] = SecondDerivatives1D(element, f, u);
            break;
            
        case 2:
            memset( bf, 0.0, sizeof(bf) );
            SecondDerivatives2D(bf, element, f, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddf[i][j] = bf[2*i+j];
                }
            }
            break;
            
        case 3:
            memset( bff, 0.0, sizeof(bff) );
            SecondDerivatives3D(bff, element, f, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddf[i][j] = bff[3*i+j];
                }
            }
            break;
            
        default:
            errorfunct("FEMNumericIntegration:globalSecondDerivativesForElement", "Element dimension not supported");
            break;
    }
    
    // ... Then add change of coordinates
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                s = s - C1[i][j][k]*df[k];
            }
            ddf[i][j] = ddf[i][j] + s;
        }
    }
    
    // Convert to contravariant base
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                for (l=0; l<dim; l++) {
                    s = s + elementMetric[3*i+k]*elementMetric[3*j+l]*ddf[k][l];
                }
            }
            cddf[i][j] = s;
        }
    }
    
    // And finally transform to global coordinates
    for (i=0; i<cdim; i++) {
        for (j=0; j<cdim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                for (l=0; l<dim; l++) {
                    s = s + dxx[i][k]*dxx[j][l]*cddf[k][l];
                }
            }
            values[3*i+j] = s;
        }
    }
}
