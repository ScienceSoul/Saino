//
//  FEMSolver.m
//  Saino
//
//  Created by Seddik hakime on 17/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMNumericIntegration.h"

#include "Utils.h"
#include "memory.h"


@implementation FEMNumericIntegration

- (id)init
{
    self = [super init];
    if (self) {
        // Initialization code here.
    }
    
    return self;
}

- (void)dealloc
{
    
}

-(BOOL)setBasis:(Element_t *)element :(Nodes_t *)nodes :(double)u :(double)v :(double)w :(BOOL)Bubbles :(int)el {
    
    int i, n, dim, cdim;
    
    n = [self Elements_Type_NumberOfNodes:el];
    dim = [self Elements_Type_dimension:el];
    cdim = [self dimension];
    
    if ([self Elements_Type_ElementCode:el] == 101) {
        basis[0] = 1.0;
        return YES;
    }
    
    for (i=0; i<n; i++) {
        basis[i] = 0.0;
    }
    
    NodalBasisFunctions(n, basis, element, u, v, w);
    
    return YES;
}

-(BOOL)setdBasisdx:(Element_t *)element :(Nodes_t *)nodes :(double)u :(double)v :(double)w :(BOOL)Bubbles :(int)el {
/*******************************************************************************/
/* Global first derivatives                                                    */
/*******************************************************************************/
    
    int i, j, k, n, q, dim, cdim;
    double **dLBasisdx;
    
    n = [self Elements_Type_NumberOfNodes:el];
    dim = [self Elements_Type_dimension:el];
    cdim = [self dimension];
    
    if ([self Elements_Type_ElementCode:el] == 101) {
        for (i=0; i<2; i++) {
            dBasisdx[0][i] = 0.0;
        }
    }
    
    dLBasisdx = doublematrix(0, n-1, 0, 2);
    
    for (i=0; i<n; i++) {
        for (j=0; j<dim; j++) {
            dLBasisdx[i][j] = 0.0;
        }
    }
    
    NodalFirstDerivatives(n, dLBasisdx, element, u, v, w);
    
    q = n;
    
    for (i=0; i<n; i++) {
        for (j=0; j<2; j++) {
            dBasisdx[i][j] = 0.0;
        }
    }
    
    [self setLtoGMap:element :nodes :u :v :w :el];
    
    for (i=0; i<q; i++) {
        for (j=0; j<cdim; j++) {
            for (k=0; k<dim; k++) {
                dBasisdx[i][j] = dBasisdx[i][j] * dLBasisdx[i][k]*[self LtoGMap:j :k];
            }
        }
    }
    
    free_dmatrix(dLBasisdx, 0, n-1, 0, 2);
    
    return YES;
    
}

-(double)dBasisdx:(int)i :(int)j {
    
    return dBasisdx[i][j];
}

-(BOOL)setddBasisddx:(Element_t *)element :(Nodes_t *)nodes :(double)u :(double)v :(double)w :(BOOL)Bubbles :(int)el {
    
    int i, j, k, n, q, dim, cdim;
    double *NodalBasis, **dLBasisdx;
    double **Values;
    
    n = [self Elements_Type_NumberOfNodes:el];
    dim = [self Elements_Type_dimension:el];
    cdim = [self dimension];
    
    NodalBasis = doublevec(0, n-1);
    dLBasisdx = doublematrix(0, n-1, 0, 2);
    Values = doublematrix(0, 2, 0, 2);
    
    for (i=0; i<n; i++) {
        NodalBasis[i] = 0.0;
    }
    
    for (i=0; i<n; i++) {
        for (j=0; j<dim; j++) {
            dLBasisdx[i][j] = 0.0;
        }
    }

    for (i=0; i<n; i++) {
        for (j=0; j<2; j++) {
            for (k=0; k<2; k++) {
                ddBasisddx[i][j][k] = 0.0;
            }
        }
    }
    
    NodalFirstDerivatives(n, dLBasisdx, element, u, v, w);
    
    [self setElementMetric:element :nodes :u :v :w :el];
    
    for (q=0; q<n; q++) {
        NodalBasis[q] = 1.0;
        [self GlobalSecondDerivatives:element :nodes :u :v :w :NodalBasis :dLBasisdx :Values :el];
        for (i=0; i<2; i++) {
            for (j=0; j<2; j++) {
                ddBasisddx[q][i][j] = Values[i][j];
            }
        }
        NodalBasis[q] = 0.0;
    }
    
    free_dvector(NodalBasis, 0, n-1);
    free_dmatrix(dLBasisdx, 0, n-1, 0, 2);
    free_dmatrix(Values, 0, 2, 0, 2);
    
    return YES;
    
}

-(void)setdx:(int)NDOFs :(Nodes_t *)nodes :(double **)dLBasisdx :(int)el{
/*******************************************************************************/
/* Partial derivatives of global coordinates with respect to local coordinates */
    
/*  Arguments:
     NDOFs              -> Number of DOFs
     
     Nodes_t *nodes     -> element nodal cooridnates 
     
     double **dLBasisdx -> Derivatives of element basis function with respect to local coordinates
     
     int el             -> Element number                                       */

/*******************************************************************************/
    int i, j, n, dim;
    double accum1, accum2, accum3;
    
    dim = [self Elements_Type_dimension:el];
    n = min((sizeof(nodes) / (sizeof(Nodes_t))), NDOFs);
    
    for (i=0; i<dim; i++) {
        accum1 = 0.0;
        accum2 = 0.0;
        accum3 = 0.0;
        for (j=0; j<n; j++) {
            accum1 = accum1 + (nodes[n].x * dLBasisdx[n][i]);
            accum2 = accum2 + (nodes[n].y * dLBasisdx[n][i]);
            accum3 = accum3 + (nodes[n].z * dLBasisdx[n][i]);
        }
        dx[0][i] = accum1;
        dx[1][i] = accum2;
        dx[2][i] = accum3;
    }
}

-(double)dx:(int)i :(int)j {
    
    return dx[i][j];
}

-(void)setCovariantMetrixTensor:(int)NDOFs :(Nodes_t *)nodes :(double **)dLBasisdx :(int)el{
/*******************************************************************************/
/*  Compute the covariant metric tensor of the element coordinate system       */

/*  Arguments:
    NDOFs              -> Number of DOFs
    
    Nodes_t *nodes     -> element nodal cooridnates 
    
    double **dLBasisdx -> Derivatives of element basis function with respect to local coordinates
 
    int el             -> Element number                                       */

/*******************************************************************************/
    
    int i, j, k, dim, cdim;
    double s;
    
    dim = [self Elements_Type_dimension:el];
    cdim = [self dimension];
    
    [self setdx:NDOFs :nodes :dLBasisdx :el];
    
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<cdim; k++) {
                s = s + ( [self dx:k :i] * [self dx:k :j]);
            }
            covariantMetrixTensor[i][j] = s;
        }
    }
}

-(double)CovariantMetrixTensor:(int)i :(int)j {
    
    return covariantMetrixTensor[i][j];
}

-(void)InvertMatrix3x3:(double **)GI :(double)detG {
    
    double s;
    
    s = 1.0 / detG;
    
    GI[0][0] = s * ( [self CovariantMetrixTensor:1 :1]*[self CovariantMetrixTensor:2 :2] - [self CovariantMetrixTensor:2 :1]*[self CovariantMetrixTensor:1 :2] );
    
    GI[1][0] = -s * ( [self CovariantMetrixTensor:1 :0]*[self CovariantMetrixTensor:2 :2] - [self CovariantMetrixTensor:2 :0]*[self CovariantMetrixTensor:1 :2] );
    
    GI[2][0] = s * ( [self CovariantMetrixTensor:1 :0]*[self CovariantMetrixTensor:2 :1] - [self CovariantMetrixTensor:2 :0]*[self CovariantMetrixTensor:1 :1] );
    
    GI[0][1] = -s * ( [self CovariantMetrixTensor:0 :1]*[self CovariantMetrixTensor:2 :2] - [self CovariantMetrixTensor:2 :1]*[self CovariantMetrixTensor:0 :2] );
    
    GI[1][1] = s * ( [self CovariantMetrixTensor:0 :0]*[self CovariantMetrixTensor:2 :2] - [self CovariantMetrixTensor:2 :0]*[self CovariantMetrixTensor:0 :2] );
    
    GI[2][1] = -s * ( [self CovariantMetrixTensor:0 :0]*[self CovariantMetrixTensor:2 :1] - [self CovariantMetrixTensor:2 :0]*[self CovariantMetrixTensor:0 :1] );
    
    GI[0][2] = s * ( [self CovariantMetrixTensor:0 :1]*[self CovariantMetrixTensor:1 :2] - [self CovariantMetrixTensor:1 :1]*[self CovariantMetrixTensor:0 :2] );
    
    GI[1][2] = -s * ( [self CovariantMetrixTensor:0 :0]*[self CovariantMetrixTensor:1 :2] - [self CovariantMetrixTensor:1 :0]*[self CovariantMetrixTensor:0 :2] );
    
    GI[2][2] = s * ( [self CovariantMetrixTensor:0 :0]*[self CovariantMetrixTensor:1 :1] - [self CovariantMetrixTensor:1 :0]*[self CovariantMetrixTensor:0 :1] );
}

-(double)detG:(Element_t *)element :(Nodes_t *)nodes :(double)u :(double)v :(double)w :(int)el{
/*******************************************************************************/
/* Compute determinant of covariant metric tensor det(J^TJ) and square root of determinant of covariant metric tensor (=sqrt(det(J^TJ))).                                                */
/*  Arguments:
     Element_t *element  -> element structure
     
     Nodes_t *nodes      -> element nodal cooridnates 
     
     double u, v, w      -> Points at which evaluate the value
 
     int el              -> Element number
     
     Function return value:
     If function failure, element is degenerate.                               */

/*******************************************************************************/
    
    int i, j, n, dim, cdim;
    double detG, **dLBasisdx;
    
    n = [self Elements_Type_NumberOfNodes:el];
    dim = [self Elements_Type_dimension:el];
    cdim = [self dimension];
    
    dLBasisdx = doublematrix(0, n-1, 0, 2);
    
    for (i=0; i<n; i++) {
        for (j=0; j<dim; j++) {
            dLBasisdx[i][j] = 0.0;
        }
    }
    
    NodalFirstDerivatives(n, dLBasisdx, element, u, v, w);
    
    [self setCovariantMetrixTensor:n :nodes :dLBasisdx :el];
    
    switch (dim) {
        // Line elements
        case 1:
            detG = [self CovariantMetrixTensor:1 :1];
            if (detG <= 0.0) {
                errorfunct("ElementMetric", "Degenerate 1D element",el);
            }
            break;
        
        // Surface elements
        case 2:
            detG = ( [self CovariantMetrixTensor:0 :0]*[self CovariantMetrixTensor:1 :1] - [self CovariantMetrixTensor:0 :1]*[self CovariantMetrixTensor:1 :0] );
            if (detG <= 0.0) {
                if (cdim < dim) {
                    warnfunct("ElementMetric", "2D element in 1D coordinate systrem?");
                }
                errorfunct("ElementMetric", "Degenerate 2D element",el);
            }
            break;

        // Volume elements    
        case 3:
            detG = [self CovariantMetrixTensor:0 :0] * ( [self CovariantMetrixTensor:1 :1]*[self CovariantMetrixTensor:2 :2] - [self CovariantMetrixTensor:1 :2]*[self CovariantMetrixTensor:2 :1] ) + [self CovariantMetrixTensor:0 :1] * ( [self CovariantMetrixTensor:1 :2]*[self CovariantMetrixTensor:2 :0] - [self CovariantMetrixTensor:1 :0]*[self CovariantMetrixTensor:2 :2] ) + [self CovariantMetrixTensor:0 :2] * ( [self CovariantMetrixTensor:1 :0]*[self CovariantMetrixTensor:2 :1] - [self CovariantMetrixTensor:1 :1]*[self CovariantMetrixTensor:2 :0] );
            if (detG <= 0.0) {
                if (cdim < dim) {
                    warnfunct("ElementMetric", "2D/3D element in 1D/2D coordinate systrem?");
                }
                errorfunct("ElementMetric", "Degenerate 3D element",el);
            }
            break;
            
        default:
            errorfunct("Method setdetJ", "Dimension not supported!!");
            break;
    }
    free_dmatrix(dLBasisdx, 0, n-1, 0, 2);
    
    return detG;
}

-(BOOL)setdetJ:(Element_t *)element :(Nodes_t *)nodes :(double)u :(double)v :(double)w :(int)el {
    
    detJ = [self detG:element :nodes :u :v :w :el];
    
    detJ = sqrt(detJ);
    
    return YES;
}

-(double)detJ {
    
    return detJ;
    
}

-(BOOL)setElementMetric:(Element_t *)element :(Nodes_t *)nodes :(double)u :(double)v :(double)w :(int)el{
/*******************************************************************************/
/* Compute contravariant metric tensor (=J^TJ)^-1 of element coordinate system */
/*******************************************************************************/
    
    int i, j, dim, cdim;
    double detG, **GI;
    
    dim = [self Elements_Type_dimension:el];
    cdim = [self dimension];
    
    detG = [self detG:element :nodes :u :v :w :el];
    
    // Convert the metric to contravariant base
    switch (dim) {
        case 1:
            elementMetric[0][0] = 1.0 / detG;
            break;
        case 2:
            elementMetric[0][0] = [self CovariantMetrixTensor:1 :1] / detG;
            elementMetric[0][1] = -[self CovariantMetrixTensor:0 :1] / detG;
            elementMetric[1][0] = -[self CovariantMetrixTensor:1 :0] / detG;
            elementMetric[1][1] = [self CovariantMetrixTensor:0 :0] / detG;
            break;
        case 3:
            GI = doublematrix(0, 2, 0, 2);
            [self InvertMatrix3x3:GI :detG];
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    elementMetric[i][j] = GI[i][j];
                }
            }
            free_dmatrix(GI, 0, 2, 0, 2);
            break;
            
        default:
            errorfunct("Method setElementMetric", "Dimension not supported!!");
            break;
    }
    
    return YES;

}

-(double)elementMetric:(int)i :(int)j {
    
    return elementMetric[i][j];
}

-(BOOL)setLtoGMap:(Element_t *)element :(Nodes_t *)nodes :(double)u :(double)v :(double)w :(int)el {
    
    int i, j, k, dim, cdim;
    double s;
    
    dim = [self Elements_Type_dimension:el];
    cdim = [self dimension];
    
    [self setElementMetric:element :nodes :u :v :w :el];
    
    for (i=0; i<cdim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                s = s + [self dx:i :k] * [self elementMetric:k :j];
            }
            LtoGMap[i][j] = s;
        }
    }
    
    return YES;
    
}

-(double)LtoGMap:(int)i :(int)j {
    
    return LtoGMap[i][j];
}

-(void)GlobalSecondDerivatives:(Element_t *)element :(Nodes_t *)nodes :(double)u :(double)v :(double)w :(double *)f :(double **)dLBasisdx :(double **)values :(int)el {
    
    int i, j, k, l, n, dim, cdim;
    
    double ***C1, ***C2, ***ddx;
    double *df;
    double **cddf, **ddf, **dxx;
    double **bf, *vecnodes;
    double accum1, accum2, accum3, accum4;
    double s;

    // Actually not quite correct... 
    if ([self Elements_Type_BasisFunctionDegree:el] <=1 ) return;
    
    n = [self Elements_Type_NumberOfNodes:el];
    dim = [self Elements_Type_dimension:el];
    cdim = [self dimension];
    
    C1 = d3tensor(0, 2, 0, 2, 0, 2);
    C2 = d3tensor(0, 2, 0, 2, 0, 2);
    ddx = d3tensor(0, 2, 0, 2, 0, 2);
    
    df = doublevec(0, 2);
    
    cddf = doublematrix(0, 2, 0, 2);
    ddf = doublematrix(0, 2, 0, 2);
    dxx = doublematrix(0, 2, 0, 2);
    
    vecnodes = doublevec(0, n-1);

    // Partial derivatives of the basis functions are given, just
    // sum for the first partial derivatives.
    for (i=0; i<2; i++) {
        df[i] = 0.0;
        for (j=0; j<2; j++) {
            dxx[i][j] = 0.0;
        }
    }
    
    switch (cdim) {
        case 1:
            for (i=0; i<dim; i++) {
                accum1 = 0.0;
                accum2 = 0.0;
                for (j=0; j<n; j++) {
                    accum1 = accum1 + nodes[n].x * dLBasisdx[n][i];
                    accum2 = accum2 + f[n] * dLBasisdx[n][i];
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
                    accum1 = accum1 + nodes[n].x * dLBasisdx[n][i];
                    accum2 = accum2 + nodes[n].y * dLBasisdx[n][i];
                    accum3 = accum3 + f[n] * dLBasisdx[n][i];
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
                    accum1 = accum1 + nodes[n].x * dLBasisdx[n][i];
                    accum2 = accum2 + nodes[n].y * dLBasisdx[n][i];
                    accum3 = accum3 + nodes[n].z * dLBasisdx[n][i];
                    accum4 = accum4 + f[n] * dLBasisdx[n][i];
                }
                dxx[0][i] = accum1;
                dxx[1][i] = accum2;
                dxx[2][i] = accum3;
                df[i] = accum4;
            }
            break;

        default:
            errorfunct("Method GlobalSecondDerivatives", "Coordinate dimension not supported!!");
            break;
    }
    
    // Get second partial derivatives with respect to local coordinates
    switch (dim) {
        case 1:
            // Line elements
            for (i=0; i<n; i++) {
                vecnodes[n] = nodes[n].x;
            }
            ddx[0][0][0] = SecondDerivatives1D(element, vecnodes, u);
            for (i=0; i<n; i++) {
                vecnodes[n] = nodes[n].y;
            }
            ddx[1][0][0] = SecondDerivatives1D(element, vecnodes, u);
            for (i=0; i<n; i++) {
                vecnodes[n] = nodes[n].z;
            }
            ddx[2][0][0] = SecondDerivatives1D(element, vecnodes, u);
            break;
            
        case 2:
            // Surface elements
            bf = doublematrix(0, 1, 0, 1);
            for (i=0; i<n; i++) {
                vecnodes[n] = nodes[n].x;
            }
            SecondDerivatives2D(bf, element, vecnodes, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddx[0][i][j] = bf[i][j];
                }
            }
            
            for (i=0; i<n; i++) {
                vecnodes[n] = nodes[n].y;
            }
            SecondDerivatives2D(bf, element, vecnodes, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddx[1][i][j] = bf[i][j];
                }
            }
            
            for (i=0; i<n; i++) {
                vecnodes[n] = nodes[n].z;
            }
            SecondDerivatives2D(bf, element, vecnodes, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddx[2][i][j] = bf[i][j];
                }
            }
            
            free_dmatrix(bf, 0, 1, 0, 1);
            break;
            
        case 3:
            // Volume elements
            bf = doublematrix(0, 2, 0, 2);
            for (i=0; i<n; i++) {
                vecnodes[n] = nodes[n].x;
            }
            SecondDerivatives3D(bf, element, vecnodes, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddx[0][i][j] = bf[i][j];
                }
            }
            
            for (i=0; i<n; i++) {
                vecnodes[n] = nodes[n].y;
            }
            SecondDerivatives3D(bf, element, vecnodes, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddx[1][i][j] = bf[i][j];
                }
            }
            
            for (i=0; i<n; i++) {
                vecnodes[n] = nodes[n].z;
            }
            SecondDerivatives3D(bf, element, vecnodes, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddx[2][i][j] = bf[i][j];
                }
            }

            free_dmatrix(bf, 0, 2, 0, 2);
            break;
            
        default:
            errorfunct("Method GlobalSecondDerivatives", "Element dimension not supported");
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
                    s = s + [self elementMetric:k :l]*C2[i][j][l];
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
            bf = doublematrix(0, 1, 0, 1);
            SecondDerivatives2D(bf, element, f, u, v);
            for (i=0; i<2; i++) {
                for (j=0; j<2; j++) {
                    ddf[i][j] = bf[1][2];
                }
            }
            
            free_dmatrix(bf, 0, 1, 0, 1);
            break;
            
        case 3:
            bf = doublematrix(0, 2, 0, 2);
            SecondDerivatives3D(bf, element, f, u, v, w);
            for (i=0; i<3; i++) {
                for (j=0; j<3; j++) {
                    ddf[i][j] = bf[i][j];
                }
            }
            
            free_dmatrix(bf, 0, 2, 0, 2);
            break;
            
        default:
            errorfunct("Method GlobalSecondDerivatives", "Element dimension not supported");
            break;
    }
    
    // ... Then add change of coordinates
    for (i=0; i<dim; i++) {
        for (j=0; j<dim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                s = s + C1[i][j][k]*df[k];
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
                    s = s + [self elementMetric:i :k]*[self elementMetric:j :l]*ddf[k][l];
                }
            }
            cddf[i][j] = s;
        }
    }
    
    // And finally transform to global coordinates
    for (i=0; i<cdim; i++) {
        for (j=0; j<cdim; j++) {
            values[i][j] = 0.0;
        }
    }
    
    for (i=0; i<cdim; i++) {
        for (j=0; j<cdim; j++) {
            s = 0.0;
            for (k=0; k<dim; k++) {
                for (l=0; l<dim; l++) {
                    s = s + dxx[i][k]*dxx[j][l]*cddf[k][l];
                }
            }
            values[i][j] = s;
        }
    }
    
    free_d3tensor(C2, 0, 2, 0, 2, 0, 2);
    free_d3tensor(C1, 0, 2, 0, 2, 0, 2);
    free_d3tensor(ddx, 0, 2, 0, 2, 0, 2);
    
    free_dvector(df, 0, 2);
    
    free_dmatrix(cddf, 0, 2, 0, 2);
    free_dmatrix(ddf, 0, 2, 0, 2);
    free_dmatrix(dxx, 0, 2, 0, 2);
    
    free_dvector(vecnodes, 0, n-1);
    
}


@end
