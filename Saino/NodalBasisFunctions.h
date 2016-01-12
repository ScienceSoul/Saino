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
#include <float.h>
#include <math.h>

#include "Constructors.h"

/****************************************************************************
 
    NodalBasisFunctions1D:
 
    Description:
 
    Arguments:
        Element_t *element -> element structure
        double u -> Point at which evaluate the value
 
    Function return value:
        double *y -> value of the quantity y = x(u)
 
****************************************************************************/
inline void NodalBasisFunctions1D(double * __nonnull y, Element_t * __nonnull element, double u) {
    
    int i, n;
    int *p = NULL;
    
    double s;
    double *coeff = NULL;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        s = 0.0;
        for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            s = s + coeff[i] * pow(u, (double)p[i]);
        }
        y[n] = s;
    }
}

/****************************************************************************
 
    NodalBasisFunction2D:
 
    Description:
 
    Arguments:
        Element_t *element -> element structure
        double u, v -> Points at which evaluate the value
 
    Function return value:
        double *y -> value of the quantity y = x(u)
 
****************************************************************************/
inline void NodalBasisFunctions2D(double * __nonnull y, Element_t * __nonnull element, double u, double v) {
    
    int i, n;
    int *p = NULL, *q = NULL;
    
    double s, s1, s2, s3;
    double *coeff = NULL;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        q = element->Type.BasisFunctions[n].q;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        // The following is the unrolled version of this loop:
        //
        //      s = 0.0;
        //      for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
        //          s = s + coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]);
        //      }
        //      y[n] = s;
        //
        // We assume a minimum of three nodes for a 2D element so we use a
        // three-way unrolling
        
        s = s1 = s2 = s3 = 0.0;
        for(i=0;i<=element->Type.BasisFunctions[n].n-3;i+=3) {
            s1 = s1 + coeff[i+0] * pow(u, (double)p[i+0]) * pow(v, (double)q[i+0]);
            s2 = s2 + coeff[i+1] * pow(u, (double)p[i+1]) * pow(v, (double)q[i+1]);
            s3 = s3 + coeff[i+2] * pow(u, (double)p[i+2]) * pow(v, (double)q[i+2]);
        }
        // Remainder loop
        for( ;i<element->Type.BasisFunctions[n].n;i++) {
            s = s + coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]);
        }
        y[n] = s + s1 + s2 + s3;
    }
}

/**********************************************************************************
 
    InterpolateInElement3D:
 
    Description: Given element structure return value of a quantity x given at
    element nodes at local coordinate point (u, v, w) inside the element.
    Element basis functions are used to compute the value. This is for 3D
    elements, and should not be called directley by the user but through the
    appropriate method in the kernel class.
 
    Arguments:
        Element_t *element -> element structure
        double *x -> Nodal values of the quantity whose value we want to know
        double u, v, w -> Points at which evaluate the value
 
    Function return value:
        double y -> value of the quantity y = x(u,v,w)
 
**********************************************************************************/
inline double InterpolateInElement3D(Element_t * __nonnull element, double * __nonnull x, double u, double v, double w) {
    
    int i, l, n;
    int *p = NULL, *q = NULL, *r = NULL;
    
    double y, s, s1, s2, s3, s4;
    double *coeff = NULL;
    
    l = element->Type.BasisFunctionDegree;
    
    if(element->Type.ElementCode == 605) {
        
        s = 0.0;
        if( w == 1) w = 1.0 - 1.0 - 12.0;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * ( (1.0-u) * (1.0-v) - w + u*v*w * s ) / 4.0;
        y = y + x[1] * ( (1.0+u) * (1.0-v) - w - u*v*w * s ) / 4.0;
        y = y + x[2] * ( (1.0+u) * (1.0+v) - w + u*v*w * s ) / 4.0;
        y = y + x[3] * ( (1.0-u) * (1.0+v) - w - u*v*w * s ) / 4.0;
        y = y + x[4] * w;
        return y;
    } else if (element->Type.ElementCode == 613) {
        
        if( w == 1) w = 1.0 - 1.0 - 12.0;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * (-u-v-1.0) * ( (1.0-u) * (1.0-v) - w + u*v*w * s ) / 4.0;
        y = y + x[1] * (u-v-1.0) * ( (1.0+u) * (1.0-v) - w - u*v*w * s ) / 4.0;
        y = y + x[2] * (u+v-1.0) * ( (1.0+u) * (1.0+v) - w + u*v*w * s ) / 4.0;
        y = y + x[3] * (-u+v-1.0) * ( (1.0-u) * (1.0+v) - w - u*v*w * s ) / 4.0;
        y = y + x[4] *w*(2.0*w-1);
        y = y + x[5] * (1.0+u-w)*(1.0-u-w)*(1.0-v-w) * s /2.0;
        y = y + x[6] * (1.0+v-w)*(1.0-v-w)*(1.0+u-w) * s /2.0;
        y = y + x[7] * (1.0+u-w)*(1.0-u-w)*(1.0+v-w) * s /2.0;
        y = y + x[8] * (1.0+v-w)*(1.0-v-w)*(1.0-u-w) * s /2.0;
        y = y + x[9] * w * (1.0-u-w) * (1.0-v-w) * s;
        y = y + x[10] * w * (1.0+u-w) * (1.0-v-w) * s;
        y = y + x[11] * w * (1.0+u-w) * (1.0+v-w) * s;
        y = y + x[12] * w * (1.0-u-w) * (1.0+v-w) * s;
        return y;
    }
    
    y = 0.0;
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        if(x[n] != 0.0) {
            
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            r = element->Type.BasisFunctions[n].r;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            // The following is the unrolled version of this loop:
            //
            //      s = 0.0;
            //      for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            //          s = s + coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
            //      }
            //      y = y + s*x[n];
            //
            // We assume a minimum of four nodes for a 3D element so we use a four-way unrolling
            
            s = s1 = s2 = s3 = s4 = 0.0;
            for(i=0;i<=element->Type.BasisFunctions[n].n-4;i+=4) {
                s1 = s1 + coeff[i+0] * pow(u, (double)p[i+0]) * pow(v, (double)q[i+0]) * pow(w, (double)r[i+0]);
                s2 = s2 + coeff[i+1] * pow(u, (double)p[i+1]) * pow(v, (double)q[i+1]) * pow(w, (double)r[i+1]);
                s3 = s3 + coeff[i+2] * pow(u, (double)p[i+2]) * pow(v, (double)q[i+2]) * pow(w, (double)r[i+2]);
                s4 = s4 + coeff[i+3] * pow(u, (double)p[i+3]) * pow(v, (double)q[i+3]) * pow(w, (double)r[i+3]);
            }
            y = y + (s1 + s2 + s3 + s4)*x[n];
            // Remainder loop
            for( ;i<element->Type.BasisFunctions[n].n;i++) {
                s = s + coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
            }
            y = y + s*x[n];
        }
    }
    
    return y;
}

/****************************************************************************
 
    NodalBasisFunctions3D:
 
    Description:
 
    Arguments:
        Element_t *element -> element structure
        double u, v, w -> Points at which evaluate the value
 
    Function return value:
        double *y -> value of the quantity y = x(u)
 
****************************************************************************/
inline void NodalBasisFunctions3D(double * __nonnull y, Element_t * __nonnull element, double u, double v, double w) {
    
    int i, n;
    int *p = NULL, *q = NULL, *r = NULL;
    
    double s, s1, s2, s3, s4;
    double *coeff = NULL;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        q = element->Type.BasisFunctions[n].q;
        r = element->Type.BasisFunctions[n].r;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        // The following is the unrolled version of this loop:
        //
        //      s = 0.0;
        //      for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
        //          s = s + coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
        //      }
        //      y[n] = s;
        //
        // We assume a minimum of four nodes for a 3D element so we use a four-way unrolling

        s = s1 = s2 = s3 = s4 = 0.0;
        for(i=0;i<=element->Type.BasisFunctions[n].n-4;i+=4) {
            s1 = s1 + coeff[i+0] * pow(u, (double)p[i+0]) * pow(v, (double)q[i+0]) * pow(w, (double)r[i+0]);
            s2 = s2 + coeff[i+1] * pow(u, (double)p[i+1]) * pow(v, (double)q[i+1]) * pow(w, (double)r[i+1]);
            s3 = s3 + coeff[i+2] * pow(u, (double)p[i+2]) * pow(v, (double)q[i+2]) * pow(w, (double)r[i+2]);
            s4 = s4 + coeff[i+3] * pow(u, (double)p[i+3]) * pow(v, (double)q[i+3]) * pow(w, (double)r[i+3]);
        }
        // Remainder loop
        for( ;i<element->Type.BasisFunctions[n].n;i++) {
            s = s + coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
        }
        y[n] = s + s1 + s2 + s3 + s4;
    }
}

inline void NodalBasisFunctions(int n, double * __nonnull Basis, Element_t * __nonnull element, double u, double v,double w) {
    int q, dim;
    double NodalBasis[n];
    
    dim = element->Type.dimension;
    
    // TODO: Add support for p-elements
    
    switch(dim) {
        case 1:
            NodalBasisFunctions1D(Basis, element, u);
            break;
            
        case 2:
            NodalBasisFunctions2D(Basis, element, u, v);
            break;
            
        case 3:
            if( element->Type.ElementCode/100 == 6) {
                
                memset( NodalBasis, 0.0, sizeof(NodalBasis) );
                for(q=0;q<n;q++) {
                    NodalBasis[q] = 1.0;
                    Basis[q] = InterpolateInElement3D(element, NodalBasis, u, v, w);
                    NodalBasis[q] = 0.0;
                }
            } else {
                NodalBasisFunctions3D(Basis, element, u, v, w);
            }
            break;
    }
}

/****************************************************************************
 
    NodalFirstDerivatives1D:
 
    Description:
 
    Arguments:
        Element_t *element -> element structure
        double u -> Point at which evaluate the value
 
    Function return value:
        double **y -> value of the quantity y = @x/@u
 
****************************************************************************/
inline void NodalFirstDerivatives1D(double * __nonnull y, Element_t * __nonnull element, double u) {
    
    int i, n;
    int *p = NULL;
    
    double s;
    double *coeff = NULL;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        s = 0.0;
        for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            if(p[i] >= 1) s = s + p[i]*coeff[i]*pow(u, (double)(p[i]-1));
        }
        y[3*n] = s;
    }
}

/****************************************************************************
 
    NodalFirstDerivatives2D:
 
    Description:
 
    Arguments:
        Element_t *element -> element structure
        double u, v -> Points at which evaluate the value
 
    Function return value:
        double **y -> value of the quantity y = @x(u,v)/@u
 
****************************************************************************/
inline void NodalFirstDerivatives2D(double * __nonnull y, Element_t * __nonnull element, double u, double v) {
    
    int i, n;
    int *p = NULL, *q = NULL;
    
    double s, s1, s2, s3, t, t1, t2, t3;
    double *coeff = NULL;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        q = element->Type.BasisFunctions[n].q;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        // The following is the unrolled version of this loop:
        //
        //       s = 0.0;
        //       t = 0.0;
        //       for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
        //          if(p[i] >= 1) s = s + p[i]*coeff[i]*pow(u, (double)(p[i]-1))*pow(v, (double)q[i]);
        //          if(q[i] >= 1) t = t + q[i]*coeff[i]*pow(u, (double)p[i])*pow(v, (double)(q[i]-1));
        //       }
        //       y[3*n] = s;
        //       y[3*n+1] = t;
        //
        // We assume a minimum of three nodes for a 2D element so we use a
        // three-way unrolling

        s = s1 = s2 = s3 = 0.0;
        t = t1 = t2 = t3 = 0.0;
        for(i=0;i<=element->Type.BasisFunctions[n].n-3;i+=3) {
            if(p[i+0] >= 1) s1 = s1 + p[i+0]*coeff[i+0]*pow(u, (double)(p[i+0]-1))*pow(v, (double)q[i+0]);
            if(p[i+1] >= 1) s2 = s2 + p[i+1]*coeff[i+1]*pow(u, (double)(p[i+1]-1))*pow(v, (double)q[i+1]);
            if(p[i+2] >= 1) s3 = s3 + p[i+2]*coeff[i+2]*pow(u, (double)(p[i+2]-1))*pow(v, (double)q[i+2]);
            
            if(q[i+0] >= 1) t1 = t1 + q[i+0]*coeff[i+0]*pow(u, (double)p[i+0])*pow(v, (double)(q[i+0]-1));
            if(q[i+1] >= 1) t2 = t2 + q[i+1]*coeff[i+1]*pow(u, (double)p[i+1])*pow(v, (double)(q[i+1]-1));
            if(q[i+2] >= 1) t3 = t3 + q[i+2]*coeff[i+2]*pow(u, (double)p[i+2])*pow(v, (double)(q[i+2]-1));
        }
        // Remainder loop
        for( ;i<element->Type.BasisFunctions[n].n;i++) {
            if(p[i] >= 1) s = s + p[i]*coeff[i]*pow(u, (double)(p[i]-1))*pow(v, (double)q[i]);
            if(q[i] >= 1) t = t + q[i]*coeff[i]*pow(u, (double)p[i])*pow(v, (double)(q[i]-1));
        }
        y[3*n] = s + s1 + s2 + s3;
        y[3*n+1] = t + t1 + t2 + t3;
    }
}


/***************************************************************************************
 
    FirstDerivativeInU3D:
 
    Description: Given element structure return value of the first partial derivative
    with respect to local coordinate u of a quantity x given at element nodes at local
    coordinate point u,v,w inside the element. Element basis functions are
    used to compute the value.
 
    Arguments:
        Element_t *element -> element structure
        double *x -> Nodal values of the quantity whose value we want to know
        double u, v, w -> Points at which evaluate the value
 
    Function return value:
        double y -> value of the quantity y = @x(u,v,w)/@u
 
***************************************************************************************/
inline double FirstDerivativeInU3D(Element_t * __nonnull element, double * __nonnull x, double u, double v, double w) {
    
    int i, l, n;
    int *p = NULL, *q = NULL, *r = NULL;
    
    double y, s, s1, s2, s3, s4;
    double *coeff = NULL;
    
    l = element->Type.BasisFunctionDegree;
    
    if(element->Type.ElementCode == 605) {
        if( w == 1.0) w = 1.0-1.0-12.0;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * ( -(1.0-v) + v*w *s ) / 4.0;
        y = y + x[1] * (  (1.0-v) - v*w * s ) / 4.0;
        y = y + x[2] * (  (1.0+v) + v*w * s ) / 4.0;
        y = y + x[3] * ( -(1.0+v) - v*w * s ) / 4.0;
        return y;
    }
    else if(element->Type.ElementCode == 613) {
        if( w == 1.0) w = 1.0-1.0-12.0;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0]  * ( -( (1.0-u) * (1.0-v) - w + u*v*w * s ) + (-u-v-1.0) * ( -(1.0-v) + v*w * s ) ) / 4.0;
        
        y = y + x[1]  * (  ( (1.0+u) * (1.0-v) - w - u*v*w * s ) + ( u-v-1.0) * (  (1.0-v) - v*w * s ) ) / 4.0;
        
        y = y + x[2]  * (  ( (1.0+u) * (1.0+v) - w + u*v*w * s ) + ( u+v-1.0) * (  (1.0+v) + v*w * s ) ) / 4.0;
        
        y = y + x[3]  * ( -( (1.0-u) * (1.0+v) - w - u*v*w * s ) + (-u+v-1.0) * ( -(1.0+v) - v*w * s ) ) / 4.0;
        
        y = y + x[4]  * 0.0;
        
        y = y + x[5]  * (  (1.0-u-w)*(1-v-w) - (1.0+u-w)*(1.0-v-w) ) * s / 2.0;
        y = y + x[6]  * (  (1.0+v-w)*(1.0-v-w) ) * s / 2.0;
        y = y + x[7]  * (  (1.0-u-w)*(1.0+v-w) - (1.0+u-w)*(1.0+v-w) ) * s / 2.0;
        y = y + x[8]  * ( -(1.0+v-w)*(1.0-v-w) ) * s / 2.0;
        
        y = y - x[9] * w * (1.0-v-w) * s;
        y = y + x[10] * w * (1.0-v-w) * s;
        y = y + x[11] * w * (1.0+v-w) * s;
        y = y - x[12] * w * (1.0+v-w) * s;
        
        return y;
    }
    
    y = 0.0;
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        if(x[n] != 0.0) {
            
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            r = element->Type.BasisFunctions[n].r;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            // The following is the unrolled version of this loop:
            //
            //      s = 0.0;
            //      for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            //          if(p[i] >= 1) s = s + p[i] * coeff[i] * pow(u, (double)(p[i]-1)) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
            //      }
            //      y = y + s*x[n];
            //
            // We assume a minimum of four nodes for a 3D element so we use a four-way unrolling

            s = s1 = s2 = s3 = s4 = 0.0;
            for(i=0;i<=element->Type.BasisFunctions[n].n-4;i+=4) {
                if(p[i+0] >= 1) s1 = s1 + p[i+0] * coeff[i+0] * pow(u, (double)(p[i+0]-1)) * pow(v, (double)q[i+0]) * pow(w, (double)r[i+0]);
                if(p[i+1] >= 1) s2 = s2 + p[i+1] * coeff[i+1] * pow(u, (double)(p[i+1]-1)) * pow(v, (double)q[i+1]) * pow(w, (double)r[i+1]);
                if(p[i+2] >= 1) s3 = s3 + p[i+2] * coeff[i+2] * pow(u, (double)(p[i+2]-1)) * pow(v, (double)q[i+2]) * pow(w, (double)r[i+2]);
                if(p[i+3] >= 1) s4 = s4 + p[i+3] * coeff[i+3] * pow(u, (double)(p[i+3]-1)) * pow(v, (double)q[i+3]) * pow(w, (double)r[i+3]);
            }
            y = y + (s1 + s2 + s3 + s4)*x[n];
            // Remainder loop
            for( ;i<element->Type.BasisFunctions[n].n;i++) {
                if(p[i] >= 1) s = s + p[i] * coeff[i] * pow(u, (double)(p[i]-1)) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
            }
            y = y + s*x[n];
        }
    }
    
    return y;
}

/***************************************************************************************
 
    FirstDerivativeInV3D:
 
    Description: Given element structure return value of the first partial derivative
    with respect to local coordinate v of a quantity x given at element nodes at local
    coordinate point u,v,w inside the element. Element basis functions are
    used to compute the value.
 
    Arguments:
        Element_t *element -> element structure
        double *x -> Nodal values of the quantity whose value we want to know
        double u, v, w -> Points at which evaluate the value
 
    Function return value:
        double y -> value of the quantity y = @x(u,v,w)/@v
 
****************************************************************************************/
inline double FirstDerivativeInV3D(Element_t * __nonnull element, double * __nonnull x, double u, double v, double w) {
    
    int i, l, n;
    int *p = NULL, *q = NULL, *r = NULL;
    
    double y, s, s1, s2, s3, s4;
    double *coeff = NULL;
    
    l = element->Type.BasisFunctionDegree;
    
    if(element->Type.ElementCode == 605) {
        if( w == 1.0) w = 1.0-1.0-12.0;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * ( -(1.0-u) + u*w * s ) / 4.0;
        y = y + x[1] * ( -(1.0+u) - u*w * s ) / 4.0;
        y = y + x[2] * (  (1.0+u) + u*w * s ) / 4.0;
        y = y + x[3] * (  (1.0-u) - u*w * s ) / 4.0;
        
        return y;
    }
    else if(element->Type.ElementCode == 613) {
        if( w == 1.0) w = 1.0-1.0-12.0;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0]  * ( -( (1.0-u) * (1.0-v) - w + u*v*w * s ) + (-u-v-1.0) * ( -(1.0-u) + u*w * s ) ) / 4.0;
        
        y = y + x[1]  * ( -( (1.0+u) * (1.0-v) - w - u*v*w * s ) + ( u-v-1.0) * ( -(1.0+u) - u*w * s ) ) / 4.0;
        
        y = y + x[2]  * (  ( (1.0+u) * (1.0+v) - w + u*v*w * s ) + ( u+v-1.0) * (  (1.0+u) + u*w * s ) ) / 4.0;
        
        y = y + x[3]  * (  ( (1.0-u) * (1.0+v) - w - u*v*w * s ) + (-u+v-1.0) * (  (1.0-u) - u*w * s ) ) / 4.0;
        
        y = y + x[4]  * 0.0;
        
        y = y - x[5]  *  (1.0+u-w)*(1.0-u-w) * s / 2.0;
        y = y + x[6]  * ( (1.0-v-w)*(1.0+u-w) - (1.0+v-w)*(1.0+u-w) ) * s / 2.0;
        y = y + x[7]  *  (1.0+u-w)*(1.0-u-w) * s / 2.0;
        y = y + x[8]  * ( (1.0-v-w)*(1.0-u-w) - (1.0+v-w)*(1.0-u-w) ) * s / 2.0;
        
        y = y - x[9] *  w * (1.0-u-w) * s;
        y = y - x[10] *  w * (1.0+u-w) * s;
        y = y + x[11] *  w * (1.0+u-w) * s;
        y = y + x[12] *  w * (1.0-u-w) * s;
        
        return y;
    }
    
    y = 0.0;
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        if(x[n] != 0.0) {
            
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            r = element->Type.BasisFunctions[n].r;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            // The following is the unrolled version of this loop:
            //
            //      s = 0.0;
            //      for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            //          if(q[i] >= 1) s = s + q[i] * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)(q[i]-1)) * pow(w, (double)r[i]);
            //      }
            //      y = y + s*x[n];
            //
            // We assume a minimum of four nodes for a 3D element so we use a four-way unrolling
            
            s = s1 = s2 = s3 = s4 = 0.0;
            for(i=0;i<=element->Type.BasisFunctions[n].n-4;i+=4) {
                if(q[i+0] >= 1) s1 = s1 + q[i+0] * coeff[i+0] * pow(u, (double)p[i+0]) * pow(v, (double)(q[i+0]-1)) * pow(w, (double)r[i+0]);
                if(q[i+1] >= 1) s2 = s2 + q[i+1] * coeff[i+1] * pow(u, (double)p[i+1]) * pow(v, (double)(q[i+1]-1)) * pow(w, (double)r[i+1]);
                if(q[i+2] >= 1) s3 = s3 + q[i+2] * coeff[i+2] * pow(u, (double)p[i+2]) * pow(v, (double)(q[i+2]-1)) * pow(w, (double)r[i+2]);
                if(q[i+3] >= 1) s4 = s4 + q[i+3] * coeff[i+3] * pow(u, (double)p[i+3]) * pow(v, (double)(q[i+3]-1)) * pow(w, (double)r[i+3]);
            }
            y = y + (s1 + s2 + s3 + s4)*x[n];
            // Remainder loop
            for( ;i<element->Type.BasisFunctions[n].n;i++) {
                if(q[i] >= 1) s = s + q[i] * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)(q[i]-1)) * pow(w, (double)r[i]);
            }
            y = y + s*x[n];
        }
    }
    
    return y;
}

/****************************************************************************************
 
    FirstDerivativeInW3D:
 
    Description: Given element structure return value of the first partial derivative
    with respect to local coordinate w of a quantity x given at element nodes at local
    coordinate point u,v,w inside the element. Element basis functions are
    used to compute the value.
 
    Arguments:
        Element_t *element -> element structure
        double *x -> Nodal values of the quantity whose value we want to know
        double u, v, w -> Points at which evaluate the value
 
    Function return value:
        double y -> value of the quantity y = @x(u,v,w)/@w
 
*****************************************************************************************/
inline double FirstDerivativeInW3D(Element_t * __nonnull element, double * __nonnull x, double u, double v, double w) {
    
    int i, l, n;
    int *p = NULL, *q = NULL, *r = NULL;
    
    double y, s, s1, s2, s3, s4;
    double *coeff = NULL;
    
    l = element->Type.BasisFunctionDegree;
    
    if(element->Type.ElementCode == 605) {
        if( w == 1.0) w = 1.0-1.0-12.0;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0] * ( -1.0 + u*v*(2.0-w) * pow(s, 2.0) ) / 4.0;
        y = y + x[1] * ( -1.0 - u*v*(2.0-w) * pow(s, 2.0) ) / 4.0;
        y = y + x[2] * ( -1.0 + u*v*(2.0-w) * pow(s, 2.0 ) ) / 4.0;
        y = y + x[3] * ( -1.0 - u*v*(2.0-w) * pow(s, 2.0) ) / 4.0;
        y = y + x[4];
        return y;
    }
    else if(element->Type.ElementCode == 613) {
        if( w == 1.0) w = 1.0-1.0-12.0;
        s = 1.0 / (1.0-w);
        
        y = 0.0;
        y = y + x[0]  * (-u-v-1.0) * ( -1.0 + u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[1]  * ( u-v-1.0) * ( -1.0 - u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[2]  * ( u+v-1.0) * ( -1.0 + u*v*pow(s, 2.0) ) / 4.0;
        y = y + x[3]  * (-u+v-1.0) * ( -1.0 - u*v*pow(s, 2.0) ) / 4.0;
        
        y = y + x[4]  * (4.0*w-1.0);
        
        y = y + x[5]  * ( ( -(1.0-u-w)*(1.0-v-w) - (1.0+u-w)*(1.0-v-w) - (1.0+u-w)*(1.0-u-w) ) * s + ( 1.0+u-w)*(1.0-u-w)*(1.0-v-w) * pow(s, 2.0) ) / 2.0;
        
        y = y + x[6]  * ( ( -(1.0-v-w)*(1.0+u-w) - (1.0+v-w)*(1.0+u-w) - (1.0+v-w)*(1.0-v-w) ) * s + ( 1.0+v-w)*(1.0-v-w)*(1.0+u-w) * pow(s, 2.0) ) / 2.0;
        
        y = y + x[7]  * ( ( -(1.0-u-w)*(1.0+v-w) - (1.0+u-w)*(1.0+v-w) - (1.0+u-w)*(1.0-u-w) ) * s + ( 1.0+u-w)*(1.0-u-w)*(1.0+v-w) * pow(s, 2.0) ) / 2.0;
        
        y = y + x[8]  * ( ( -(1.0-v-w)*(1.0-u-w) - (1.0+v-w)*(1.0-u-w) - (1.0+v-w)*(1.0-v-w) ) * s + ( 1.0+v-w)*(1.0-v-w)*(1.0-u-w) * pow(s, 2.0) ) / 2.0;
        
        y = y + x[9] * ( ( (1.0-u-w) * (1.0-v-w) - w * (1.0-v-w) - w * (1.0-u-w) ) * s  + w * (1.0-u-w) * (1.0-v-w) * pow(s, 2.0) );
        
        y = y + x[10] * ( ( (1.0+u-w) * (1.0-v-w) - w * (1.0-v-w) - w * (1.0+u-w) ) * s  + w * (1.0+u-w) * (1.0-v-w) * pow(s, 2.0) );
        
        y = y + x[11] * ( ( (1.0+u-w) * (1.0+v-w) - w * (1.0+v-w) - w * (1.0+u-w) ) * s  + w * (1.0+u-w) * (1.0+v-w) * pow(s, 2.0) );
        
        y = y + x[12] * ( ( (1.0-u-w) * (1.0+v-w) - w * (1.0+v-w) - w * (1.0-u-w) ) * s  + w * (1.0-u-w) * (1.0+v-w) * pow(s, 2.0) );
        
        return y;
    }
    
    y = 0.0;
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        if(x[n] != 0.0) {
            
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            r = element->Type.BasisFunctions[n].r;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            // The following is the unrolled version of this loop:
            //
            //      s = 0.0;
            //      for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            //          if(r[i] >= 1) s = s + r[i] * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)(r[i]-1));
            //      }
            //      y = y + s*x[n];
            //
            // We assume a minimum of four nodes for a 3D element so we use a four-way unrolling

            s = s1 = s2 = s3 = s4 = 0.0;
            for(i=0;i<=element->Type.BasisFunctions[n].n-4;i+=4) {
                if(r[i+0] >= 1) s1 = s1 + r[i+0] * coeff[i+0] * pow(u, (double)p[i+0]) * pow(v, (double)q[i+0]) * pow(w, (double)(r[i+0]-1));
                if(r[i+1] >= 1) s2 = s2 + r[i+1] * coeff[i+1] * pow(u, (double)p[i+1]) * pow(v, (double)q[i+1]) * pow(w, (double)(r[i+1]-1));
                if(r[i+2] >= 1) s3 = s3 + r[i+2] * coeff[i+2] * pow(u, (double)p[i+2]) * pow(v, (double)q[i+2]) * pow(w, (double)(r[i+2]-1));
                if(r[i+3] >= 1) s4 = s4 + r[i+3] * coeff[i+3] * pow(u, (double)p[i+3]) * pow(v, (double)q[i+3]) * pow(w, (double)(r[i+3]-1));
            }
            y = y + (s1 + s2 + s3 + s4)*x[n];
            // Remainder loop
            for( ;i<element->Type.BasisFunctions[n].n;i++) {
                if(r[i] >= 1) s = s + r[i] * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)(r[i]-1));
            }
            y = y + s*x[n];
        }
    }
    
    return y;
}

/****************************************************************************
 
    NodalFirstDerivatives3D:
 
    Description:
 
    Arguments:
        Element_t *element -> element structure
        double u, v, w -> Points at which evaluate the value
 
    Function return value:
        double **y -> value of the quantity y =  @x(u,v,w)/@u
 
****************************************************************************/
inline void NodalFirstDerivatives3D(double * __nonnull y, Element_t * __nonnull element, double u, double v, double w) {
    
    int i, n;
    int *p = NULL, *q = NULL, *r = NULL;
    
    double s, s1, s2, s3, s4;
    double t, t1, t2, t3, t4;
    double z, z1, z2, z3, z4;
    double *coeff = NULL;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        q = element->Type.BasisFunctions[n].q;
        r = element->Type.BasisFunctions[n].r;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        // The following is the unrolled version of this loop:
        //
        //      s = 0.0;
        //      t = 0.0;
        //      z = 0.0;
        //
        //      for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
        //         if(p[i] >= 1) s = s + p[i]*coeff[i]*pow(u, (double)(p[i]-1))*pow(v, (double)q[i])*pow(w, (double)r[i]);
        //         if(q[i] >= 1) t = t + q[i]*coeff[i]*pow(u, (double)p[i])*pow(v, (double)(q[i]-1))*pow(w, (double)r[i]);
        //         if(r[i] >= 1) z = z + r[i]*coeff[i]*pow(u, (double)p[i])*pow(v, (double)q[i])*pow(w, (double)(r[i]-1));
        //      }
        //      y[3*n] = s;
        //      y[3*n+1] = t;
        //      y[3*n+2] = z;
        //
        // We assume a minimum of four nodes for a 3D element so we use a four-way unrolling
        
        s = s1 = s2 = s3 = s4 = 0.0;
        t = t1 = t2 = t3 = t4 = 0.0;
        z = z1 = z2 = z3 = z4 = 0.0;
        for(i=0;i<=element->Type.BasisFunctions[n].n-4;i+=4) {
            if(p[i+0] >= 1) s1 = s1 + p[i+0]*coeff[i+0]*pow(u, (double)(p[i+0]-1))*pow(v, (double)q[i+0])*pow(w, (double)r[i+0]);
            if(p[i+1] >= 1) s2 = s2 + p[i+1]*coeff[i+1]*pow(u, (double)(p[i+1]-1))*pow(v, (double)q[i+1])*pow(w, (double)r[i+1]);
            if(p[i+2] >= 1) s3 = s3 + p[i+2]*coeff[i+2]*pow(u, (double)(p[i+2]-1))*pow(v, (double)q[i+2])*pow(w, (double)r[i+2]);
            if(p[i+3] >= 1) s4 = s4 + p[i+3]*coeff[i+3]*pow(u, (double)(p[i+3]-1))*pow(v, (double)q[i+3])*pow(w, (double)r[i+3]);

            if(q[i+0] >= 1) t1 = t1 + q[i+0]*coeff[i+0]*pow(u, (double)p[i+0])*pow(v, (double)(q[i+0]-1))*pow(w, (double)r[i+0]);
            if(q[i+1] >= 1) t2 = t2 + q[i+1]*coeff[i+1]*pow(u, (double)p[i+1])*pow(v, (double)(q[i+1]-1))*pow(w, (double)r[i+1]);
            if(q[i+2] >= 1) t3 = t3 + q[i+2]*coeff[i+2]*pow(u, (double)p[i+2])*pow(v, (double)(q[i+2]-1))*pow(w, (double)r[i+2]);
            if(q[i+3] >= 1) t4 = t4 + q[i+3]*coeff[i+3]*pow(u, (double)p[i+3])*pow(v, (double)(q[i+3]-1))*pow(w, (double)r[i+3]);
            
            if(r[i+0] >= 1) z1 = z1 + r[i+0]*coeff[i+0]*pow(u, (double)p[i+0])*pow(v, (double)q[i+0])*pow(w, (double)(r[i+0]-1));
            if(r[i+1] >= 1) z2 = z2 + r[i+1]*coeff[i+1]*pow(u, (double)p[i+1])*pow(v, (double)q[i+1])*pow(w, (double)(r[i+1]-1));
            if(r[i+2] >= 1) z3 = z3 + r[i+2]*coeff[i+2]*pow(u, (double)p[i+2])*pow(v, (double)q[i+2])*pow(w, (double)(r[i+2]-1));
            if(r[i+3] >= 1) z4 = z4 + r[i+3]*coeff[i+3]*pow(u, (double)p[i+3])*pow(v, (double)q[i+3])*pow(w, (double)(r[i+3]-1));

        }
        for( ;i<element->Type.BasisFunctions[n].n;i++) {
            if(p[i] >= 1) s = s + p[i]*coeff[i]*pow(u, (double)(p[i]-1))*pow(v, (double)q[i])*pow(w, (double)r[i]);
            if(q[i] >= 1) t = t + q[i]*coeff[i]*pow(u, (double)p[i])*pow(v, (double)(q[i]-1))*pow(w, (double)r[i]);
            if(r[i] >= 1) z = z + r[i]*coeff[i]*pow(u, (double)p[i])*pow(v, (double)q[i])*pow(w, (double)(r[i]-1));
        }
        y[3*n] = s + s1 + s2 + s3 + s4;
        y[3*n+1] = t + t1 + t2 + t3 + t4;
        y[3*n+2] = z + z1 + z2 + z3 + z4;
    }
}

inline void NodalFirstDerivatives(int n, double * __nonnull dLBasisdx, Element_t * __nonnull element, double u, double v, double w) {
    
    int q, dim;
    double NodalBasis[n];
    
    dim = element->Type.dimension;
    
    switch(dim) {
        case 1:
            NodalFirstDerivatives1D(dLBasisdx, element, u);
            break;
            
        case 2:
            NodalFirstDerivatives2D(dLBasisdx, element, u, v);
            break;
        case 3:
            if(element->Type.ElementCode / 100 == 6) {
                
                memset( NodalBasis, 0.0, sizeof(NodalBasis) );
                for(q=0;q<n;q++) {
                    NodalBasis[q] = 1.0;
                    dLBasisdx[3*q] = FirstDerivativeInU3D(element, NodalBasis, u, v, w);
                    dLBasisdx[3*q+1] = FirstDerivativeInV3D(element, NodalBasis, u, v, w);
                    dLBasisdx[3*q+2] = FirstDerivativeInW3D(element, NodalBasis, u, v, w);
                    NodalBasis[q] = 0.0;
                }
            }
            else {
                NodalFirstDerivatives3D(dLBasisdx, element, u, v, w);
            }
            break;
    }
}

/*******************************************************************************************************
 
    SecondDerivatives1D:
 
    Description: Given element structure return value of the second partial derivative with
    respect to local coordinate of a quantity x given at element nodes at local
    coordinate point u inside the element. Element basis functions are used to
    compute the value.
 
    Arguments:
        Element_t *element -> element structure
        double *nodes      -> Nodal values of the quantity whose partial derivative we want to know
        double u           -> Point at which evaluate the partial derivative
 
    Function return value:
        double y -> value of the quantity y = @x/@u
 
*******************************************************************************************************/
inline double SecondDerivatives1D(Element_t * __nonnull element, double * __nonnull nodes, double u) {
    
    double y;
    
    int i, n;
    int *p = NULL;
    double s;
    double *coeff = NULL;
    
    y = 0.0;
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        if (nodes[n] != 0.0) {
            p = element->Type.BasisFunctions[n].p;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            s = 0.0;
            for (i=0; i<element->Type.BasisFunctions[n].n; i++) {
                if (p[i] >= 2) s = s + p[i] * (p[i]-1) * coeff[i] * pow(u, (double)(p[i]-2));
            }
            y = y + s *nodes[n];
        }
    }
    
    return y;
}

/*********************************************************************************************************
 
    SecondDerivatives2D:
 
    Description: Given element structure return value of the second partial derivatives
    with respect to local coordinates of a quantity x given at element nodes at local
    coordinate point u,v inside the element. Element basis functions are used to
    compute the value.
 
    Arguments:
        double **ddx       -> Return matrix of second partial derivatives
        Element_t *element -> element structure
        double *nodes      -> Nodal values of the quantity whose partial derivative we want to know
        double u, v        -> Points at which evaluate the partial derivative
 
    Function return value:
        double y -> value of the quantity s = @^2x(u,v)/@v^2
 
*********************************************************************************************************/
inline void SecondDerivatives2D(double * __nonnull ddx, Element_t * __nonnull element, double * __nonnull nodes, double u, double v) {
    
    int i, n;
    int *p = NULL, *q = NULL;
    
    double s;
    double *coeff = NULL;
    
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        
        if (nodes[n] != 0.0) {
            
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            // @^2x/@u^2
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (p[i] >= 2) s = s + p[i] * (p[i]-1) * coeff[i] * pow(u, (double)(p[i]-2)) * pow(v, (double)q[i]);
            }
            ddx[0] = ddx[0] + s*nodes[n];
            
            // @^2x/@u@v
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (p[i] >= 1 && q[i] >= 1) s = s + p[i] * q[i] * coeff[i] * pow(u, (double)(p[i]-1)) * pow(v, (double)(q[i]-1));
            }
            ddx[1] = ddx[1] + s*nodes[n];
            
            // @^2x/@v^2
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (q[i] >= 2) s = s + q[i] * (q[i]-1) * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)(q[i]-2));
            }
            ddx[3] = ddx[3] + s*nodes[n];
            
        }
    }
    
    ddx[2] = ddx[1];
}

/*****************************************************************************************************
 
    SecondDerivatives3D:
 
    Description: Given element structure return value of the second partial derivatives
    with respect to local coordinates of i quantity x given at element nodes at local
    coordinate point u,v, w inside the element. Element basis functions are used
    to compute the value.
 
    Arguments:
        double **ddx       -> Return matrix of second partial derivatives
        Element_t *element -> element structure
        double *nodes      -> Nodal values of the quantity whose partial derivative we want to know
        double u, v, w      -> Points at which evaluate the partial derivative
 
    Function return value:
        double y -> value of the quantity s = @^2x(u,v)/@v^2
 
*****************************************************************************************************/
inline void SecondDerivatives3D(double * __nonnull ddx, Element_t * __nonnull element, double * __nonnull nodes, double u, double v, double w) {
    
    int i, n;
    int *p = NULL, *q = NULL, *r = NULL;
    
    double s;
    double *coeff = NULL;
    
    for (n=0; n<element->Type.NumberOfNodes; n++) {
        
        if (nodes[n] != 0.0) {
            
            p = element->Type.BasisFunctions[n].p;
            q = element->Type.BasisFunctions[n].q;
            r = element->Type.BasisFunctions[n].r;
            coeff = element->Type.BasisFunctions[n].coeff;
            
            // @^2x/@u^2
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (p[i] >= 2) s = s + p[i] * (p[i]-1) * coeff[i] * pow(u, (double)(p[i]-2)) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
            }
            ddx[0] = ddx[0] + s*nodes[n];
            
            // @^2x/@u@v
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (p[i] >= 1 && q[i] >= 1) s = s + p[i] * q[i] * coeff[i] * pow(u, (double)(p[i]-1)) * pow(v, (double)(q[i]-1)) * pow(w, (double)r[i]);
            }
            ddx[1] = ddx[1] + s*nodes[n];
            
            // @^2x/@u@w
            s = 0.0;
            for (i=1; i<element->Type.NumberOfNodes; i++) {
                if (p[i] >= 1 && r[i] >= 1) s = s + p[i] * r[i] * coeff[i] * pow(u, (double)(p[i]-1)) * pow(v, (double)q[i]) * pow(w, (double)(r[i]-1));
            }
            ddx[2] = ddx[2] + s*nodes[n];
            
            // @^2x/@v^2
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (q[i] >= 2) s = s + q[i] * (q[i]-1) * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)(q[i]-2)) * pow(w, (double)r[i]);
            }
            ddx[4] = ddx[4] + s*nodes[n];
            
            // @^2x/@v@w
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (q[i] >= 1 && r[i]>= 1) s = s + q[i] * r[i] * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)(q[i]-1)) * pow(w, (double)(r[i]-1));
            }
            ddx[5] = ddx[5] + s*nodes[n];
            
            // @^2x/@w^2
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (r[i] >= 2) s = s + r[i] * (r[i]-1) * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)(r[i]-2));
            }
            ddx[8] = ddx[8] + s*nodes[n];
            
        }
    }
    
    ddx[3] = ddx[1];
    ddx[6] = ddx[2];
    ddx[7] = ddx[5];
}

