//
//  NodalBasisFunction.c
//  Saino
//
//  Created by Seddik hakime on 30/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#include <float.h>
#include <math.h>

#include "NodalBasisFunctions.h"


void NodalBasisFunctions1D(double *y, Element_t *element, double u) {
/****************************************************************************
 NodalBasisFunctions1D:
 
 Description:           
 
 Arguments:
    Element_t *element -> element structure
    double u -> Point at which evaluate the value
 
 Function return value:
    double *y -> value of the quantity y = x(u)
****************************************************************************/
    
    int i, n;
    int *p;
    
    double s;
    double *coeff;
    
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


void NodalFirstDerivatives1D(double **y, Element_t *element, double u) {
/****************************************************************************
 NodalFirstDerivatives1D:
 
 Description:           
 
 Arguments:
    Element_t *element -> element structure
    double u -> Point at which evaluate the value
 
 Function return value:
    double **y -> value of the quantity y = @x/@u
****************************************************************************/
    
    int i, n;
    int *p;
    
    double s;
    double *coeff;

    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        s = 0.0;
        for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            if(p[i] >= 1) s = s + p[i]*coeff[i]*pow(u, (double)(p[i]-1));
        }
        y[n][0] = s;
    }
}


void NodalBasisFunctions2D(double *y, Element_t *element, double u, double v) {
/****************************************************************************
 NodalBasisFunction2D:
 
 Description:           
 
 Arguments:
    Element_t *element -> element structure
    double u, v -> Points at which evaluate the value
 
 Function return value:
    double *y -> value of the quantity y = x(u)
****************************************************************************/

    int i, n;
    int *p, *q;
    
    double s;
    double *coeff;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        q = element->Type.BasisFunctions[n].q;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        s = 0.0;
        for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            s = s + coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]);
        }
        y[n] = s;
    }
}

void NodalFirstDerivatives2D(double **y, Element_t *element, double u, double v) {
/****************************************************************************
 NodalFirstDerivatives2D:
 
 Description:           
 
 Arguments:
    Element_t *element -> element structure
    double u, v -> Points at which evaluate the value
 
 Function return value:
    double **y -> value of the quantity y = @x(u,v)/@u
****************************************************************************/

    
    int i, n;
    int *p, *q;
    
    double s, t;
    double *coeff;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        q = element->Type.BasisFunctions[n].q;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        s = 0.0;
        t = 0.0;
        for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            if(p[i] >= 1) s = s + p[i]*coeff[i]*pow(u, (double)(p[i]-1))*pow(v, (double)q[i]);
            if(q[i] >= 1) t = t + q[i]*coeff[i]*pow(u, (double)p[i])*pow(v, (double)(q[i]-1));
        }
        y[n][0] = s;
        y[n][1] = t;
    }
}


double InterpolateInElement3D(Element_t *element, double *x, double u, double v, double w) {
/****************************************************************************
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
****************************************************************************/

    int i, l, n;
    int *p, *q, *r;
    
    double y, s;
    double *coeff;
    
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
            
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                s = s + coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
            }
            y = y + s*x[n];
        }
    }
    
    return y;
}


void NodalBasisFunctions3D(double *y, Element_t *element, double u, double v, double w) {
/****************************************************************************
 NodalBasisFunctions3D:
 
 Description:           
 
 Arguments:
    Element_t *element -> element structure
    double u, v, w -> Points at which evaluate the value
 
 Function return value:
    double *y -> value of the quantity y = x(u)
****************************************************************************/

    int i, n;
    int *p, *q, *r;
    
    double s;
    double *coeff;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        q = element->Type.BasisFunctions[n].q;
        r = element->Type.BasisFunctions[n].r;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        s = 0.0;
        for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            s = s + coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
        }
        y[n] = s;
    }
}


double FirstDerivativeInU3D(Element_t *element, double *x, double u, double v, double w) {
/****************************************************************************
 InterpolateInElement3D:
 
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
****************************************************************************/
    
    int i, l, n;
    int *p, *q, *r;
    
    double y, s;
    double *coeff;
    
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
            
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if(p[i] >= 1) s = s + p[i] * coeff[i] * pow(u, (double)(p[i]-1)) * pow(v, (double)q[i]) * pow(w, (double)r[i]);
            }
            y = y + s*x[n];
        }
    }
    
    return y;
}


double FirstDerivativeInV3D(Element_t *element, double *x, double u, double v, double w) {
/****************************************************************************
 InterpolateInElement3D:
 
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
****************************************************************************/
    
    int i, l, n;
    int *p, *q, *r;
    
    double y, s;
    double *coeff;
    
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
            
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if(q[i] >= 1) s = s + q[i] * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)(q[i]-1)) * pow(w, (double)r[i]);
            }
            y = y + s*x[n];
        }
    }
    
    return y;
}


double FirstDerivativeInW3D(Element_t *element, double *x, double u, double v, double w) {
/****************************************************************************
 InterpolateInElement3D:
 
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
****************************************************************************/
    
    int i, l, n;
    int *p, *q, *r;
    
    double y, s;
    double *coeff;
    
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
            
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if(r[i] >= 1) s = s + r[i] * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)(r[i]-1));
            }
            y = y + s*x[n];
        }
    }
    
    return y;
}


void NodalFirstDerivatives3D(double **y, Element_t *element, double u, double v, double w) {
/****************************************************************************
 NodalBasisFunctions3D:
 
 Description:           
 
 Arguments:
    Element_t *element -> element structure
    double u, v, w -> Points at which evaluate the value
 
 Function return value:
    double **y -> value of the quantity y =  @x(u,v,w)/@u
****************************************************************************/

    int i, n;
    int *p, *q, *r;
    
    double s, t, z;
    double *coeff;
    
    for(n=0;n<element->Type.NumberOfNodes;n++) {
        
        p = element->Type.BasisFunctions[n].p;
        q = element->Type.BasisFunctions[n].q;
        r = element->Type.BasisFunctions[n].r;
        coeff = element->Type.BasisFunctions[n].coeff;
        
        s = 0.0;
        t = 0.0;
        z = 0.0;
        
        for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
            if(p[i] >= 1) s = s + p[i]*coeff[i]*pow(u, (double)(p[i]-1))*pow(v, (double)q[i])*pow(w, (double)r[i]);
            if(q[i] >= 1) t = t + q[i]*coeff[i]*pow(u, (double)p[i])*pow(v, (double)(q[i]-1))*pow(w, (double)r[i]);
            if(r[i] >= 1) z = z + r[i]*coeff[i]*pow(u, (double)p[i])*pow(v, (double)q[i])*pow(w, (double)(r[i]-1));
        }
        y[n][0] = s;
        y[n][1] = t;
        y[n][2] = z;
    }
}


void NodalBasisFunctions(int n, double *Basis, Element_t *element, double u, double v,double w)
{
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

void NodalFirstDerivatives(int n, double **dLBasisdx, Element_t *element, double u, double v, double w) {
    
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
                    dLBasisdx[q][0] = FirstDerivativeInU3D(element, NodalBasis, u, v, w);
                    dLBasisdx[q][1] = FirstDerivativeInV3D(element, NodalBasis, u, v, w);
                    dLBasisdx[q][2] = FirstDerivativeInW3D(element, NodalBasis, u, v, w);
                    NodalBasis[q] = 0.0;
                }
            }
            else {
                NodalFirstDerivatives3D(dLBasisdx, element, u, v, w);
            }
            break;
    }
}

double SecondDerivatives1D(Element_t* element, double *nodes, double u) {
/****************************************************************************
 SecondDerivatives1D:
 
 Description:
    Given element structure return value of the second partial derivative with
    respect to local coordinate of a quantity x given at element nodes at local
    coordinate point u inside the element. Element basis functions are used to
    compute the value. 
 
 Arguments:
 Element_t *element -> element structure
 double *nodes      -> Nodal values of the quantity whose partial derivative we want to know
 double u           -> Point at which evaluate the partial derivative
 
 Function return value:
 double y -> value of the quantity y = @x/@u
****************************************************************************/

    
    double y;
    
    int i, n;
    int *p;
    double s;
    double *coeff;
    
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

void SecondDerivatives2D(double **ddx, Element_t* element, double *nodes, double u, double v) {
/****************************************************************************
 SecondDerivatives2D:
 
 Description:
 Given element structure return value of the second partial derivatives with
 respect to local coordinates of a quantity x given at element nodes at local
 coordinate point u,v inside the element. Element basis functions are used to
 compute the value. 

 Arguments:
 double **ddx       -> Return matrix of second partial derivatives
 Element_t *element -> element structure
 double *nodes      -> Nodal values of the quantity whose partial derivative we want to know
 double u, v        -> Points at which evaluate the partial derivative
 
 Function return value:
 double y -> value of the quantity s = @^2x(u,v)/@v^2
****************************************************************************/
    
    int i, n;
    int *p, *q;
    
    double s;
    double *coeff;
    
    memset( *ddx, 0.0, (2*2)*sizeof(double) );
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
            ddx[0][0] = ddx[0][0] + s*nodes[n];
                        
            // @^2x/@u@v
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (p[i] >= 1 && q[i] >= 1) s = s + p[i] * q[i] * coeff[i] * pow(u, (double)(p[i]-1)) * pow(v, (double)(q[i]-1));
            }
            ddx[0][1] = ddx[0][1] + s*nodes[n];
            
            // @^2x/@v^2
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (q[i] >= 2) s = s + q[i] * (q[i]-1) * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)(q[i]-2));
            }
            ddx[1][1] = ddx[1][1] + s*nodes[n];
            
        }   
    }
    
    ddx[1][0] = ddx[0][1];
}

void SecondDerivatives3D(double **ddx, Element_t* element, double *nodes, double u, double v, double w) {
/****************************************************************************
 SecondDerivatives3D:
 
 Description:
 Given element structure return value of the second partial derivatives with
 respect to local coordinates of i quantity x given at element nodes at local
 coordinate point u,v, w inside the element. Element basis functions are used 
 to compute the value. 
 
 Arguments:
 double **ddx       -> Return matrix of second partial derivatives
 Element_t *element -> element structure
 double *nodes      -> Nodal values of the quantity whose partial derivative we want to know
 double u, v, w      -> Points at which evaluate the partial derivative
 
 Function return value:
 double y -> value of the quantity s = @^2x(u,v)/@v^2
****************************************************************************/
    
    int i, n;
    int *p, *q, *r;
    
    double s;
    double *coeff;
    
    memset( *ddx, 0.0, (3*3)*sizeof(double) );
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
            ddx[0][0] = ddx[0][0] + s*nodes[n];
                        
            // @^2x/@u@v
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (p[i] >= 1 && q[i] >= 1) s = s + p[i] * q[i] * coeff[i] * pow(u, (double)(p[i]-1)) * pow(v, (double)(q[i]-1)) * pow(w, (double)r[i]);
            }
            ddx[0][1] = ddx[0][1] + s*nodes[n];
                        
            // @^2x/@u@w
            s = 0.0;
            for (i=1; i<element->Type.NumberOfNodes; i++) {
                if (p[i] >= 1 && r[i] >= 1) s = s + p[i] * r[i] * coeff[i] * pow(u, (double)(p[i]-1)) * pow(v, (double)q[i]) * pow(w, (double)(r[i]-1));
            }
            ddx[0][2] = ddx[0][2] + s*nodes[n];
                        
            // @^2x/@v^2
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (q[i] >= 2) s = s + q[i] * (q[i]-1) * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)(q[i]-2)) * pow(w, (double)r[i]);
            }
            ddx[1][1] = ddx[1][1] + s*nodes[n];

            // @^2x/@v@w
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (q[i] >= 1 && r[i]>= 1) s = s + q[i] * r[i] * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)(q[i]-1)) * pow(w, (double)(r[i]-1));
            }
            ddx[1][2] = ddx[1][2] + s*nodes[n];
                        
            // @^2x/@w^2
            s = 0.0;
            for(i=0;i<element->Type.BasisFunctions[n].n;i++) {
                if (r[i] >= 2) s = s + r[i] * (r[i]-1) * coeff[i] * pow(u, (double)p[i]) * pow(v, (double)q[i]) * pow(w, (double)(r[i]-2));
            }
            ddx[2][2] = ddx[2][2] + s*nodes[n];
        
        }
    }
    
    ddx[1][0] = ddx[0][1];
    ddx[2][0] = ddx[0][2];
    ddx[2][1] = ddx[1][2];
}



























