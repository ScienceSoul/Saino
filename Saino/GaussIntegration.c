//
//  GaussIntegration.c
//  Saino
//
//  Created by Seddik hakime on 16/06/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#include <Accelerate/Accelerate.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>

#include "GaussIntegration.h"

static bool GInit = false;

/***************************************************************************/
/* Storage for 1D Gauss points and the weights. The values are computed on */
/* the fly (i.e., GaussQuadraturePoints1D). These vales are used for quads */
/* and bricks as well.                                                     */
/***************************************************************************/
static double AllPoints[MAXN][MAXN], AllWeights[MAXN][MAXN];

/***************************************************************************/
/* Triangle - 1 point rule; exact integration of x^py^q, p+q=<=1           */
/***************************************************************************/
static double UTriangle1P[1] = {0.3333333333333333E+00};
static double VTriangle1P[1] = {0.3333333333333333E+00};
static double STriangle1P[1] = {1.0000000000000000E+00};

/***************************************************************************/
/* Triangle - 3 point rule; exact integration of x^py^q, p+q<=2            */
/***************************************************************************/
static double UTriangle3P[3] = {
    0.16666666666666667E+00, 0.66666666666666667E+00, 0.16666666666666667E+00
};

static double VTriangle3P[3] = {
    0.16666666666666667E+00, 0.16666666666666667E+00, 0.66666666666666667E+00
};

static double STriangle3P[3] = {
    0.33333333333333333E+00, 0.33333333333333333E+00, 0.33333333333333333E+00
};

/***************************************************************************/
/* Triangle - 4 point rule; exact integration of x^py^q, p+q<=3            */
/***************************************************************************/
static double UTriangle4P[4] = {
    0.33333333333333333E+00, 0.2000000000000000E+00, 
    0.60000000000000000E+00, 0.20000000000000000E+00
};

static double VTriangle4P[4] = {
    0.33333333333333333E+00, 0.2000000000000000E+00, 
    0.20000000000000000E+00, 0.60000000000000000E+00
};

static double STriangle4P[4] = {
    -0.56250000000000000E+00, 0.52083333333333333E+00, 
    0.52083333333333333E+00, 0.52083333333333333E+00
};

/***************************************************************************/
/* Triangle - 6 point rule; exact integration of x^py^q, p+q,=4            */
/***************************************************************************/
static double UTriangle6P[6] = {
    0.091576213509771E+00, 0.816847572980459E+00, 0.091576213509771E+00,
    0.445948490915965E+00, 0.108103018168070E+00, 0.445948490915965E+00
};

static double VTriangle6P[6] = {
    0.091576213509771E+00, 0.091576213509771E+00, 0.816847572980459E+00,
    0.445948490915965E+00, 0.445948490915965E+00, 0.108103018168070E+00
};

static double STriangle6P[6] = {
    0.109951743655322E+00, 0.109951743655322E+00, 0.109951743655322E+00,
    0.223381589678011E+00, 0.223381589678011E+00, 0.223381589678011E+00
};

/***************************************************************************/
/* Triangle - 7 point rule; exact integration of x^py^q, p+q,=5            */
/***************************************************************************/
static double UTriangle7P[7] = {
    0.333333333333333E+00, 0.101286507323456E+00, 0.797426985353087E+00,
    0.101286507323456E+00, 0.470142064105115E+00, 0.059715871789770E+00,
    0.470142064105115E+00
};

static double VTriangle7P[7] = {
    0.333333333333333E+00, 0.101286507323456E+00, 0.101286507323456E+00,
    0.797426985353087E+00, 0.470142064105115E+00, 0.470142064105115E+00,
    0.059715871789770E+00
};

static double STriangle7P[7] = {
    0.225000000000000E+00, 0.125939180544827E+00, 0.125939180544827E+00,
    0.125939180544827E+00, 0.132394152788506E+00, 0.132394152788506E+00,
    0.132394152788506E+00
};

/***************************************************************************/
/* Triangle - 11 point rule; exact integration of x^py^q, p+q<=6           */
/***************************************************************************/
static double Utriangle11P[11] = {
    0.3019427231413448E-01, 0.5298143569082113E-01, 0.4972454892773975E-01,
    0.7697772693248785E-01, 0.7008117469890058E+00, 0.5597774797709894E+00,
    0.5428972301980696E+00, 0.3437947421925572E+00, 0.2356669356664465E+00,
    0.8672623210691472E+00, 0.2151020995173866E+00
};

static double VTriangle11P[11] = {
    0.2559891985673773E+00, 0.1748087863744473E-01, 0.6330812033358987E+00,
    0.8588528075577063E+00, 0.2708520519075563E+00, 0.1870602768957014E-01,
    0.2027008533579804E+00, 0.5718576583152437E+00, 0.1000777578531811E+00,
    0.4654861310422605E-01, 0.3929681357810497E+00
};

static double STriangle11P[11] = {
    0.3375321205342688E-01, 0.1148426034648707E-01, 0.4197958777582435E-01,
    0.3098130358202468E-01, 0.2925899761167147E-01, 0.2778515729102349E-01,
    0.8323049608963519E-01, 0.6825761580824108E-01, 0.6357334991651026E-01,
    0.2649352562792455E-01, 0.8320249389723097E-01
};

/***************************************************************************/
/* Triangle - 12 point; exact integration of x^py^q, p+q,=7                */
/***************************************************************************/
static double UTriangle12P[12] = {
    0.6232720494911090E+00, 0.3215024938520235E+00, 0.5522545665686063E-01,
    0.2777161669760336E+00, 0.5158423343535236E+00, 0.2064414986704435E+00,
    0.3432430294535058E-01, 0.3047265008682535E+00, 0.6609491961864082E+00,
    0.6238226509441210E-01, 0.8700998678316921E+00, 0.6751786707389329E-01
};

static double VTriangle12P[12] = {
    0.3215024938520235E+00, 0.5522545665686063E-01, 0.6232720494911090E+00,
    0.5158423343535236E+00, 0.2064414986704435E+00, 0.2777161669760336E+00,
    0.3047265008682535E+00, 0.6609491961864082E+00, 0.3432430294535058E-01,
    0.8700998678316921E+00, 0.6751786707389329E-01, 0.6238226509441210E-01
};

static double STriangle12P[12] = {
    0.4388140871440586E-01, 0.4388140871440586E-01, 0.4388140871440587E-01,
    0.6749318700971417E-01, 0.6749318700971417E-01, 0.6749318700971417E-01,
    0.2877504278510970E-01, 0.2877504278510970E-01, 0.2877504278510969E-01,
    0.2651702815743698E-01, 0.2651702815743698E-01, 0.2651702815743698E-01
};

/***************************************************************************/
/* Triangle - 17 point rule; exact integration of x^py^q, p+q<=8           */
/***************************************************************************/
static double UTriangle17P[17] = {
    0.2292423642627924E+00, 0.4951220175479885E-01, 0.3655948407066446E+00,
    0.4364350639589269E+00, 0.1596405673569602E+00, 0.9336507149305228E+00,
    0.5219569066777245E+00, 0.7110782758797098E+00, 0.5288509041694864E+00,
    0.1396967677642513E-01, 0.4205421906708996E-01, 0.4651359156686354E-01,
    0.1975981349257204E+00, 0.7836841874017514E+00, 0.4232808751402256E-01,
    0.4557097415216423E+00, 0.2358934246935281E+00
};

static double VTriangle17P[17] = {
    0.5117407211006358E+00, 0.7589103637479163E+00, 0.1529647481767193E+00,
    0.3151398735074337E-01, 0.5117868393288316E-01, 0.1964516824106966E-01,
    0.2347490459725670E+00, 0.4908682577187765E-01, 0.4382237537321878E+00,
    0.3300210677033395E-01, 0.2088758614636060E+00, 0.9208929246654702E+00,
    0.2742740954674795E+00, 0.1654585179097472E+00, 0.4930011699833554E+00,
    0.4080804967846944E+00, 0.7127872162741824E+00
};

static double STriangle17P[17] = {
    0.5956595662857148E-01, 0.2813390230006461E-01, 0.3500735477096827E-01,
    0.2438077450393263E-01, 0.2843374448051010E-01, 0.7822856634218779E-02, 
    0.5179111341004783E-01, 0.3134229539096806E-01, 0.2454951584925144E-01,
    0.5371382557647114E-02, 0.2571565514768072E-01, 0.1045933340802507E-01,
    0.4937780841212319E-01, 0.2824772362317584E-01, 0.3218881684015661E-01,
    0.2522089247693226E-01, 0.3239087356572598E-01
};

/***************************************************************************/
/* Triangle - 20 point rule; exact integration of x^py^q, p+q<=9           */
/***************************************************************************/
static double UTriangle20P[20] = {
    0.2469118866487856E-01, 0.3348782965514246E-00, 0.4162560937597861E-00,
    0.1832492889417431E-00, 0.2183952668281443E-00, 0.4523362527443628E-01,
    0.4872975112073226E-00, 0.7470127381316580E-00, 0.7390287107520658E-00,
    0.3452260444515281E-01, 0.4946745467572288E-00, 0.3747439678780460E-01,
    0.2257524791528391E-00, 0.9107964437563798E-00, 0.4254445629399445E-00,
    0.1332215072275240E-00, 0.5002480151788234E-00, 0.4411517722238269E-01,
    0.1858526744057914E-00, 0.6300024376672695E-00
};

static double VTriangle20P[20] = {
    0.4783451248176442E-00, 0.3373844236168988E-00, 0.1244378463254732E-00, 
    0.6365569723648120E-00, 0.3899759363237886E-01, 0.9093437140096456E-00, 
    0.1968266037596590E-01, 0.2191311129347709E-00, 0.3833588560240875E-01, 
    0.7389795063475102E-00, 0.4800989285800525E-00, 0.2175137165318176E-00, 
    0.7404716820879975E-00, 0.4413531509926682E-01, 0.4431292142978816E-00,
    0.4440953593652837E-00, 0.1430831401051367E-00, 0.4392970158411878E-01,
    0.1973209364545017E-00, 0.1979381059170009E-00
};

static double STriangle20P[20] = {
    0.1776913091122958E-01, 0.4667544936904065E-01, 0.2965283331432967E-01,
    0.3880447634997608E-01, 0.2251511457011248E-01, 0.1314162394636178E-01, 
    0.1560341736610505E-01, 0.1967065434689744E-01, 0.2247962849501080E-01, 
    0.2087108394969067E-01, 0.1787661200700672E-01, 0.2147695865607915E-01,
    0.2040998247303970E-01, 0.1270342300533680E-01, 0.3688713099356314E-01,
    0.3813199811535777E-01, 0.1508642325812160E-01, 0.1238422287692121E-01,
    0.3995072336992735E-01, 0.3790911262589247E-01
};

/***************************************************************************/
/* Tetrahedron - 1 point rule; exact integration of x^py^qz^r, p+q+r<=1    */
/***************************************************************************/
static double UTetra1P[1] = {0.25E-00};
static double Vtetra1P[1] = {0.25E-00};
static double WTetra1P[1] = {0.25E-00};
static double STetra1P[1] = {1.00E-00};

/***************************************************************************/
/* Tetrahedron - 4 point rule; exact integration of x^py^qz^r, p+q+r<=2    */
/***************************************************************************/
static double UTetra4P[4] = {
    0.1757281246520584E+00, 0.2445310270213291E+00,
    0.5556470949048655E+00, 0.0240937534217468E+00
};

static double VTetra4P[4] = {
    0.5656137776620919E+00, 0.0501800797762026E+00,
    0.1487681308666864E+00, 0.2354380116950194E+00
};

static double WTetra4P[4] = {
    0.2180665126782654E+00, 0.5635595064952189E+00,
    0.0350112499848832E+00, 0.1833627308416330E+00
};

static double STetra4P[4] = {
    0.2500000000000000E+00, 0.2500000000000000E+00,
    0.2500000000000000E+00, 0.2500000000000000E+00
};

/***************************************************************************/
/* Tetrahedron - 5 point rule; exact integration of x^py^qz^r, p+q+r<=3    */
/***************************************************************************/
static double UTetra5P[5] = {
    0.25000000000000000E+00, 0.50000000000000000E+00,
    0.16666666666666667E+00, 0.16666666666666667E+00,
    0.16666666666666667E+00
};

static double VTetra5P[5] = {
    0.25000000000000000E+00, 0.16666666666666667E+00,
    0.50000000000000000E+00, 0.16666666666666667E+00,
    0.16666666666666667E+00
};

static double WTetra5P[5] = {
    0.25000000000000000E+00, 0.16666666666666667E+00,
    0.16666666666666667E+00, 0.50000000000000000E+00,
    0.16666666666666667E+00
};

static double STetra5P[5] = {
   -0.80000000000000000E+00, 0.45000000000000000E+00,
    0.45000000000000000E+00, 0.45000000000000000E+00,
    0.45000000000000000E+00
};

/***************************************************************************/
/* Tetrahedron - 11 point rule; exact integration of x^py^qz^r, p+q+r<=4   */
/***************************************************************************/
static double UTetra11P[11] = {
    0.3247902050850455E+00, 0.4381969657060433E+00, 0.8992592373310454E-01,
    0.1092714936292849E+00, 0.3389119319942253E-01, 0.5332363613904868E-01,
    0.1935618747806815E+00, 0.4016250624424964E-01, 0.3878132182319405E+00,
    0.7321489692875428E+00, 0.8066342495294049E-01
};

static double VTetra11P[11] = {
    0.4573830181783998E+00, 0.9635325047480842E-01, 0.3499588148445295E+00,
    0.1228957438582778E+00, 0.4736224692062527E-01, 0.4450376952468180E+00,
    0.2165626476982170E+00, 0.8033385922433729E+00, 0.7030897281814283E-01,
    0.1097836536360084E+00, 0.1018859284267242E-01
};

static double WTetra11P[11] = {
    0.1116787541193331E+00, 0.6966288385119494E-01, 0.5810783971325720E-01,
    0.3424607753785182E+00, 0.7831772466208499E+00, 0.3688112094344830E+00,
    0.5872345323698884E+00, 0.6178518963560731E-01, 0.4077342860913465E+00,
    0.9607290317342082E-01, 0.8343823045787845E-01
};

static double STetra11P[11] = {
    0.1677896627448221E+00, 0.1128697325878004E+00, 0.1026246621329828E+00,
    0.1583002576888426E+00, 0.3847841737508437E-01, 0.1061709382037234E+00,
    0.5458124994014422E-01, 0.3684475128738168E-01, 0.1239234851349682E+00,
    0.6832098141300300E-01, 0.3009586149124714E-01
};

void DerivPoly(int n, double *Q, double *P) {
    
    int i;
    
    for(i=0;i<n;i++) {
        Q[i] = P[i] * (n-(i+1)+1);
    }
}

double EvalPoly(int n, double *P, double x) {
    int i;
    double s;
    
    s = 0.0;
    for(i=0;i<n+1;i++) {
        s = s * x + P[i];
    }
    
    return s;
}

void RefineRoots(int n, double *P, double *Q, double *Points) {
    
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

/**************************************************************************************
    Function to compute gaussian integration points and weights in [-1,1] as roots
    of Legendre polynomials
**************************************************************************************/
void GaussQuadraturePoints1D(int n) {
   
    int i, j, k, np, info;
    int arg1, arg2, arg3, arg4;
    char ch1, ch2;
    double Points[n], Weights[n];
    double A[n/2][n/2], *AT, s, Work[8*n];
    double P[n+1], Q[n], P0[n], P1[n+1];
    double sum;
        
    // One point is trivial
    if (n <=1) {
        AllPoints[0][0] = 0.0;
        AllWeights[0][0] = 2.0;
        return;
    }
    
    memset( Points, 0.0, sizeof(Points) );
    memset( Weights, 0.0, sizeof(Weights) );
    
    // Compute coefficients of n:th Legendre polynomial from the recurrence:
    // (i+p)P_{i+1}(x) = (2i+1)*x*P_i{x} - i*P_{i-1}(x), P_{0} = 1; P_{1} = x;
    // Caveat: Computed coefficients inaccurate for n > ~15
    
    memset( P, 0.0, sizeof(P) );
    P0[0] = 1.0;
    P1[0] = 1.0;
    P1[1] = 0.0;
    
    for (i=1;i<=n-1;i++) {
        for (j=0;j<i+1;j++) {
            P[j] = (2*i+1) * P1[j] / (i+1);
        }
        k = 0;
        for (j=2;j<i+2;j++) {
            P[j] = P[j] - i*P0[k] / (i+1);
            k++;
        }
        for(j=0;j<i+1;j++) {
            P0[j] = P1[j];
        }
        for(j=0;j<i+2;j++) {
            P1[j] = P[j];
        }
    }
    
    // Odd n implicates zero as one of the roots...
    np = n - n % 2;
    
    // Variable substitution: y=x^2
    np = np / 2;
    for(i=0;i<np+1;i++) {
        P[i] = P[2*i];
    }
    
    // Solve the roots of the polynomial by forming a matrix whose
    // characteristic is the n:th Legendre polynomial and solving 
    // for the eigenvalues
    memset( *A, 0.0, (n/2*n/2)*sizeof(double) );
    for(i=0;i<np-1;i++) {
        A[i][i+1] = 1.0;
    }
    for(i=0;i<np;i++) {
        A[np-1][i] = -P[(np-1)+2-(i+1)] / P[0];
    }

    ch1 = 'N';
    ch2 = 'N';
    arg1 = n/2;
    arg2 = 1;
    arg3 = 1;
    arg4 = 8*n;
    
    // Transform the matrix for LAPACK, column-major order
    AT = doublevec(0, (n/2*n/2)-1);
    for (i=0; i<n/2; i++) {
        for (j=0; j<n/2; j++) {
            AT[j+(n/2)*i] = A[j][i];
        }
    }
    dgeev_(&ch1, &ch2, &np, AT, &arg1, Points, P0, Work, &arg2, Work, &arg3, Work, &arg4, &info);
    
    // Backsubstitute from y=x^2
    for(i=0;i<np+1;i++) {
        Q[i] = P[i];
    }
    
    memset( P, 0.0, sizeof(P) );
    for(i=0;i<np+1;i++) {
        P[2*i] = Q[i];
    }
    
    for(i=0;i<np;i++) {
        Q[i] = Points[i];
    }
    for(i=0;i<np;i++) {
        Points[2*(i+1)-2] = +sqrt(Q[i]);
        Points[2*(i+1)-1] = -sqrt(Q[i]);
    }
    if ((n % 2) == 1) Points[n-1] = 0.0;
    
    DerivPoly(n,Q,P);
    RefineRoots(n,P,Q,Points);
    
    // Check for roots
    for(i=0;i<n;i++) {
        
        s = EvalPoly(n,P,Points[i]);
        if( fabs(s) > 1.0e-12 ) {
            warnfunct("GaussQuadraturePoints1D", "-------------------------------------------------------");
            warnfunct("GaussQuadraturePoints1D", "Computed integration point");
            warnfunct("GaussQuadraturePoints1D", "seems to be inaccurate");
            warnfunct("GaussQuadraturePoints1D", "Points req.:", n);
            warnfunct("GaussQuadraturePoints1D", "Residual:");
            printf("%f\n", s);
            warnfunct("GaussQuadraturePoints1D", "Point: +-: ");
            printf("%f\n", sqrt(Points[i]));
             warnfunct("GaussQuadraturePoints1D", "-------------------------------------------------------");
        }
    }
    
    // Finally, the integration weights computed by
    // W_i = 2/( (1-x_i^2)*Q(x_i)^2 ), x_i is the i:th root, and Q(x) = dP(x) / dx
    DerivPoly(n,Q,P);
    
    for(i=0;i<n;i++) {
        
        s = EvalPoly(n-1, Q, Points[i]);
        Weights[i] = 2.0 / ((1.0-pow(Points[i], 2.0))*pow(s,2.0));
    }
    
    
    // Make really sure the weights add up:
    sum = 0.0;
    for(i=0;i<n;i++) {
        sum = sum + Weights[i];
    }
    for(i=0;i<n;i++) {
        Weights[i] = 2.0 * Weights[i] / sum;
    }
    
    //Copy data so that it's available later
    for(i=0;i<n;i++) {
        
        AllPoints[i][n-1] = Points[i];
        AllWeights[i][n-1] = Weights[i];
    }
    
    free_dvector(AT, 0, (n/2*n/2)-1);
}

void GaussQuadratureInit(GaussIntegrationPoints *pt) {

    if (GInit == false) {
        GInit = true;
        for(int n=1;n<=MAXN;n++) {
            GaussQuadraturePoints1D(n);
        }
        // The caller which gets the pointer to the structure of Gauss integration points
        // is responsible for releasing this memory by calling GaussQuadratureDeallocation()
        pt->u = doublevec(0, MAX_INTEGRATION_POINTS-1);
        pt->v = doublevec(0, MAX_INTEGRATION_POINTS-1);
        pt->w = doublevec(0, MAX_INTEGRATION_POINTS-1);
        pt->s = doublevec(0, MAX_INTEGRATION_POINTS-1);
        
        if (pt->u == NULL || pt->v == NULL || pt->w == NULL || pt->s == NULL) {
            errorfunct("GaussQuadratureInit", "Memory allocation error.");
        }
    }
}

void GaussQuadratureDeallocation(GaussIntegrationPoints *pt) {
    if (pt->u != NULL) {
        free_dvector(pt->u, 0, MAX_INTEGRATION_POINTS-1);
        pt->u = NULL;
    }
    
    if (pt->v != NULL) {
        free_dvector(pt->v, 0, MAX_INTEGRATION_POINTS-1);
        pt->v = NULL;
    }
    
    if (pt->w != NULL) {
        free_dvector(pt->w, 0, MAX_INTEGRATION_POINTS-1);
        pt->w = NULL;
    }
    
    if (pt->s != NULL) {
        free_dvector(pt->s, 0, MAX_INTEGRATION_POINTS-1);
        pt->s = NULL;
    }
    free(pt);
    GInit = false;
}

void GaussQuadrature0D(int n, GaussIntegrationPoints *pt) {
    
    GaussQuadratureInit(pt);
    
    pt->n = 1;
    pt->u[0] = 0.0;
    pt->v[0] = 0.0;
    pt->w[0] = 0.0;
    pt->s[0] = 1.0;
}

/****************************************************************************
 Description: 
    Return gaussian integration points for 1D line element            
                                                                           
 Arguments:
    int n -> Number of points in the requested rule      
    GaussIntegrationPoints *pt -> integration point structure
****************************************************************************/
void GaussQuadrature1D(int n, GaussIntegrationPoints *pt) {
    
    int i;
    
    GaussQuadratureInit(pt);
    
    if(n < 1 || n > MAXN) {
        pt->n = 0;
        warnfunct("GaussQuadrature1D", "Invalid number of points: ", n);
        return;
    }
    
    pt->n = n;
    for(i=0;i<n;i++) {
        pt->u[i] = AllPoints[i][n-1];
        pt->v[i] = 0.0;
        pt->w[i] = 0.0;
        pt->s[i] = AllWeights[i][n-1];
    }
}

/****************************************************************************
 Description: 
    Return gaussian integration points for 2D triangle element         
                                                                            
 Arguments:
    int n -> Number of points in the requested rule      
    GaussIntegrationPoints *pt -> integration point structure
****************************************************************************/
void GaussQuadratureTriangle(int n, GaussIntegrationPoints *pt) {
    
    int i;
    double buffer;
    
    GaussQuadratureInit(pt);

    switch(n) {
        case 1:
            pt->u[0] = UTriangle1P[0];
            pt->v[0] = VTriangle1P[0];
            pt->s[0] = STriangle1P[0] / 2.0;
            pt->n = 1;
            break;
        case 3:
            for(i=0;i<n;i++) {
                pt->u[i] = UTriangle3P[i];
                pt->v[i] = VTriangle3P[i];
                pt->s[i] = STriangle3P[i] / 2.0;
            }
            pt->n = 3;
            break;
        case 4:
            for(i=0;i<n;i++) {
                pt->u[i] = UTriangle4P[i];
                pt->v[i] = VTriangle4P[i];
                pt->s[i] = STriangle4P[i] / 2.0;
            }
            pt->n = 4;
            break;
        case 6:
            for(i=0;i<n;i++) {
                pt->u[i] = UTriangle6P[i];
                pt->v[i] = VTriangle6P[i];
                pt->s[i] = STriangle6P[i] / 2.0;
            }
            pt->n = 6;
            break;
        case 7:
            for(i=0;i<n;i++) {
                pt->u[i] = UTriangle7P[i];
                pt->v[i] = VTriangle7P[i];
                pt->s[i] = STriangle7P[i] / 2.0;
            }
            pt->n = 7;
            break;
        case 11:
            for(i=0;i<n;i++) {
                pt->u[i] = Utriangle11P[i];
                pt->v[i] = VTriangle11P[i];
                pt->s[i] = STriangle11P[i];
            }
            pt->n = 11;
            break;
        case 12:
            for(i=0;i<n;i++) {
                pt->u[i] = UTriangle12P[i];
                pt->v[i] = VTriangle12P[i];
                pt->s[i] = STriangle12P[i];
            }
            pt->n = 12;
            break;
        case 17:
            for(i=0;i<n;i++) {
                pt->u[i] = UTriangle17P[i];
                pt->v[i] = VTriangle17P[i];
                pt->s[i] = STriangle17P[i];
            }
            pt->n = 17;
            break;
        case 20:
            for(i=0;i<n;i++) {
                pt->u[i] = UTriangle20P[i];
                pt->v[i] = VTriangle20P[i];
                pt->s[i] = STriangle20P[i];
            }
            pt->n = 20;
            break;
        default:
          
            GaussQuadratureQuad(n, pt);
            for(i=0;i<pt->n;i++) {
                pt->v[i] = (pt->v[i] + 1.0)/2.0;
                pt->u[i] = (pt->u[i] + 1.0) / 2.0 * (1.0-pt->v[i]);
                pt->s[i] = pt->s[i] * ( 1.0 - pt->v[i] );
            }
            buffer = 0;
            for(i=0;i<pt->n;i++) {
                buffer = buffer + pt->s[i];
            }
            for(i=0;i<pt->n;i++) {
                pt->s[i] = 0.5 * pt->s[i] / buffer;
            }
            break;
    }
    
    for(i=0;i<n;i++) {
        pt->w[i] = 0.0;
    }
}

/****************************************************************************
 Description: 
    Return gaussian integration points for 3D tetra element            
                                                                            
 Arguments:
    int n -> Number of points in the requested rule      
    GaussIntegrationPoints *pt -> integration point structure
****************************************************************************/
void GaussQuadratureTetra(int n, GaussIntegrationPoints *pt) {
    
    int i;
    double ScaleFactor;
    
    GaussQuadratureInit(pt);
    
    switch(n) {
        case 1:
            pt->u[0] = UTetra1P[0];
            pt->v[0] = Vtetra1P[0];
            pt->w[0] = WTetra1P[0];
            pt->s[0] = STetra1P[0] / 6.0;
            pt->n = 1;
            break;
        case 4:
            for(i=0;i<n;i++) {
                pt->u[i] = UTetra4P[i];
                pt->v[i] = VTetra4P[i];
                pt->w[i] = WTetra4P[i];
                pt->s[i] = STetra4P[i] / 6.0;
            }
            pt->n = 4;
            break;
        case 5:
            for(i=0;i<n;i++) {
                pt->u[i] = UTetra5P[i];
                pt->v[i] = VTetra5P[i];
                pt->w[i] = WTetra5P[i];
                pt->s[i] = STetra5P[i] / 6.0;
            }
            pt->n = 5;
            break;
        case 11:
            for(i=0;i<n;i++) {
                pt->u[i] = UTetra11P[i];
                pt->v[i] = VTetra11P[i];
                pt->w[i] = WTetra11P[i];
                pt->s[i] = STetra11P[i] / 6.0;
            }
            pt->n = 11;
            break;
        default:
            
            GaussQuadratureBrick(n, pt);
            for(i=0;i<pt->n;i++) {
                
                ScaleFactor = 0.5;
                pt->u[i] = ( pt->u[i] + 1.0 ) * ScaleFactor;
                pt->v[i] = ( pt->v[i] + 1.0 ) * ScaleFactor;
                pt->w[i] = ( pt->w[i] + 1.0 ) * ScaleFactor;
                pt->s[i] = pt->s[i] * pow(ScaleFactor, 3.0);
                
                ScaleFactor = 1.0 - pt->w[i];
                pt->u[i] = pt->u[i] * ScaleFactor;
                pt->v[i] = pt->v[i] * ScaleFactor;
                pt->s[i] = pt->s[i] * pow(ScaleFactor, 2.0);
                
                ScaleFactor = 1.0 - pt->v[i] / ScaleFactor;
                pt->u[i] = pt->u[i] * ScaleFactor;
                pt->s[i] = pt->s[i] * ScaleFactor;
            }
            break;
    }
}

/****************************************************************************
 Description: 
    Return gaussian integration points for 3D prism element            
                                                                            
 Arguments:
    int n -> Number of points in the requested rule      
    GaussIntegrationPoints *pt -> integration point structure
****************************************************************************/
void GaussQuadraturePyramid(int np, GaussIntegrationPoints *pt) {
    
    int i, j, k, n, t;
    
    GaussQuadratureInit(pt);
    
    n = pow((double)np, 1.0/3.0) + 0.5;
    
    if( n < 1 || n > MAXN ) {
        pt->n = 0;
        warnfunct("GaussQuadraturePyramid", "Invalid number of points: ", n);
        return;
    }
    
    t = 0;
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            for(k=0;k<n;k++) {
                pt->u[t] = AllPoints[k][n-1];
                pt->v[t] = AllPoints[j][n-1];
                pt->w[t] = AllPoints[i][n-1];
                pt->s[t] = AllWeights[i][n-1] * AllWeights[j][n-1] * AllWeights[k][n-1];
                t++;
            }
        }
    }
    pt->n = t;
    
    for(t=0;t<pt->n;t++) {
        pt->w[t] = (pt->w[t] + 1.0) / 2.0;
        pt->u[t] = pt->u[t] * (1.0-pt->w[t]);
        pt->v[t] = pt->v[t] * (1.0-pt->w[t]);
        pt->s[t] = pt->s[t] * pow((1.0-pt->w[t]), 2.0) / 2.0;
    }
}

void GaussQuadratureWedge(int np, GaussIntegrationPoints *pt) {
    
    int i, j, k, n, t;
    
    GaussQuadratureInit(pt);
    
    n = pow( (double)np, 1.0/3.0) + 0.5;
    
    if ( n < 1 || n > MAXN) {
        pt->n = 0;
        warnfunct("GaussQuadratureWedge", "Invalid number of points: ", n);
        return;
    }
    
    t = 0;
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            for(k=0;k<n;k++) {
                pt->u[t] = AllPoints[k][n-1];
                pt->v[t] = AllPoints[j][n-1];
                pt->w[t] = AllPoints[i][n-1];
                pt->s[t] = AllWeights[i][n-1] * AllWeights[j][n-1] * AllWeights[k][n-1];
                t++;
            }
        }
    }
    pt->n = t;
    
    for(i=0;i<pt->n;i++) {
        pt->v[i] = ( pt->v[i] + 1.0 )/2.0;
        pt->u[i] = ( pt->u[i] + 1.0 )/2.0 * (1.0-pt->v[i]);
        pt->s[i] = pt->s[i] * (1.0-pt->v[i])/4.0;
    }
}

void GaussQuadratureQuad(int np, GaussIntegrationPoints *pt) {
    
    int i, j, n, t;
    
    GaussQuadratureInit(pt);
    
    n = sqrt( (double)np ) + 0.5;
    
    if( n < 1 || n > MAXN ) {
        pt->n = 0.0;
        warnfunct("GaussQuadratureQuad", "Invalid number of points: ", n);
        return;
    }
    
    t = 0;
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            pt->u[t] = AllPoints[j][n-1];
            pt->v[t] = AllPoints[i][n-1];
            pt->s[t] = AllWeights[i][n-1] * AllWeights[j][n-1];
            t++;
        }
    }
    pt->n = t;
}

void GaussQuadratureBrick(int np, GaussIntegrationPoints *pt) {
    
    int i, j, k, n, t;
    
    GaussQuadratureInit(pt);
    
    n = pow( (double)np, 1.0/3.0) + 0.5;
    
    if ( n < 1 || n > MAXN) {
        pt->n = 0;
        warnfunct("GaussQuadratureBrick", "Invalid number of points: ", n);
        return;
    }
    
    t = 0;
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            for(k=0;k<n;k++) {
                pt->u[t] = AllPoints[k][n-1];
                pt->v[t] = AllPoints[j][n-1];
                pt->w[t] = AllPoints[i][n-1];
                pt->s[t] = AllWeights[i][n-1] * AllWeights[j][n-1] * AllWeights[k][n-1];
                t++;
            }
        }
    }
    pt->n = t;
}

GaussIntegrationPoints* GaussQuadrature(Element_t *element, int *np, int *relOrder) {
    
    int n, eldim, p1d;
    GaussIntegrationPoints *IntegStuff;
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
                    errorfunct("GaussIntegrationPoints", "Number of integration points must remain positive!");
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
                warnfunct("GaussQuadrature", "RelOrder can only be {-1, 0, 1}!");
            }
        }
    } else {
        if (pElement == true) {
            n = element->Pdefs->GaussPoints;
        } else {
            n = element->Type.GaussPoints;
        }
    }
    
    IntegStuff = (GaussIntegrationPoints*)malloc(sizeof(GaussIntegrationPoints));
    IntegStuff->u = NULL;
    IntegStuff->v = NULL;
    IntegStuff->w = NULL;
    IntegStuff->s = NULL;
    
    switch (element->Type.ElementCode / 100) {
            
        case 1:
            GaussQuadrature0D(n, IntegStuff);
            break;
            
        case 2:
            GaussQuadrature1D(n, IntegStuff);
            break;
            
        case 3:
            // TODO: add support for p element
            GaussQuadratureTriangle(n, IntegStuff);
            break;
          
        case 4:
            GaussQuadratureQuad(n, IntegStuff);
            break;
            
        case 5:
             // TODO: add support for p element
            GaussQuadratureTetra(n, IntegStuff);
            break;
            
        case 6:
             // TODO: add support for p element
            GaussQuadraturePyramid(n, IntegStuff);
            break;
            
        case 7:
             // TODO: add support for p element
            GaussQuadratureWedge(n, IntegStuff);
            break;
            
        case 8:
            GaussQuadratureBrick(n, IntegStuff);
            break;
            
    }
    
    return IntegStuff;
}