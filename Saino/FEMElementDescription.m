//
//  FEMElementDescription.m
//  Saino
//
//  Created by Hakime Seddik on 05/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMElementDescription.h"

#import "Numerics.h"
#import "Utils.h"
#import "GaussIntegration.h"

static double AEPS = 10.0 * DBL_EPSILON;

@implementation FEMElementDescription

- (id)init
{
    int i;
    
    self = [super init];
    if (self) {
        // Initialization code here.
       
        point = intmatrix(0, 0, 0, 0);
        line = intmatrix(0, 0, 0, 1);
        triangle = intmatrix(0, 2, 0, 1);
        quad = intmatrix(0, 3, 0, 1);
        tetra = intmatrix(0, 5, 0, 1);
        prism = intmatrix(0, 7, 0, 1);
        wedge = intmatrix(0, 8, 0, 1);
        brick = intmatrix(0, 11, 0, 1);
        
        for (i=0; i<8; i++) {
            initialized[i] = NO;
        }
        
    }
    
    return self;
}

- (void)dealloc
{
    free_imatrix(point, 0, 0, 0, 0);
    free_imatrix(line, 0, 0, 0, 1);
    free_imatrix(triangle, 0, 2, 0, 1);
    free_imatrix(quad, 0, 3, 0, 1);
    free_imatrix(tetra, 0, 5, 0, 1);
    free_imatrix(prism, 0, 7, 0, 1);
    free_imatrix(wedge, 0, 8, 0, 1);
    free_imatrix(brick, 0, 11, 0, 1);
}


-(int **)getEdgeMap:(int)elementFamily {
    
    int **edgeMag;
    
    switch (elementFamily) {
        case 1:
            edgeMag = point;
        case 2:
            edgeMag = line;
            break;
        case 3:
            edgeMag = triangle;
            break;
        case 4:
            edgeMag = quad;
            break;
        case 5:
            edgeMag = tetra;
            break;
        case 6:
            edgeMag = prism;
            break;
        case 7:
            edgeMag = wedge;
            break;
        case 8:
            edgeMag = brick;
            break;
    }
    
    if (initialized[elementFamily-1] == NO) {
        initialized[elementFamily-1] = YES;
        switch (elementFamily) {
            case 1:
                edgeMag[0][0] = 0;
                break;
            case 2:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                break;
            case 3:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 2;
                edgeMag[2][1] = 0;
                break;
            case 4:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 2;
                edgeMag[2][1] = 3;
                
                edgeMag[3][0] = 3;
                edgeMag[3][1] = 0;
                break;
            case 5:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][0] = 2;
                
                edgeMag[2][0] = 2;
                edgeMag[2][1] = 0;
                
                edgeMag[3][0] = 0;
                edgeMag[3][1] = 3;
                
                edgeMag[4][0] = 1;
                edgeMag[4][1] = 3;
                
                edgeMag[5][0] = 2;
                edgeMag[5][1] = 3;
                break;
            case 6:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 3;
                edgeMag[2][1] = 2;
                
                edgeMag[3][0] = 0;
                edgeMag[3][1] = 3;
                
                edgeMag[4][0] = 0;
                edgeMag[4][1] = 4;
                
                edgeMag[5][0] = 1;
                edgeMag[5][1] = 4;
                
                edgeMag[6][0] = 2;
                edgeMag[6][1] = 4;
                
                edgeMag[7][0] = 3;
                edgeMag[7][1] = 4;
                break;
            case 7:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 2;
                edgeMag[2][1] = 0;
                
                edgeMag[3][0] = 3;
                edgeMag[3][1] = 4;
                
                edgeMag[4][0] = 4;
                edgeMag[4][1] = 5;
                
                edgeMag[5][0] = 5;
                edgeMag[5][1] = 3;
                
                edgeMag[6][0] = 0;
                edgeMag[6][1] = 3;
                
                edgeMag[7][0] = 1;
                edgeMag[7][1] = 4;
                
                edgeMag[8][0] = 2;
                edgeMag[8][1] = 5;
                break;
            case 8:
                edgeMag[0][0] = 0;
                edgeMag[0][1] = 1;
                
                edgeMag[1][0] = 1;
                edgeMag[1][1] = 2;
                
                edgeMag[2][0] = 3;
                edgeMag[2][1] = 2;
                
                edgeMag[3][0] = 0;
                edgeMag[3][1] = 3;
                
                edgeMag[4][0] = 4;
                edgeMag[4][1] = 5;
                
                edgeMag[5][0] = 5;
                edgeMag[5][1] = 6;
                
                edgeMag[6][0] = 7;
                edgeMag[6][1] = 6;
                
                edgeMag[7][0] = 4;
                edgeMag[7][1] = 7;
                
                edgeMag[8][0] = 0;
                edgeMag[8][1] = 4;
                
                edgeMag[9][0] = 1;
                edgeMag[9][1] = 5;
                
                edgeMag[10][0] = 2;
                edgeMag[10][1] = 6;
                
                edgeMag[11][0] = 3;
                edgeMag[11][1] = 7;
                break;
        }
    }
    
    return edgeMag;
    
}

-(double)elementDiameter:(Element_t *)element: (Nodes_t *)nodes {
/*****************************************************************************************
    Retrive element diameter parameter for stabilization
 
    Element_t *element  ->  element
    Nodes_t *nodes      ->  nodal coodinate arrays of the element
 
*****************************************************************************************/
    
    int i, j, k, family;
    double x0, y0, z0, hk, a, s, cx, cy, cz;
    double j11, j12, j13, j21, j22, j23, g11, g12, g22;
    int **edgeMap, size;
    
    family = element->Type.ElementCode / 100;
    
    switch (family) {
        case 1:
            hk = 0.0;
            break;
        case 3: // Triangular element
            j11 = nodes[1].x - nodes[0].x;
            j12 = nodes[1].y - nodes[0].y;
            j13 = nodes[1].z - nodes[0].z;
            j21 = nodes[2].x - nodes[0].x;
            j22 = nodes[2].y - nodes[0].y;
            j23 = nodes[2].z - nodes[0].z;
            g11 = pow(j11, 2.0) + pow(j12, 2.0) + pow(j13, 2.0);
            g12 = j11*j21 + j12+j22 + j13*j23;
            g22 = pow(j21, 2.0) + pow(j22, 2.0) + pow(j23, 2.0);
            a = sqrt(g11*g22 - pow(g12, 2.0)) / 2.0;
            
            cx = ( nodes[0].x + nodes[1].x + nodes[2].x ) / 3.0;
            cy = ( nodes[0].y + nodes[1].y + nodes[2].y ) / 3.0;
            cz = ( nodes[0].z + nodes[1].z + nodes[2].z ) / 3.0;
            
            s = pow((nodes[0].x-cx), 2.0) + pow((nodes[0].y-cy), 2.0) + pow((nodes[0].z-cz), 2.0);
            s = s + pow((nodes[1].x-cx), 2.0) + pow((nodes[1].y-cy), 2.0) + pow((nodes[1].z-cz), 2.0);
            s = s + pow((nodes[2].x-cx), 2.0) + pow((nodes[2].y-cy), 2.0) + pow((nodes[2].z-cz), 2.0);
            
            hk = 16.0 * a * a / (3.0 * s);
            break;
        case 4: // Quadrilateral element
            cx = pow((nodes[1].x-nodes[0].x), 2.0) + pow((nodes[1].y-nodes[0].y), 2.0) + pow((nodes[1].z-nodes[0].z), 2.0);
            cy = pow((nodes[3].x-nodes[0].x), 2.0) + pow((nodes[3].y-nodes[0].y), 2.0) + pow((nodes[3].z-nodes[0].z), 2.0);
            hk = 2.0 * cx * cy / (cx + cy);
        default:
            edgeMap = [self getEdgeMap:family];
            hk = HUGE_VAL;
            switch (family) {
                case 2:
                    size = 1;
                    break;
                case 3:
                    size = 3;
                    break;
                case 4:
                    size = 4;
                    break;
                case 5:
                    size = 6;
                    break;
                case 6:
                    size = 8;
                    break;
                case 7:
                    size = 9;
                    break;
                case 8:
                    size = 12;
                    break;
            }
            for (i=0; i<size; i++) {
                j = edgeMap[i][0];
                k = edgeMap[i][1];
                x0 = nodes[j].x - nodes[k].x;
                y0 = nodes[j].y - nodes[k].y;
                z0 = nodes[j].z - nodes[k].z;
                hk = min( hk, (pow(x0, 2.0)+pow(y0, 2.0)+pow(z0, 2.0)) );
            }
            break;
    }
    edgeMap = NULL;
    return sqrt(hk);
}

-(void)computeStabilizationParameter:(Element_t *)element: (Nodes_t *)nodes: (FEMMesh *)mesh: (int)n: (double)mk: (double *)hk {
/*******************************************************************************************************
    Compute convection diffusion equation stabilization parameter for each and every element of the 
    model by solving the largest eigenvalue of
        
        Lu = \lambda Gu,
        L = (\nabla^2 u,\nabla^2 w), G = (\nabla u,\nabla w)
 
*******************************************************************************************************/
 
    int i, j, k, p, q, t, dim;
    double *eigr, s, *ddp, *ddq, ***dNodalBasisdx;
    double u, v, w, **l, **g, *work, sum1, sum2;
    double *buffer1, *buffer2;
    GaussIntegrationPoints integCompound;
    FEMNumericIntegration *numericIntegration;
    BOOL stat;
    
    // Used for lapack routine
    int itype, order, lda, ldb, lwork, info;
    char *jobz, *uplo;
    
    eigr = doublevec(0, n-1);
    ddp = doublevec(0, 2);
    ddq = doublevec(0, 2);
    dNodalBasisdx = d3tensor(0, n-1, 0, n-1, 0, 2);
    l = doublematrix(0, (n-1)-1, 0, (n-1)-1);
    g = doublematrix(0, (n-1)-1, 0, (n-1)-1);
    work = doublevec(0, (16*n)-1);
    
    numericIntegration = [[FEMNumericIntegration alloc] init];
    [numericIntegration allocation:mesh];
    
    if (element->Type.BasisFunctionDegree <= 1) {
        switch (element->Type.ElementCode) {
            case 202:
            case 303:
            case 404:
            case 504:
            case 605:
            case 706:
                mk = 1.0 / 3.0;
                break;
            case 808:
                mk = 1.0 / 6.0;
        }
        if (hk != NULL) *hk = [self elementDiameter:element :nodes];
        return;
    }
    
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            for (k=0; k<3; k++) {
                dNodalBasisdx[i][j][k] = 0.0;
            }
        }
    }
    for (p=0; p<n; p++) {
        u = element->Type.NodeU[p];
        v = element->Type.NodeV[p];
        w = element->Type.NodeW[p];
        stat = [numericIntegration setBasisFirstDerivativeForElement:element 
                                                        elementNodes:nodes 
                                                              inMesh:mesh 
                                                firstEvaluationPoint:u 
                                               secondEvaluationPoint:v 
                                                thirdEvaluationPoint:w 
                                                         withBubbles:NO 
                                                         basisDegree:NULL];
        for (i=0; i<n; i++) {
            for (j=0; j<3; j++) {
                dNodalBasisdx[i][p][j] = [numericIntegration basisFirstDerivative:i :j];
            }
        }
    }
    
    dim = [mesh dimension];
    integCompound = GaussQuadrature(element);
    
    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1; j++) {
            l[i][j] = 0.0;
            g[i][j] = 0.0;
        }
    }
    for (t =0; t<integCompound.n; t++) {
        u = integCompound.u[t];
        v = integCompound.v[t];
        w = integCompound.w[t];
        
        stat = [numericIntegration setMetricDeterminantForElement:element 
                                                     elementNodes:nodes inMesh:mesh 
                                             firstEvaluationPoint:u 
                                            secondEvaluationPoint:v 
                                             thirdEvaluationPoint:w];
        
        s = [numericIntegration metricDeterminant] * integCompound.s[t];
        
        for (p=1; p<n; p++) {
            for (q=1; q<n; q++) {
               memset( ddp, 0.0, (3*sizeof(ddp)) );
               memset( ddq, 0.0, (3*sizeof(ddq)) );
               for (i=0; i<dim; i++) {
                   g[p-1][q-1] = g[p-1][q-1] + s * [numericIntegration basisFirstDerivative:p :i] * [numericIntegration basisFirstDerivative:q :i];
                   sum1 = 0.0;
                   sum2 = 0.0;
                   for (j=0; j<n; j++) {
                       sum1 = sum1 + dNodalBasisdx[p][j][i] * [numericIntegration basisFirstDerivative:j :i];
                       sum2= sum2 + dNodalBasisdx[q][j][i] * [numericIntegration basisFirstDerivative:j :i];
                   }
                   ddp[i] = ddp[i] + sum1;
                   ddq[i] = ddq[i] + sum2;
               }
                sum1 = 0.0;
                sum2 = 0.0;
               for (i=0; i<3; i++) {
                   sum1 = sum1 + ddp[i];
                   sum2 = sum2 + ddq[i];
               }
                l[p-1][q-1] = l[p-1][q-1] + s * sum1 * sum2; 
            }
        }
    }
    stat = YES;
    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1; j++) {
            if (l[i][j] < AEPS) {
                continue;
            } else {
                stat = NO;
                break;
            }
        }
    }
    if (stat == YES) {
        mk = 1.0 / 3.0;
        if (hk != NULL) {
            *hk = [self elementDiameter:element :nodes];
        }
        return;
    }
    
    // To do's: Here we serialize the arrays l[][] and g[][] because we need to do that for using 
    // them in dsygv_(). Later optimize so that we work initially on serialzed arrays so that 
    // we don't need to do a copy.
    
    buffer1 = doublevec(0, ((n-1)*(n-1))-1);
    buffer2 = doublevec(0, ((n-1)*(n-1))-1);
    
    k = 0;
    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1; j++) {
            buffer1[k] = l[i][j];
            k++;
        }
    }
    k=0;
    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1; j++) {
            buffer2[k] = g[i][j];
            k++;
        }
    }
    
    itype = 1;
    jobz = "N";
    uplo = "U";
    order = n-1;
    lda = n-1;
    ldb = n-1;
    lwork = 12*n;
    dsygv_(&itype, jobz, uplo, &order, buffer1, &lda, buffer2, &ldb, eigr, work, &lwork, &info);
    if (info < 0 || info > 0) {
        warnfunct("computeStabilizationParameter", "Error in lapack routine dsygv. Error code:");
        printf("%d\n", info);
        errorfunct("computeStabilizationParameter", "Program terminating now...");
    }
    mk = eigr[n-2];
    
    free_dvector(buffer1, 0, ((n-1)*(n-1))-1);
    free_dvector(buffer2, 0, ((n-1)*(n-1))-1);
    
    if (mk < 10.0 * AEPS) {
        mk = 1.0 / 3.0;
        if (hk != NULL) {
            *hk = [self elementDiameter:element :nodes];
        }
        return;
    }
    
    if (hk != NULL) {
        *hk = sqrt( 2.0 / (mk * element->Type.StabilizationMK) );
        mk = min( 1.0/3.0, element->Type.StabilizationMK );
    } else {
        switch (element->Type.ElementCode / 100) {
            case 2:
            case 4:
            case 8:
                mk = 4.0 * mk;
                break;
        }
        mk = min( 1.0/3.0, 2.0/mk );
    }

    free_dvector(eigr, 0, n-1);
    free_dvector(ddp, 0, 2);
    free_dvector(ddq, 0, 2);
    free_d3tensor(dNodalBasisdx, 0, n-1, 0, n-1, 0, 2);
    free_dmatrix(l, 0, (n-1)-1, 0, (n-1)-1);
    free_dmatrix(g, 0, (n-1)-1, 0, (n-1)-1);
    free_dvector(work, 0, (16*n)-1);
    
    [numericIntegration deallocation:mesh];
    
}

-(ElementType_t *)getElementType:(int)code inMesh:(FEMMesh *)mesh stabilization:(BOOL *)computeStab {
    
    int i;
    
    Element_t *elm;
    ElementType_t *element;
    Nodes_t *nodes;
    
    element = elementTypeList;
    
    while (element != NULL) {
        if (code == element->ElementCode) break;
        element = element->NextElementType;
    }
    
    if (element == NULL) {
        errorfunct("getElementType", "Element type code not found:");
        printf("%d\n", code);
        errorfunct("getElementType", "Ignoring element.");
        return NULL;
    }
    
    if (computeStab != NULL) {
        if (*computeStab == NO) return element;
    }
    
    if (element->StabilizationMK == 0.0) {
        elm = (Element_t*) malloc( sizeof(Element_t));
        elm->Type = *element;
        elm->BDOFs = 0;
        elm->DGDOFs = 0;
        elm->Pdefs = NULL;
        elm->DGIndexes = NULL;
        elm->FaceIndexes = NULL;
        elm->BubbleIndexes = NULL;
        nodes = (Nodes_t*)malloc(sizeof(Nodes_t) * element->NumberOfNodes );
        for (i=0; i<element->NumberOfNodes; i++) {
            nodes[i].x = element->NodeU[i];
            nodes[i].y = element->NodeV[i];
            nodes[i].z = element->NodeW[i];
        }
        [self computeStabilizationParameter:elm :nodes :mesh :element->NumberOfNodes :element->StabilizationMK :NULL];
        free(elm);
        free(nodes);
    }
    
    *element = elm->Type;
    return element;
    
}
    

@end
