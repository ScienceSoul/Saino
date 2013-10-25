#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// Implementation for three nodes elements
struct BasisFunctionsTiangle {
    
    int p[3], q[3], r[3];
    double coeff[3];
    
};

/***************************************************************************/
/* Triangle - 3 point rule; exact integration of x^py^q, p+q<=2            */
/***************************************************************************/
constant double UTriangle3P[3] = {
    0.16666666666666667E+00, 0.66666666666666667E+00, 0.16666666666666667E+00
};

constant double VTriangle3P[3] = {
    0.16666666666666667E+00, 0.16666666666666667E+00, 0.66666666666666667E+00
};

constant double STriangle3P[3] = {
    0.33333333333333333E+00, 0.33333333333333333E+00, 0.33333333333333333E+00
};

#define CL_DBL_MIN          2.225073858507201383090e-308

void computeBasisFunctions(struct BasisFunctionsTiangle basisFunctions[]) {
    
    int i, j, k, l, basisFunctionDegree = 0, pivot[3], upow, vpow, wpow;
    int maxDeg, maxDeg3, maxDeg2;
    float basisTerms[3], nodeU[3], nodeV[3];
    float u, v, w;
    float a[3][3], s, swap;
    
    maxDeg = 4;
    maxDeg3 = pow((float)maxDeg, 3);
    maxDeg2 = pow((float)maxDeg, 2);
    
    basisTerms[0] = 1.0;
    basisTerms[1] = 2.0;
    basisTerms[2] = 5.0;
    
    nodeU[0] = 0.0; nodeU[1] = 1.0; nodeU[2] = 0.0;
    nodeV[0] = 0.0; nodeV[1] = 0.0; nodeV[2] = 1.0;
    
    for (i=0; i<3; i++) {
        u = nodeU[i];
        v = nodeV[i];
        for (j=0; j<3; j++) {
            k = basisTerms[j] - 1;
            vpow = k / maxDeg;
            upow = k & (maxDeg - 1);
            
            if (upow == 0) {
                a[i][j] = 1.0;
            } else {
                a[i][j] = pow(u, upow);
            }
            
            if (vpow != 0) {
                a[i][j] = a[i][j] * pow(v, vpow);
            }
            
            basisFunctionDegree = max(basisFunctionDegree, upow);
            basisFunctionDegree = max(basisFunctionDegree, vpow);
        }
    }
    
    // Invert matrix
    
    //LU decomposition
    
    for (i=0; i<3; i++) {
        j = i;
        for (k=i+1; k<3; k++) {
            if (fabs(a[i][k]) > fabs(a[i][j])) j = k;
        }
        
        if (fabs(a[i][j]) == 0.0) {
            return;
        }
        
        pivot[i] = j;
        
        if (j != i) {
            for (k=0; k<=i; k++) {
                swap = a[k][j];
                a[k][j] = a[k][i];
                a[k][i] = swap;
            }
        }
        
        for (k=i+1; k<3; k++) {
            a[i][k] = a[i][k] / a[i][i];
        }
        
        for (k=i+1; k<3; k++) {
            if (j != i) {
                swap = a[k][i];
                a[k][i] = a[k][j];
                a[k][j] = swap;
            }
            
            for (l=i+1; l<3; l++) {
                a[k][l] = a[k][l] - a[k][i] * a[i][l];
            }
        }
    }
    pivot[3-1] = 3-1;
    if (fabs(a[3-1][3-1]) == 0.0) {
        // matrix is (at least almost) singular
        // If we get here that's not good, but we really shoudn't as we are not taking any question
    }
    
    for (i=0; i<3; i++) {
        if (fabs(a[i][i]) == 0.0) {
            return;
        }
        a[i][i] = 1.0 / a[i][i];
    }
    
    for (i=3-2; i>=0; i--) {
        for (j=3-1; j>=i+1; j--) {
            s = -a[i][j];
            for (k=i+1; k<=j-1; k++) {
                s = s - a[i][k]*a[k][j];
            }
            a[i][j] = s;
        }
    }
    
    for (i=3-2; i>=0; i--) {
        for (j=3-1; j>=i+1; j--) {
            s = 0.0;
            for (k=i+1; k<=j; k++) {
                s = s - a[j][k]*a[k][i];
            }
            a[j][i] = a[i][i]*s;
        }
    }
    
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            s = 0.0;
            for (k=max(i,j); k<3; k++) {
                if (k != i) {
                    s = s + a[i][k]*a[k][j];
                } else {
                    s = s + a[k][j];
                }
            }
            a[i][j] = s;
        }
    }
    
    for (i=3-1; i>=0; i--) {
        if (pivot[i] != i) {
            for (j=0; j<3; j++) {
                s = a[i][j];
                a[i][j] = a[pivot[i]][j];
                a[pivot[i]][j] = s;
            }
        }
    }
    
    upow = 0;
    vpow = 0;
    wpow = 0;
    
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            k = basisTerms[j] - 1;
            vpow = k / maxDeg;
            upow = k & (maxDeg - 1);
            basisFunctions[i].p[j] = upow;
            basisFunctions[i].q[j] = vpow;
            basisFunctions[i].r[j] = wpow;
            basisFunctions[i].coeff[j] = a[j][i];
        }
    }
}

void computeBasis(double basis[], struct BasisFunctionsTiangle basisFunctions[], double u, double v) {
    
    double s;
    
    // Basis 2D three nodes triangle element only
    for(int n=0; n<3; n++) {
        s = 0.0;
        for(int i=0; i<3; i++) {
            s = s + basisFunctions[n].coeff[i] * pow((float)u, basisFunctions[n].p[i]) * pow((float)v, basisFunctions[n].q[i]);
        }
        basis[n] = s;
    }
}

void NodalFirstDerivatives2D(float y[],struct BasisFunctionsTiangle basisFunctions[], double u, double v) {
    
    int i, n;
    float s, t;
    
    for(n=0;n<3;n++) {
        s = 0.0f;
        t = 0.0f;
        for(i=0;i<3;i++) {
            if(basisFunctions[n].p[i] >= 1) s = s + basisFunctions[n].p[i] * basisFunctions[n].coeff[i] * pow((float)u, (basisFunctions[n].p[i]-1)) * pow((float)v, basisFunctions[n].q[i]);
            if(basisFunctions[n].q[i] >= 1) t = t + basisFunctions[n].q[i] * basisFunctions[n].coeff[i] * pow((float)u, basisFunctions[n].p[i]) * pow((float)v, (basisFunctions[n].q[i]-1));
        }
        y[2*n] = s;
        y[2*n+1] = t;
    }
}

void setDxForElement(float dx[], struct BasisFunctionsTiangle basisFunctions[],  double u, double v, int globalID, int numberElementDofs, int kernelDof, __global int *elementPermutationStore, __global double *nodesInfo) {
    
    int i, j, k, n, dim;
    float accum1, accum2, accum3;
    float dLBasisdx[6];
    
    NodalFirstDerivatives2D(dLBasisdx, basisFunctions, u, v);
    
    for (i=0; i<2; i++) { // Two dimensional element
        accum1 = 0.0f;
        accum2 = 0.0f;
        accum3 = 0.0f;
        for (j=0; j<3; j++) { // Three nodes element
            k =  elementPermutationStore[(globalID*numberElementDofs)+j];
            accum1 = accum1 + (nodesInfo[kernelDof*k+14] * dLBasisdx[2*j+i]);
            accum2 = accum2 + (nodesInfo[kernelDof*k+15] * dLBasisdx[2*j+i]);
            accum3 = accum3 + (nodesInfo[kernelDof*k+16] * dLBasisdx[2*j+i]);
        }
        dx[2*0+i] = accum1;
        dx[2*1+i] = accum2;
        dx[2*2+i] = accum3;
    }
}

// Function returns sqrt(detG) and also computes dBasisdx
double setBasisFirstDerivativeMetricDeterminant (double dBasisdx[], struct BasisFunctionsTiangle basisFunctions[], double u, double v, int globalID, int numberElementDofs, int kernelDof, __global int *elementPermutationStore, __global double *nodesInfo) {
    
    int i, j, k;
    float covariantMetricTensor[2][2], detG, dx[4], elementMetric[2][2], ltoGMap[2][2], s;
    float dLBasisdx[6];
    
    setDxForElement(dx, basisFunctions, u, v, globalID, numberElementDofs, kernelDof, elementPermutationStore, nodesInfo);
    
    for (i=0; i<2; i++) { // Two dimensional element
        for (j=0; j<2; j++) { // Two dimensional element
            s = 0.0f;
            for (k=0; k<2; k++) { // Two dimensional mesh
                s = s + ( dx[2*k+i] * dx[2*k+j]);
            }
            covariantMetricTensor[i][j] = s;
        }
    }
    
    // Surface element
    detG = (covariantMetricTensor[0][0]*covariantMetricTensor[1][1] - covariantMetricTensor[0][1]*covariantMetricTensor[1][0]);
    if (detG <= FLT_MIN) {
        // Don't do anything in the OpenCL kernel but dangerous is something goes wrong
    };
    
    NodalFirstDerivatives2D(dLBasisdx, basisFunctions, u, v);
    
    // Two dimensional element
    elementMetric[0][0] = covariantMetricTensor[1][1] / detG;
    elementMetric[0][1] = -covariantMetricTensor[0][1] / detG;
    elementMetric[1][0] = -covariantMetricTensor[1][0] / detG;
    elementMetric[1][1] = covariantMetricTensor[0][0] / detG;
    
    for (i=0; i<2; i++) { // Two dimensional mesh
        for (j=0; j<2; j++) { // Two dimensional element
            s = 0.0f;
            for (k=0; k<2; k++) { // Two dimensional element
                s = s + dx[2*i+k] * elementMetric[k][j];
            }
            ltoGMap[i][j] = s;
        }
    }
    
    for (i=0; i<6; i++) {
        dBasisdx[i] = 0.0;
    }
    for (i=0; i<3; i++) { // Three nodes element
        for (j=0; j<2; j++) { // Two dimensional mesh
            for (k=0; k<2; k++) { // Two dimensional element
                dBasisdx[2*i+j] = dBasisdx[2*i+j] + dLBasisdx[2*i+k]*ltoGMap[j][k];
            }
        }
    }
    
    return sqrt(detG);
}


__kernel void heatSolutionAssembly(__global int *diag, __global int *rows, __global int *cols, __global double *values, __global double *rhs, __global int *colorMapping, __global double *nodesInfo, __global int *elementPermutationStore, int positionInColorMapping, int dimension, int numberOfNodes, int numberElementDofs, int nBasis, int kernelDof, int varDofs) {
    
    int workItemID = get_global_id(0);
    int i, j, k, l, c, col, n, p, q, row, t, globalID;
    struct BasisFunctionsTiangle basisFunctions[3];
    double basis[3], dBasisdx[6];
    double u[3], v[3], s[3];
    double a, c0, c1, ct, force, load, m, rho, ss, sum;
    double c2[3][3];
    double massMatrix[nBasis][nBasis], stiffMatrix[nBasis][nBasis], forceVector[nBasis];
    float detJ;
    
    // Triangle three gauss points
    for(i=0;i<3;i++) {
        u[i] = UTriangle3P[i];
        v[i] = VTriangle3P[i];
        s[i] = STriangle3P[i] / 2.0;
    }
    n = 3;
    
    for (i=0; i<nBasis; i++) {
        for (j=0; j<nBasis; j++) {
            massMatrix[i][j] = 0.0;
            stiffMatrix[i][j] = 0.0;
        }
        forceVector[i] = 0.0;
    }
    
    globalID = colorMapping[positionInColorMapping+workItemID];
    
    computeBasisFunctions(basisFunctions);
    
    for (t=0; t<n; t++) {
        
        computeBasis(basis, basisFunctions, u[t], v[t]);
        
        detJ = setBasisFirstDerivativeMetricDeterminant(dBasisdx, basisFunctions, u[t], v[t], globalID, numberElementDofs, kernelDof, elementPermutationStore, nodesInfo);
        ss = detJ * s[t];
        
        // Coefficient of the convection and time derivative terms at the integration point
        c0 = 0.0;
        c1 = 0.0;
        ct = 0.0;
        for (i=0; i<numberOfNodes; i++) {
            j =  elementPermutationStore[(globalID*numberElementDofs)+i];
            c0 = c0 + nodesInfo[kernelDof*j+13]*basis[i];
            c1 = c1 + nodesInfo[kernelDof*j+12]*basis[i];
            ct = ct + nodesInfo[kernelDof*j]*basis[i];
        }
        
        // Coefficient of the diffusion term & its derivatives at the integration point
        rho = 0.0;
        for (i=0; i<numberOfNodes; i++) {
            j = elementPermutationStore[(globalID*numberElementDofs)+i];
            rho = rho + nodesInfo[kernelDof*j+10]*basis[i];
        }
        
        q = 1;
        for (i=0; i<dimension; i++) {
            for (j=0; j<dimension; j++) {
                sum = 0.0;
                for (k=0; k<numberOfNodes; k++) {
                    p = elementPermutationStore[(globalID*numberElementDofs)+k];
                    sum = sum + nodesInfo[kernelDof*p+q]*basis[k];
                }
                q++;
                c2[i][j] = sum;
            }
            // Jump further in the list if actual dimension is smaller than maximum dimension
            if (dimension == 2) {
                q = q + 1;
            } else if (dimension == 1) {
                q = q + 2;
            }
        }
        
        // Loop over basis functions of both unknowns and weights
        for (p=0; p<nBasis; p++) {
            for (q=0; q<nBasis; q++) {
                // The diffusive-convective equation without stabilization
                m = ct * basis[q] * basis[p];
                a = c0 * basis[q] * basis[p];
                
                // The diffusion term
                for (i=0; i<dimension; i++) {
                    for (j=0; j<dimension; j++) {
                        a = a + c2[i][j] * dBasisdx[2*q+i] * dBasisdx[2*p+j];
                    }
                }
                stiffMatrix[p][q] = stiffMatrix[p][q] + ss * a;
                massMatrix[p][q] = massMatrix[p][q] + ss * m;
            }
        }
        
        force = 0.0;
        for (i=0; i<numberOfNodes; i++) {
            j = elementPermutationStore[(globalID*numberElementDofs)+i];
            force = force + nodesInfo[kernelDof*j+11]*basis[i];
        }
        
        for (p=0; p<nBasis; p++) {
            load = force * basis[p];
            forceVector[p] = forceVector[p] + ss * load;
        }
    }
    
    if (varDofs == 1) {
        for (i=0; i<numberElementDofs; i++) {
            row = elementPermutationStore[(globalID*numberElementDofs)+i];
            if (row < 0) continue;
            for (j=0; j<numberElementDofs; j++) {
                col = elementPermutationStore[(globalID*numberElementDofs)+j];
                if (col < 0) continue;
                if (col >= row) {
                    for (c=diag[row]; c<=rows[row+1]-1; c++) {
                        if (cols[c] == col) {
                            values[c] = values[c] + stiffMatrix[i][j];
                            break;
                        }
                    }
                } else {
                    for (c=rows[row]; c<=diag[row]-1; c++) {
                        if (cols[c] == col) {
                            values[c] = values[c] + stiffMatrix[i][j];
                            break;
                        }
                    }
                }
            }
        }
    } else {
        for (i=0; i<numberElementDofs; i++) {
            for (k=1; k<=varDofs; k++) {
                if (elementPermutationStore[(globalID*numberElementDofs)+i] < 0) continue;
                row = varDofs * (elementPermutationStore[(globalID*numberElementDofs)+i]+1) - k;
                for (j=0; j<numberElementDofs; j++) {
                    for (l=1; l<=varDofs; l++) {
                        if (elementPermutationStore[(globalID*numberElementDofs)+j] < 0) continue;
                        col = varDofs * (elementPermutationStore[(globalID*numberElementDofs)+j]+1) - l;
                        if (col >= row) {
                            for (c=diag[row]; c<=rows[row+1]-1; c++) {
                                if (cols[c] == col) {
                                    values[c] = values[c] + stiffMatrix[varDofs*(i+1)-k][varDofs*(j+1)-l];
                                    break;
                                }
                            }
                        } else {
                            for (c=rows[row]; c<=diag[row]-1; c++) {
                                if (cols[c] == col) {
                                    values[c] = values[c] + stiffMatrix[varDofs*(i+1)-k][varDofs*(j+1)-l];
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    for (i=0; i<numberElementDofs; i++) {
        if (elementPermutationStore[(globalID*numberElementDofs)+i] >= 0) {
            for (j=0; j<varDofs; j++) {
                k = varDofs * elementPermutationStore[(globalID*numberElementDofs)+i] + j;
                rhs[k] = rhs[k] + forceVector[varDofs*i+j];
            }
        }
    }
}