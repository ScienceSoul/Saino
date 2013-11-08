#pragma OPENCL EXTENSION cl_khr_fp64 : enable


constant float STriangle3P[3] = {
    0.333333333f, 0.333333333f, 0.333333333f
    //0.33333333333333333E+00, 0.33333333333333333E+00, 0.33333333333333333E+00
};

__kernel void heatSolutionAssembly(__global double *values, __global double *rhs, __global int *diag, __global int *rows, __global int *cols, __global int *colorMapping, __global int *elementNodeIndexesStore, __global int *permutation, __global double *nodesX, __global double *nodesY, __global double *nodesZ, int positionInColorMapping, int dimension, int numberOfNodes, int numberElementDofs, int nBasis, int varDofs) {
    
    int workItemID = get_global_id(0);
    int i, j, k, l, c, col, n, p, q, row, t, globalID, pivot[3], upow, vpow, wpow, temp;
    int basisTerms[3];
    float basis[3];
    double dBasisdx[3][2], stiffMatrix[3][3], aa;
    float c2[3][3];
    float c0, c1, ct, m, s, force, load, rho, ss, sum, tt, massMatrix[3][3], forceVector[3];
    float accum1, accum2, accum3;
    float detJ;
    float covariantMetricTensor[2][2], dx[3][2], elementMetric[2][2], ltoGMap[2][2];
    float dLBasisdx[3][2];
    // Basis functions implementation for three nodes elements
    int basisFunctionsP[3][3];
    int basisFunctionsQ[3][3];
    int basisFunctionsR[3][3];
    float basisFunctionsCoeff[3][3];
    float maxDeg, maxDeg3, maxDeg2;
    float nodeU[3], nodeV[3];
    float u, v;
    float a[3][3], swap;
    /***************************************************************************/
    /* Triangle - 3 point rule; exact integration of x^py^q, p+q<=2            */
    /***************************************************************************/
    float UTriangle3P[3] = {
        0.166666666f, 0.666666666f, 0.166666666f
        //0.16666666666666667E+00, 0.66666666666666667E+00, 0.16666666666666667E+00
    };
    
    float VTriangle3P[3] = {
        0.166666666f, 0.166666666f, 0.666666666f
        //0.16666666666666667E+00, 0.16666666666666667E+00, 0.66666666666666667E+00
    };
    
    n = 3;
    
    for (i=0; i<nBasis; i++) {
        for (j=0; j<nBasis; j++) {
            massMatrix[i][j] = 0.0f;
            stiffMatrix[i][j] = 0.0f;
        }
        forceVector[i] = 0.0f;
    }
    
    globalID = colorMapping[positionInColorMapping+workItemID];
    
    maxDeg = 4.0f;
    temp = 4;
    maxDeg3 = pow(maxDeg, 3);
    maxDeg2 = pow(maxDeg, 2);
    
    basisTerms[0] = 1;
    basisTerms[1] = 2;
    basisTerms[2] = 5;
    
    nodeU[0] = 0.0f; nodeU[1] = 1.0f; nodeU[2] = 0.0f;
    nodeV[0] = 0.0f; nodeV[1] = 0.0f; nodeV[2] = 1.0f;
    
    
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            basisFunctionsP[i][j] = 0;
            basisFunctionsQ[i][j] = 0;
            basisFunctionsR[i][j] = 0;
            basisFunctionsCoeff[i][j] = 0.0f;
        }
    }
    
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            a[i][j] = 0.0f;
        }
    }
    
    for (i=0; i<3; i++) {
        u = nodeU[i];
        v = nodeV[i];
        for (j=0; j<3; j++) {
            k = basisTerms[j] - 1;
            vpow = k / temp;
            upow = k & (temp - 1);
            
            if (upow == 0) {
                a[i][j] = 1.0f;
            } else {
                a[i][j] = pow(u, upow);
            }
            
            if (vpow != 0) {
                a[i][j] = a[i][j] * pow(v, vpow);
            }
        }
    }
    
    // Invert matrix
    
    //LU decomposition
    
    for (i=0; i<3; i++) {
        j = i;
        for (k=i+1; k<3; k++) {
            if (fabs(a[i][k]) > fabs(a[i][j])) j = k;
        }
        
        // matrix is (at least almost) singular
        // If we get here that's not good, but we really shoudn't as we are not taking any question
        //if (fabs(a[i][j]) == 0.0f) {
        //}
        
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
    // matrix is (at least almost) singular
    // If we get here that's not good, but we really shoudn't as we are not taking any question
    //if (fabs(a[3-1][3-1]) == 0.0f) {
    //}
    
    for (i=0; i<3; i++) {
        // matrix is (at least almost) singular
        // If we get here that's not good, but we really shoudn't as we are not taking any question
        //if (fabs(a[i][i]) == 0.0f) {
        //}
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
            vpow = k / temp;
            upow = k & (temp - 1);
            basisFunctionsP[i][j] = upow;
            basisFunctionsQ[i][j] = vpow;
            basisFunctionsR[i][j] = wpow;
            basisFunctionsCoeff[i][j] = a[j][i];
        }
    }
    
    for (t=0; t<n; t++) {
        
        for(i=0; i<3; i++) {
            basis[i] = 0.0f;
        }
        
        // Basis 2D three nodes triangle element only
        for(n=0; n<3; n++) {
            s = 0.0f;
            for(i=0; i<3; i++) {
                s = s + basisFunctionsCoeff[n][i] * pow(UTriangle3P[t], basisFunctionsP[n][i]) * pow(VTriangle3P[t], basisFunctionsQ[n][i]);
            }
            basis[n] = s;
        }
        
        // dLBasisDx
        for(n=0;n<3;n++) {
            s = 0.0;
            tt = 0.0;
            for(i=0;i<3;i++) {
                if(basisFunctionsP[n][i] >= 1) s = s + basisFunctionsP[n][i] * basisFunctionsCoeff[n][i] * pow(UTriangle3P[t], (basisFunctionsP[n][i]-1)) * pow(VTriangle3P[t], basisFunctionsQ[n][i]);
                if(basisFunctionsQ[n][i] >= 1) tt = tt + basisFunctionsQ[n][i] * basisFunctionsCoeff[n][i] * pow(UTriangle3P[t], basisFunctionsP[n][i]) * pow(VTriangle3P[t], (basisFunctionsQ[n][i]-1));
            }
            dLBasisdx[n][0] = s;
            dLBasisdx[n][1] = tt;
        }
        
        for (i=0; i<2; i++) { // Two dimensional element
            accum1 = 0.0;
            accum2 = 0.0;
            accum3 = 0.0;
            for (j=0; j<3; j++) { // Three nodes element
                k =  permutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]];
                accum1 = accum1 + (nodesX[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] * dLBasisdx[j][i]);
                accum2 = accum2 + (nodesY[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] * dLBasisdx[j][i]);
                accum3 = accum3 + (nodesZ[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] * dLBasisdx[j][i]);
            }
            dx[0][i] = accum1;
            dx[1][i] = accum2;
            dx[2][i] = accum3;
        }
        
        for (i=0; i<2; i++) { // Two dimensional element
            for (j=0; j<2; j++) { // Two dimensional element
                ss = 0.0;
                for (k=0; k<2; k++) { // Two dimensional mesh
                    ss = ss + ( dx[k][i] * dx[k][j]);
                }
                covariantMetricTensor[i][j] = ss;
            }
        }
        
        // Surface element
        detJ = (covariantMetricTensor[0][0]*covariantMetricTensor[1][1] - covariantMetricTensor[0][1]*covariantMetricTensor[1][0]);
        if (detJ <= FLT_MIN) {
            // Don't do anything in the OpenCL kernel but dangerous is something goes wrong
        };
        
        // Two dimensional element
        elementMetric[0][0] = covariantMetricTensor[1][1] / detJ;
        elementMetric[0][1] = -covariantMetricTensor[0][1] / detJ;
        elementMetric[1][0] = -covariantMetricTensor[1][0] / detJ;
        elementMetric[1][1] = covariantMetricTensor[0][0] / detJ;
        
        for (i=0; i<2; i++) { // Two dimensional mesh
            for (j=0; j<2; j++) { // Two dimensional element
                ss = 0.0;
                for (k=0; k<2; k++) { // Two dimensional element
                    ss = ss + dx[i][k] * elementMetric[k][j];
                }
                ltoGMap[i][j] = ss;
            }
        }
        
        for (i=0; i<3; i++) {
            for (j=0; j<2; j++) {
                dBasisdx[i][j] = 0.0;
            }
        }
        
        for (i=0; i<3; i++) { // Three nodes element
            for (j=0; j<2; j++) { // Two dimensional mesh
                for (k=0; k<2; k++) { // Two dimensional element
                    dBasisdx[i][j] = dBasisdx[i][j] + dLBasisdx[i][k]*ltoGMap[j][k];
                }
            }
        }

        detJ = sqrt(detJ);
        s = detJ * STriangle3P[t] / 2.0;
        
        // Coefficient of the convection and time derivative terms at the integration point
        c0 = 0.0;
        c1 = 0.0;
        ct = 0.0;
        for (i=0; i<numberOfNodes; i++) {
            c0 = c0 + 0.0f*basis[i];
            c1 = c1 + 0.0f*basis[i];
            ct = ct + 0.0*basis[i];
        }
        
        // Coefficient of the diffusion term & its derivatives at the integration point
        rho = 0.0;
        for (i=0; i<numberOfNodes; i++) {
            rho = rho + 1.0*basis[i];
        }
        
        for (i=0; i<dimension; i++) {
            for (j=0; j<dimension; j++) {
                sum = 0.0;
                if ((i == 0 && j == 0) || (i == 1 && j == 1)) {
                    for (k=0; k<numberOfNodes; k++) {
                        sum = sum + 1.0*basis[k];
                    }
                }
                c2[i][j] = sum;
            }
        }
        
        // Loop over basis functions of both unknowns and weights
        for (p=0; p<nBasis; p++) {
            for (q=0; q<nBasis; q++) {
                // The diffusive-convective equation without stabilization
                m = ct * basis[q] * basis[p];
                aa = c0 * basis[q] * basis[p];
                
                // The diffusion term
                for (i=0; i<dimension; i++) {
                    for (j=0; j<dimension; j++) {
                        aa = aa + c2[i][j] * dBasisdx[q][i] * dBasisdx[p][j];
                    }
                }
                stiffMatrix[p][q] = stiffMatrix[p][q] + s * aa;
                massMatrix[p][q] = massMatrix[p][q] + s * m;
            }
        }
        
        force = 0.0;
        for (i=0; i<numberOfNodes; i++) {
            force = force + 1.0*basis[i];
        }
        
        for (p=0; p<nBasis; p++) {
            load = force * basis[p];
            forceVector[p] = forceVector[p] + s * load;
        }
    }

    if (varDofs == 1) {
        for (i=0; i<numberElementDofs; i++) {
            row = permutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]];
            if (row < 0) continue;
            for (j=0; j<numberElementDofs; j++) {
                col = permutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]];
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
                if (permutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]] < 0) continue;
                row = varDofs * (permutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]]+1) - k;
                for (j=0; j<numberElementDofs; j++) {
                    for (l=1; l<=varDofs; l++) {
                        if (permutation[elementNodeIndexesStore[(globalID*numberElementDofs)]+j] < 0) continue;
                        col = varDofs * (permutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]+1) - l;
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
        if (permutation[elementNodeIndexesStore[(globalID*numberElementDofs)]+i] >= 0) {
            for (j=0; j<varDofs; j++) {
                k = varDofs * permutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]] + j;
                rhs[k] = rhs[k] + forceVector[varDofs*i+j];
            }
        }
    }
}