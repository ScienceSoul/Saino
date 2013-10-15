#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// No support for convection yet
kernel void heatSolutionAssembly(global int *diag, global int *rows, global int *cols, global double *values, global double *rhs, global double *nodesInfo, global int *elementPermutationStore, global double *elementBasis, global double *elementdBasis, int dimension, int numberOfNodes, int numberElementDofs, int numberOfIntegrationPoints, int nBasis, int kernelDof, int varDofs) {
    
    int globalID = get_global_id(0);
    int i, j, k, l, c, col, p, q, row, t;
    double a, c0, c1, ct, force, load, m, rho, sum;
    double c2[3][3];
    double massMatrix[nBasis][nBasis], stiffMatrix[nBasis][nBasis], forceVector[nBasis];
    
    for (i=0; i<nBasis; i++) {
        for (j=0; j<nBasis; j++) {
            massMatrix[i][j] = 0.0;
            stiffMatrix[i][j] = 0.0;
        }
        forceVector[i] = 0.0;
    }
    
    for (t=0; t<numberOfIntegrationPoints; t++) {
        
        // Coefficient of the convection and time derivative terms at the integration point
        c0 = 0.0;
        c1 = 0.0;
        ct = 0.0;
        for (i=0; i<numberOfNodes; i++) {
            j =  elementPermutationStore[(globalID*numberElementDofs)+i];
            c0 = c0 + nodesInfo[kernelDof*j+25]*elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+i];
            c1 = c1 + nodesInfo[kernelDof*j+21]*elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+i];
            ct = ct + nodesInfo[kernelDof*j]*elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+i];
        }
        
        // Coefficient of the diffusion term & its derivatives at the integration point
        rho = 0.0;
        for (i=0; i<numberOfNodes; i++) {
            j = elementPermutationStore[(globalID*numberElementDofs)+i];
            rho = rho + nodesInfo[kernelDof*j+11]*elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+i];
        }
        
        q = 1;
        for (i=0; i<dimension; i++) {
            for (j=0; j<dimension; j++) {
                sum = 0.0;
                for (k=0; k<numberOfNodes; k++) {
                    p = elementPermutationStore[(globalID*numberElementDofs)+k];
                    sum = sum + nodesInfo[kernelDof*p+q]*elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+k];
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
                m = ct * elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+q] * elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+p];
                a = c0 * elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+q] * elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+p];
                
                // The diffusion term
                for (i=0; i<dimension; i++) {
                    for (j=0; j<dimension; j++) {
                        a = a + c2[i][j] * elementdBasis[(globalID*((numberOfNodes*3)*numberOfNodes))+(t*(numberOfNodes*3))+(q*3)+i] * elementdBasis[(globalID*((numberOfNodes*3)*numberOfNodes))+(t*(numberOfNodes*3))+(p*3)+j];
                    }
                }
                stiffMatrix[p][q] = stiffMatrix[p][q] + elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+3] * a;
                massMatrix[p][q] = massMatrix[p][q] + elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+3] * m;
            }
        }
        
        force = 0.0;
        for (i=0; i<numberOfNodes; i++) {
            j = elementPermutationStore[(globalID*numberElementDofs)+i];
            force = force + nodesInfo[kernelDof*j+13]*elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+i];
        }
        
        for (p=0; p<nBasis; p++) {
            load = force * elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+p];
            forceVector[p] = forceVector[p] + elementBasis[(globalID*((numberOfNodes+1)*numberOfNodes))+(t*(numberOfNodes+1))+3] * load;
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