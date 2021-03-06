//===----------------------------------------------------------------------===//
//  NavierStokesAssemblyKernel.cl
//  OpenCL kernel source
//
//  Copyright © 2016 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===/

#ifdef KERNEL_FP_32
#define REAL float
#define ZERO 0.0f
#define ONE 1.0f
#define TWO 2.0f
#define FOUR 4.0f
#define EIGHT 8.0f
#define IP_U 0.57735027f
#define IP_V 0.57735027f
#define IP_W 0.57735027f
#define IP_S 1.0f
#define viscoExponent (1.0f/3.0f)
#endif

#ifdef KERNEL_FP_64
#define REAL double
#define ZERO 0.0
#define ONE 1.0
#define TWO 2.0
#define FOUR 4.0
#define EIGHT 8.0
#define IP_U 0.5773502691896257
#define IP_V 0.5773502691896257
#define IP_W 0.5773502691896257
#define IP_S 1.0
#define viscoExponent (1.0/3.0)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

#ifdef MESH_DIMENSION_2
#define CDIM 2
#endif

#ifdef MESH_DIMENSION_3
#define CDIM 3
#endif

#ifdef ELEMENT_DIMENSION_2
#define DIM 2
#endif

#ifdef ELEMENT_DIMENSION_3
#define DIM 3
#endif

#ifdef MODEL_DIMENSION_2
#define MDIM 2
#endif

#ifdef MODEL_DIMENSION_3
#define MDIM 3
#endif

#ifdef UNIT_OFFSET
#define OFFSET 1
#else
#define OFFSET 0
#endif

#ifdef LOOP_UNROLL_1
#define UNROLL_DEPTH 3
#endif

#ifdef LOOP_UNROLL_2
#define UNROLL_DEPTH 4
#endif

inline void basis3D(__global REAL *BasisFunctions, REAL basis[], REAL u, REAL v, REAL w, int numberOfNodes);
inline REAL dBasisdx3D(__global REAL *BasisFunctions,
                __global int *elementNodeIndexesStore,
                __global REAL *nodesX,
                __global REAL *nodesY,
                __global REAL *nodesZ,
                REAL dBasisdx[][3],
                REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int globalID);

inline void basis3D(__global REAL *BasisFunctions, REAL basis[], REAL u, REAL v, REAL w, int numberOfNodes) {
    
    REAL s;
    int stride = numberOfNodes + OFFSET;
    
    for(int n=0; n<numberOfNodes+OFFSET; n++) {
        s = ZERO;
        for(int i=0; i<numberOfNodes+OFFSET; i++) {
            s = s + BasisFunctions[(i+stride*3)+(n*(4*stride))] * pow(u, BasisFunctions[i+(n*(4*stride))])
                * pow(v, BasisFunctions[(i+stride)+(n*(4*stride))])
                * pow(w, BasisFunctions[(i+stride*2)+(n*(4*stride))]);
        }
        basis[n] = s;
    }
}

inline REAL dBasisdx3D(__global REAL *BasisFunctions,
                __global int *elementNodeIndexesStore,
                __global REAL *nodesX,
                __global REAL *nodesY,
                __global REAL *nodesZ,
                REAL dBasisdx[][3],
                REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int globalID) {
    
    int stride = numberOfNodes + OFFSET;
    REAL dLBasisdx[24] = {};
    REAL dx[9] = {};
    REAL covariantMetricTensor[9] = {};
    REAL elementMetric[9] = {};
    REAL GI[9] = {};
    REAL ltoGMap[9] = {};
    REAL accum1, accum2, accum3, detJ, s, t, z;
    
    // Nodal derivatives for 3D element
    for(int n=0; n<numberOfNodes+OFFSET; n++) {
        s = ZERO;
        t = ZERO;
        z = ZERO;
        for(int i=0; i<numberOfNodes+OFFSET; i++) {
            if(BasisFunctions[i+(n*(4*stride))] >= ONE) s = s
                + BasisFunctions[i+(n*(4*stride))]*BasisFunctions[(i+stride*3)+(n*(4*stride))]
                * pow(u, (BasisFunctions[i+(n*(4*stride))]-1))
                * pow(v, BasisFunctions[(i+stride)+(n*(4*stride))])
                * pow(w, BasisFunctions[(i+stride*2)+(n*(4*stride))]);
            if(BasisFunctions[(i+stride)+(n*(4*stride))] >= ONE) t = t
                + BasisFunctions[(i+stride)+(n*(4*stride))]
                * BasisFunctions[(i+stride*3)+(n*(4*stride))]
                * pow(u, BasisFunctions[i+(n*(4*stride))])
                * pow(v, (BasisFunctions[(i+stride)+(n*(4*stride))]-1))
                * pow(w, BasisFunctions[(i+stride*2)+(n*(4*stride))]);
            if(BasisFunctions[(i+stride*2)+(n*(4*stride))] >= ONE) z = z
                + BasisFunctions[(i+stride*2)+(n*(4*stride))]
                * BasisFunctions[(i+stride*3)+(n*(4*stride))]
                * pow(u, BasisFunctions[i+(n*(4*stride))])
                * pow(v, BasisFunctions[(i+stride)+(n*(4*stride))])
                * pow(w, (BasisFunctions[(i+stride*2)+(n*(4*stride))]-1));
        }
        dLBasisdx[3*n]   = s;
        dLBasisdx[3*n+1] = t;
        dLBasisdx[3*n+2] = z;
    }
    
    for (int i=0; i<DIM; i++) {
        accum1 = ZERO;
        accum2 = ZERO;
        accum3 = ZERO;
        for (int j=0; j<numberOfNodes; j++) {
            accum1 = accum1 + (nodesX[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] * dLBasisdx[3*j+i]);
            accum2 = accum2 + (nodesY[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] * dLBasisdx[3*j+i]);
            accum3 = accum3 + (nodesZ[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] * dLBasisdx[3*j+i]);
        }
        dx[i]   = accum1;
        dx[3+i] = accum2;
        dx[6+i] = accum3;
    }
    
    for (int i=0; i<DIM; i++) {
        for (int j=0; j<DIM; j++) {
            s = ZERO;
            for (int k=0; k<CDIM; k++) {
                s = s + ( dx[k*3+i] * dx[k*3+j]);
            }
            covariantMetricTensor[3*i+j] = s;
        }
    }
    
     // Volume elements
    detJ = ZERO;
    detJ = covariantMetricTensor[0] * ( covariantMetricTensor[4]*covariantMetricTensor[8] - covariantMetricTensor[5]*covariantMetricTensor[7] )
    + covariantMetricTensor[1] * ( covariantMetricTensor[5]*covariantMetricTensor[6] - covariantMetricTensor[3]*covariantMetricTensor[8] )
    + covariantMetricTensor[2] * ( covariantMetricTensor[3]*covariantMetricTensor[7] - covariantMetricTensor[4]*covariantMetricTensor[6] );
#ifdef KERNEL_FP_32
    if (detJ <= FLT_MIN) {
        return ZERO; // Something is going wrong (probably degenerated element)
    };
#elif KERNEL_FP_64
    if (detJ <= DBL_MIN) {
        return ZERO;
    };
#endif
    
    // Convert the metric to contravariant base
    // First invert matrix
    s = ONE / detJ;
    GI[0] =  s * ( covariantMetricTensor[4]*covariantMetricTensor[8] - covariantMetricTensor[7]*covariantMetricTensor[5] );
    GI[3] = -s * ( covariantMetricTensor[3]*covariantMetricTensor[8] - covariantMetricTensor[6]*covariantMetricTensor[5] );
    GI[6] =  s * ( covariantMetricTensor[3]*covariantMetricTensor[7] - covariantMetricTensor[6]*covariantMetricTensor[4] );
    
    GI[1] = -s * ( covariantMetricTensor[1]*covariantMetricTensor[8] - covariantMetricTensor[7]*covariantMetricTensor[2] );
    GI[4] =  s * ( covariantMetricTensor[0]*covariantMetricTensor[8] - covariantMetricTensor[6]*covariantMetricTensor[2] );
    GI[7] = -s * ( covariantMetricTensor[0]*covariantMetricTensor[7] - covariantMetricTensor[6]*covariantMetricTensor[1] );
    
    GI[2] =  s * ( covariantMetricTensor[1]*covariantMetricTensor[5] - covariantMetricTensor[4]*covariantMetricTensor[2] );
    GI[5] = -s * ( covariantMetricTensor[0]*covariantMetricTensor[5] - covariantMetricTensor[3]*covariantMetricTensor[2] );
    GI[8] =  s * ( covariantMetricTensor[0]*covariantMetricTensor[4] - covariantMetricTensor[3]*covariantMetricTensor[1] );
    
    // Only for 3D element
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            elementMetric[3*i+j] = GI[3*i+j];
        }
    }
    
    // Map local to global
    for (int i=0; i<CDIM; i++) {
        for (int j=0; j<DIM; j++) {
            s = ZERO;
            for (int k=0; k<DIM; k++) {
                s = s + dx[3*i+k] * elementMetric[k*3+j];
            }
            ltoGMap[3*i+j] = s;
        }
    }
    
    // dBasisdx
    for(int i=0; i<8; i++) {
        for (int j=0; j<3; j++) {
            dBasisdx[i][j] = ZERO;
        }
    }
    for (int i=0; i<numberOfNodes; i++) {
        for (int j=0; j<CDIM; j++) {
            for (int k=0; k<DIM; k++) {
                dBasisdx[i][j] = dBasisdx[i][j] + dLBasisdx[3*i+k]*ltoGMap[3*j+k];
            }
        }
    }
    return sqrt(detJ);
}

__kernel void NavierStokesCompose(__global REAL *BasisFunctions,
                                  __global REAL *values,
                                  __global REAL *rhs,
                                  __global int *diag,
                                  __global int *rows,
                                  __global int *cols,
                                  __global int *colorMapping,
                                  __global int *elementNodeIndexesStore,
                                  __global REAL *nodesX,
                                  __global REAL *nodesY,
                                  __global REAL *nodesZ,
                                  __global REAL *varSolution,
                                  __global int *varPermutation,
                                  __global char *newtonLinearization,      // Just a pointer to an integer value,
                                                                           // 1 for true, 0 for false
                                  REAL density,
                                  REAL viscosity,
                                  REAL load,                               // load = Flow bodyforce 2 or 3
                                  REAL hk,
                                  REAL mk,
                                  int positionInColorMapping,
                                  int numberOfNodes,
                                  int numberElementDofs,
                                  int nBasis,
                                  int varDofs) {
    
    // DIM=3, the only model dimension we support (given by MDIM) so that
    // c = MDIM = dim = model.dimension (see the original implementation)
    // Model dimension and mesh dimension are the same
    
    int k;
    int stride = numberOfNodes + OFFSET;
    int col, row;
    REAL basis[8] = {};
    REAL dBasisdx[8][3] = {};
    REAL dNodalBasisdx[8][8][3] = {};
    REAL dmudx[3] = {};
    REAL force[4] = {};
    REAL grad[3][3] = {};
    REAL gradT[3][3] = {};
    REAL strain[3][3] = {};
    REAL su[8][4][4] = {};
    REAL sw[8][4][4] = {};
    REAL cc[4][4];
    REAL stiff[32][32];
    REAL mass[32][32];
    REAL forceVector[32];
    REAL jacM[32][32];
    REAL c2, c3, delta, detJ, mu, muder0, muder, rho, u, v, w, s, ss, sum, tau;
    
    // Definitions for 8 nodes octahedron
    REAL NodeU[8] = {-ONE, ONE, ONE, -ONE, -ONE, ONE, ONE, -ONE};
    REAL NodeV[8] = {-ONE, -ONE, ONE, ONE, -ONE, -ONE, ONE, ONE};
    REAL NodeW[8] = {-ONE, -ONE, -ONE, -ONE, ONE, ONE, ONE, ONE};
    
    // Gauss points for 8 nodes octahedron
    REAL ip_u[8] = {IP_U, -IP_U, IP_U, -IP_U, IP_U, -IP_U, IP_U, -IP_U};
    REAL ip_v[8] = {IP_V, IP_V, -IP_V, -IP_V, IP_V, IP_V, -IP_V, -IP_V};
    REAL ip_w[8] = {IP_W, IP_W, IP_W, IP_W, -IP_W, -IP_W, -IP_W, -IP_W};
    REAL ip_s[8] = {IP_S, IP_S, IP_S, IP_S, IP_S, IP_S, IP_S, IP_S};
    
    REAL criticalShearRate=1.0e-10;
    bool viscNewtonLin;
    
    int workItemID = get_global_id(0);
    int globalID = colorMapping[positionInColorMapping+workItemID];
    
    if (globalID < 0) return;
    
    for(int i=0; i<varDofs*numberOfNodes; i++) {
        for(int j=0; j<varDofs*numberOfNodes; j++) {
            stiff[i][j] = ZERO;
        }
    }
    for(int i=0; i<varDofs*numberOfNodes; i++) {
        for(int j=0; j<varDofs*numberOfNodes; j++) {
            mass[i][j] = ZERO;
        }
    }
    for(int i=0; i<varDofs*numberOfNodes; i++) {
        forceVector[i] = ZERO;
    }

    for(int i=0; i<varDofs*numberOfNodes; i++) {
        for(int j=0; j<varDofs*numberOfNodes; j++) {
            jacM[i][j] = ZERO;
        }
    }
    
    // Always use stabilization, bubbles not supported
    for (int p=0; p<numberOfNodes; p++) {
        u = NodeU[p];
        v = NodeV[p];
        w = NodeW[p];
        detJ = dBasisdx3D(BasisFunctions, elementNodeIndexesStore, nodesX, nodesY, nodesZ, dBasisdx,
                   u, v, w, numberOfNodes, numberElementDofs, globalID);
        for (int i=0; i<numberOfNodes; i++) {
            for (int j=0; j<3; j++) {
                dNodalBasisdx[i][p][j] = dBasisdx[i][j];
            }
        }
    }
    
    // Now we start integrating
    // Assume that number of Gauss points is equal to the number of nodes
    for (int t=0; t<numberOfNodes; t++) {
        u = ip_u[t];
        v = ip_v[t];
        w = ip_w[t];
        
        basis3D(BasisFunctions, basis, u, v, w, numberOfNodes);
        detJ = dBasisdx3D(BasisFunctions, elementNodeIndexesStore, nodesX, nodesY, nodesZ, dBasisdx,
                   u, v, w, numberOfNodes, numberElementDofs, globalID);
        s = detJ * ip_s[t];
        
        rho = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
            rho = rho + density * basis[i];
        }
        
        for (int i=0; i<3; i++) {
            grad[0][i] = ZERO;
            grad[1][i] = ZERO;
            grad[2][i] = ZERO;
            for (int j=0; j<numberOfNodes; j++) {
                grad[0][i] = grad[0][i] + varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]] * dBasisdx[j][i];
                grad[1][i] = grad[1][i] + varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]+1] * dBasisdx[j][i];
                grad[2][i] = grad[2][i] + varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]+2] * dBasisdx[j][i];
            }
        }
        
        // Force at integration point
        // Flow bodyforce 3 for MDIM=3
        force[2] = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
            force[2] = force[2] + load * basis[i];
        }

        // Always assume isotropic material
        // --------------------------------
        
        mu = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
            mu = mu + viscosity * basis[i];
        }
        // Viscous Non-Newtonian
        // The viscosity exponent for ice
        c2 = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
            c2 = c2 + viscoExponent * basis[i];
        }
        
        // Second invariant for Cartesian coordinate
        ss = ZERO;
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                ss = ss + (grad[i][j] + grad[j][i]) * (grad[i][j] + grad[j][i]);
            }
        }
        ss = ss / TWO;
        
        muder0 = mu * (c2-ONE)/TWO * pow(ss, (c2-ONE)/TWO - ONE);
        
        // Critical Shear Rate
        c3 = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
            c3 = c3 + criticalShearRate * basis[i];
        }
        if (ss < pow(c3, TWO)) {
            ss = pow(c3, TWO);
            muder0 = ZERO;
        }
        
        // The effective viscosity
        mu = mu * pow(ss, (c2-ONE)/TWO);
        
        viscNewtonLin = (*newtonLinearization && muder0 != ZERO) ? true : false;
        if (viscNewtonLin) {
            // Transpose
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    gradT[j][i] = grad[i][j];
                }
            }
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    strain[i][j] = (grad[i][j] + gradT[i][j]) / TWO;
                }
            }
        }
        
        // Always do the stabilization
        // ---------------------------
        
        for (int i=0; i<3; i++) {
            dmudx[i] = ZERO;
            for (int j=0; j<numberOfNodes; j++) {
                dmudx[i] = dmudx[i] + viscosity * dBasisdx[j][i];
            }
        }
        // Stabilization parameters Tau and Delta
        delta = ZERO;
        tau = mk * pow(hk, TWO) / (EIGHT * mu);
        
        for(int i=0; i<numberOfNodes; i++) {
            for(int j=0; j<4; j++) {
                for(k=0; k<4; k++) {
                    su[i][j][k] = ZERO;
                    sw[i][j][k] = ZERO;
                }
            }
        }
        for (int p=0; p<numberOfNodes; p++) {
            for (int i=0; i<MDIM; i++) {
                su[p][i][MDIM] = su[p][i][MDIM] + dBasisdx[p][i];
                for (int j=0; j<MDIM; j++) {
                    su[p][i][i] = su[p][i][i] - dmudx[j] * dBasisdx[p][j];
                    su[p][i][j] = su[p][i][j] - dmudx[j] * dBasisdx[p][i];
                    sum = ZERO;
                    for (k=0; k<numberOfNodes; k++) {
                        sum = sum + dNodalBasisdx[p][k][j] * dBasisdx[k][j];
                    }
                    su[p][i][i] = su[p][i][i] - mu * sum;
                    sum = ZERO;
                    for (k=0; k<numberOfNodes; k++) {
                        sum = sum + dNodalBasisdx[p][k][i] * dBasisdx[k][j];
                    }
                    su[p][i][j] = su[p][i][j] - mu * sum;
                }
                sw[p][MDIM][i] = sw[p][MDIM][i] + dBasisdx[p][i];
                for (int j=0; j<MDIM; j++) {
                    sw[p][i][i] = sw[p][i][i] - dmudx[j] * dBasisdx[p][j];
                    sw[p][j][i] = sw[p][j][i] - dmudx[j] * dBasisdx[p][i];
                    sum = ZERO;
                    for (k=0; k<numberOfNodes; k++) {
                        sum = sum + dNodalBasisdx[p][k][j] * dBasisdx[k][j];
                    }
                    sw[p][i][i] = sw[p][i][i] - mu * sum;
                    sum = ZERO;
                    for (k=0; k<numberOfNodes; k++) {
                        sum = sum + dNodalBasisdx[p][k][i] * dBasisdx[k][j];
                    }
                    sw[p][j][i] = sw[p][j][i] - mu * sum;
                }
            }
        }
        
        // Loop over basis functions (of both unknowns and weights)
        for (int p=0; p<nBasis; p++) {
            for (int q=0; q<nBasis; q++) {
                // Mass matrix
                // Momentum equations
                for (int i=0; i<MDIM; i++) {
                    mass[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = mass[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * rho * basis[q] * basis[p];
                }
                if (viscNewtonLin) {
                    for (int i=0; i<MDIM; i++) {
                        sum = 0.0;
                        for (k=0; k<3; k++) {
                            sum = sum + strain[i][k] * dBasisdx[q][k];
                        }
                        muder = muder0 * FOUR * sum;
                        for (int j=0; j<MDIM; j++) {
                            sum = 0.0;
                            for (k=0; k<3; k++) {
                                sum = sum + strain[j][k] * dBasisdx[p][k];
                            }
                            jacM[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = jacM[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * TWO * muder * sum;
                        }
                    }
                }
                
                for (int i=0; i<MDIM; i++) {
                    for (int j=0; j<MDIM; j++) {
                        // Isotropic and no Laplace discretization
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * mu * dBasisdx[q][j] * dBasisdx[p][j];
                        
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * mu * dBasisdx[q][i] * dBasisdx[p][j];
                    }
                    // Pressure term
                    stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] - s * basis[q] * dBasisdx[p][i];
                    
                    // Continuity equation
                    stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] + s * dBasisdx[q][i] * basis[p];
                }
                
                // Add stabilization
                for (int i=0; i<MDIM; i++) {
                    for (int j=0; j<MDIM+1; j++) {
                        mass[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = mass[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * tau * rho * basis[q] * sw[p][j][i];
                    }
                    for (int j=0; j<MDIM; j++) {
                        stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * delta * dBasisdx[q][i] * dBasisdx[p][j];
                    }
                }
                for (int i=0; i<MDIM+1; i++) {
                    for (int j=0; j<MDIM+1; j++) {
                        sum = ZERO;
                        for (k=0; k<MDIM; k++) {
                            sum = sum + sw[p][i][k]*su[q][k][j];
                        }
                        cc[i][j] = sum;
                    }
                }
                for (int i=0; i<MDIM+1; i++) {
                    for (int j=0; j<MDIM+1; j++) {
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * tau * cc[i][j];
                    }
                }
            }
        }
        
        // The right hand side...
        for (int p=0; p<nBasis; p++) {
            for (int i=0; i<MDIM+1; i++) {
                forceVector[((MDIM+1)*p)+i] = forceVector[((MDIM+1)*p)+i] + s * rho * force[i] * basis[p];
            }
            // Add stabilization
            for (int i=0; i<MDIM; i++) {
                for (int j=0; j<MDIM+1; j++) {
                    forceVector[((MDIM+1)*p)+j] = forceVector[((MDIM+1)*p)+j] + s * tau * rho * force[i] * sw[p][j][i];
                }
            }
        }
    }
    
    if (viscNewtonLin) {
        REAL sol[32];
        for(int i=0; i<varDofs*numberOfNodes; i++) {
            sol[i] = ZERO;
        }
        
        int kk = 0;
        for (int i=0; i<(MDIM+1)*numberOfNodes; i+=(MDIM+1)) {
            sol[i] = varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+kk]]];
            kk++;
        }
        kk = 0;
        for (int i=1; i<(MDIM+1)*numberOfNodes; i+=(MDIM+1)) {
            sol[i] = varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+kk]]+1];
            kk++;
        }
        kk = 0;
        for (int i=2; i<(MDIM+1)*numberOfNodes; i+=(MDIM+1)) {
            sol[i] = varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+kk]]+2];
            kk++;
        }
        
        int p = (MDIM+1) * nBasis;
        for (int i=0; i<p; i++) {
            for (int j=0; j<p; j++) {
                stiff[i][j] = stiff[i][j] + jacM[i][j];
            }
        }
        REAL yy[32];
        for (int i=0; i<p; i++) {
            yy[i] = ZERO;
        }
        for (int i=0; i<p; i++) {
            sum = ZERO;
            for (int j=0; j<p; j++) {
                sum = sum + jacM[i][j] * sol[j];
            }
            yy[i] = sum;
        }
        for (int i=0; i<p; i++) {
            forceVector[i] = forceVector[i] + yy[i];
        }
    }
    
    // The contribution of the local matrix to the global matrix
    // Only for DOF > 1
    for (int i=0; i<numberElementDofs; i++) {
        for (k=1; k<=varDofs; k++) {
            if (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]] < 0) continue;
            row = varDofs * (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]]+1) - k;
            for (int j=0; j<numberElementDofs; j++) {
                for (int l=1; l<=varDofs; l++) {
                    if (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)]+j] < 0) continue;
                    col = varDofs * (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]+1) - l;
                    if (col >= row) {
                        for (int c=diag[row]; c<=rows[row+1]-1; c++) {
                            if (cols[c] == col) {
                                values[c] = values[c] + stiff[varDofs*(i+1)-k][varDofs*(j+1)-l];
                                break;
                            }
                        }
                    } else {
                        for (int c=rows[row]; c<=diag[row]-1; c++) {
                            if (cols[c] == col) {
                                values[c] = values[c] + stiff[varDofs*(i+1)-k][varDofs*(j+1)-l];
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // The right-hand side
    for (int i=0; i<numberElementDofs; i++) {
        if (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)]+i] >= 0) {
            for (int j=0; j<varDofs; j++) {
                k = varDofs * varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]] + j;
                rhs[k] = rhs[k] + forceVector[varDofs*i+j];
            }
        }
    }
}
