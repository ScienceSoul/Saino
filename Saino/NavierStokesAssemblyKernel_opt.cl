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
#define COEFF 0.125f
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
#define COEFF 0.125
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

typedef struct _basis_functions {
    int p[8];
    int q[8];
    int r[8];
    REAL coeff[8];
} _basis_functions;

typedef struct _private_basis_functions {
    REAL coeff[8];
} _private_basis_functions;

typedef struct _element_info {
    REAL nodes_x[8];
    REAL nodes_y[8];
    REAL nodes_z[8];
    REAL vx[8];
    REAL vy[8];
    REAL vz[8];
#ifdef KERNEL_BASIS_DBASISDX
    REAL dBasisdx[8][8][3];
    REAL detJ[8];
#endif
    REAL padding;
} _element_info;

typedef struct {
    REAL basis[8][8];
} _element_basis;

typedef struct {
    REAL dBasisdx[8][8][3];
    REAL detJ[8];
} _element_dBasisdx;

typedef struct {
    REAL nodes_x[8];
    REAL nodes_y[8];
    REAL nodes_z[8];
    REAL vx[8];
    REAL vy[8];
    REAL vz[8];
    int perm[8];
} _nodal_data;

typedef struct {
    int indexes[1024];
} _non_zero;

typedef struct {
    REAL stiff[32][32];
    REAL force[32];
} _stiff_force;

#ifdef GLOBAL_BASIS_FUNCTIONS
    inline void basis3D(__global _basis_functions *basis_functions, REAL basis[], REAL u, REAL v, REAL w, int numberOfNodes);
#else
    inline void basis3D(_private_basis_functions *basis_functions, int p[], int q[], int r[], REAL basis[], REAL u, REAL v, REAL w, int numberOfNodes);
#endif

#ifdef USE_GPU_LOCAL_MEM
#ifdef GLOBAL_BASIS_FUNCTIONS
    inline REAL dBasisdx3D(__global _basis_functions *basis_functions,
                       __local _element_info *element_info,
                       REAL dBasisdx[][3],
                       REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int localID);
#else
    inline REAL dBasisdx3D(_private_basis_functions *basis_functions,
                       int p[], int q[], int r[],
                       __local _element_info *element_info,
                       REAL dBasisdx[][3],
                       REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int localID);
#endif
#else
#ifdef GLOBAL_BASIS_FUNCTIONS
    inline REAL dBasisdx3D(__global _basis_functions *basis_functions
#ifdef NODAL_DATA
                           , __global _nodal_data *nodal_data,
#else
                           , __global int *elementNodeIndexesStore,
                           __global REAL *nodesX,
                           __global REAL *nodesY,
                           __global REAL *nodesZ,
#endif
                           REAL dBasisdx[][3],
                           REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int globalID);
#else
inline REAL dBasisdx3D(_private_basis_functions *basis_functions,
                        int p[], int q[], int r[]
#ifdef NODAL_DATA
                       , __global _nodal_data *nodal_data,
#else
                       , __global int *elementNodeIndexesStore,
                       __global REAL *nodesX,
                       __global REAL *nodesY,
                       __global REAL *nodesZ,
#endif
                       REAL dBasisdx[][3],
                       REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int globalID);
#endif
#endif


#ifdef GLOBAL_BASIS_FUNCTIONS
    inline void basis3D(__global _basis_functions *basis_functions, REAL basis[], REAL u, REAL v, REAL w, int numberOfNodes) {
#else
    inline void basis3D(_private_basis_functions *basis_functions, int p[], int q[], int r[], REAL basis[], REAL u, REAL v, REAL w, int numberOfNodes) {
#endif

    REAL s;
    for(int n=0; n<numberOfNodes; n++) {
        s = ZERO;
        for(int i=0; i<numberOfNodes; i++) {
#ifdef GLOBAL_BASIS_FUNCTIONS
            s = s + basis_functions[n].coeff[i] * pow(u, basis_functions[n].p[i])
            * pow(v, basis_functions[n].q[i])
            * pow(w, basis_functions[n].r[i]);
#else
            s = s + basis_functions[n].coeff[i] * pow(u, p[i])
            * pow(v, q[i])
            * pow(w, r[i]);
#endif
        }
        basis[n] = s;
    }
}

#ifdef USE_GPU_LOCAL_MEM
#ifdef GLOBAL_BASIS_FUNCTIONS
        inline REAL dBasisdx3D(__global _basis_functions *basis_functions,
                               __local _element_info *element_info,
                               REAL dBasisdx[][3],
                               REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int localID) {
#else
        inline REAL dBasisdx3D(_private_basis_functions *basis_functions,
                               int p[], int q[], int r[],
                               __local _element_info *element_info,
                               REAL dBasisdx[][3],
                               REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int localID) {
#endif
#else
#ifdef GLOBAL_BASIS_FUNCTIONS
        inline REAL dBasisdx3D(__global _basis_functions *basis_functions
#ifdef NODAL_DATA
                               , __global _nodal_data *nodal_data,
#else
                               , __global int *elementNodeIndexesStore,
                               __global REAL *nodesX,
                               __global REAL *nodesY,
                               __global REAL *nodesZ,
#endif
                               REAL dBasisdx[][3],
                               REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int globalID) {
#else
        inline REAL dBasisdx3D(_private_basis_functions *basis_functions,
                                int p[], int q[], int r[]
#ifdef NODAL_DATA
                               , __global _nodal_data *nodal_data,
#else
                                , __global int *elementNodeIndexesStore,
                                __global REAL *nodesX,
                               __global REAL *nodesY,
                               __global REAL *nodesZ,
#endif
                                REAL dBasisdx[][3],
                                REAL u, REAL v, REAL w, int numberOfNodes, int numberElementDofs, int globalID) {
#endif
#endif


    REAL dLBasisdx[24] = {};
    REAL dx[9] = {};
    REAL covariantMetricTensor[9] = {};
    REAL elementMetric[9] = {};
    REAL GI[9] = {};
    REAL ltoGMap[9] = {};
    REAL accum1, accum2, accum3, detJ, s, t, z;
    
    // Nodal derivatives for 3D element
    for(int n=0; n<numberOfNodes; n++) {
        s = ZERO;
        t = ZERO;
        z = ZERO;
        for(int i=0; i<numberOfNodes; i++) {
#ifdef GLOBAL_BASIS_FUNCTIONS
            if(basis_functions[n].p[i] >= ONE) s = s
                + basis_functions[n].p[i] * basis_functions[n].coeff[i]
                * pow(u, basis_functions[n].p[i]-1)
                * pow(v, basis_functions[n].q[i])
                * pow(w, basis_functions[n].r[i]);
#else
            if(p[i] >= ONE) s = s
                + p[i] * basis_functions[n].coeff[i]
                * pow(u, p[i]-1)
                * pow(v, q[i])
                * pow(w, r[i]);
#endif
#ifdef GLOBAL_BASIS_FUNCTIONS
            if(basis_functions[n].q[i] >= ONE) t = t
                + basis_functions[n].q[i] * basis_functions[n].coeff[i]
                * pow(u, basis_functions[n].p[i])
                * pow(v, basis_functions[n].q[i]-1)
                * pow(w, basis_functions[n].r[i]);
#else
            if(q[i] >= ONE) t = t
                + q[i] * basis_functions[n].coeff[i]
                * pow(u, p[i])
                * pow(v, q[i]-1)
                * pow(w, r[i]);
#endif
#ifdef GLOBAL_BASIS_FUNCTIONS
            if(basis_functions[n].r[i] >= ONE) z = z
                + basis_functions[n].r[i] * basis_functions[n].coeff[i]
                * pow(u, basis_functions[n].p[i])
                * pow(v, basis_functions[n].q[i])
                * pow(w, basis_functions[n].r[i]-1);
#else
            if(r[i] >= ONE) z = z
                + r[i] * basis_functions[n].coeff[i]
                * pow(u, p[i])
                * pow(v, q[i])
                * pow(w, r[i]-1);
#endif
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
#ifdef USE_GPU_LOCAL_MEM
            accum1 = accum1 + (element_info[localID].nodes_x[j] * dLBasisdx[3*j+i]);
            accum2 = accum2 + (element_info[localID].nodes_y[j] * dLBasisdx[3*j+i]);
            accum3 = accum3 + (element_info[localID].nodes_z[j] * dLBasisdx[3*j+i]);
#else
    #ifdef NODAL_DATA
            accum1 = accum1 + (nodal_data[globalID].nodes_x[j] * dLBasisdx[3*j+i]);
            accum2 = accum2 + (nodal_data[globalID].nodes_y[j] * dLBasisdx[3*j+i]);
            accum3 = accum3 + (nodal_data[globalID].nodes_z[j] * dLBasisdx[3*j+i]);
    #else
            accum1 = accum1 + (nodesX[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] * dLBasisdx[3*j+i]);
            accum2 = accum2 + (nodesY[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] * dLBasisdx[3*j+i]);
            accum3 = accum3 + (nodesZ[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] * dLBasisdx[3*j+i]);
    #endif
#endif
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
#ifdef DEBUG
#ifdef KERNEL_FP_32
    if (detJ <= FLT_MIN) {
        return ZERO; // Something is going wrong (probably degenerated element)
    };
#elif KERNEL_FP_64
    if (detJ <= DBL_MIN) {
        return ZERO;
    };
#endif
#endif
    
    // Convert the metric to contravariant base
    // First invert matrix
    s = native_divide(ONE,detJ);
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
    return native_sqrt(detJ);
}

#ifdef KERNEL_BASIS_DBASISDX

__kernel void ComputeBasisDBasisdx(__global _element_basis *global_basis,
                                   __global _element_dBasisdx *global_dBasisdx
#ifdef NODAL_DATA
                                   , __global _nodal_data *nodal_data,
#else
                                   , __global REAL *nodesX,
                                   __global REAL *nodesY,
                                   __global REAL *nodesZ,
                                   __global int *elementNodeIndexesStore,
#endif
                                   int numberOfNodes,
                                   int numberElementDofs,
                                   int numberOfElements
#ifdef GLOBAL_BASIS_FUNCTIONS
                                   , __global _basis_functions *basis_functions
#endif
#ifdef USE_GPU_LOCAL_MEM
                                   , __local _element_info *element_info
#endif
                                   ) {
    
    REAL basis[8] = {};
    REAL dBasisdx[8][3] = {};
    REAL dLBasisdx[24] = {};
    REAL dx[9] = {};
    REAL covariantMetricTensor[9] = {};
    REAL elementMetric[9] = {};
    REAL GI[9] = {};
    REAL ltoGMap[9] = {};
    REAL accum1, accum2, accum3, detJ, s, t, u, v, w, z;

    // Gauss points for 8 nodes octahedron element
    REAL ip_u[8] = {IP_U, -IP_U, IP_U, -IP_U, IP_U, -IP_U, IP_U, -IP_U};
    REAL ip_v[8] = {IP_V, IP_V, -IP_V, -IP_V, IP_V, IP_V, -IP_V, -IP_V};
    REAL ip_w[8] = {IP_W, IP_W, IP_W, IP_W, -IP_W, -IP_W, -IP_W, -IP_W};
    REAL ip_s[8] = {IP_S, IP_S, IP_S, IP_S, IP_S, IP_S, IP_S, IP_S};
    
    int global_id = get_global_id(0);
    int local_id = get_local_id(0);
    
#ifndef GLOBAL_BASIS_FUNCTIONS
    int p[8], q[8], r[8];
    p[0] = 0; p[1] = 1; p[2] = 0; p[3] = 1;
    p[4] = 0; p[5] = 1; p[6] = 0; p[7] = 1;
    
    q[0] = 0; q[1] = 0; q[2] = 1; q[3] = 1;
    q[4] = 0; q[5] = 0; q[6] = 1; q[7] = 1;
    
    r[0] = 0; r[1] = 0; r[2] = 0; r[3] = 0;
    r[4] = 1; r[5] = 1; r[6] = 1; r[7] = 1;
    
    _private_basis_functions basis_functions[8];
    basis_functions[0].coeff[0] =  COEFF; basis_functions[0].coeff[1] = -COEFF;
    basis_functions[0].coeff[2] = -COEFF; basis_functions[0].coeff[3] =  COEFF;
    basis_functions[0].coeff[4] = -COEFF; basis_functions[0].coeff[5] =  COEFF;
    basis_functions[0].coeff[6] =  COEFF; basis_functions[0].coeff[7] = -COEFF;
    
    basis_functions[1].coeff[0] =  COEFF; basis_functions[1].coeff[1] =  COEFF;
    basis_functions[1].coeff[2] = -COEFF; basis_functions[1].coeff[3] = -COEFF;
    basis_functions[1].coeff[4] = -COEFF; basis_functions[1].coeff[5] = -COEFF;
    basis_functions[1].coeff[6] =  COEFF; basis_functions[1].coeff[7] =  COEFF;
    
    basis_functions[2].coeff[0] =  COEFF; basis_functions[2].coeff[1] =  COEFF;
    basis_functions[2].coeff[2] =  COEFF; basis_functions[2].coeff[3] =  COEFF;
    basis_functions[2].coeff[4] = -COEFF; basis_functions[2].coeff[5] = -COEFF;
    basis_functions[2].coeff[6] = -COEFF; basis_functions[2].coeff[7] = -COEFF;
    
    basis_functions[3].coeff[0] =  COEFF; basis_functions[3].coeff[1] = -COEFF;
    basis_functions[3].coeff[2] =  COEFF; basis_functions[3].coeff[3] = -COEFF;
    basis_functions[3].coeff[4] = -COEFF; basis_functions[3].coeff[5] =  COEFF;
    basis_functions[3].coeff[6] = -COEFF; basis_functions[3].coeff[7] =  COEFF;
    
    basis_functions[4].coeff[0] =  COEFF; basis_functions[4].coeff[1] = -COEFF;
    basis_functions[4].coeff[2] = -COEFF; basis_functions[4].coeff[3] =  COEFF;
    basis_functions[4].coeff[4] =  COEFF; basis_functions[4].coeff[5] = -COEFF;
    basis_functions[4].coeff[6] = -COEFF; basis_functions[4].coeff[7] =  COEFF;
    
    basis_functions[5].coeff[0] =  COEFF; basis_functions[5].coeff[1] =  COEFF;
    basis_functions[5].coeff[2] = -COEFF; basis_functions[5].coeff[3] = -COEFF;
    basis_functions[5].coeff[4] =  COEFF; basis_functions[5].coeff[5] =  COEFF;
    basis_functions[5].coeff[6] = -COEFF; basis_functions[5].coeff[7] = -COEFF;
    
    basis_functions[6].coeff[0] =  COEFF; basis_functions[6].coeff[1] =  COEFF;
    basis_functions[6].coeff[2] =  COEFF; basis_functions[6].coeff[3] =  COEFF;
    basis_functions[6].coeff[4] =  COEFF; basis_functions[6].coeff[5] =  COEFF;
    basis_functions[6].coeff[6] =  COEFF; basis_functions[6].coeff[7] =  COEFF;
    
    basis_functions[7].coeff[0] =  COEFF; basis_functions[7].coeff[1] = -COEFF;
    basis_functions[7].coeff[2] =  COEFF; basis_functions[7].coeff[3] = -COEFF;
    basis_functions[7].coeff[4] =  COEFF; basis_functions[7].coeff[5] = -COEFF;
    basis_functions[7].coeff[6] =  COEFF; basis_functions[7].coeff[7] = -COEFF;
#endif
    
#if defined (USE_GPU_LOCAL_MEM) || defined (ENABLE_WORK_GROUPS)
    if (global_id >= numberOfElements) return;
#endif
 
#ifdef USE_GPU_LOCAL_MEM
    #ifdef NODAL_DATA
    for (int i=0; i<numberOfNodes; i++) {
        element_info[local_id].nodes_x[i] = nodal_data[global_id].nodes_x[i];
        element_info[local_id].nodes_y[i] = nodal_data[global_id].nodes_y[i];
        element_info[local_id].nodes_z[i] = nodal_data[global_id].nodes_z[i];
    }
    #else
    for (int i=0; i<numberOfNodes; i++) {
        element_info[local_id].nodes_x[i] = nodesX[elementNodeIndexesStore[(global_id*numberElementDofs)+i]];
        element_info[local_id].nodes_y[i] = nodesY[elementNodeIndexesStore[(global_id*numberElementDofs)+i]];
        element_info[local_id].nodes_z[i] = nodesZ[elementNodeIndexesStore[(global_id*numberElementDofs)+i]];
    }
    #endif
    
    barrier(CLK_LOCAL_MEM_FENCE);
#endif

    for (int t=0; t<numberOfNodes; t++) {
        u = ip_u[t];
        v = ip_v[t];
        w = ip_w[t];
#ifdef USE_GPU_LOCAL_MEM
#ifdef GLOBAL_BASIS_FUNCTIONS
        detJ = dBasisdx3D(basis_functions, element_info, dBasisdx, u, v, w, numberOfNodes, numberElementDofs, local_id);
#else
        detJ = dBasisdx3D(basis_functions, p, q, r,
                          element_info, dBasisdx, u, v, w, numberOfNodes, numberElementDofs, local_id);
#endif
#else
#ifdef GLOBAL_BASIS_FUNCTIONS
    #ifdef NODAL_DATA
        detJ = dBasisdx3D(basis_functions, nodal_data, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, global_id);
    #else
        detJ = dBasisdx3D(basis_functions, elementNodeIndexesStore, nodesX, nodesY, nodesZ, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, global_id);
    #endif
#else
    #ifdef NODAL_DATA
        detJ = dBasisdx3D(basis_functions, p, q, r, nodal_data, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, global_id);
    #else
        detJ = dBasisdx3D(basis_functions, p, q, r, elementNodeIndexesStore, nodesX, nodesY, nodesZ, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, global_id);
    #endif
#endif
#endif

#ifdef USE_GPU_LOCAL_MEM
        for (int i=0; i<numberOfNodes; i++) {
            for (int j=0; j<3; j++) {
                element_info[local_id].dBasisdx[t][i][j] = dBasisdx[i][j];
            }
        }
        element_info[local_id].detJ[t] = detJ;
#else
        for (int i=0; i<numberOfNodes; i++) {
            for (int j=0; j<3; j++) {
                global_dBasisdx[global_id].dBasisdx[t][i][j] = dBasisdx[i][j];
            }
        }
        global_dBasisdx[global_id].detJ[t] = detJ;

#endif
    }
    
#ifdef USE_GPU_LOCAL_MEM
    barrier(CLK_LOCAL_MEM_FENCE);
    
    for (int t=0; t<numberOfNodes; t++) {
        for (int i=0; i<numberOfNodes; i++) {
            for (int j=0; j<3; j++) {
                global_dBasisdx[global_id].dBasisdx[t][i][j] = element_info[local_id].dBasisdx[t][i][j];
            }
        }
        global_dBasisdx[global_id].detJ[t] =  element_info[local_id].detJ[t];
    }
#endif
    
    if (global_id == 0) {
        for (int t=0; t<numberOfNodes; t++) {
            u = ip_u[t];
            v = ip_v[t];
            w = ip_w[t];
#ifdef GLOBAL_BASIS_FUNCTIONS
            basis3D(basis_functions, basis, u, v, w, numberOfNodes);
#else
            basis3D(basis_functions, p, q, r, basis, u, v, w, numberOfNodes);
#endif
            for (int i=0; i<numberOfNodes; i++) {
                global_basis->basis[t][i] = basis[i];
            }
        }
    }
}

#endif

__kernel void AssemblyByColoringStiffForceCompute(
#ifdef COLORING_ASSEMBLY
                                  
#ifdef GLOBAL_STIFF_FORCE
                                  __global _stiff_force *stiff_force,
#endif
                                  __global REAL *values,
                                  __global REAL *rhs
#ifdef PRECOMPUTE_NZ
                                  ,__global _non_zero *non_zeros,
#else
                                  ,__global int *diag,
                                  __global int *rows,
                                  __global int *cols,
#endif
                                  __global int *colorMapping
                                  
#endif
                                  
#ifdef NONZEROS_ASSEMBLY
                               __global REAL *element_stiffs,
                               __global REAL *element_forces
#endif
                                  
#ifdef NODAL_DATA
                                  , __global _nodal_data *nodal_data,
#else
                                  , __global int *elementNodeIndexesStore,
                                  __global REAL *nodesX,
                                  __global REAL *nodesY,
                                  __global REAL *nodesZ,
                                  __global REAL *varSolution,
                                  __global int *varPermutation,
#endif
                                  __global char *newtonLinearization,      // Just a pointer to an integer value, 1 for true, 0 for false
                                  REAL density,
                                  REAL viscosity,
                                  REAL load,                               // load = Flow bodyforce 2 or 3
                                  REAL hk,
                                  REAL mk,
#ifdef COLORING_ASSEMBLY
                                  int positionInColorMapping,
#endif
                                                  
#ifdef NONZEROS_ASSEMBLY
                                  int numberOfElements,
#endif
                                  int numberOfNodes,
                                  int numberElementDofs,
                                  int nBasis,
                                  int varDofs
#ifdef GLOBAL_BASIS_FUNCTIONS
                                  , __global _basis_functions *basis_functions
#endif
#ifdef KERNEL_BASIS_DBASISDX
                                  , __global _element_basis *global_basis,
                                  __global _element_dBasisdx *global_dBasisdx
#endif
#ifdef USE_GPU_LOCAL_MEM
                                  , __local _element_info *element_info
#endif
                                  ) {
    
    // DIM=3, the only model dimension we support (given by MDIM) so that
    // c = MDIM = dim = model.dimension (see the original implementation)
    // Model dimension and mesh dimension are the same
    
    int k;
    
#ifdef COLORING_ASSEMBLY
    int col, row, start, end;
#endif
    
#ifdef KERNEL_BASIS_DBASISDX
    REAL basis[8][8];
#else
    REAL basis[8] = {};
#endif
    REAL dBasisdx[8][3] = {};
    REAL dNodalBasisdx[8][8][3] = {};
    REAL dmudx[3] = {};
    REAL force[4] = {};
    REAL grad[3][3] = {};
#ifdef NEWTONLINEAR
    REAL gradT[3][3] = {};
    REAL strain[3][3] = {};
#endif
    REAL su[8][4][4] = {};
    REAL sw[8][4][4] = {};
    REAL cc[4][4];
#ifndef GLOBAL_STIFF_FORCE
    REAL stiff[32][32];
    REAL forceVector[32];
#endif
#ifdef TRANSIENT
    REAL mass[32][32];
#endif
#ifdef NEWTONLINEAR
    REAL jacM[32][32];
#endif
    REAL c2, c3, delta, detJ, mu, muder0, muder, rho, u, v, w, s, ss, sum, tau;
    
    REAL criticalShearRate=1.0e-10;
    bool viscNewtonLin;
    
    // Definitions for 8 nodes octahedron element
    REAL NodeU[8] = {-ONE, ONE, ONE, -ONE, -ONE, ONE, ONE, -ONE};
    REAL NodeV[8] = {-ONE, -ONE, ONE, ONE, -ONE, -ONE, ONE, ONE};
    REAL NodeW[8] = {-ONE, -ONE, -ONE, -ONE, ONE, ONE, ONE, ONE};
    
    // Gauss points for 8 nodes octahedron element
    REAL ip_u[8] = {IP_U, -IP_U, IP_U, -IP_U, IP_U, -IP_U, IP_U, -IP_U};
    REAL ip_v[8] = {IP_V, IP_V, -IP_V, -IP_V, IP_V, IP_V, -IP_V, -IP_V};
    REAL ip_w[8] = {IP_W, IP_W, IP_W, IP_W, -IP_W, -IP_W, -IP_W, -IP_W};
    REAL ip_s[8] = {IP_S, IP_S, IP_S, IP_S, IP_S, IP_S, IP_S, IP_S};
    
#ifndef GLOBAL_BASIS_FUNCTIONS
    int basis_p[8], basis_q[8], basis_r[8];
    basis_p[0] = 0; basis_p[1] = 1; basis_p[2] = 0; basis_p[3] = 1;
    basis_p[4] = 0; basis_p[5] = 1; basis_p[6] = 0; basis_p[7] = 1;
    
    basis_q[0] = 0; basis_q[1] = 0; basis_q[2] = 1; basis_q[3] = 1;
    basis_q[4] = 0; basis_q[5] = 0; basis_q[6] = 1; basis_q[7] = 1;
    
    basis_r[0] = 0; basis_r[1] = 0; basis_r[2] = 0; basis_r[3] = 0;
    basis_r[4] = 1; basis_r[5] = 1; basis_r[6] = 1; basis_r[7] = 1;
    
    _private_basis_functions basis_functions[8];
    basis_functions[0].coeff[0] =  COEFF; basis_functions[0].coeff[1] = -COEFF;
    basis_functions[0].coeff[2] = -COEFF; basis_functions[0].coeff[3] =  COEFF;
    basis_functions[0].coeff[4] = -COEFF; basis_functions[0].coeff[5] =  COEFF;
    basis_functions[0].coeff[6] =  COEFF; basis_functions[0].coeff[7] = -COEFF;
    
    basis_functions[1].coeff[0] =  COEFF; basis_functions[1].coeff[1] =  COEFF;
    basis_functions[1].coeff[2] = -COEFF; basis_functions[1].coeff[3] = -COEFF;
    basis_functions[1].coeff[4] = -COEFF; basis_functions[1].coeff[5] = -COEFF;
    basis_functions[1].coeff[6] =  COEFF; basis_functions[1].coeff[7] =  COEFF;
    
    basis_functions[2].coeff[0] =  COEFF; basis_functions[2].coeff[1] =  COEFF;
    basis_functions[2].coeff[2] =  COEFF; basis_functions[2].coeff[3] =  COEFF;
    basis_functions[2].coeff[4] = -COEFF; basis_functions[2].coeff[5] = -COEFF;
    basis_functions[2].coeff[6] = -COEFF; basis_functions[2].coeff[7] = -COEFF;
    
    basis_functions[3].coeff[0] =  COEFF; basis_functions[3].coeff[1] = -COEFF;
    basis_functions[3].coeff[2] =  COEFF; basis_functions[3].coeff[3] = -COEFF;
    basis_functions[3].coeff[4] = -COEFF; basis_functions[3].coeff[5] =  COEFF;
    basis_functions[3].coeff[6] = -COEFF; basis_functions[3].coeff[7] =  COEFF;
    
    basis_functions[4].coeff[0] =  COEFF; basis_functions[4].coeff[1] = -COEFF;
    basis_functions[4].coeff[2] = -COEFF; basis_functions[4].coeff[3] =  COEFF;
    basis_functions[4].coeff[4] =  COEFF; basis_functions[4].coeff[5] = -COEFF;
    basis_functions[4].coeff[6] = -COEFF; basis_functions[4].coeff[7] =  COEFF;
    
    basis_functions[5].coeff[0] =  COEFF; basis_functions[5].coeff[1] =  COEFF;
    basis_functions[5].coeff[2] = -COEFF; basis_functions[5].coeff[3] = -COEFF;
    basis_functions[5].coeff[4] =  COEFF; basis_functions[5].coeff[5] =  COEFF;
    basis_functions[5].coeff[6] = -COEFF; basis_functions[5].coeff[7] = -COEFF;
    
    basis_functions[6].coeff[0] =  COEFF; basis_functions[6].coeff[1] =  COEFF;
    basis_functions[6].coeff[2] =  COEFF; basis_functions[6].coeff[3] =  COEFF;
    basis_functions[6].coeff[4] =  COEFF; basis_functions[6].coeff[5] =  COEFF;
    basis_functions[6].coeff[6] =  COEFF; basis_functions[6].coeff[7] =  COEFF;
    
    basis_functions[7].coeff[0] =  COEFF; basis_functions[7].coeff[1] = -COEFF;
    basis_functions[7].coeff[2] =  COEFF; basis_functions[7].coeff[3] = -COEFF;
    basis_functions[7].coeff[4] =  COEFF; basis_functions[7].coeff[5] = -COEFF;
    basis_functions[7].coeff[6] =  COEFF; basis_functions[7].coeff[7] = -COEFF;
#endif
    
#ifdef COLORING_ASSEMBLY
    int workItemID = get_global_id(0);
    int globalID = colorMapping[positionInColorMapping+workItemID];
#endif

#ifdef NONZEROS_ASSEMBLY
    int globalID = get_global_id(0);
#endif
    
    //int local_size = get_local_size(0);
    int local_id = get_local_id(0);
    
#ifdef COLORING_ASSEMBLY
    #if defined (USE_GPU_LOCAL_MEM) || defined (ENABLE_WORK_GROUPS)
        if (globalID < 0) return;
    #endif
#endif
    
#ifdef NONZEROS_ASSEMBLY
    #if defined (USE_GPU_LOCAL_MEM) || defined (ENABLE_WORK_GROUPS)
        if (globalID >= numberOfElements) return;
    #endif
#endif
    
#ifdef KERNEL_BASIS_DBASISDX
    for(int i=0; i<numberOfNodes; i++) {
        for(int j=0; j<numberOfNodes; j++) {
            basis[i][j] = global_basis->basis[i][j];
        }
    }
#endif
    
#ifdef USE_GPU_LOCAL_MEM
    #ifdef NODAL_DATA
    for (int i=0; i<numberOfNodes; i++) {
        element_info[local_id].nodes_x[i] = nodal_data[globalID].nodes_x[i];
        element_info[local_id].nodes_y[i] = nodal_data[globalID].nodes_y[i];
        element_info[local_id].nodes_z[i] = nodal_data[globalID].nodes_z[i];
        
        element_info[local_id].vx[i] = nodal_data[globalID].vx[i];
        element_info[local_id].vy[i] = nodal_data[globalID].vy[i];
        element_info[local_id].vz[i] = nodal_data[globalID].vz[i];
    }
    #else
    for (int i=0; i<numberOfNodes; i++) {
        element_info[local_id].nodes_x[i] = nodesX[elementNodeIndexesStore[(globalID*numberElementDofs)+i]];
        element_info[local_id].nodes_y[i] = nodesY[elementNodeIndexesStore[(globalID*numberElementDofs)+i]];
        element_info[local_id].nodes_z[i] = nodesZ[elementNodeIndexesStore[(globalID*numberElementDofs)+i]];
        
        element_info[local_id].vx[i] = varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]]];
        element_info[local_id].vy[i] = varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]]+1];
        element_info[local_id].vz[i] = varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]]+2];
    }
    #endif
    
#ifdef KERNEL_BASIS_DBASISDX
    for(int i=0; i<numberOfNodes; i++) {
        for(int j=0; j<numberOfNodes; j++) {
            for(k=0; k<3; k++) {
                element_info[local_id].dBasisdx[i][j][k] = global_dBasisdx[globalID].dBasisdx[i][j][k];
            }
        }
        element_info[local_id].detJ[i] = global_dBasisdx[globalID].detJ[i];
    }
#endif
    
    barrier(CLK_LOCAL_MEM_FENCE);
#endif

    for(int i=0; i<varDofs*numberOfNodes; i++) {
        for(int j=0; j<varDofs*numberOfNodes; j++) {
#ifdef GLOBAL_STIFF_FORCE
            stiff_force[globalID].stiff[i][j] = ZERO;
#else
            stiff[i][j] = ZERO;
#endif
        }
    }
    
#ifdef TRANSIENT
    for(int i=0; i<varDofs*numberOfNodes; i++) {
        for(int j=0; j<varDofs*numberOfNodes; j++) {
            mass[i][j] = ZERO;
        }
    }
#endif
    
    for(int i=0; i<varDofs*numberOfNodes; i++) {
#ifdef GLOBAL_STIFF_FORCE
        stiff_force[globalID].force[i] = ZERO;
#else
        forceVector[i] = ZERO;
#endif
    }
    
#ifdef NEWTONLINEAR
    for(int i=0; i<varDofs*numberOfNodes; i++) {
        for(int j=0; j<varDofs*numberOfNodes; j++) {
            jacM[i][j] = ZERO;
        }
    }
#endif
    
    // Always use stabilization, bubbles not supported
    for (int p=0; p<numberOfNodes; p++) {
        u = NodeU[p];
        v = NodeV[p];
        w = NodeW[p];
#ifdef USE_GPU_LOCAL_MEM
#ifdef GLOBAL_BASIS_FUNCTIONS
        detJ = dBasisdx3D(basis_functions, element_info, dBasisdx, u, v, w, numberOfNodes, numberElementDofs, local_id);
#else
        detJ = dBasisdx3D(basis_functions, basis_p, basis_q, basis_r,
                          element_info, dBasisdx, u, v, w, numberOfNodes, numberElementDofs, local_id);
#endif
#else
#ifdef GLOBAL_BASIS_FUNCTIONS
    #ifdef NODAL_DATA
        detJ = dBasisdx3D(basis_functions, nodal_data, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, globalID);
    #else
        detJ = dBasisdx3D(basis_functions, elementNodeIndexesStore, nodesX, nodesY, nodesZ, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, globalID);
    #endif
#else
    #ifdef NODAL_DATA
        detJ = dBasisdx3D(basis_functions, basis_p, basis_q, basis_r, nodal_data, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, globalID);
    #else
        detJ = dBasisdx3D(basis_functions, basis_p, basis_q, basis_r, elementNodeIndexesStore, nodesX, nodesY, nodesZ, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, globalID);
    #endif
#endif
#endif
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

#ifndef KERNEL_BASIS_DBASISDX
#ifdef GLOBAL_BASIS_FUNCTIONS
        basis3D(basis_functions, basis, u, v, w, numberOfNodes);
#else
        basis3D(basis_functions, basis_p, basis_q, basis_r, basis, u, v, w, numberOfNodes);
#endif
#endif

#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
        s = element_info[local_id].detJ[t] * ip_s[t];
    #else
        s = global_dBasisdx[globalID].detJ[t] * ip_s[t];
    #endif
#else
    #ifdef USE_GPU_LOCAL_MEM
      #ifdef GLOBAL_BASIS_FUNCTIONS
        detJ = dBasisdx3D(basis_functions, element_info, dBasisdx, u, v, w, numberOfNodes, numberElementDofs, local_id);
      #else
        detJ = dBasisdx3D(basis_functions, basis_p, basis_q, basis_r,
                          element_info, dBasisdx, u, v, w, numberOfNodes, numberElementDofs, local_id);
      #endif
    #else
      #ifdef GLOBAL_BASIS_FUNCTIONS
        #ifdef NODAL_DATA
            detJ = dBasisdx3D(basis_functions, nodal_data, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, globalID);
        #else
            detJ = dBasisdx3D(basis_functions, elementNodeIndexesStore, nodesX, nodesY, nodesZ, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, globalID);
        #endif
      #else
        #ifdef NODAL_DATA
            detJ = dBasisdx3D(basis_functions, basis_p, basis_q, basis_r, nodal_data, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, globalID);
        #else
            detJ = dBasisdx3D(basis_functions, basis_p, basis_q, basis_r, elementNodeIndexesStore, nodesX, nodesY, nodesZ, dBasisdx,
                          u, v, w, numberOfNodes, numberElementDofs, globalID);
        #endif
      #endif
    #endif
        s = detJ * ip_s[t];
#endif
        
        rho = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
#ifdef KERNEL_BASIS_DBASISDX
            rho = rho + density * basis[t][i];
#else
            rho = rho + density * basis[i];
#endif
        }
        
        for (int i=0; i<3; i++) {
            grad[0][i] = ZERO;
            grad[1][i] = ZERO;
            grad[2][i] = ZERO;
            for (int j=0; j<numberOfNodes; j++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                grad[0][i] = grad[0][i] + element_info[local_id].vx[j] * element_info[local_id].dBasisdx[t][j][i];
                grad[1][i] = grad[1][i] + element_info[local_id].vy[j] * element_info[local_id].dBasisdx[t][j][i];
                grad[2][i] = grad[2][i] + element_info[local_id].vz[j] * element_info[local_id].dBasisdx[t][j][i];
    #else
        #ifdef NODAL_DATA
                grad[0][i] = grad[0][i] + nodal_data[globalID].vx[j] * global_dBasisdx[globalID].dBasisdx[t][j][i];
                grad[1][i] = grad[1][i] + nodal_data[globalID].vy[j] * global_dBasisdx[globalID].dBasisdx[t][j][i];
                grad[2][i] = grad[2][i] + nodal_data[globalID].vz[j] * global_dBasisdx[globalID].dBasisdx[t][j][i];

        #else
                grad[0][i] = grad[0][i] + varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]]   * global_dBasisdx[globalID].dBasisdx[t][j][i];
                grad[1][i] = grad[1][i] + varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]+1] * global_dBasisdx[globalID].dBasisdx[t][j][i];
                grad[2][i] = grad[2][i] + varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]+2] * global_dBasisdx[globalID].dBasisdx[t][j][i];
        #endif
    #endif
#else
    #ifdef USE_GPU_LOCAL_MEM
                grad[0][i] = grad[0][i] + element_info[local_id].vx[j] * dBasisdx[j][i];
                grad[1][i] = grad[1][i] + element_info[local_id].vy[j] * dBasisdx[j][i];
                grad[2][i] = grad[2][i] + element_info[local_id].vz[j] * dBasisdx[j][i];
    #else
        #ifdef NODAL_DATA
                grad[0][i] = grad[0][i] + nodal_data[globalID].vx[j] * dBasisdx[j][i];
                grad[1][i] = grad[1][i] + nodal_data[globalID].vy[j] * dBasisdx[j][i];
                grad[2][i] = grad[2][i] + nodal_data[globalID].vz[j] * dBasisdx[j][i];
        #else
                grad[0][i] = grad[0][i] + varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]]   * dBasisdx[j][i];
                grad[1][i] = grad[1][i] + varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]+1] * dBasisdx[j][i];
                grad[2][i] = grad[2][i] + varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]]+2] * dBasisdx[j][i];
        #endif
    #endif
#endif
            }
        }
        
        // Force at integration point
        // Flow bodyforce 3 for MDIM=3
        force[2] = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
#ifdef KERNEL_BASIS_DBASISDX
            force[2] = force[2] + load * basis[t][i];
#else
            force[2] = force[2] + load * basis[i];
#endif
        }
        
        // Always assume isotropic material
        // --------------------------------
        
        mu = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
#ifdef KERNEL_BASIS_DBASISDX
            mu = mu + viscosity * basis[t][i];
#else
            mu = mu + viscosity * basis[i];
#endif
        }
        // Viscous Non-Newtonian
        // The viscosity exponent for ice
        c2 = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
#ifdef KERNEL_BASIS_DBASISDX
            c2 = c2 + viscoExponent * basis[t][i];
#else
            c2 = c2 + viscoExponent * basis[i];
#endif
        }
        
        // Second invariant for Cartesian coordinate
        ss = ZERO;
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                ss = ss + (grad[i][j] + grad[j][i]) * (grad[i][j] + grad[j][i]);
            }
        }
        ss = native_divide(ss,TWO);
        
        muder0 = mu * (c2-ONE)/TWO * pow(ss, (c2-ONE)/TWO - ONE);
        
        // Critical Shear Rate
        c3 = ZERO;
        for (int i=0; i<numberOfNodes; i++) {
#ifdef KERNEL_BASIS_DBASISDX
            c3 = c3 + criticalShearRate * basis[t][i];
#else
            c3 = c3 + criticalShearRate * basis[i];
#endif
        }
        
        ss = (ss < pow(c3, TWO)) ? pow(c3, TWO) : ss;
        muder0 = (ss < pow(c3, TWO)) ? ZERO : muder0;
        
        // The effective viscosity
        mu = mu * pow(ss, (c2-ONE)/TWO);
#ifdef NEWTONLINEAR
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
                    strain[i][j] = native_divide(grad[i][j] + gradT[i][j],TWO);
                }
            }
        }
#endif
        
        // Always do the stabilization
        // ---------------------------
        
        for (int i=0; i<3; i++) {
            dmudx[i] = ZERO;
            for (int j=0; j<numberOfNodes; j++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                dmudx[i] = dmudx[i] + viscosity * element_info[local_id].dBasisdx[t][j][i];
    #else
                dmudx[i] = dmudx[i] + viscosity * global_dBasisdx[globalID].dBasisdx[t][j][i];
    #endif
#else
                dmudx[i] = dmudx[i] + viscosity * dBasisdx[j][i];
#endif
            }
        }
        // Stabilization parameters Tau and Delta
        delta = ZERO;
        tau = mk * native_divide(pow(hk, TWO),EIGHT * mu);
        
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
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                su[p][i][MDIM] = su[p][i][MDIM] + element_info[local_id].dBasisdx[t][p][i];
    #else
                su[p][i][MDIM] = su[p][i][MDIM] + global_dBasisdx[globalID].dBasisdx[t][p][i];
    #endif
#else
                su[p][i][MDIM] = su[p][i][MDIM] + dBasisdx[p][i];
#endif
                for (int j=0; j<MDIM; j++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                    su[p][i][i] = su[p][i][i] - dmudx[j] * element_info[local_id].dBasisdx[t][p][j];
                    su[p][i][j] = su[p][i][j] - dmudx[j] * element_info[local_id].dBasisdx[t][p][i];
    #else
                    su[p][i][i] = su[p][i][i] - dmudx[j] * global_dBasisdx[globalID].dBasisdx[t][p][j];
                    su[p][i][j] = su[p][i][j] - dmudx[j] * global_dBasisdx[globalID].dBasisdx[t][p][i];
    #endif
#else
                    su[p][i][i] = su[p][i][i] - dmudx[j] * dBasisdx[p][j];
                    su[p][i][j] = su[p][i][j] - dmudx[j] * dBasisdx[p][i];
#endif
                    sum = ZERO;
                    for (k=0; k<numberOfNodes; k++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                        sum = sum + dNodalBasisdx[p][k][j] * element_info[local_id].dBasisdx[t][k][j];
    #else
                        sum = sum + dNodalBasisdx[p][k][j] * global_dBasisdx[globalID].dBasisdx[t][k][j];
    #endif
#else
                        sum = sum + dNodalBasisdx[p][k][j] * dBasisdx[k][j];
#endif
                    }
                    su[p][i][i] = su[p][i][i] - mu * sum;
                    sum = ZERO;
                    for (k=0; k<numberOfNodes; k++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                        sum = sum + dNodalBasisdx[p][k][i] * element_info[local_id].dBasisdx[t][k][j];
    #else
                        sum = sum + dNodalBasisdx[p][k][i] * global_dBasisdx[globalID].dBasisdx[t][k][j];
    #endif
#else
                        sum = sum + dNodalBasisdx[p][k][i] * dBasisdx[k][j];
#endif
                    }
                    su[p][i][j] = su[p][i][j] - mu * sum;
                }
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                sw[p][MDIM][i] = sw[p][MDIM][i] + element_info[local_id].dBasisdx[t][p][i];
    #else
                sw[p][MDIM][i] = sw[p][MDIM][i] + global_dBasisdx[globalID].dBasisdx[t][p][i];
    #endif
#else
                sw[p][MDIM][i] = sw[p][MDIM][i] + dBasisdx[p][i];
#endif
                for (int j=0; j<MDIM; j++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                    sw[p][i][i] = sw[p][i][i] - dmudx[j] * element_info[local_id].dBasisdx[t][p][j];
                    sw[p][j][i] = sw[p][j][i] - dmudx[j] * element_info[local_id].dBasisdx[t][p][i];
    #else
                    sw[p][i][i] = sw[p][i][i] - dmudx[j] * global_dBasisdx[globalID].dBasisdx[t][p][j];
                    sw[p][j][i] = sw[p][j][i] - dmudx[j] * global_dBasisdx[globalID].dBasisdx[t][p][i];
    #endif
#else
                    sw[p][i][i] = sw[p][i][i] - dmudx[j] * dBasisdx[p][j];
                    sw[p][j][i] = sw[p][j][i] - dmudx[j] * dBasisdx[p][i];
#endif
                    sum = ZERO;
                    for (k=0; k<numberOfNodes; k++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                        sum = sum + dNodalBasisdx[p][k][j] * element_info[local_id].dBasisdx[t][k][j];
    #else
                        sum = sum + dNodalBasisdx[p][k][j] * global_dBasisdx[globalID].dBasisdx[t][k][j];
    #endif
#else
                        sum = sum + dNodalBasisdx[p][k][j] * dBasisdx[k][j];
#endif
                    }
                    sw[p][i][i] = sw[p][i][i] - mu * sum;
                    sum = ZERO;
                    for (k=0; k<numberOfNodes; k++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                        sum = sum + dNodalBasisdx[p][k][i] * element_info[local_id].dBasisdx[t][k][j];
    #else
                        sum = sum + dNodalBasisdx[p][k][i] * global_dBasisdx[globalID].dBasisdx[t][k][j];
    #endif
#else
                        sum = sum + dNodalBasisdx[p][k][i] * dBasisdx[k][j];
#endif
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
#ifdef TRANSIENT
#ifdef KERNEL_BASIS_DBASISDX
                    mass[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = mass[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * rho * basis[t][q] * basis[t][p];
#else
                    mass[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = mass[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * rho * basis[q] * basis[p];
#endif
#endif
                }
#ifdef NEWTONLINEAR
                if (viscNewtonLin) {
                    for (int i=0; i<MDIM; i++) {
                        sum = 0.0;
                        for (k=0; k<3; k++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                            sum = sum + strain[i][k] * element_info[local_id].dBasisdx[t][q][k];
    #else
                            sum = sum + strain[i][k] * global_dBasisdx[globalID].dBasisdx[t][q][k];
    #endif
#else
                            sum = sum + strain[i][k] * dBasisdx[q][k];
#endif
                        }
                        muder = muder0 * FOUR * sum;
                        for (int j=0; j<MDIM; j++) {
                            sum = 0.0;
                            for (k=0; k<3; k++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
                                sum = sum + strain[j][k] * element_info[local_id].dBasisdx[t][p][k];
    #else
                                sum = sum + strain[j][k] * global_dBasisdx[globalID].dBasisdx[t][p][k];
    #endif
#else
                                sum = sum + strain[j][k] * dBasisdx[p][k];
#endif
                            }
                            jacM[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = jacM[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * TWO * muder * sum;
                        }
                    }
                }
#endif
                
                for (int i=0; i<MDIM; i++) {
                    for (int j=0; j<MDIM; j++) {
                        // Isotropic and no Laplace discretization
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
        #ifdef GLOBAL_STIFF_FORCE
                        stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * mu * element_info[local_id].dBasisdx[t][q][j] * element_info[local_id].dBasisdx[t][p][j];
                        stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * mu * element_info[local_id].dBasisdx[t][q][i] * element_info[local_id].dBasisdx[t][p][j];
        #else
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * mu * element_info[local_id].dBasisdx[t][q][j] * element_info[local_id].dBasisdx[t][p][j];
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * mu * element_info[local_id].dBasisdx[t][q][i] * element_info[local_id].dBasisdx[t][p][j];
        #endif
    #else
        #ifdef GLOBAL_STIFF_FORCE
                        stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * mu * global_dBasisdx[globalID].dBasisdx[t][q][j] * global_dBasisdx[globalID].dBasisdx[t][p][j];
                        stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * mu * global_dBasisdx[globalID].dBasisdx[t][q][i] * global_dBasisdx[globalID].dBasisdx[t][p][j];
        #else
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * mu * global_dBasisdx[globalID].dBasisdx[t][q][j] * global_dBasisdx[globalID].dBasisdx[t][p][j];
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * mu * global_dBasisdx[globalID].dBasisdx[t][q][i] * global_dBasisdx[globalID].dBasisdx[t][p][j];
        #endif
    #endif
#else
    #ifdef GLOBAL_STIFF_FORCE
                        stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * mu * dBasisdx[q][j] * dBasisdx[p][j];
                        stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * mu * dBasisdx[q][i] * dBasisdx[p][j];
    #else
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+i] + s * mu * dBasisdx[q][j] * dBasisdx[p][j];
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * mu * dBasisdx[q][i] * dBasisdx[p][j];
    #endif
#endif
                    }
                    // Pressure term
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
        #ifdef GLOBAL_STIFF_FORCE
                    stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] - s * basis[t][q] * element_info[local_id].dBasisdx[t][p][i];
        #else
                    stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] - s * basis[t][q] * element_info[local_id].dBasisdx[t][p][i];
        #endif
    #else
        #ifdef GLOBAL_STIFF_FORCE
                    stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] - s * basis[t][q] * global_dBasisdx[globalID].dBasisdx[t][p][i];
        #else
                    stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] - s * basis[t][q] * global_dBasisdx[globalID].dBasisdx[t][p][i];
        #endif
    #endif
#else
    #ifdef GLOBAL_STIFF_FORCE
                    stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] - s * basis[q] * dBasisdx[p][i];
    #else
                    stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+MDIM] - s * basis[q] * dBasisdx[p][i];
    #endif
#endif
                    
                    // Continuity equation
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
        #ifdef GLOBAL_STIFF_FORCE
                    stiff_force[globalID].stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] = stiff_force[globalID].stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] + s * element_info[local_id].dBasisdx[t][q][i] * basis[t][p];
        #else
                    stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] + s * element_info[local_id].dBasisdx[t][q][i] * basis[t][p];
        #endif
    #else
        #ifdef GLOBAL_STIFF_FORCE
                    stiff_force[globalID].stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] = stiff_force[globalID].stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] + s * global_dBasisdx[globalID].dBasisdx[t][q][i] * basis[t][p];
        #else
                    stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] + s * global_dBasisdx[globalID].dBasisdx[t][q][i] * basis[t][p];
        #endif
    #endif
#else
    #ifdef GLOBAL_STIFF_FORCE
                    stiff_force[globalID].stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] = stiff_force[globalID].stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] + s * dBasisdx[q][i] * basis[p];
    #else
                    stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+MDIM][((MDIM+1)*q)+i] + s * dBasisdx[q][i] * basis[p];
    #endif
#endif
                }
                
                // Add stabilization
                for (int i=0; i<MDIM; i++) {
                    for (int j=0; j<MDIM+1; j++) {
#ifdef TRANSIENT
#ifdef KERNEL_BASIS_DBASISDX
                        mass[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = mass[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * tau * rho * basis[t][q] * sw[p][j][i];
#else
                        mass[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = mass[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * tau * rho * basis[q] * sw[p][j][i];
#endif
#endif
                    }
                    for (int j=0; j<MDIM; j++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef USE_GPU_LOCAL_MEM
        #ifdef GLOBAL_STIFF_FORCE
                        stiff_force[globalID].stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = stiff_force[globalID].stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * delta * element_info[local_id].dBasisdx[t][q][i] * element_info[local_id].dBasisdx[t][p][j];
        #else
                        stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * delta * element_info[local_id].dBasisdx[t][q][i] * element_info[local_id].dBasisdx[t][p][j];
        #endif
    #else
        #ifdef GLOBAL_STIFF_FORCE
                        stiff_force[globalID].stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = stiff_force[globalID].stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * delta * global_dBasisdx[globalID].dBasisdx[t][q][i] * global_dBasisdx[globalID].dBasisdx[t][p][j];
        #else
                        stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * delta * global_dBasisdx[globalID].dBasisdx[t][q][i] * global_dBasisdx[globalID].dBasisdx[t][p][j];
        #endif
    #endif
#else
    #ifdef GLOBAL_STIFF_FORCE
                        stiff_force[globalID].stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = stiff_force[globalID].stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * delta * dBasisdx[q][i] * dBasisdx[p][j];
    #else
                        stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] = stiff[((MDIM+1)*p)+j][((MDIM+1)*q)+i] + s * delta * dBasisdx[q][i] * dBasisdx[p][j];
    #endif
#endif
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
#ifdef GLOBAL_STIFF_FORCE
                        stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff_force[globalID].stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * tau * cc[i][j];
#else
                        stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] = stiff[((MDIM+1)*p)+i][((MDIM+1)*q)+j] + s * tau * cc[i][j];
#endif
                    }
                }
            }
        }
        
        // The right hand side...
        for (int p=0; p<nBasis; p++) {
            for (int i=0; i<MDIM+1; i++) {
#ifdef KERNEL_BASIS_DBASISDX
    #ifdef GLOBAL_STIFF_FORCE
                stiff_force[globalID].force[((MDIM+1)*p)+i] = stiff_force[globalID].force[((MDIM+1)*p)+i] + s * rho * force[i] * basis[t][p];
    #else
                forceVector[((MDIM+1)*p)+i] = forceVector[((MDIM+1)*p)+i] + s * rho * force[i] * basis[t][p];
    #endif
#else
    #ifdef GLOBAL_STIFF_FORCE
                stiff_force[globalID].force[((MDIM+1)*p)+i] = stiff_force[globalID].force[((MDIM+1)*p)+i] + s * rho * force[i] * basis[p];
    #else
                forceVector[((MDIM+1)*p)+i] = forceVector[((MDIM+1)*p)+i] + s * rho * force[i] * basis[p];
    #endif
#endif
            }
            // Add stabilization
            for (int i=0; i<MDIM; i++) {
                for (int j=0; j<MDIM+1; j++) {
#ifdef GLOBAL_STIFF_FORCE
                    stiff_force[globalID].force[((MDIM+1)*p)+j] = stiff_force[globalID].force[((MDIM+1)*p)+j] + s * tau * rho * force[i] * sw[p][j][i];
#else
                    forceVector[((MDIM+1)*p)+j] = forceVector[((MDIM+1)*p)+j] + s * tau * rho * force[i] * sw[p][j][i];
#endif
                }
            }
        }
    }

#ifdef NEWTONLINEAR
    if (viscNewtonLin) {
        REAL sol[32];
        for(int i=0; i<varDofs*numberOfNodes; i++) {
            sol[i] = ZERO;
        }
        
        int kk = 0;
        for (int i=0; i<(MDIM+1)*numberOfNodes; i+=(MDIM+1)) {
#ifdef USE_GPU_LOCAL_MEM
            sol[i] = element_info[local_id].vx[kk];
#else
    #ifdef NODAL_DATA
            sol[i] = nodal_data[globalID].vx[kk];
    #else
            sol[i] = varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+kk]]];
    #endif
#endif
            kk++;
        }
        kk = 0;
        for (int i=1; i<(MDIM+1)*numberOfNodes; i+=(MDIM+1)) {
#ifdef USE_GPU_LOCAL_MEM
            sol[i] = element_info[local_id].vy[kk];
#else
    #ifdef NODAL_DATA
            sol[i] = nodal_data[globalID].vy[kk];
    #else
            sol[i] = varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+kk]]+1];
    #endif
#endif
            kk++;
        }
        kk = 0;
        for (int i=2; i<(MDIM+1)*numberOfNodes; i+=(MDIM+1)) {
#ifdef USE_GPU_LOCAL_MEM
            sol[i] = element_info[local_id].vz[kk];
#else
    #ifdef NODAL_DATA
            sol[i] = nodal_data[globalID].vz[kk];
    #else
            sol[i] = varSolution[varDofs*varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+kk]]+2];
    #endif
#endif
            kk++;
        }
        
        int p = (MDIM+1) * nBasis;
        for (int i=0; i<p; i++) {
            for (int j=0; j<p; j++) {
#ifdef GLOBAL_STIFF_FORCE
                stiff_force[globalID].stiff[i][j] = stiff_force[globalID].stiff[i][j] + jacM[i][j];
#else
                stiff[i][j] = stiff[i][j] + jacM[i][j];
#endif
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
#ifdef GLOBAL_STIFF_FORCE
            stiff_force[globalID].force[i] = stiff_force[globalID].force[i] + yy[i];
#else
            forceVector[i] = forceVector[i] + yy[i];
#endif
        }
    }
#endif
    
#ifdef COLORING_ASSEMBLY // The coloring method does the final global and vector contribution here
    
    // The contribution of the local matrix to the global matrix
    // Only for DOF > 1
#ifdef PRECOMPUTE_NZ
    int c;
    col = 0;
    for (int i=0; i<numberElementDofs; i++) {
        for (k=1; k<=varDofs; k++) {
            for (int j=0; j<numberElementDofs; j++) {
                for (int l=1; l<=varDofs; l++) {
                    c = non_zeros[globalID].indexes[col];
#ifdef GLOBAL_STIFF_FORCE
                    values[c] = values[c] + stiff_force[globalID].stiff[varDofs*(i+1)-k][varDofs*(j+1)-l];
#else
                    values[c] = values[c] + stiff[varDofs*(i+1)-k][varDofs*(j+1)-l];
#endif
                    col++;
                }
            }
        }
    }
    
#else
    for (int i=0; i<numberElementDofs; i++) {
        for (k=1; k<=varDofs; k++) {
#ifdef CHECK_NEG_PERM
    #ifdef NODAL_DATA
            if (nodal_data[globalID].perm[i] < 0) continue;
    #else
            if (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]] < 0) continue;
    #endif
#endif

#ifdef NODAL_DATA
            row = varDofs * (nodal_data[globalID].perm[i] + 1) - k;
#else
            row = varDofs * (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]] + 1) - k;
#endif
            for (int j=0; j<numberElementDofs; j++) {
                for (int l=1; l<=varDofs; l++) {
#ifdef CHECK_NEG_PERM
    #ifdef NODAL_DATA
                    if (nodal_data[globalID].perm[j] < 0) continue;
    #else
                    if (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)]+j] < 0) continue;
    #endif
#endif
                    
#ifdef NODAL_DATA
                    col = varDofs * (nodal_data[globalID].perm[j] + 1) - l;
#else
                    col = varDofs * (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+j]] + 1) - l;
#endif
                    start = (col >= row) ? diag[row] : rows[row];
                    end = (col >= row) ? rows[row+1]-1 : diag[row]-1;
                    for (int c=start; c<=end; c++) {
                        if (cols[c] == col) {
#ifdef GLOBAL_STIFF_FORCE
                            values[c] = values[c] + stiff_force[globalID].stiff[varDofs*(i+1)-k][varDofs*(j+1)-l];
#else
                            values[c] = values[c] + stiff[varDofs*(i+1)-k][varDofs*(j+1)-l];
#endif
                            break;
                        }
                    }
                }
            }
        }
    }
#endif
    
    // The right-hand side
    for (int i=0; i<numberElementDofs; i++) {
#ifdef CHECK_NEG_PERM
    #ifdef NODAL_DATA
        if (nodal_data[globalID].perm[i] >= 0) {
    #else
        if (varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)]+i] >= 0) {
    #endif
#endif
            for (int j=0; j<varDofs; j++) {
#ifdef NODAL_DATA
                k = varDofs * nodal_data[globalID].perm[i] + j;
#else
                k = varDofs * varPermutation[elementNodeIndexesStore[(globalID*numberElementDofs)+i]] + j;
#endif
#ifdef GLOBAL_STIFF_FORCE
                rhs[k] = rhs[k] + stiff_force[globalID].force[varDofs*i+j];
#else
                rhs[k] = rhs[k] + forceVector[varDofs*i+j];
#endif
            }
#ifdef CHECK_NEG_PERM
        }
#endif
    }
        
#endif // End of the coloring method specific code
        
#ifdef NONZEROS_ASSEMBLY // If we do the assembly using the nonzeros method, we only need to store all computed local matrices and vectors
    int indx = globalID * ( (numberElementDofs*varDofs)*(numberElementDofs*varDofs) );
    for(int i=0; i<numberElementDofs*varDofs; i++) {
        for(int j=0; j<numberElementDofs*varDofs; j++) {
            element_stiffs[indx] = stiff[i][j];
            indx++;
        }
    }
    indx = globalID * (numberElementDofs*varDofs);
    for(int i=0; i<numberElementDofs*varDofs; i++) {
        element_forces[indx] = forceVector[i];
        indx++;
    }
        
#endif
        
}

__kernel void AssemblyGlobalMatrixByNonZeros(__global REAL *element_stiffs, __global int *reduction, __global REAL *values, int blockSize, int nzPerThread) {
    
    int i;
    
    int groupID = get_group_id(0);
    int localID = get_local_id(0);

    if (localID >= blockSize) return;
    
    int endBlock = blockSize * (nzPerThread*9);
    int id = (groupID * endBlock) + localID;
    
    REAL t = ZERO;
    int count = 0;
    while (count < endBlock) {
        i = reduction[id];
        
        if (i > 0) {
            t = t + element_stiffs[i-1];
        } else if (i < 0) {
            values[-i-1] = t;
            t = ZERO;
        }
        id = id + blockSize;
        count = count + blockSize;
    }
}

__kernel void AssemblyGlobalVectorByNonZeros(__global REAL *element_forces, __global int *reduction, __global REAL *rhs, int blockSize, int nzPerThread) {
    
    int i;
    
    int groupID = get_group_id(0);
    int localID = get_local_id(0);

    if (localID >= blockSize) return;
    
    int endBlock = blockSize * (nzPerThread*9);
    int id = (groupID * endBlock) + localID;
    
    REAL t = ZERO;
    int count = 0;
    while (count < endBlock) {
        i = reduction[id];
        
        if (i > 0) {
            t = t + element_forces[i-1];
        } else if (i < 0) {
            rhs[-i-1] = t;
            t = ZERO;
        }
        id = id + blockSize;
        count = count + blockSize;
    }
}
