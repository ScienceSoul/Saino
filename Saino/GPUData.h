//
//  GPUData.h
//  Saino
//
//  Created by Hakime Seddik on 26/04/2016.
//  Copyright © 2016 Institute of Low Temperature Science. All rights reserved.
//

#include <OpenCL/OpenCL.h>

// GPU data structure for ice flow simulations
typedef struct {
    
    double * __nullable basisFunctions_dp;
    double * __nullable nodesX_dp;
    double * __nullable nodesY_dp;
    double * __nullable nodesZ_dp;
    double * __nullable nodesvx_dp;
    double * __nullable nodesvy_dp;
    double * __nullable nodesvz_dp;
    double * __nullable varSol_dp;
    double * __nullable matValues_dp;
    double * __nullable matRHS_dp;
    double density_dp;
    double viscosity_dp;
    double gravity_dp;
    double hk_dp;
    double mk_dp;
    
    float * __nullable basisFunctions_sp;
    float * __nullable nodesX_sp;
    float * __nullable nodesY_sp;
    float * __nullable nodesZ_sp;
    float * __nullable nodesvx_sp;
    float * __nullable nodesvy_sp;
    float * __nullable nodesvz_sp;
    float * __nullable varSol_sp;
    float * __nullable matValues_sp;
    float * __nullable matRHS_sp;
    float density_sp;
    float viscosity_sp;
    float gravity_sp;
    float hk_sp;
    float mk_sp;
    
    // Memory locations to which the real variables
    // point according to the precision used
    void * __nullable basisFunctions_v;
    void * __nullable nodesX_v;
    void * __nullable nodesY_v;
    void * __nullable nodesZ_v;
    void * __nullable nodesvx_v;
    void * __nullable nodesvy_v;
    void * __nullable nodesvz_v;
    void * __nullable varSol_v;
    void * __nullable matValues_v;
    void * __nullable matRHS_v;
    void * __nullable density_v;
    void * __nullable viscosity_v;
    void * __nullable gravity_v;
    void * __nullable hk_v;
    void * __nullable mk_v;
    
    int * __nullable perm;
    
} ice_flow_gpu;

typedef struct {
    int p[8];
    int q[8];
    int r[8];
    float coeff[8];
} basis_functions_f;

typedef struct {
    int p[8];
    int q[8];
    int r[8];
    double coeff[8];
} basis_functions_d;

typedef struct {
    float nodes_x[8];
    float nodes_y[8];
    float nodes_z[8];
    float vx[8];
    float vy[8];
    float vz[8];
    float padding;
} element_info_f;

typedef struct {
    float nodes_x[8];
    float nodes_y[8];
    float nodes_z[8];
    float vx[8];
    float vy[8];
    float vz[8];
    float dBasisdx[8][8][3];
    float detJ[8];
    float padding;
} element_info_extended_f;

typedef struct {
    double nodes_x[8];
    double nodes_y[8];
    double nodes_z[8];
    double vx[8];
    double vy[8];
    double vz[8];
    double padding;
} element_info_d;

typedef struct {
    double nodes_x[8];
    double nodes_y[8];
    double nodes_z[8];
    double vx[8];
    double vy[8];
    double vz[8];
    double dBasisdx[8][8][3];
    double detJ[8];
    double padding;
} element_info_extended_d;

typedef struct {
    float basis[8][8];
} element_basis_f;

typedef struct {
    double basis[8][8];
} element_basis_d;

typedef struct {
    float dBasisdx[8][8][3];
    float detJ[8];
} element_dBasisdx_f;

typedef struct {
    double dBasisdx[8][8][3];
    double detJ[8];
} element_dBasisdx_d;

typedef struct {
    float nodes_x[8];
    float nodes_y[8];
    float nodes_z[8];
    float vx[8];
    float vy[8];
    float vz[8];
    int perm[8];
} nodal_data_f;

typedef struct {
    double nodes_x[8];
    double nodes_y[8];
    double nodes_z[8];
    double vx[8];
    double vy[8];
    double vz[8];
    int perm[8];
} nodal_data_d;

typedef struct {
    int indexes[1024];
} non_zero;

typedef struct {
    float stiff[32][32];
    float force[32];
} stiff_force_f;

typedef struct {
    double stiff[32][32];
    double force[32];
} stiff_force_d;

typedef struct {
    double * __nullable nodesX;
    double * __nullable nodesY;
    double * __nullable nodesZ;
    double * __nullable nodesvx;
    double * __nullable nodesvy;
    double * __nullable nodesvz;
    int * __nullable perm;
} bucket;

typedef struct {
    bucket * __nullable buckets;
    int numberOfBuckets;
    
    ice_flow_gpu * __nullable data;
    int * __nullable perm;
    
    cl_mem __nullable _nodal_nodes_x;
    cl_mem __nullable _nodal_nodes_y;
    cl_mem __nullable _nodal_nodes_z;
    cl_mem __nullable _nodal_vx;
    cl_mem __nullable _nodal_vy;
    cl_mem __nullable _nodal_vz;
    cl_mem __nullable _nodal_perm;
} colors_data;

void init_gpu_data(ice_flow_gpu * __nonnull data);
void __attribute__((overloadable)) init_basis_functions(basis_functions_f * __nonnull data);
void __attribute__((overloadable)) init_basis_functions(basis_functions_d * __nonnull data);

void __attribute__((overloadable)) init_nodal_data(nodal_data_f * __nonnull data, int size);
void __attribute__((overloadable)) init_nodal_data(nodal_data_d * __nonnull data, int size);

void init_nz_indexes(non_zero * __nonnull data, int size);

void __attribute__((overloadable)) init_stiff_force(stiff_force_f * __nonnull data, int size);
void __attribute__((overloadable)) init_stiff_force(stiff_force_d * __nonnull data, int size);

void init_color_data(colors_data * __nonnull data, int size);

