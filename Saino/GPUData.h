//===----------------------------------------------------------------------===//
//  GPUData.h
//  Saino
//
//  Created by Hakime Seddik on 26/04/2016.
//  Copyright © 2016 ScienceSoul. All rights reserved.
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
//===----------------------------------------------------------------------===//

#include <OpenCL/OpenCL.h>

// GPU data structure for ice flow simulations
typedef struct {
    
    double * _Nullable basisFunctions_dp;
    double * _Nullable nodesX_dp;
    double * _Nullable nodesY_dp;
    double * _Nullable nodesZ_dp;
    double * _Nullable nodesvx_dp;
    double * _Nullable nodesvy_dp;
    double * _Nullable nodesvz_dp;
    double * _Nullable varSol_dp;
    double * _Nullable matValues_dp;
    double * _Nullable matRHS_dp;
    double density_dp;
    double viscosity_dp;
    double gravity_dp;
    double hk_dp;
    double mk_dp;
    
    float * _Nullable basisFunctions_sp;
    float * _Nullable nodesX_sp;
    float * _Nullable nodesY_sp;
    float * _Nullable nodesZ_sp;
    float * _Nullable nodesvx_sp;
    float * _Nullable nodesvy_sp;
    float * _Nullable nodesvz_sp;
    float * _Nullable varSol_sp;
    float * _Nullable matValues_sp;
    float * _Nullable matRHS_sp;
    float density_sp;
    float viscosity_sp;
    float gravity_sp;
    float hk_sp;
    float mk_sp;
    
    // Memory locations to which the real variables
    // point according to the precision used
    void * _Nullable basisFunctions_v;
    void * _Nullable nodesX_v;
    void * _Nullable nodesY_v;
    void * _Nullable nodesZ_v;
    void * _Nullable nodesvx_v;
    void * _Nullable nodesvy_v;
    void * _Nullable nodesvz_v;
    void * _Nullable varSol_v;
    void * _Nullable matValues_v;
    void * _Nullable matRHS_v;
    void * _Nullable density_v;
    void * _Nullable viscosity_v;
    void * _Nullable gravity_v;
    void * _Nullable hk_v;
    void * _Nullable mk_v;
    
    int * _Nullable perm;
    
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
    double * _Nullable nodesX;
    double * _Nullable nodesY;
    double * _Nullable nodesZ;
    double * _Nullable nodesvx;
    double * _Nullable nodesvy;
    double * _Nullable nodesvz;
    int * _Nullable perm;
} bucket;

typedef struct {
    bucket * _Nullable buckets;
    int numberOfBuckets;
    
    ice_flow_gpu * _Nullable data;
    int * _Nullable perm;
    
    cl_mem _Nullable _nodal_nodes_x;
    cl_mem _Nullable _nodal_nodes_y;
    cl_mem _Nullable _nodal_nodes_z;
    cl_mem _Nullable _nodal_vx;
    cl_mem _Nullable _nodal_vy;
    cl_mem _Nullable _nodal_vz;
    cl_mem _Nullable _nodal_perm;
} colors_data;

void init_gpu_data(ice_flow_gpu * _Nonnull data);
void __attribute__((overloadable)) init_basis_functions(basis_functions_f * _Nonnull data);
void __attribute__((overloadable)) init_basis_functions(basis_functions_d * _Nonnull data);

void __attribute__((overloadable)) init_nodal_data(nodal_data_f * _Nonnull data, int size);
void __attribute__((overloadable)) init_nodal_data(nodal_data_d * _Nonnull data, int size);

void init_nz_indexes(non_zero * _Nonnull data, int size);

void __attribute__((overloadable)) init_stiff_force(stiff_force_f * _Nonnull data, int size);
void __attribute__((overloadable)) init_stiff_force(stiff_force_d * _Nonnull data, int size);

void init_color_data(colors_data * _Nonnull data, int size);

