//===----------------------------------------------------------------------===//
//  GPUData.c
//  Saino
//
//  Created by Hakime Seddik on 02/06/2016.
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

#include <stdio.h>
#include <string.h>
#include "GPUData.h"

void init_gpu_data(ice_flow_gpu * _Nonnull data) {
    
    *data = (ice_flow_gpu){
        .basisFunctions_dp=NULL,
        .nodesX_dp=NULL,
        .nodesY_dp=NULL,
        .nodesZ_dp=NULL,
        .nodesvx_dp=NULL,
        .nodesvy_dp=NULL,
        .nodesvz_dp=NULL,
        .varSol_dp=NULL,
        .matValues_dp=NULL,
        .matRHS_dp=NULL,
        
        .basisFunctions_sp=NULL,
        .nodesX_sp=NULL,
        .nodesY_sp=NULL,
        .nodesZ_sp=NULL,
        .nodesvx_sp=NULL,
        .nodesvy_sp=NULL,
        .nodesvz_sp=NULL,
        .varSol_sp=NULL,
        .matValues_sp=NULL,
        .matRHS_sp=NULL,
        
        .basisFunctions_v=NULL,
        .nodesX_v=NULL,
        .nodesY_v=NULL,
        .nodesZ_v=NULL,
        .nodesvx_v=NULL,
        .nodesvy_v=NULL,
        .nodesvz_v=NULL,
        .varSol_v=NULL,
        .matValues_v=NULL,
        .matRHS_v=NULL,
        .density_v = NULL,
        .viscosity_v=NULL,
        .gravity_v=NULL,
        .hk_v=NULL,
        .mk_v=NULL,
        
        .perm=NULL,
        
        .density_dp = 0.0,
        .viscosity_dp=0.0,
        .gravity_dp=0.0,
        .hk_dp=0.0,
        .mk_dp=0.0,
        
        .density_sp = 0.0f,
        .viscosity_sp=0.0f,
        .gravity_sp=0.0f,
        .hk_sp=0.0f,
        .mk_sp=0.0f,
    };
}

void __attribute__((overloadable)) init_basis_functions(basis_functions_f * _Nonnull data) {
    
    *data = (basis_functions_f) {
        .p = {0, 0, 0, 0, 0, 0, 0, 0},
        .q = {0, 0, 0, 0, 0, 0, 0, 0},
        .r = {0, 0, 0, 0, 0, 0, 0, 0},
        .coeff = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
    };
}

void __attribute__((overloadable)) init_basis_functions(basis_functions_d * _Nonnull data) {
    
    *data = (basis_functions_d) {
        .p = {0, 0, 0, 0, 0, 0, 0, 0},
        .q = {0, 0, 0, 0, 0, 0, 0, 0},
        .r = {0, 0, 0, 0, 0, 0, 0, 0},
        .coeff = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    };
}

void __attribute__((overloadable)) init_nodal_data(nodal_data_f * _Nonnull data, int size) {
    
    for (int i=0; i<size; i++) {
        memset(data[i].nodes_x, 0.0f, sizeof(data[i].nodes_x) );
        memset(data[i].nodes_y, 0.0f, sizeof(data[i].nodes_y) );
        memset(data[i].nodes_z, 0.0f, sizeof(data[i].nodes_z) );
        
        memset(data[i].vx, 0.0f, sizeof(data[i].vx) );
        memset(data[i].vy, 0.0f, sizeof(data[i].vy) );
        memset(data[i].vz, 0.0f, sizeof(data[i].vz) );
        
        memset(data[i].perm, -1, sizeof(data[i].perm) );
    }
}

void __attribute__((overloadable)) init_nodal_data(nodal_data_d * _Nonnull data, int size) {
    
    for (int i=0; i<size; i++) {
        memset(data[i].nodes_x, 0.0, sizeof(data[i].nodes_x) );
        memset(data[i].nodes_y, 0.0, sizeof(data[i].nodes_y) );
        memset(data[i].nodes_z, 0.0, sizeof(data[i].nodes_z) );
        
        memset(data[i].vx, 0.0, sizeof(data[i].vx) );
        memset(data[i].vy, 0.0, sizeof(data[i].vy) );
        memset(data[i].vz, 0.0, sizeof(data[i].vz) );
        
        memset(data[i].perm, -1, sizeof(data[i].perm) );
    }
}

void init_nz_indexes(non_zero * _Nonnull data, int size) {
    
    for (int i=0; i<size; i++) {
        memset(data[i].indexes, -1, sizeof(data[i].indexes) );
    }
}

void __attribute__((overloadable)) init_stiff_force(stiff_force_f * _Nonnull data, int size) {
    for (int i=0; i<size; i++) {
        memset(*(data[i].stiff), 0.0f, sizeof(float)*(32*32));
        memset(data[i].force, 0.0f, sizeof(data[i].force));
    }
}

void __attribute__((overloadable)) init_stiff_force(stiff_force_d * _Nonnull data, int size) {
    for (int i=0; i<size; i++) {
        memset(*(data[i].stiff), 0.0, sizeof(double)*(32*32));
        memset(data[i].force, 0.0, sizeof(data[i].force));
    }
}

void init_color_data(colors_data * _Nonnull data, int size) {
    for (int i=0; i<size; i++) {
        data[i].buckets = NULL;
        data[i].numberOfBuckets = 0;
        data[i].data = NULL;
        data[i].perm = NULL;
        
        data[i]._nodal_nodes_x = NULL;
        data[i]._nodal_nodes_y = NULL;
        data[i]._nodal_nodes_z = NULL;
        data[i]._nodal_vx = NULL;
        data[i]._nodal_vy = NULL;
        data[i]._nodal_vz = NULL;
        data[i]._nodal_perm = NULL;
    }
}

