//
//  GPUData.c
//  Saino
//
//  Created by Hakime Seddik on 02/06/2016.
//  Copyright © 2016 Institute of Low Temperature Science. All rights reserved.
//

#include <stdio.h>
#include <string.h>

#include "GPUData.h"

void init_gpu_data(ice_flow_gpu * __nonnull data) {
    
    *data = (ice_flow_gpu){
        .basisFunctions_dp=NULL,
        .nodesX_dp=NULL,
        .nodesY_dp=NULL,
        .nodesZ_dp=NULL,
        .varSol_dp=NULL,
        .matValues_dp=NULL,
        .matRHS_dp=NULL,
        
        .basisFunctions_sp=NULL,
        .nodesX_sp=NULL,
        .nodesY_sp=NULL,
        .nodesZ_sp=NULL,
        .varSol_sp=NULL,
        .matValues_sp=NULL,
        .matRHS_sp=NULL,
        
        .basisFunctions_v=NULL,
        .nodesX_v=NULL,
        .nodesY_v=NULL,
        .nodesZ_v=NULL,
        .varSol_v=NULL,
        .matValues_v=NULL,
        .matRHS_v=NULL,
        .density_v = NULL,
        .viscosity_v=NULL,
        .gravity_v=NULL,
        .hk_v=NULL,
        .mk_v=NULL,
        
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

void __attribute__((overloadable)) init_basis_functions(basis_functions_f * __nonnull data) {
    
    *data = (basis_functions_f) {
        .p = {0, 0, 0, 0, 0, 0, 0, 0},
        .q = {0, 0, 0, 0, 0, 0, 0, 0},
        .r = {0, 0, 0, 0, 0, 0, 0, 0},
        .coeff = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
    };
}

void __attribute__((overloadable)) init_basis_functions(basis_functions_d * __nonnull data) {
    
    *data = (basis_functions_d) {
        .p = {0, 0, 0, 0, 0, 0, 0, 0},
        .q = {0, 0, 0, 0, 0, 0, 0, 0},
        .r = {0, 0, 0, 0, 0, 0, 0, 0},
        .coeff = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    };
}

void __attribute__((overloadable)) init_nodal_data(nodal_data_f * __nonnull data, int size) {
    
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

void __attribute__((overloadable)) init_nodal_data(nodal_data_d * __nonnull data, int size) {
    
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


