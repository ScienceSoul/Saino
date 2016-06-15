//
//  GPUData.c
//  Saino
//
//  Created by Hakime Seddik on 02/06/2016.
//  Copyright © 2016 Institute of Low Temperature Science. All rights reserved.
//

#include <stdio.h>

#include "GPUData.h"

void init_gpu_data(ice_flow_gpu * __nonnull data) {
    
    *data = (ice_flow_gpu){
        .basisFunctions_dp=NULL,
        .nodesUVW_dp=NULL,
        .gaussPoints_dp=NULL,
        .nodesX_dp=NULL,
        .nodesY_dp=NULL,
        .nodesZ_dp=NULL,
        .varSol_dp=NULL,
        .matValues_dp=NULL,
        .matRHS_dp=NULL,
        
        .basisFunctions_sp=NULL,
        .nodesUVW_sp=NULL,
        .gaussPoints_sp=NULL,
        .nodesX_sp=NULL,
        .nodesY_sp=NULL,
        .nodesZ_sp=NULL,
        .varSol_sp=NULL,
        .matValues_sp=NULL,
        .matRHS_sp=NULL,
        
        .basisFunctions_v=NULL,
        .nodesUVW_v=NULL,
        .gaussPoints_v=NULL,
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

