//
//  GPUData.h
//  Saino
//
//  Created by Hakime Seddik on 26/04/2016.
//  Copyright © 2016 Institute of Low Temperature Science. All rights reserved.
//
    
// GPU data structure for ice flow simulations
typedef struct {
    
    double * __nullable basisFunctions_dp;
    double * __nullable nodesUVW_dp;
    double * __nullable gaussPoints_dp;
    double * __nullable nodesX_dp;
    double * __nullable nodesY_dp;
    double * __nullable nodesZ_dp;
    double * __nullable varSol_dp;
    double * __nullable matValues_dp;
    double * __nullable matRHS_dp;
    double density_dp;
    double viscosity_dp;
    double gravity_dp;
    double hk_dp;
    double mk_dp;
    
    float * __nullable basisFunctions_sp;
    float * __nullable nodesUVW_sp;
    float * __nullable gaussPoints_sp;
    float * __nullable nodesX_sp;
    float * __nullable nodesY_sp;
    float * __nullable nodesZ_sp;
    float * __nullable varSol_sp;
    float * __nullable matValues_sp;
    float * __nullable matRHS_sp;
    float density_sp;
    float viscosity_sp;
    float gravity_sp;
    float hk_sp;
    float mk_sp;
    
    // Memory locations to which the above variable
    // point according to the type used
    void * __nullable basisFunctions_v;
    void * __nullable nodesUVW_v;
    void * __nullable gaussPoints_v;
    void * __nullable nodesX_v;
    void * __nullable nodesY_v;
    void * __nullable nodesZ_v;
    void * __nullable varSol_v;
    void * __nullable matValues_v;
    void * __nullable matRHS_v;
    void * __nullable density_v;
    void * __nullable viscosity_v;
    void * __nullable gravity_v;
    void * __nullable hk_v;
    void * __nullable mk_v;
    
} ice_flow_gpu;

void init_gpu_data(ice_flow_gpu * __nonnull data);