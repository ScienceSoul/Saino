//
//  GPUData.h
//  Saino
//
//  Created by Hakime Seddik on 26/04/2016.
//  Copyright © 2016 Institute of Low Temperature Science. All rights reserved.
//
    
// GPU data structure for ice flow simulations
typedef struct {
    double * __nullable basisFunctions;
    double * __nullable nodesUVW;
    double * __nullable gaussPoints;
    double density;
    double viscosity;
    double gravity;
    double hk;
    double mk;
} ice_flow_gpu;