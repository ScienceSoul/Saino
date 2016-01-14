//
//  FEMDifferentials.m
//  Saino
//
//  Created by Seddik hakime on 18/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMDifferentials.h"
#import "FEMListUtilities.h"
#import "FEMUtilities.h"
#import "FEMCoordinateSystems.h"
#import "FEMLinearAlgebra.h"

@interface FEMDifferentials ()
-(void)FEMDifferentials_computeLorentzMagnetic:(double * __nonnull)b dhdx:(double[][3])dhdx permeability:(double)mu sqrtMetric:(double)sqrtMetric metric:(double[][3])metric symbols:(double[][3][3])symb lorentzForce:(double * __nonnull)lorentzForce model:(FEMModel * __nonnull)model;
-(double)FEMDifferentials_computeMagneticHeat:(double[3])b dhdx:(double[][3])dhdx permeability:(double)mu sqrtMetric:(double)sqrtMetric metric:(double[][3])metric symbols:(double[][3][3])symb model:(FEMModel * __nonnull)model;
@end

@implementation FEMDifferentials

#pragma mark Private methods

-(void)FEMDifferentials_computeLorentzMagnetic:(double * __nonnull)b dhdx:(double[][3])dhdx permeability:(double)mu sqrtMetric:(double)sqrtMetric metric:(double[][3])metric symbols:(double[][3][3])symb lorentzForce:(double * __nonnull)lorentzForce model:(FEMModel * __nonnull)model {
    
    int i, j, k, l, m;
    double bc[3], jc[3], ji[3], perm[3][3][3], r, s;
    
    if (model.coordinates == cartesian) {
        ji[0] = dhdx[2][1] - dhdx[1][2];
        ji[1] = dhdx[0][2] - dhdx[2][0];
        ji[2] = dhdx[1][0] - dhdx[0][1];
        lorentzForce[0] = ji[1]*b[2] - ji[2]*b[1];
        lorentzForce[1] = ji[2]*b[0] - ji[0]*b[2];
        lorentzForce[2] = ji[0]*b[1] - ji[1]*b[0];
        return;
    }
    
    r = sqrtMetric;
    
    if (model.coordinates == cylindric_symmetric) {
        ji[0] = -dhdx[2][1];
        ji[1] = dhdx[2][0];
        if (r > 1.0e-10) {
            ji[1] = ji[1] + b[2] / (r*mu);
        } else {
            ji[1] = ji[1] + ji[1];
        }
        ji[2] = dhdx[0][1] - dhdx[1][0];
        lorentzForce[0] = ji[2]*b[1] - ji[1]*b[2];
        lorentzForce[1] = ji[0]*b[2] - ji[2]*b[0];
        
        // You might want to use SI units for the azimuthal component,
        // if you compute Lorentz force at nodal points and symmetry axis,
        // otherwise you divide by zero.
        
        if (r > 1.0e-10) {
            lorentzForce[2] = (ji[1]*b[0] - ji[0]*b[1]) / r;
        } else {
            lorentzForce[2] = 0.0;
        }
        return;
    }
    
    memset( **perm, 0.0, (3*3*3)*sizeof(double) );
    perm[0][1][2] = -1.0 / sqrtMetric;
    perm[0][2][1] =  1.0 / sqrtMetric;
    perm[1][0][2] =  1.0 / sqrtMetric;
    perm[1][2][0] = -1.0 / sqrtMetric;
    perm[2][0][1] = -1.0 / sqrtMetric;
    perm[2][1][0] =  1.0 / sqrtMetric;
    
    memset( bc, 0.0, sizeof(bc) );
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            bc[i] = bc[i] + metric[i][j]*b[j];
        }
    }

    memset( ji, 0.0, sizeof(ji) );
    for (i=0; i<3; i++) {
        s = 0.0;
        for (j=0; j<3; j++) {
            for (k=0; k<3; k++) {
                if (perm[i][j][k] != 0.0) {
                    for (l=0; l<3; l++) {
                        s = s + perm[i][j][k]*metric[j][l]*dhdx[l][k];
                    }
                    for (m=0; m<3; m++) {
                        s = s + perm[i][j][k]*metric[j][l]*symb[k][m][l]*b[m]/mu;
                    }
                }
            }
        }
        ji[i] = s;
    }
    
    memset( jc, 0.0, sizeof(jc) );
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            jc[i] = jc[i] + metric[i][j]*ji[j];
        }
    }
    
    for (i=0; i<3; i++) {
        s = 0.0;
        for (j=0; j<3; j++) {
            for (k=0; k<3; k++) {
                if (perm[i][j][k] != 0) {
                    s = s + perm[i][j][k] * jc[k] * bc[j];
                }
            }
        }
        lorentzForce[i] = s;
    }
}

-(double)FEMDifferentials_computeMagneticHeat:(double[3])b dhdx:(double[][3])dhdx permeability:(double)mu sqrtMetric:(double)sqrtMetric metric:(double[][3])metric symbols:(double[][3][3])symb model:(FEMModel * __nonnull)model {
    
    int i, j, k, l, m;
    double bc[3], jc[3], jh, ji[3], perm[3][3][3], r, s;
    
    if (model.coordinates == cartesian) {
        ji[0] = dhdx[2][1] - dhdx[1][2];
        ji[1] = dhdx[0][2] - dhdx[2][0];
        ji[2] = dhdx[1][0] - dhdx[0][1];
        jh = ji[0]*ji[0] + ji[1]*ji[1] + ji[2]*ji[2];
        return jh;
    }
    
    if (model.coordinates == cylindric_symmetric) {
        r = sqrtMetric;
        ji[0] = -dhdx[2][1];
        ji[1] = b[2]/(r*mu) + dhdx[2][0];
        ji[2] = dhdx[0][1] - dhdx[1][0];
        jh = ji[0]*ji[0] + ji[1]*ji[1] + ji[2]*ji[2];
        return jh;
    }
    
    memset( **perm, 0.0, (3*3*3)*sizeof(double) );
    perm[0][1][2] = -1.0 / sqrtMetric;
    perm[0][2][1] =  1.0 / sqrtMetric;
    perm[1][0][2] =  1.0 / sqrtMetric;
    perm[1][2][0] = -1.0 / sqrtMetric;
    perm[2][0][1] = -1.0 / sqrtMetric;
    perm[2][1][0] =  1.0 / sqrtMetric;
    
    memset( bc, 0.0, sizeof(bc) );
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            bc[i] = bc[i] + metric[i][j]*b[j];
        }
    }
    
    memset( ji, 0.0, sizeof(ji) );
    for (i=0; i<3; i++) {
        s = 0.0;
        for (j=0; j<3; j++) {
            for (k=0; k<3; k++) {
                if (perm[i][j][k] != 0.0) {
                    for (l=0; l<3; l++) {
                        s = s + perm[i][j][k]*metric[j][l]*dhdx[l][k];
                    }
                    for (m=0; m<3; m++) {
                        s = s + perm[i][j][k]*metric[j][l]*symb[k][m][l]*b[m]/mu;
                    }
                }
            }
        }
        ji[i] = s;
    }
    
    memset( jc, 0.0, sizeof(jc) );
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            jc[i] = jc[i] + metric[i][j]*ji[j];
        }
    }
    
    jh = 0.0;
    for (i=0; i<3; i++) {
        jh = jh + ji[i] * jc[i];
    }
    return jh;
}

#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
    }
    
    return self;
}


/*************************************************************************************************************
 
    Computes the Lorentz force resulting from a magnetic field at integration point (u,v,w).
 
**************************************************************************************************************/
-(void)lorentzForceElement:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int)n integrationU:(double)u integrationV:(double)v integrationW:(double)w lorentzForce:(double * __nonnull)lorentzForce mesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration coordinateSystems:(FEMCoordinateSystems * __nonnull)coordinateSystems listUtilities:(FEMListUtilities * __nonnull)listUtilities utilities:(FEMUtilities * __nonnull)utilities {
    
    int i, j, k;
    double mu, sqrtElementMetric, sqrtMetric, x, y, z;
    double B[3], dHdx[3][3], dSymb[3][3][3][3], ExtMx[n], ExtMy[n], ExtMz[n], metric[3][3], symb[3][3][3];
    double *basis = NULL, **basisFirstDerivative = NULL;
    FEMVariable *mx, *my, *mz;
    FEMVariable *mfx, *mfy, *mfz;
    FEMMaterial *materialAtID = nil;
    variableArraysContainer *mxContainers = NULL, *myContainers = NULL, *mzContainers = NULL;
    variableArraysContainer *mfxContainers = NULL, *mfyContainers = NULL, *mfzContainers = NULL;
    listBuffer permeability = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer extMx = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer extMy = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer extMz = { NULL, NULL, NULL, NULL, 0, 0, 0};
    BOOL any, found, stat;
    
    mx = [utilities getVariableFrom:model.variables model:model name:@"magnetic field 1" onlySearch:NULL maskName:nil info:&found];
    my = [utilities getVariableFrom:model.variables model:model name:@"magnetic field 2" onlySearch:NULL maskName:nil info:&found];
    mz = [utilities getVariableFrom:model.variables model:model name:@"magnetic field 3" onlySearch:NULL maskName:nil info:&found];
    if (mx == nil) return;
    
    mxContainers = mx.getContainers;
    myContainers = my.getContainers;
    mzContainers = mz.getContainers;
    
    basis = integration.basis;
    basisFirstDerivative = integration.basisFirstDerivative;
    
    any = NO;
    for (i=0; i<mxContainers->sizePerm; i++) {
        if (mxContainers->Perm[i] < 0) {
            any = YES;
            break;
        }
    }
    if (any == YES) return;
    
    k = [(model.bodies)[element->BodyID-1][@"material"] intValue];
    if (k < 1) k = 1;
    if (k > model.numberOfMaterials) k = model.numberOfMaterials;
    materialAtID = (model.materials)[k-1];
    
    found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"magnetic permeability" numberOfNodes:n indexes:element->NodeIndexes buffer:&permeability minValue:NULL maxValue:NULL];
    found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"applied magnetic field 1" numberOfNodes:n indexes:element->NodeIndexes buffer:&extMx minValue:NULL maxValue:NULL];
    found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"applied magnetic field 2" numberOfNodes:n indexes:element->NodeIndexes buffer:&extMy minValue:NULL maxValue:NULL];
    found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"applied magnetic field 3" numberOfNodes:n indexes:element->NodeIndexes buffer:&extMz minValue:NULL maxValue:NULL];
    
    mfx = [utilities getVariableFrom:model.variables model:model name:@"magnetic flux density 1" onlySearch:NULL maskName:nil info:&found];
    mfy = [utilities getVariableFrom:model.variables model:model name:@"magnetic flux density 2" onlySearch:NULL maskName:nil info:&found];
    mfz = [utilities getVariableFrom:model.variables model:model name:@"magnetic flux density 3" onlySearch:NULL maskName:nil info:&found];
    if (mfx != nil) {
        mfxContainers = mfx.getContainers;
        mfyContainers = mfy.getContainers;
        mfzContainers = mfz.getContainers;
        for (i=0; i<n; i++) {
            ExtMx[i] = extMx.vector[i] + mfxContainers->Values[mfxContainers->Perm[element->NodeIndexes[i]]];
            ExtMy[i] = extMy.vector[i] + mfyContainers->Values[mfyContainers->Perm[element->NodeIndexes[i]]];
            ExtMz[i] = extMz.vector[i] + mfzContainers->Values[mfzContainers->Perm[element->NodeIndexes[i]]];
        }
    }
    
    // Get element info
    stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
    stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
    stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
    sqrtElementMetric = integration.metricDeterminant;

    memset(B, 0.0, sizeof(B) );
    for (i=0; i<n; i++) {
        B[0] = B[0] + basis[i]*mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
        B[1] = B[1] + basis[i]*myContainers->Values[myContainers->Perm[element->NodeIndexes[i]]];
        B[2] = B[2] + basis[i]*mzContainers->Values[mzContainers->Perm[element->NodeIndexes[i]]];
    }
    
    B[0] = B[0] + cblas_ddot(n, ExtMx, 1, basis, 1);
    B[1] = B[1] + cblas_ddot(n, ExtMy, 1, basis, 1);
    B[2] = B[2] + cblas_ddot(n, ExtMz, 1, basis, 1);
    
    memset( *dHdx, 0.0, (3*3)*sizeof(double) );
    for (i=0; i<3; i++) {
        for (j=0; j<n; j++) {
            dHdx[0][i] = dHdx[0][i] + basisFirstDerivative[j][i]*mxContainers->Values[mxContainers->Perm[element->NodeIndexes[j]]]/permeability.vector[j];
            dHdx[1][i] = dHdx[1][i] + basisFirstDerivative[j][i]*myContainers->Values[myContainers->Perm[element->NodeIndexes[j]]]/permeability.vector[j];
            dHdx[2][i] = dHdx[2][i] + basisFirstDerivative[j][i]*mzContainers->Values[mzContainers->Perm[element->NodeIndexes[j]]]/permeability.vector[j];
        }
    }
    
    // Get coordinate system info
    x = cblas_ddot(n, nodes->x, 1, basis, 1);
    y = cblas_ddot(n, nodes->y, 1, basis, 1);
    z = cblas_ddot(n, nodes->z, 1, basis, 1);
    [coordinateSystems coordinateSystemInfoModel:model metric:metric sqrtMetric:&sqrtMetric symbols:symb dSymbols:dSymb coordX:x coordY:y coordZ:z];
    if (model.coordinates != cartesian) {
        FEMLinearAlgebra *linearAlgebra = [[FEMLinearAlgebra alloc] init];
        double **matrix = doublematrix(0, 2, 0, 2);
        memcpy(*matrix, *metric, (3*3)*sizeof(double));
        [linearAlgebra invertMatrix:matrix ofSize:3];
        memcpy(*metric, *matrix, (3*3)*sizeof(double));
        free_dmatrix(matrix, 0, 2, 0, 2);
    }
    
    mu = cblas_ddot(n, permeability.vector, 1, basis, 1);
    [self FEMDifferentials_computeLorentzMagnetic:B dhdx:dHdx permeability:mu sqrtMetric:sqrtMetric metric:metric symbols:symb lorentzForce:lorentzForce model:model];
    
    if (permeability.vector != NULL) {
        free_dvector(permeability.vector, 0, permeability.m-1);
        permeability.vector = NULL;
    }
    if (extMx.vector != NULL) {
        free_dvector(extMx.vector, 0, extMx.m-1);
        extMx.vector = NULL;
    }
    if (extMy.vector != NULL) {
        free_dvector(extMy.vector, 0, extMy.m-1);
        extMy.vector = NULL;
    }
    if (extMz.vector != NULL) {
        free_dvector(extMz.vector, 0, extMz.m-1);
        extMz.vector = NULL;
    }
}


/*************************************************************************************************
    Compute the Joule heating at integration point (u, v, w) given the appropriate electrostatic
    or magnetic field that indicates the current through a conductor.
*************************************************************************************************/
-(double)jouleHeatElement:(Element_t * __nonnull)element nodes:(Nodes_t * __nonnull)nodes numberOfNodes:(int)n integrationU:(double)u integrationV:(double)v integrationW:(double)w mesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration listUtilities:(FEMListUtilities * __nonnull)listUtilities {
 
    int i, j, k, bf_id, jouleNode;
    static int prevElementBodyID = -1;
    double b[3], dhdx[3][3], dSymb[3][3][3][3], elcond, jouleheat, metric[3][3], sqrtElementMetric, sqrtMetric=0.0, sum, symb[3][3][3], x, y, z;
    BOOL all, found, stat;
    static BOOL jouleHeat = NO;
    FEMBodyForce *bodyForceAtID;
    FEMMaterial *materialAtID;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    variableArraysContainer *mxContainers = NULL;
    
    jouleheat = 0.0;
    
    if (element->BodyID-1 != prevElementBodyID) {
        prevElementBodyID = element->BodyID-1;
        if ((model.bodies)[element->BodyID-1][@"body force"] != nil) {
            bf_id = [(model.bodies)[element->BodyID-1][@"body force"] intValue];
        } else {
            return jouleheat;
        }
        bodyForceAtID = (model.bodyForces)[bf_id-1];
        jouleHeat = [listUtilities listGetLogical:model inArray:bodyForceAtID.valuesList forVariable:@"joule heat" info:&found];
    }
    if (jouleHeat == NO) return jouleheat;
    
    jouleNode = 0;
    if (jouleNode == 0) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *mx = [utilities getVariableFrom:model.variables model:model name:@"joule field" onlySearch:NULL maskName:nil info:&found];
        if (found == YES) {
            mxContainers = mx.getContainers;
            all = YES;
            for (i=0; i<n; i++) {
                if (mxContainers->Perm[element->NodeIndexes[i]] < 0) {
                    all = NO;
                    break;
                }
            }
            if (all == YES) jouleNode = 2;
        }
    }
       
    if (jouleNode == 0) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *mx = [utilities getVariableFrom:model.variables model:model name:@"magnetic field 1" onlySearch:NULL maskName:nil info:&found];
        if (found == YES) {
            mxContainers = mx.getContainers;
            all = YES;
            for (i=0; i<n; i++) {
                if (mxContainers->Perm[element->NodeIndexes[i]] < 0) {
                    all = NO;
                    break;
                }
            }
            if (all == YES) jouleNode = 1;

        }
    }
    
    if (jouleNode == 0) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *mx = [utilities getVariableFrom:model.variables model:model name:@"potential" onlySearch:NULL maskName:nil info:&found];
        if (found == YES) {
            mxContainers = mx.getContainers;
            all = YES;
            for (i=0; i<n; i++) {
                if (mxContainers->Perm[element->NodeIndexes[i]] < 0) {
                    all = NO;
                    break;
                }
            }
            if (all == YES) jouleNode = 3;
        }
    }
    
    if (jouleNode == 0) return jouleheat;
    
    // Get element info
    stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
    stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
    stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
    sqrtElementMetric = integration.metricDeterminant;
    
    if (model.coordinates != cartesian) {
        x = cblas_ddot(n, nodes->x, 1, integration.basis, 1);
        y = cblas_ddot(n, nodes->y, 1, integration.basis, 1);
        z = cblas_ddot(n, nodes->z, 1, integration.basis, 1);
        FEMCoordinateSystems *coordinateSystems = [[FEMCoordinateSystems alloc] init];
        [coordinateSystems coordinateSystemInfoModel:model metric:metric sqrtMetric:&sqrtMetric symbols:symb dSymbols:dSymb coordX:x coordY:y coordZ:z];
        FEMLinearAlgebra *linearAlgebra = [[FEMLinearAlgebra alloc] init];
        double **matrix = doublematrix(0, 2, 0, 2);
        memcpy(*matrix, *metric, (3*3)*sizeof(double));
        [linearAlgebra invertMatrix:matrix ofSize:3];
        memcpy(*metric, *matrix, (3*3)*sizeof(double));
        free_dmatrix(matrix, 0, 2, 0, 2);
    }
    
    // All models need electric conductivity
    k = [(model.bodies)[element->BodyID-1][@"material"] intValue];
    materialAtID = (model.materials)[k-1];
    
    found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"electrical conductivity" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
    if (found == YES) {
        NSLog(@"FEMDifferentials:jouleHeatElement: use electrical conductivity instead of electric.\n");
    } else {
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"electric conductivity" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
    }
    elcond = 0.0;
    if (found == YES) {
        elcond = cblas_ddot(n, buffer.vector, 1, integration.basis, 1);
    }
    
    // Joule heating takes place only when there is ohmic resistance
    if (elcond < DBL_MIN) return jouleheat;
    
    // Magnetic induction equation
    if (jouleheat == 1) {
        listBuffer permeability = { NULL, NULL, NULL, NULL, 0, 0, 0};
        listBuffer extmx = { NULL, NULL, NULL, NULL, 0, 0, 0};
        listBuffer extmy = { NULL, NULL, NULL, NULL, 0, 0, 0};
        listBuffer extmz = { NULL, NULL, NULL, NULL, 0, 0, 0};
        
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"magnetic permeability" numberOfNodes:n indexes:element->NodeIndexes buffer:&permeability minValue:NULL maxValue:NULL];
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        FEMVariable *mx = [utilities getVariableFrom:model.variables model:model name:@"magnetic field 1" onlySearch:NULL maskName:nil info:&found];
        FEMVariable *my = [utilities getVariableFrom:model.variables model:model name:@"magnetic field 2" onlySearch:NULL maskName:nil info:&found];
        FEMVariable *mz = [utilities getVariableFrom:model.variables model:model name:@"magnetic field 3" onlySearch:NULL maskName:nil info:&found];
        variableArraysContainer *mxContainers = mx.getContainers;
        variableArraysContainer *myContainers = my.getContainers;
        variableArraysContainer *mzContainers = mz.getContainers;
        
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"applied magnetic field 1" numberOfNodes:n indexes:element->NodeIndexes buffer:&extmx minValue:NULL maxValue:NULL];
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"applied magnetic field 2" numberOfNodes:n indexes:element->NodeIndexes buffer:&extmy minValue:NULL maxValue:NULL];
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"applied magnetic field 2" numberOfNodes:n indexes:element->NodeIndexes buffer:&extmz minValue:NULL maxValue:NULL];
        
        FEMVariable *mfx = [utilities getVariableFrom:model.variables model:model name:@"magnetic flux density 1" onlySearch:NULL maskName:nil info:&found];
        FEMVariable *mfy = [utilities getVariableFrom:model.variables model:model name:@"magnetic flux density 2" onlySearch:NULL maskName:nil info:&found];
        FEMVariable *mfz = [utilities getVariableFrom:model.variables model:model name:@"magnetic flux density 3" onlySearch:NULL maskName:nil info:&found];
        variableArraysContainer *mfxContainers = mfx.getContainers;
        variableArraysContainer *mfyContainers = mfy.getContainers;
        variableArraysContainer *mfzContainers = mfz.getContainers;
        if (mfx != nil) {
            for (i=0; i<n; i++) {
                extmx.vector[i] = extmx.vector[i] + mfxContainers->Values[mfxContainers->Perm[element->NodeIndexes[i]]];
                extmy.vector[i] = extmy.vector[i] + mfyContainers->Values[mfyContainers->Perm[element->NodeIndexes[i]]];
                extmz.vector[i] = extmz.vector[i] + mfzContainers->Values[mfzContainers->Perm[element->NodeIndexes[i]]];
            }
        }
        
        b[0] = 0.0;
        b[1] = 0.0;
        b[2] = 0.0;
        for (i=0; i<n; i++) {
            b[0] = b[0] + integration.basis[i]*mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
            b[1] = b[1] + integration.basis[i]*myContainers->Values[myContainers->Perm[element->NodeIndexes[i]]];
            b[2] = b[2] + integration.basis[i]*mzContainers->Values[mzContainers->Perm[element->NodeIndexes[i]]];
        }
        
        b[0] = b[0] + cblas_ddot(n, extmx.vector, 1, integration.basis, 1);
        b[1] = b[1] + cblas_ddot(n, extmy.vector, 1, integration.basis, 1);
        b[2] = b[2] + cblas_ddot(n, extmz.vector, 1, integration.basis, 1);
        
        double mu = cblas_ddot(n, permeability.vector, 1, integration.basis, 1);
        memset( *dhdx, 0.0, (3*3)*sizeof(double) );
        for (i=0; i<3; i++) {
            for (j=0; j<n; j++) {
                dhdx[0][i] = dhdx[0][i] + (integration.basisFirstDerivative[j][i] * mxContainers->Values[mxContainers->Perm[element->NodeIndexes[j]]]/permeability.vector[j]);
            }
            for (j=0; j<n; j++) {
                dhdx[1][i] = dhdx[1][i] + (integration.basisFirstDerivative[j][i] * myContainers->Values[myContainers->Perm[element->NodeIndexes[j]]]/permeability.vector[j]);
            }
            for (j=0; j<n; j++) {
                dhdx[2][i] = dhdx[2][i] + (integration.basisFirstDerivative[j][i] * mzContainers->Values[mzContainers->Perm[element->NodeIndexes[j]]]/permeability.vector[j]);
            }
        }
        
        jouleheat = [self FEMDifferentials_computeMagneticHeat:b dhdx:dhdx permeability:mu sqrtMetric:sqrtMetric metric:metric symbols:symb model:model] / elcond;
        
        if (permeability.vector != NULL) {
            free_dvector(permeability.vector, 0, permeability.m-1);
            permeability.vector = NULL;
        }
        if (extmx.vector != NULL) {
            free_dvector(extmx.vector, 0, extmx.m-1);
            extmx.vector = NULL;
        }
        if (extmy.vector != NULL) {
            free_dvector(extmy.vector, 0, extmy.m-1);
            extmy.vector = NULL;
        }
        if (extmz.vector != NULL) {
            free_dvector(extmz.vector, 0, extmz.m-1);
            extmz.vector = NULL;
        }
    }
    
    // Axisymmetric vector potential
    if (jouleNode == 2) {
        sum = 0.0;
        for (i=0; i<n; i++) {
            sum = sum + integration.basis[i] * mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
        }
        jouleheat = elcond * sum;
    }
    
    // Static current conduction
    if (jouleNode == 3) {
        // The electric field at integration point
        b[0] = 0.0;
        b[1] = 0.0;
        b[2] = 0.0;
        for (i=0; i<n; i++) {
            b[0] = b[0] + integration.basisFirstDerivative[i][0] * mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
            b[1] = b[1] + integration.basisFirstDerivative[i][1] * mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
            b[2] = b[2] + integration.basisFirstDerivative[i][2] * mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
        }
        jouleheat = elcond * cblas_ddot(3, b, 1, b, 1);
    }
    
    if (buffer.vector != NULL) {
        free_dvector(buffer.vector, 0, buffer.m-1);
        buffer.vector = NULL;
    }
    
    return jouleheat;
}


@end
