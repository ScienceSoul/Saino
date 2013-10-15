//
//  FEMDifferentials.m
//  Saino
//
//  Created by Seddik hakime on 18/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMDifferentials.h"
#import "FEMBodyForce.h"
#import "FEMListUtilities.h"
#import "FEMUtilities.h"
#import "FEMNumericIntegration.h"
#import "FEMCoordinateSystems.h"
#import "FEMLinearAlgebra.h"
#import "FEMMaterial.h"

@interface FEMDifferentials ()
-(double)FEMDifferentials_computeMagneticHeat:(double[3])b dhdx:(double[][3])dhdx permeability:(double)mu sqrtMetric:(double)sqrtMetric metric:(double[][3])metric symbols:(double[][3][3])symb model:(FEMModel *)model;
@end

@implementation FEMDifferentials

#pragma mark Private methods

-(double)FEMDifferentials_computeMagneticHeat:(double[3])b dhdx:(double[][3])dhdx permeability:(double)mu sqrtMetric:(double)sqrtMetric metric:(double[][3])metric symbols:(double[][3][3])symb model:(FEMModel *)model {
    
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

/*************************************************************************************************
    Compute the Joule heating at integration point (u, v, w) given the appropriate electrostatic
    or magnetic field that indicates the current through a conductor
*************************************************************************************************/
-(double)jouleHeatElement:(Element_t *)element nodes:(Nodes_t *)nodes numberOfNodes:(int)n integrationU:(double)u integrationV:(double)v integrationW:(double)w mesh:(FEMMesh *)mesh model:(FEMModel *)model {
 
    int i, j, k, bf_id, jouleNode;
    double b[3], dhdx[3][3], dSymb[3][3][3][3], elcond, jouleheat, metric[3][3], sqrtElementMetric, sqrtMetric, sum, symb[3][3][3], x, y, z;
    BOOL all, found, stat;
    FEMBodyForce *bodyForceAtID;
    FEMMaterial *materialAtID;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    variableArraysContainer *mxContainers = NULL;
    
    jouleheat = 0.0;
    
    if ((model.bodies)[element->BodyID-1][@"body force"] != nil) {
        bf_id = [(model.bodies)[element->BodyID-1][@"body force"] intValue];
    } else {
        return jouleheat;
    }

    bodyForceAtID = (model.bodyForces)[bf_id-1];
    FEMListUtilities *listUtilities = [[FEMListUtilities alloc] init];
    if ([listUtilities listGetLogical:model inArray:bodyForceAtID.valuesList forVariable:@"joule heat" info:&found] == NO) return jouleheat;
    
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
    FEMNumericIntegration *integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) errorfunct("FEMDifferentials:jouleHeatElement", "Allocation error in FEMNumericIntegration!");
    stat = [integration setBasisForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
    stat = [integration setBasisFirstDerivativeForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
    stat = [integration setMetricDeterminantForElement:element elementNodes:nodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
    sqrtElementMetric = integration.metricDeterminant;
    
    if (model.coordinates != cartesian) {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        for (i=0; i<n; i++) {
            x = x + nodes->x[i]*integration.basis[i];
            y = y + nodes->y[i]*integration.basis[i];
            z = z + nodes->z[i]*integration.basis[i];
        }
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
        NSLog(@"FEMDifferentials:jouleHeatElement: use electrical conductivity instead of electric\n");
    } else {
        found = [listUtilities listGetReal:model inArray:materialAtID.valuesList forVariable:@"electric conductivity" numberOfNodes:n indexes:element->NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
    }
    elcond = 0.0;
    if (found == YES) {
        for (i=0; i<n; i++) {
            elcond = elcond + buffer.vector[i]*integration.basis[i];
        }
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
        
        sum = 0.0;
        for (i=0; i<n; i++) {
            sum = sum + integration.basis[i]*extmx.vector[i];
        }
        b[0] = b[0] + sum;
        sum = 0.0;
        for (i=0; i<n; i++) {
            sum = sum + integration.basis[i]*extmy.vector[i];
        }
        b[1] = b[1] + sum;
        sum = 0.0;
        for (i=0; i<n; i++) {
            sum = sum + integration.basis[i]*extmz.vector[i];
        }
        b[2] = b[2] + sum;
        
        double mu = 0.0;
        for (i=0; i<n; i++) {
            mu = mu + integration.basis[i]*permeability.vector[i];
        }
        for (i=0; i<3; i++) {
            sum = 0.0;
            for (j=0; j<n; j++) {
                sum = sum + (integration.basisFirstDerivative[j][i] * mxContainers->Values[mxContainers->Perm[element->NodeIndexes[j]]]/permeability.vector[j]);
            }
            dhdx[0][i] = sum;
            sum = 0.0;
            for (j=0; j<n; j++) {
                sum = sum + (integration.basisFirstDerivative[j][i] * myContainers->Values[myContainers->Perm[element->NodeIndexes[j]]]/permeability.vector[j]);
            }
            dhdx[1][i] = sum;
            sum = 0.0;
            for (j=0; j<n; j++) {
                sum = sum + (integration.basisFirstDerivative[j][i] * mzContainers->Values[mzContainers->Perm[element->NodeIndexes[j]]]/permeability.vector[j]);
            }
            dhdx[2][i] = sum;
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
            sum = sum + integration.basis[i]*mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
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
            b[0] = b[0] + integration.basisFirstDerivative[i][0]*mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
            b[1] = b[1] + integration.basisFirstDerivative[i][1]*mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
            b[2] = b[2] + integration.basisFirstDerivative[i][2]*mxContainers->Values[mxContainers->Perm[element->NodeIndexes[i]]];
        }
        sum = 0.0;
        for (i=0; i<3; i++) {
            sum = sum + b[i]*b[i];
        }
        jouleheat = elcond * sum;
    }
    
    if (buffer.vector != NULL) {
        free_dvector(buffer.vector, 0, buffer.m-1);
        buffer.vector = NULL;
    }
    [integration deallocation:mesh];
    
    return jouleheat;
}


@end
