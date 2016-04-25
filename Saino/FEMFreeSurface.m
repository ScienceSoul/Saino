//
//  FEMFreeSurface.m
//  Saino
//
//  Created by Seddik hakime on 29/05/2014.
//  Copyright (c) 2014 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMFreeSurface.h"

#import "FEMCore.h"
#import "FEMSolution.h"
#import "FEMListUtilities.h"
#import "FEMUtilities.h"
#import "FEMElementDescription.h"
#import "FEMElementUtils.h"
#import "GaussIntegration.h"
#import "Utils.h"

@interface FEMFreeSurface ()
-(void)FEMFreeSurface_poissonSolution:(double * __nonnull)poissonSolution xCoord:(double * __nonnull)xCoord yCoord:(double * __nonnull)yCoord zCoord:(double * __nonnull)zCoord moved:(int)moved model:(FEMModel * __nonnull)model core:(FEMCore * __nonnull)core integration:(FEMNumericIntegration * __nonnull)integration listUtilities:(FEMListUtilities * __nonnull)listUtilities;
@end

@implementation FEMFreeSurface

#pragma mark Private methods...

-(void)FEMFreeSurface_poissonSolution:(double * __nonnull)poissonSolution xCoord:(double * __nonnull)xCoord yCoord:(double * __nonnull)yCoord zCoord:(double * __nonnull)zCoord moved:(int)moved model:(FEMModel * __nonnull)model core:(FEMCore * __nonnull)core integration:(FEMNumericIntegration * __nonnull)integration listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    int i, j, k, n, p, q, t, dim;
    int *cPerm;
    double a, detJ, s, u, v, w;
    double *forceVector, **localMatrix;
    double *basis = NULL, **basisFirstDerivative = NULL;
    BOOL found, stat;
    FEMMatrix *cMatrix = nil;
    Element_t *elements = NULL;
    Nodes_t *nodes = NULL;
    GaussIntegrationPoints *IP = NULL;
    
    basis = integration.basis;
    basisFirstDerivative = integration.basisFirstDerivative;
    
    FEMSolution *solution = [[FEMSolution alloc] init];
    solution.mesh = (FEMMesh *)model.mesh;
    [solution.solutionInfo setObject:@"cgs" forKey:@"linear system iterative method"];
    [solution.solutionInfo setObject:@500 forKey:@"linear system maximum iterations"];
    [solution.solutionInfo setObject:@1.0e-09 forKey:@"linear system convergence tolerance"];
    [solution.solutionInfo setObject:@"ilu0" forKey:@"linear system preconditioning"];
    [solution.solutionInfo setObject:@1 forKey:@"linear system residual output"];
    
    forceVector = doublevec(0, model.numberOfNodes-1);
    cPerm = intvec(0, model.numberOfNodes-1);
    
    FEMElementUtils *elementUtils = [[FEMElementUtils alloc] init];
    cMatrix = [elementUtils createMatrixInModel:model forSolution:solution mesh:solution.mesh dofs:1 permutation:cPerm sizeOfPermutation:model.numberOfNodes matrixFormat:MATRIX_CRS optimizeBandwidth:NO equationName:nil discontinuousGalerkinSolution:NULL globalBubbles:NULL nodalDofsOnly:NULL projectorDofs:NULL];
    if (cMatrix == nil) fatal("FEMFreeSurface_poissonSolveModel", "Can't create matrix.");
    
    localMatrix = doublematrix(0, model.maxElementNodes-1, 0, model.maxElementNodes-1);
    nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    nodes->x = doublevec(0, model.maxElementNodes-1);
    nodes->y = doublevec(0, model.maxElementNodes-1);
    nodes->z = doublevec(0, model.maxElementNodes-1);
    
    solution.timeOrder = 0;
    
    FEMMatrixCRS *crsMatrix = [[FEMMatrixCRS alloc] init];
    [crsMatrix zeroMatrix:cMatrix];
    
    memset(forceVector, 0.0, model.numberOfNodes*sizeof(double) );
    elements = model.getElements;
    for (i=0; i<model.numberOfBulkElements; i++) {
        
        dim = elements[i].Type.dimension;
        n = elements[i].Type.NumberOfNodes;
        for (j=0; j<n; j++) {
            nodes->x[j] = xCoord[elements[i].NodeIndexes[j]];
            nodes->y[j] = yCoord[elements[i].NodeIndexes[j]];
            nodes->z[j] = zCoord[elements[i].NodeIndexes[j]];
        }
        
        IP = GaussQuadrature(&elements[i], NULL, NULL);
        memset(*localMatrix, 0.0, (model.maxElementNodes*model.maxElementNodes)*sizeof(double) );
        for (t=0; t<IP->n; t++) {
            
            u = IP->u[t];
            v = IP->v[t];
            w = IP->w[t];
            
            // Basis funciton values and derivatives at the integration point
            stat = [integration setBasisForElement:&elements[i] elementNodes:nodes inMesh:solution.mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
            stat = [integration setBasisFirstDerivativeForElement:&elements[i] elementNodes:nodes inMesh:solution.mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
            stat = [integration setMetricDeterminantForElement:&elements[i] elementNodes:nodes inMesh:solution.mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
            detJ = integration.metricDeterminant;
            s = detJ * IP->s[t];
            
            for (p=0; p<n; p++) {
                for (q=0; q<n; q++) {
                    a = basisFirstDerivative[p][moved] * basisFirstDerivative[q][moved];
                    localMatrix[p][q] = localMatrix[p][q] + s * a;
                }
            }
        }
        [crsMatrix glueLocalMatrix:localMatrix inMatrix:cMatrix numberOfNodes:n dofs:1 indexes:elements[i].NodeIndexes];
    }
    
    for (t=model.numberOfBulkElements; t<model.numberOfBulkElements+model.numberOfBoundaryElements; t++) {
        
        n = elements[t].Type.NumberOfNodes;
        for (FEMBoundaryCondition *boundary in model.boundaryConditions) {
            if (elements[t].BoundaryInfo->Constraint == boundary.tag) {
                if ([listUtilities listGetLogical:model inArray:boundary.valuesList forVariable:@"free moving" info:&found] == YES) {
                    for (j=0; j<n; j++) {
                        k = elements[t].NodeIndexes[j];
                        forceVector[k] = poissonSolution[k];
                        [crsMatrix zeroRowInMatrix:cMatrix numberOfRows:k];
                        [crsMatrix setElementInMatrix:cMatrix row:k col:k value:1.0];
                    }
                }
            }
        }
    }
    
    [core iterativeSolveMatrix:cMatrix result:poissonSolution rhs:forceVector dimensions:NULL solution:solution];
    
    free_dvector(forceVector, 0, model.numberOfNodes-1);
    free_ivector(cPerm, 0, model.numberOfNodes-1);
    free_dmatrix(localMatrix, 0, model.maxElementNodes-1, 0, model.maxElementNodes-1);
    
    free_dvector(nodes->x, 0, model.maxElementNodes-1);
    free_dvector(nodes->y, 0, model.maxElementNodes-1);
    free_dvector(nodes->z, 0, model.maxElementNodes-1);
    free(nodes);
    
    [cMatrix deallocation];
    [solution deallocation];
    cMatrix = nil;
}

#pragma mark Public methods...

-(void)moveBoundaryModel:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration relax:(double)relax {
    
    int i, ii, j, k, m, n, t, fiter, which;
    double dxdu=0.0, dydu=0.0, dzdu=0.0, dxdv=0.0, dydv=0.0, dzdv=0.0, feps, r, s, sum, u, v, ux, uy, uz, x, x1, y1, y, z;
    double **averagedNormals, dLBasisdx[16][2], *nodalBasis, nrm[3], *xCoord, *yCoord, *zCoord;
    BOOL found, *visited, *turned, xMoved, yMoved, zMoved;
    NSArray *bc = nil;
    FEMSolution *solution;
    FEMMesh *mesh = nil;
    FEMVariable *velocity1 = nil, *velocity2 = nil, *velocity3 = nil;
    Element_t *boundary = NULL, *elements = NULL;
    Nodes_t *nodes = NULL, *boundaryNodes = NULL, *elementNodes;
    variableArraysContainer *velocity1Containers = NULL, *velocity2Containers = NULL, *velocity3Containers = NULL;
    
    FEMCore *core = [FEMCore sharedCore];
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    FEMElementDescription *elementDescription = [[FEMElementDescription alloc] init];
    
    solution = (FEMSolution *)model.solution;
    elements = model.getElements;
    nodes = model.getNodes;
    
    if ((solution.solutionInfo)[@"free surface convergence tolerance"] != nil) {
        feps = [(solution.solutionInfo)[@"free surface convergence tolerance"] doubleValue];
    } else feps = 1.0e-12;

    if ((solution.solutionInfo)[@"free surface maximum iterations"] != nil) {
        fiter = [(solution.solutionInfo)[@"free surface maximum iterations"] intValue];
    } else fiter = 100;
    
    velocity1 = [utilities getVariableFrom:model.variables model:model name:@"velocity 1" onlySearch:NULL maskName:nil info:&found];
    velocity2 = [utilities getVariableFrom:model.variables model:model name:@"velocity 2" onlySearch:NULL maskName:nil info:&found];
    velocity3 = [utilities getVariableFrom:model.variables model:model name:@"velocity 3" onlySearch:NULL maskName:nil info:&found];
    if (velocity1 != nil) velocity1Containers = velocity1.getContainers;
    if (velocity2 != nil) velocity2Containers = velocity2.getContainers;
    if (velocity3 != nil) velocity3Containers = velocity3.getContainers;
    
    visited = (BOOL*)malloc(sizeof(BOOL) * model.numberOfNodes );
    turned = (BOOL*)malloc(sizeof(BOOL) * model.numberOfNodes );
    memset(visited, 0.0, model.numberOfNodes*sizeof(BOOL) );
    memset(turned, 0.0, model.numberOfNodes*sizeof(BOOL) );
    
    xCoord = doublevec(0, model.numberOfNodes-1);
    xMoved = NO;
    memcpy(xCoord, nodes->x, model.numberOfNodes*sizeof(double));
    
    yCoord = doublevec(0, model.numberOfNodes-1);
    yMoved = NO;
    memcpy(yCoord, nodes->y, model.numberOfNodes*sizeof(double));
    
    zCoord = doublevec(0, model.numberOfNodes-1);
    zMoved = NO;
    memcpy(zCoord, nodes->z, model.numberOfNodes*sizeof(double));
    
    averagedNormals = doublematrix(0, model.numberOfBoundaryConditions-1, 0, 2);
    memset(*averagedNormals, 0.0, (model.numberOfBoundaryConditions*3)*sizeof(double) );
    
    nodalBasis = doublevec(0, 15);
    memset(nodalBasis, 0.0, 16*sizeof(double) );
    
    ux = 0.0; uy = 0.0; uz = 0.0;
    // Check normal direction first in a separate loop, becuase within the node moving loop,
    // the parent elements might be a little funny....!
    mesh = solution.mesh;
    for (t=0; t<mesh.numberOfBoundaryElements; t++) {
        boundary = [core getBoundaryElement:solution atIndex:t];
        if ([core isActiveBoundaryElement:boundary inSolution:solution model:model] == NO) continue;
        
        bc = [core getBoundaryCondition:model forElement:boundary];
        if (bc == nil) continue;
        if ([listUtilities listGetLogical:model inArray:bc forVariable:@"free surface" info:&found] == NO) continue;
        
        n = boundary->Type.NumberOfNodes;
        if (boundaryNodes == NULL) {
            boundaryNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
            initNodes(boundaryNodes);
        }
        [core getNodes:solution model:model inElement:boundary resultNodes:boundaryNodes numberOfNodes:NULL mesh:nil];
        
        // Go through boundary element nodes
        for (i=0; i<n; i++) {
            k = boundary->NodeIndexes[i];
            if (boundary->Type.dimension == 1) {
                if (i == 0) {
                    j = 1;
                } else {
                    j = 0;
                }
                j = boundary->NodeIndexes[j];
                x1 = xCoord[j] - xCoord[k];
                y1 = yCoord[j] - yCoord[k];
                
                ux = velocity1Containers->Values[velocity1Containers->Perm[k]];
                uy = velocity2Containers->Values[velocity2Containers->Perm[k]];
                
                if (ux * x1 + uy * y1 > 0) continue;
            }
            
            // Should not move the same node twice, so check if already done
            if (visited[k] == NO) {
                visited[k] = YES;
                
                // 2D case, compute normal
                if (boundary->Type.dimension == 1) {
                    u = boundary->Type.NodeU[i];
                    
                    // Basis function derivatives with respect to local coordinates
                    memset(nodalBasis, 0.0, n*sizeof(double) );
                    for (j=0; j<n; j++) {
                        nodalBasis[j] = 1.0;
                        dLBasisdx[j][0] = [elementDescription firstDerivative1DInElement:boundary nodalValues:nodalBasis evalutationPoint:u];
                        nodalBasis[j] = 0.0;
                    }
                    memset(nrm, 0.0, sizeof(nrm) );
                    for (j=0; j<n; j++) {
                         nrm[0] =  nrm[0] + boundaryNodes->y[j] * dLBasisdx[j][0];
                    }
                    nrm[0] = -nrm[0];
                    for (j=0; j<n; j++) {
                        nrm[1] = nrm[1] + boundaryNodes->x[j] * dLBasisdx[j][0];
                    }
                    nrm[2] = 0.0;
                } else {
                    // 3D case, compute normal
                    u = boundary->Type.NodeU[i];
                    v = boundary->Type.NodeV[i];
                    
                    memset(nodalBasis, 0.0, n*sizeof(double) );
                    for (j=0; j<n; j++) {
                        nodalBasis[j] = 1.0;
                        dLBasisdx[j][0] = [elementDescription firstDerivativeU2DInElement:boundary nodalValues:nodalBasis evaluatedAt:u andAt:v];
                        dLBasisdx[j][1] = [elementDescription firstDerivativeV2DInElement:boundary nodalValues:nodalBasis evaluatedAt:u andAt:v];
                        nodalBasis[j] = 0.0;
                    }
                    
                    dxdu = 0.0; dydu = 0.0; dzdu = 0.0;
                    dxdv = 0.0; dydv = 0.0; dzdv = 0.0;
                    for (j=0; j<n; j++) {
                        dxdu = dxdu + boundaryNodes->x[j] * dLBasisdx[j][0];
                        dydu = dydu + boundaryNodes->y[j] * dLBasisdx[j][0];
                        dzdu = dzdu + boundaryNodes->z[j] * dLBasisdx[j][0];
                        
                        dxdv = dxdv + boundaryNodes->x[j] * dLBasisdx[j][1];
                        dydv = dydv + boundaryNodes->y[j] * dLBasisdx[j][1];
                        dzdv = dzdv + boundaryNodes->z[j] * dLBasisdx[j][1];
                    }
                    nrm[0] = dydu * dzdv - dydv * dzdu;
                    nrm[1] = dxdv * dzdu - dxdu * dzdv;
                    nrm[2] = dxdu * dydv - dxdv * dydu;
                }
                
                // Turn the normal to point outwards, or torwards less dense material
                x = nodes->x[k];
                y = nodes->y[k];
                z = nodes->z[k];
                [elementDescription checkNormalDirectionInBDElement:boundary forNormals:nrm mesh:mesh x:x y:y z:z turn:&turned[k]];
                
                m = boundary->BoundaryInfo->Constraint;
                vDSP_svesqD(nrm, 1, &sum, 3);
                r = sqrt(sum);
                for (j=0; j<3; j++) {
                    averagedNormals[m-1][j] = averagedNormals[m-1][j] + fabs(nrm[j]) / r;
                }
            }
        }
    }
    free_dvector(boundaryNodes->x, 0, boundaryNodes->numberOfNodes-1);
    free_dvector(boundaryNodes->y, 0, boundaryNodes->numberOfNodes-1);
    free_dvector(boundaryNodes->z, 0, boundaryNodes->numberOfNodes-1);
    free(boundaryNodes);
    
    // Iterate until convergence
    for (ii=1; ii<=fiter; ii++) {
        s = 0.0;
        memset(visited, 0.0, model.numberOfNodes*sizeof(BOOL) );
        
        for (t=0; t<mesh.numberOfBoundaryElements; t++) {
            boundary = [core getBoundaryElement:solution atIndex:t];
            if ([core isActiveBoundaryElement:boundary inSolution:solution model:model] == NO) continue;
            
            bc = [core getBoundaryCondition:model forElement:boundary];
            if (bc == nil) continue;
            if ([listUtilities listGetLogical:model inArray:bc forVariable:@"free surface" info:&found] == NO) continue;
            
            i = model.dimension;
            int minv = 1;
            which = [listUtilities listGetInteger:model inArray:bc forVariable:@"free coordinate" info:&found minValue:&minv maxValue:&i];
            if (found == NO) {
                i = boundary->BoundaryInfo->Constraint;
                for (j=0; j<3; j++) {
                    nrm[j] = averagedNormals[i-1][j];
                }
                which = 1;
                for (i=1; i<3; i++) {
                    if (nrm[i] > nrm[which-1]) which = i+1;
                }
            }
            
            n = boundary->Type.NumberOfNodes;
            if (boundaryNodes == NULL) {
                boundaryNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
                initNodes(boundaryNodes);
            }
            [core getNodes:solution model:model inElement:boundary resultNodes:boundaryNodes numberOfNodes:NULL mesh:nil];
            
            for (i=0; i<n; i++) {
                k = boundary->NodeIndexes[i];
                if (boundary->Type.dimension == 1) {
                    if (i == 0) {
                        j = 1;
                    } else {
                        j = 0;
                    }
                    j = boundary->NodeIndexes[j];
                    x1 = xCoord[j] - xCoord[k];
                    y1 = yCoord[j] - yCoord[k];
                    
                    ux = velocity1Containers->Values[velocity1Containers->Perm[k]];
                    uy = velocity2Containers->Values[velocity2Containers->Perm[k]];
                    
                    if (ux * x1 + uy * y1 > 0) continue;
                }
                
                // Should not move the same node twice, so check if already done
                if (visited[k] == NO) {
                    visited[k] = YES;
                    
                    // 2D case, compute normal
                    if (boundary->Type.dimension == 1) {
                        u = boundary->Type.NodeU[i];
                        
                        // Basis function derivatives with respect to local coordinates
                        memset(nodalBasis, 0.0, n*sizeof(double) );
                        for (j=0; j<n; j++) {
                            nodalBasis[j] = 1.0;
                            dLBasisdx[j][0] = [elementDescription firstDerivative1DInElement:boundary nodalValues:nodalBasis evalutationPoint:u];
                            nodalBasis[j] = 0.0;
                        }
                        memset(nrm, 0.0, sizeof(nrm) );
                        for (j=0; j<n; j++) {
                            nrm[0] =  nrm[0] + boundaryNodes->y[j] * dLBasisdx[j][0];
                        }
                        nrm[0] = -nrm[0];
                        for (j=0; j<n; j++) {
                            nrm[1] = nrm[1] + boundaryNodes->x[j] * dLBasisdx[j][0];
                        }
                        nrm[2] = 0.0;
                    } else {
                        // 3D case, compute normal
                        u = boundary->Type.NodeU[i];
                        v = boundary->Type.NodeV[i];
                        
                        memset(nodalBasis, 0.0, n*sizeof(double) );
                        for (j=0; j<n; j++) {
                            nodalBasis[j] = 1.0;
                            dLBasisdx[j][0] = [elementDescription firstDerivativeU2DInElement:boundary nodalValues:nodalBasis evaluatedAt:u andAt:v];
                            dLBasisdx[j][1] = [elementDescription firstDerivativeV2DInElement:boundary nodalValues:nodalBasis evaluatedAt:u andAt:v];
                            nodalBasis[j] = 0.0;
                        }
                        
                        dxdu = 0.0; dydu = 0.0; dzdu = 0.0;
                        dxdv = 0.0; dydv = 0.0; dzdv = 0.0;
                        for (j=0; j<n; j++) {
                            dxdu = dxdu + boundaryNodes->x[j] * dLBasisdx[j][0];
                            dydu = dydu + boundaryNodes->y[j] * dLBasisdx[j][0];
                            dzdu = dzdu + boundaryNodes->z[j] * dLBasisdx[j][0];
                            
                            dxdv = dxdv + boundaryNodes->x[j] * dLBasisdx[j][1];
                            dydv = dydv + boundaryNodes->y[j] * dLBasisdx[j][1];
                            dzdv = dzdv + boundaryNodes->z[j] * dLBasisdx[j][1];
                        }
                        nrm[0] = dydu * dzdv - dydv * dzdu;
                        nrm[1] = dxdv * dzdu - dxdu * dzdv;
                        nrm[2] = dxdu * dydv - dxdv * dydu;
                    }
                    
                    // Turn the normal to point outwards or towards less dense material
                    if (turned[k] == YES) {
                        for (j=0; j<3; j++) {
                            nrm[j] = -nrm[j];
                        }
                    }
                    
                    // Now then, lets move the node so that u.n will be reduced
                    if (boundary->Type.dimension == 1) {
                        // TODO: this want handle the three node line
                        // 2D case, move the nodes...
                        r = ux * nrm[0] + uy * nrm[1];
                        if (which == 2) {
                            if (fabs(ux) > AEPS) {
                                if (turned[k] == YES) {
                                    nodes->y[k] = nodes->y[k] - r / (ux * dLBasisdx[i][0]);
                                } else {
                                    nodes->y[k] = nodes->y[k] + r / (ux * dLBasisdx[i][0]);
                                }
                                yMoved = YES;
                            }
                        } else {
                            if (fabs(uy) > AEPS) {
                                if (turned[k] == YES) {
                                    nodes->x[k] = nodes->x[k] + r /(uy * dLBasisdx[i][0]);
                                } else {
                                    nodes->x[k] = nodes->x[k] - r / (uy * dLBasisdx[i][0]);
                                }
                                xMoved = YES;
                            }
                        }
                    } else {
                        // 3D case, move the nodes...
                        // TODO: this is just guesswork, no testing done...
                        ux = velocity1Containers->Values[velocity1Containers->Perm[k]];
                        uy = velocity2Containers->Values[velocity2Containers->Perm[k]];
                        uz = velocity3Containers->Values[velocity3Containers->Perm[k]];
                        
                        r = ux * nrm[0] + uy * nrm[1] + uz * nrm[2];
                        if (which == 1) {
                            if (fabs(uy) > AEPS || fabs(uz) > AEPS) {
                                nodes->x[k] = nodes->x[k] + r / ( (dzdu * dLBasisdx[i][1] - dzdv * dLBasisdx[i][0]) * uy +
                                              (dydv * dLBasisdx[i][0] - dydu * dLBasisdx[i][1]) * uz );
                                xMoved = YES;
                            }
                        } else if (which == 2) {
                            if (fabs(ux) > AEPS || fabs(uz) > AEPS) {
                                nodes->y[k] = nodes->y[k] + r / ( (dzdv * dLBasisdx[i][0] - dzdu * dLBasisdx[i][1]) * ux +
                                              (dxdu * dLBasisdx[i][1] - dxdv * dLBasisdx[i][0]) * uz );
                                yMoved = YES;
                            }
                        } else {
                            if (fabs(ux) > AEPS || fabs(uy) > AEPS) {
                                nodes->z[k] = nodes->z[k] + r / ( (dydu * dLBasisdx[i][1] - dydv * dLBasisdx[i][0]) * ux +
                                              (dxdv * dLBasisdx[i][0] - dxdu * dLBasisdx[i][1]) * uy );
                                zMoved = YES;
                            }
                        }
                    }
                    vDSP_svesqD(nrm, 1, &sum, 3);
                    s = s + pow(( (ux * nrm[0] + uy * nrm[1] + uz * nrm[2]) / sqrt(sum) ), 2.0);
                }
            }
        }
        
        s = sqrt(s);
        fprintf(stdout, "FEMFreeSurface:moveBoundaryModel: iter: %d, free surface residual: %f.\n", ii, s);
        if (s < feps) {
            free_dvector(boundaryNodes->x, 0, boundaryNodes->numberOfNodes-1);
            free_dvector(boundaryNodes->y, 0, boundaryNodes->numberOfNodes-1);
            free_dvector(boundaryNodes->z, 0, boundaryNodes->numberOfNodes-1);
            free(boundaryNodes);
            break;
        }
    }
    
    if (xMoved == YES) {
        [self FEMFreeSurface_poissonSolution:nodes->x xCoord:xCoord yCoord:yCoord zCoord:zCoord moved:0 model:model core:core integration:integration listUtilities:listUtilities];
        if (relax != 1.0) {
            for (i=0; i<model.numberOfNodes; i++) {
                nodes->x[i] = relax * nodes->x[i] + (1.0 - relax) * xCoord[i];
            }
        }
    }
    
    if (yMoved == YES) {
        [self FEMFreeSurface_poissonSolution:nodes->y xCoord:xCoord yCoord:yCoord zCoord:zCoord moved:1 model:model core:core integration:integration listUtilities:listUtilities];
        if (relax != 1.0) {
            for (i=0; i<model.numberOfNodes; i++) {
                nodes->y[i] = relax * nodes->y[i] + (1.0 - relax) * yCoord[i];
            }
        }
    }
    
    if (zMoved == YES) {
        [self FEMFreeSurface_poissonSolution:nodes->z xCoord:xCoord yCoord:yCoord zCoord:zCoord moved:2 model:model core:core integration:integration listUtilities:listUtilities];
        if (relax != 1.0) {
            for (i=0; i<model.numberOfNodes; i++) {
                nodes->z[i] = relax * nodes->z[i] + (1.0 - relax) * zCoord[i];
            }
        }
    }
    
    elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    elementNodes->x = doublevec(0, model.numberOfNodes-1);
    elementNodes->y = doublevec(0, model.numberOfNodes-1);
    elementNodes->z = doublevec(0, model.numberOfNodes-1);
    for (i=0; i<model.numberOfBulkElements; i++) {
        n = elements[i].Type.NumberOfNodes;
        for (j=0; j<n; j++) {
            elementNodes->x[j] = nodes->x[elements[i].NodeIndexes[j]];
            elementNodes->y[j] = nodes->y[elements[i].NodeIndexes[j]];
            elementNodes->z[j] = nodes->z[elements[i].NodeIndexes[j]];
        }
        [elementDescription computeStabilizationParameterInElement:&elements[i] nodes:elementNodes mesh:mesh numberOfNodes:n mk:&elements[i].StabilizationMK hk:&elements[i].hK];
    }
    free_dvector(elementNodes->x, 0, model.numberOfNodes-1);
    free_dvector(elementNodes->y, 0, model.numberOfNodes-1);
    free_dvector(elementNodes->z, 0, model.numberOfNodes-1);
    free(elementNodes);
    
    free(visited);
    free(turned);
    free_dvector(xCoord, 0, model.numberOfNodes-1);
    free_dvector(yCoord, 0, model.numberOfNodes-1);
    free_dvector(zCoord, 0, model.numberOfNodes-1);
    free_dmatrix(averagedNormals, 0, model.numberOfBoundaryConditions-1, 0, 2);
    free_dvector(nodalBasis, 0, 15);
}

@end
