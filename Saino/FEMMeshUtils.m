//
//  FEMMeshUtils.m
//  Saino
//
//  Created by Seddik hakime on 28/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMMeshUtils.h"
#import "FEMCore.h"
#import "FEMElementUtils.h"
#import "FEMListUtilities.h"
#import "FEMPElementMaps.h"
#import "FEMElementDescription.h"
#import "FEMUtilities.h"
#import "FEMNumericIntegration.h"
#import "FEMInterpolation.h"
#import "FEMListMatrix.h"
#import "FEMMatrixCRS.h"
#import "FEMElementUtils.h"
#import "FEMLinearAlgebra.h"
#import "FEMInterpolateMeshToMesh.h"
#import "FEMPost.h"
#import "Utils.h"
#import "TimeProfile.h"
#import "GaussIntegration.h"

@interface FEMMeshUtils ()
-(void)FEMMeshUtils_assignConstraints:(FEMMesh * __nonnull)mesh;
-(void)FEMMeshUtils_fixFaceEdges:(FEMMesh * __nonnull)mesh;
-(Element_t * __nullable)FEMMeshUtils_getEntityForElement:(Element_t * __nonnull)element edge:(int)number inMesh:(FEMMesh * __nonnull)mesh ;
-(void)FEMMeshUtils_unitSegmentDivisionTable:(double * __nonnull)w levels:(int)n model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(void)FEMMeshUtils_coordinateTransformationNodalType:(NSString * __nonnull)type vector:(double * __nonnull)vector solution:(FEMSolution * __nonnull)solution;
-(void)FEMMeshUtils_markHaloNodesMesh:(FEMMesh * __nonnull)mesh haloNodes:(BOOL * __nullable)haloNodes sizeHaloNodes:(int * __nullable)size;
-(Element_t * __nullable)FEMMeshUtils_findFaceParent:(Element_t * __nonnull)parent element:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model;
-(FEMMatrix * __nullable)FEMMeshUtils_weightedProjectorDiscontinousMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model boundary:(int)bc listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(FEMMatrix * __nonnull)FEMMeshUtils_nodalProjectorDiscontinuousMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model boundary:(int)bc listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(BOOL)FEMMeshUtils_createInterfaceMeshesModel:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh masterBoundary:(int)mbd targetBoundary:(int)trgt bMesh1:(FEMMesh * __nonnull)bMesh1 bMesh2:(FEMMesh * __nonnull)bMesh2;
-(void)FEMMeshUtils_mapInterfaceCoordinateMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 BCParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(void)FEMMeshUtils_cylinderFitMesh:(FEMMesh * __nonnull)pMesh pParams:(FEMBoundaryCondition * __nonnull)pParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(void)FEMMeshUtils_rotationInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams cylindrical:(BOOL)cylindrical radius:(double)radius fullCircle:(BOOL * __nonnull)fullCircle model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(void)FEMMeshUtils_radiusInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(void)FEMMeshUtils_flatInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(void)FEMMeshUtils_planeInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(void)FEMMeshUtils_checkInterfaceAngleMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 angles:(double * __nonnull)angles gotAngles:(BOOL * __nonnull)gotAngles;
-(void)FEMMeshUtils_overlayInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities;
-(FEMMatrix * __nonnull)FEMMeshUtils_nodalProjectorMesh2:(FEMMesh * __nonnull)bMesh2 mesh1:(FEMMesh * __nonnull)bMesh1 useQuadrantTree:(BOOL)useQuadrantTree repeating:(BOOL)repeating antiRepeating:(BOOL)antiRepeating model:(FEMModel * __nonnull)model;
-(void)FEMMeshUtils_setProjectorRowSum:(FEMMatrix * __nonnull)projector;
@end

@implementation FEMMeshUtils

#pragma mark Private methods

-(void)FEMMeshUtils_assignConstraints:(FEMMesh * __nonnull)mesh {
    
    int i, j, k, l, n, nfound, count=0;
    int *faceIndexes;
    Element_t *elements = NULL, *element = NULL, *boundary, *faces = NULL, *face = NULL;
    
    elements = mesh.getElements;
    
    for (i=0; i<mesh.numberOfBoundaryElements; i++) {
            
        boundary = &elements[mesh.numberOfBulkElements];
        
        element = boundary->BoundaryInfo->Left;
        if (element == NULL) element = boundary->BoundaryInfo->Right;
        if (element == NULL) continue;
        
        switch (boundary->Type.dimension) {
            case 1:
                faces = mesh.getEdges;
                faceIndexes = element->EdgeIndexes;
                count = element->sizeEdgeIndexes;
                break;
            case 2:
                faces = mesh.getFaces;
                faceIndexes = element->FaceIndexes;
                count = element->sizeFaceIndexes;
                break;
            default:
                faces = NULL;
                faceIndexes = NULL;
                break;
        }
        
        if (faces == NULL || faceIndexes == NULL) continue;

        for (j=0; j<count; j++) {
            face = &faces[faceIndexes[j]];
            
            n = boundary->Type.NumberOfNodes;
            if (n == 0) continue;
            nfound = 0;
            for (k=0; k<n; k++) {
                for (l=0; l<n; l++) {
                    if (boundary->NodeIndexes[k] == face->NodeIndexes[l]) nfound++;
                }
            }
            if (nfound == n) {
                face->BoundaryInfo = boundary->BoundaryInfo;
                break;
            }
        }
    }
}

-(void)FEMMeshUtils_fixFaceEdges:(FEMMesh * __nonnull)mesh {
    
    int i, j, k, l, n, swap;
    int edgeind[4], i1[2], i2[2];
    Element_t *edges = NULL, *faces = NULL;
    BOOL condition;
    
    faces = mesh.getFaces;
    edges = mesh.getEdges;
    
    for (i=0; i<mesh.numberOfFaces; i++) {
        n = faces[i].Type.NumberOfEdges;
        for (l=0; l<n; l++) {
            edgeind[l] = faces[i].EdgeIndexes[l];
        }
        for (j=0; j<n; j++) {
            for (l=0; l<2; l++) {
                i1[l] = edges[edgeind[j]].NodeIndexes[l];
            }
            if (i1[0] > i1[1]) {
                swap = i1[0];
                i1[0] = i1[1];
                i1[1] = swap;
            }
            for (k=0; k<n; k++) {
                i2[0] = k;
                i2[1] = k+1;
                if (i2[1] > n-1) i2[1] = 0;
                for (l=0; l<2; l++) {
                    i2[l] = faces[i].NodeIndexes[i2[l]];
                }
                if (i2[0] > i2[1]) {
                    swap = i2[0];
                    i2[0] = i2[1];
                    i2[1] = swap;
                }
                condition = YES;
                for (l=0; l<2; l++) {
                    if (i1[l] != i2[l]) {
                        condition = NO;
                        break;
                    }
                }
                if (condition == YES) {
                    faces[i].EdgeIndexes[k] = edgeind[k];
                    break;
                }
            }
        }
    }
}

-(Element_t * __nullable)FEMMeshUtils_getEntityForElement:(Element_t * __nonnull)element edge:(int)number inMesh:(FEMMesh * __nonnull)mesh {
    
    Element_t *edges = NULL, *faces = NULL, *entity = NULL;
    
    edges = mesh.getEdges;
    faces = mesh.getFaces;
    
    entity = NULL;
    // Switch by element dimension
    switch (element->Type.dimension) {
        case 2:
            entity = &edges[element->EdgeIndexes[number]];
            break;
        case 3:
            entity = &faces[element->FaceIndexes[number]];
            break;
        default:
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_getEntityForElement: unsupported dimension.");
            return entity;
            break;
    }
    return entity;
}

/*************************************************************************************
 
    Create node distribution for a unit segment x \in [0,1] with n elements
    i.e. n+1 nodes. There are different options for the type of distribution.
    1) Even distribution
    2) Geometric distribution
    3) Arbitrary distribution determined by a functional dependence
    Note that the 3rd algorithm involves iterative solution of the nodal
    positions and is therefore not bullet-proof.
 
*************************************************************************************/
-(void)FEMMeshUtils_unitSegmentDivisionTable:(double * __nonnull)w levels:(int)n model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    BOOL found, gotRatio;
    
    // Linear distribution and initial guess for the generic case
    
    // Geometric division
    double q = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"extruded mesh ratio" info:&gotRatio minValue:NULL maxValue:NULL];
    if (gotRatio) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_unitSegmentDivisionTable: creating geometric division.");
        
        double h1 = (1.0 - pow(q, (1.0/n))) / (1.0 - q);
        w[0] = 0.0;
        double hn = h1;
        for (int i=1; i<=n-1; i++) {
            w[i] = w[i-1] + hn;
            hn = hn * pow(q, (1.0/n));
        }
        w[n] = 1.0;
    }
    // Generic division given by a function
    else if ([listUtilities listCheckPresentVariable:@"extruded mesh density" inArray:model.simulation.valuesList] == YES) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_unitSegmentDivisionTable: creating functional division.");
        
        // Initial guess is an even distribution
        for (int i=0; i<=n; i++) {
            w[i] = i / (1.0 * n);
        }
        double *wold = doublevec(0, n);
        double *h = doublevec(0, n-1);
        memcpy(wold, w, (n+1)*sizeof(double));
        
        // Parameters that determine the accuracy of the iteration
        int maxiter = 10000;
        double err_eps = 1.0e-6;
        
        // Iterate to have a density distribution
        int iter;
        double err, minhn, xn;
        for (iter=1; iter<=maxiter; iter++) {
            
            minhn = HUGE_VAL;
            memcpy(wold, w, (n+1)*sizeof(double));
            
            // Compute the point in the local mesh xn \in [0,1]
            // and get the mesh parameter for that element from
            // external function
            for (int i=1; i<=n; i++) {
                xn = (w[i] + w[i-1]) / 2.0;
                minhn = min(minhn, w[i]-w[i-1]);
                h[i-1] = [listUtilities listGetValueParameter:model inArray:model.simulation.valuesList forVariable:@"extruded mesh density" value:xn info:&found minValue:NULL maxValue:NULL];
                if (h[i-1] < DBL_EPSILON) fatal("FEMMeshUtils:FEMMeshUtils_unitSegmentDivisionTable", "Given value for h[i] was negative.");
            }
            
            // Utilize symmetric Gaus-Seidel to compute the new positions w(i)
            // from a weighted mean of the desired elemental densities h(i).
            // Note that something more clever could be applied here.
            // This is was just a first implementation....
            for (int i=1; i<=n-1; i++) {
                w[i] = (w[i-1]*h[i]+w[i+1]*h[i-1]) / (h[i-1]+h[i]);
            }
            for (int i=n-1; i>=1; i--) {
                w[i] = (w[i-1]*h[i]+w[i+1]*h[i-1]) / (h[i-1]+h[i]);
            }
            
            // If the maximum error is small compared to the minimum element size then exit
            double max = -HUGE_VAL;
            for (int i=0; i<=n; i++) {
                if (fabs(w[i]-wold[i])>max) {
                    max = fabs(w[i]-wold[i]);
                }
            }
            err = max/minhn;
            if (err < err_eps) {
                fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_unitSegmentDivisionTable: convergence obtained in %d iterations.\n", iter);
                break;
            }
        }
        if (iter > maxiter) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_unitSegmentDivisionTable: no convergence obtained for the unit mesh division.\n");
        }
        free_dvector(wold, 0, n);
        free_dvector(h, 0, n-1);
    } else { // Uniform division
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_unitSegmentDivisionTable: creating linear division.\n");
        for (int i=0; i<=n; i++) {
            w[i] = i/(1.0 * n);
        }
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_unitSegmentDivisionTable: mesh division ready.\n");
        for (int i=0; i<=n; i++) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_unitSegmentDivisionTable: w(%d): %f\n", i, w[i]);
        }
    }
}

-(void)FEMMeshUtils_coordinateTransformationNodalType:(NSString * __nonnull)type vector:(double * __nonnull)vector solution:(FEMSolution * __nonnull)solution {
    
    double coeff=0.0, rtmp[3];
    static BOOL visited = NO;
    
    if (visited == NO) {
        if ([(solution.solutionInfo)[@"angles in degrees"] boolValue] == YES) {
            coeff = 180.0 / M_PI;
        } else {
            coeff = 1.0;
        }
        visited = YES;
    }
    
    if ([type isEqualToString:@"cartesian to cylindrical"] == YES) {
        rtmp[0] = sqrt(pow(vector[0], 2.0) + pow(vector[1], 2.0));
        rtmp[1] = coeff * atan2(vector[1], vector[0]);
        rtmp[2] = vector[2];
    } else if ([type isEqualToString:@"cylindrical to cartesiean"] == YES) {
        rtmp[0] = cos(vector[1] / coeff) * vector[0];
        rtmp[1] = sin(vector[1] / coeff) * vector[0];
        rtmp[2] = vector[2];
    } else {
        fprintf(stderr, "FEMMeshUtils:FEMMeshUtils_coordinateTransformationNodalType: unknown transformation: %s.\n", [type UTF8String]);
        fatal("FEMMeshUtils:FEMMeshUtils_coordinateTransformationNodalType");
    }
    memcpy(vector, rtmp, sizeof(rtmp));
}

/************************************************************
 
    Method corresponds to Elmer from git on October 27 2015
 
************************************************************/
-(void)FEMMeshUtils_markHaloNodesMesh:(FEMMesh * __nonnull)mesh haloNodes:(BOOL * __nullable)haloNodes sizeHaloNodes:(int * __nullable)size {
    
    BOOL foundHaloNodes;
    
    // Check whether we need to skip some elements and nodes on the halo boundary.
    // We don't want to create additional nodes on the nodes that are on the halo only
    // since they would just create further need for new halo...
    foundHaloNodes = NO;
    
    //TODO: add support for parallem run.
    // The rest of the method implementation is only needed for parrallel runs
}

/************************************************************
 
    Method corresponds to Elmer from git on October 27 2015
 
************************************************************/
-(Element_t * __nullable)FEMMeshUtils_findFaceParent:(Element_t * __nonnull)parent element:(Element_t * __nonnull)element model:(FEMModel * __nonnull)model {
    
    int n;
    Element_t *ptr = NULL;
    
    FEMMesh *mesh = (FEMMesh *)model.mesh;
    Element_t *faces = mesh.getFaces;
    for (int i=0; i<parent->Type.NumberOfFaces; i++) {
        ptr = &faces[parent->FaceIndexes[i]];
        n = 0;
        for (int j=0; j<ptr->Type.NumberOfNodes; j++) {
            for (int k=0; k<element->Type.NumberOfNodes; k++) {
                if (ptr->NodeIndexes[j] == element->NodeIndexes[k]) n++;
            }
        }
        if (n == ptr->Type.NumberOfNodes) break;
    }
    
    return ptr;
}

/**********************************************************************************************
 
    Create a Galerkin projector related to discontinous interface. This uses the information 
    stored when the discontinuous interface was first coined. This enables simple one-to-one 
    mapping. Integration weight is used for the nodel projector to allow physical jump 
    conditions. For the edge dofs there is no such jumps and hence the projector uses
    weights of one.
 
    Method corresponds to Elmer from git on October 27 2015

**********************************************************************************************/
-(FEMMatrix * __nullable)FEMMeshUtils_weightedProjectorDiscontinousMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model boundary:(int)bc listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    int actSides, e1, e2, e12, indp=0, i1, i2, j1, j2, n, nn, **newMap, **oldMap, parentFound, parentMissing, posSides;
    double coeff, detJ, x, u, v, val, w, weight;
    double *basis = NULL, **basisFirstDerivative = NULL;
    FEMMatrix *projector;
    FEMBoundaryCondition *bcParams;
    FEMSolution *solution;
    bool *edgeDone = NULL, *haloNodes = NULL;
    BOOL all, allLeft, allRight, axisSym, checkHaloNodes=NO, doEdges, doNodes, found, localConstraints=NO, nodalJump=NO, noHalo=NO, setDiag,
         setDiagEdges, stat;
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: creating projector for discontinous boundary %d.\n", bc);
    
    if (mesh.isDiscontinuousMesh == NO) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: discontinuous mesh not created.\n");
        return nil;
    }
    
    int j = 0;
    for (FEMBoundaryCondition *boundary in model.boundaryConditions) {
        if ([listUtilities listGetLogical:model inArray:boundary.valuesList forVariable:@"discontinuous boundary" info:&found] == YES) j++;
    }
    if (j > 1) fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: only one BC (not %d) for discontinuous boundary.\n", j);
    
    bcParams = (model.boundaryConditions)[bc];
    double scale = [listUtilities listGetConstReal:model inArray:bcParams.valuesList forVariable:@"mortar bc scaling" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) scale = -1.0;
    
    nodalJump = [listUtilities listCheckPrefix:@"mortar bc coefficient" inArray:bcParams.valuesList];
    if (nodalJump == NO) [listUtilities listCheckPrefix:@"mortar bc resistivity" inArray:bcParams.valuesList];
    
    // Take the full weight when creating the constraints since the value will
    // not be communicated
    if (model.solution != nil) solution = (FEMSolution *)model.solution;
    if (solution != nil) {
        if ((solution.solutionInfo)[@"partition local projector"] != nil) {
            localConstraints = [(solution.solutionInfo)[@"partition local projector"] boolValue];
        } else {
            if ((solution.solutionInfo)[@"partition local constraints"] != nil) {
                localConstraints = [(solution.solutionInfo)[@"partition local constraints"] boolValue];
            }
        }
    }
    
    // Don't consider halo when creating discontinuity
    if (solution != nil) {
        if((solution.solutionInfo)[@"projector no halo"] != nil) {
            noHalo = [(solution.solutionInfo)[@"projector no halo"] boolValue];
        }
    }
    
    // Don't consider single halo nodes when creating discontinuity
    if (solution != nil) {
        if ((solution.solutionInfo)[@"projector no halo nodes"] != nil) {
            checkHaloNodes = [(solution.solutionInfo)[@"projector no halo nodes"] boolValue];
        }
    }
    if (checkHaloNodes == YES) {
        // TODO: add the call to FEMMeshUtils_markHaloNodesMesh:haloNodes:sizeHaloNodes: for parallel run
    }
    
    doEdges = NO;
    if (solution != nil) {
        if([(solution.solutionInfo)[@"projector skip edges"] boolValue] == YES) {
            doEdges = NO;
        }
    } else if ([listUtilities listGetLogical:model inArray:bcParams.valuesList forVariable:@"projector skip edges" info:&found] == YES) {
        doEdges = NO;
    } else {
        doEdges = (mesh.numberOfEdges > 0) ? YES : NO;
    }
    if (doEdges == YES && mesh.numberOfEdges == 0) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: edge basis requested but mesh has no edges.\n");
        doEdges = NO;
    }
    
    doNodes = NO;
    if (solution != nil) {
        if ([(solution.solutionInfo)[@"projector skip nodes"] boolValue] == YES) {
            doNodes = NO;
        }
    } else if ([listUtilities listGetLogical:model inArray:bcParams.valuesList forVariable:@"projector skip nodes" info:&found] == YES) {
        doNodes = NO;
    } else {
        doNodes = (mesh.numberOfNodes > 0) ? YES : NO;
    }
    
    // Should the projector be diagonal or mass matrix type
    setDiag = [listUtilities listGetLogical:model inArray:bcParams.valuesList forVariable:@"mortar bc diag" info:&found];
    if (found == NO) setDiag = [listUtilities listGetLogical:model inArray:bcParams.valuesList forVariable:@"use biorthogonal basis" info:&found];
    
    // If we want to eliminate the constraints, we have to hava a biorthogonal basis
    if (found == NO) {
        if (solution != nil) {
            if ((solution.solutionInfo)[@"eliminate linear constraints"] != nil) {
                setDiag = [(solution.solutionInfo)[@"eliminate linear constraints"] boolValue];
            }
        }
        if (setDiag == YES) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: setting > use biorthogonal basis < to YES to enable elimination.\n");
        }
    }
    
    setDiagEdges = [listUtilities listGetLogical:model inArray:bcParams.valuesList forVariable:@"mortar bc diag edges" info:&found];
    if (found == NO) setDiagEdges = setDiag;
    
    // Integration weights should follow the metrics if we want physical nodal jumps
    axisSym = NO;
    if (model.coordinates == axis_symmetric || model.coordinates == cylindric_symmetric) {
        if (nodalJump == YES) {
            axisSym = YES;
        } else if (solution != nil) {
            if ((solution.solutionInfo)[@"projector metrics"] != nil) {
                axisSym = [(solution.solutionInfo)[@"projector metrics"] boolValue];
            }
        }
        if (axisSym == YES) fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: projector will be weighted for axi-symmetry.\n");
    }
    
    Nodes_t *elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    elementNodes->x = doublevec(0, mesh.maxElementDofs-1);
    elementNodes->y = doublevec(0, mesh.maxElementDofs-1);
    elementNodes->z = doublevec(0, mesh.maxElementDofs-1);
    
    int *indexes = intvec(0, mesh.maxElementDofs-1);
    int *disContIndexes = intvec(0, mesh.maxElementDofs-1);
    
    double **wBasis = doublematrix(0, mesh.maxElementDofs-1, 0, 2);
    double **wBasis2 = doublematrix(0, mesh.maxElementDofs-1, 0, 2);
    double **rotWBasis = doublematrix(0, mesh.maxElementDofs-1, 0, 2);
    
    memset(indexes, -1, mesh.maxElementDofs*sizeof(int) );
    memset(disContIndexes, -1, mesh.maxElementDofs*sizeof(int) );
    
    int *nodePerm = mesh.getDiscontinousPerm;
    int noOrigNodes = mesh.numberOfNodes;
    int noDiscontNodes = 0;
    for (int i=0; i<noOrigNodes; i++) {
        if (nodePerm[i] >= 0) noDiscontNodes++;
    }
    
    int indpoffset;
    if (doNodes == YES) {
        indpoffset = noDiscontNodes;
    } else {
        indpoffset = 0;
    }
    int *invPerm = NULL;
    int invPermSize = indpoffset;
    
    FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
    
    // Compute the number of potential edges. This mimics the loop that really creates the
    // projector below
    Element_t *elements = mesh.getElements, *left = NULL, *newFace = NULL, *oldFace = NULL, *swap = NULL, *right = NULL;
    if (doEdges == YES) {
        edgeDone = boolvec(0, mesh.numberOfEdges-1);
        memset(edgeDone, 0, n*sizeof(bool));
        indp = indpoffset;
        
        for (int t=0; t<mesh.numberOfBoundaryElements; t++) {
            
            if (elements[mesh.numberOfBulkElements+t].BoundaryInfo->Constraint != bcParams.tag) continue;
            left = elements[mesh.numberOfBulkElements+t].BoundaryInfo->Left;
            right = elements[mesh.numberOfBulkElements+t].BoundaryInfo->Right;
            
            if (left == NULL || right == NULL) continue;
            
            actSides = 0;
            if (left != NULL) {
                // TODO: add support for parallel run
            }
            if (right != NULL) {
                // TODO: add support for parallel run
            }
            if (noHalo == YES && actSides == 0) continue;
            
            // Consistently choose the face with the old edges
            allLeft = YES;
            for (int i=0; i<left->Type.NumberOfNodes; i++) {
                if (left->NodeIndexes[i] >= noOrigNodes) {
                    allLeft = NO;
                    break;
                }
            }
            allRight = YES;
            for (int i=0; i<right->Type.NumberOfNodes; i++) {
                if (right->NodeIndexes[i] >= noOrigNodes) {
                    allRight = NO;
                    break;
                }
            }
            if (allLeft == YES) {
                oldFace = left;
            } else if (allRight == YES) {
                oldFace = right;
            } else {
                fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: neither face is purely old.\n");
                continue;
            }
            
            oldMap = [elementDescription getEdgeMap:oldFace->Type.ElementCode/100];
            for (int i=0; i<oldFace->Type.NumberOfEdges; i++) {
                e1 = oldFace->EdgeIndexes[i];
                if (edgeDone[e1] == true) continue;
                
                i1 = oldFace->NodeIndexes[oldMap[i][0]];
                i2 = oldFace->NodeIndexes[oldMap[i][1]];
                
                // i1 and i2 were already checked to be "old" nodes
                if (nodePerm[i1] < 0) continue;
                if (nodePerm[i2] < 0) continue;
                
                indp++;
                edgeDone[e1] = true;
            }
        }
        invPermSize = indp;
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: size of invPerm estimated to be %d.\n", invPermSize);
        free_bvector(edgeDone, 0, mesh.numberOfEdges-1);
    }
    
    GaussIntegrationPoints *IP = NULL;
    
    FEMNumericIntegration *integration = [[FEMNumericIntegration alloc] init];
    if ([integration allocation:mesh] == NO) fatal("FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh", "Allocation error in FEMNumericIntegration.");
    
    FEMListMatrix *listmatrix = [[FEMListMatrix alloc] init];
    
    basis = integration.basis;
    basisFirstDerivative = integration.basisFirstDerivative;
    
    Nodes_t *nodes = mesh.getNodes;

    // Ok, nothing to do just go and tidy things up
    if (invPermSize == 0) goto jump;
    
    // Create a list of matrix that allows for unspecified entries in the matrix
    // structure to be introduced
    projector = [[FEMMatrix alloc] init];
    projector.format = MATRIX_LIST;
    projector.projectorType = PROJECTOR_TYPE_GALERKIN;
    projector.projectorBC = bc;
    
    // Create the inverse permutation needed when the projector matrix is added
    // to the global matrix
    matrixArraysContainer *matrixContainers = projector.getContainers;
    matrixContainers->InvPerm = intvec(0, invPermSize-1);
    matrixContainers->sizeInvPerm = invPermSize;
    memset(matrixContainers->InvPerm, -1, invPermSize*sizeof(int) );
    
    // Projector for the nodal dofs
    if (doNodes == YES) {
        parentMissing = 0;
        parentFound = 0;
        for (int t=0; t<mesh.numberOfBoundaryElements; t++) {
            n = elements[mesh.numberOfBulkElements+t].Type.NumberOfNodes;
            memcpy(indexes, elements[mesh.numberOfBulkElements+t].NodeIndexes, n*sizeof(int));
            
            if (elements[mesh.numberOfBulkElements+t].BoundaryInfo->Constraint != bcParams.tag) continue;
            
            left = elements[mesh.numberOfBulkElements+t].BoundaryInfo->Left;
            right = elements[mesh.numberOfBulkElements+t].BoundaryInfo->Right;
            
            posSides = 0;
            actSides = 0;
            if (left != NULL) {
                posSides++;
                // TODO: add support for parallel run
            }
            if (right != NULL) {
                posSides++;
                // TODO: add support for parallel run
            }
            
            if (noHalo == YES && actSides == 0) continue;
            
            if (localConstraints == YES) {
                coeff = 1.0;
            } else {
                coeff = (1.0 * actSides) / posSides;
            }
            if (fabs(coeff) < DBL_MIN) continue;
            
            parentFound++;
            
            for (int i=0; i<n; i++) {
                elementNodes->x[i] = nodes->x[indexes[i]];
                elementNodes->y[i] = nodes->y[indexes[i]];
                elementNodes->z[i] = nodes->z[indexes[i]];
            }
            
            all = YES;
            for (int i=0; i<n; i++) {
                if (nodePerm[indexes[i]] >= 0) {
                    all = NO;
                    break;
                }
            }
            if (all == YES) continue;
            
            if (checkHaloNodes == YES) {
                all = YES;
                for (int i=0; i<n; i++) {
                    if (haloNodes[indexes[i]] != true) {
                        all = NO;
                    }
                }
                if (all == YES) continue;
            }
            
            // Get the indexes on the other side of the discontinuous boundary
            for (int i=0; i<n; i++) {
                j = nodePerm[indexes[i]];
                if (j < 0) {
                    disContIndexes[i] = indexes[i];
                } else {
                    disContIndexes[i] = j + noOrigNodes;
                }
            }
            
            IP = GaussQuadrature(&elements[mesh.numberOfBulkElements+t], NULL, NULL);
            for (j=0; j<IP->n; j++) {
                u = IP->u[j];
                v = IP->v[j];
                w = IP->w[j];
                
                stat = [integration setBasisForElement:&elements[mesh.numberOfBulkElements+t] elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                stat = [integration setMetricDeterminantForElement:&elements[mesh.numberOfBulkElements+t] elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
                detJ = integration.metricDeterminant;
                
                weight = coeff * detJ * IP->s[j];
                if (axisSym == YES) {
                    x = cblas_ddot(n, integration.basis, 1, elementNodes->x, 1);
                    weight = weight * x;
                }
                for (int p=0; p<n; p++) {
                    indp = nodePerm[indexes[p]];
                    if (indp < 0) continue;
                    if (checkHaloNodes == YES) {
                        if (haloNodes[indexes[p]] == true) continue;
                    }
                    val = weight * basis[p];
                    
                    // Only set for the nodes are really used
                    invPerm[indp] = indexes[p];
                    
                    if (setDiag == YES) {
                        [listmatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:indp andIndex:indexes[p] value:val setValue:NULL];
                        [listmatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:indp andIndex:disContIndexes[p] value:scale*val setValue:NULL];
                    } else {
                        for (int q=0; q<n; q++) {
                            indp = nodePerm[indexes[p]];
                            if (indp < 0) continue;
                            
                            if (checkHaloNodes == YES) {
                                if (haloNodes[indexes[p]] == true) continue;
                            }
                            [listmatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:indp andIndex:indexes[q] value:basis[q]*val setValue:NULL];
                            [listmatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:indp andIndex:disContIndexes[q] value:scale*basis[q]*val setValue:NULL];
                        }
                    }
                }
            }
        }
        if (parentMissing > 0) {
            //TODO: add support for parallel run
        }
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: created projector for %d discontinuous nodes.\n", noDiscontNodes);
    }
    
    // Create the projector also for edge dofs if they exist and are
    // requested
    if (doEdges == YES) {
        FEMPElementMaps *elementMaps = [[FEMPElementMaps alloc] init];
        parentMissing = 0;
        parentFound = 0;
        n = mesh.numberOfNodes;
        
        val = 1.0;
        scale = 1.0;
        
        indp = indpoffset;
        int *eqInd = intvec(0, mesh.numberOfEdges-1);
        memset(eqInd, -1, mesh.numberOfEdges*sizeof(int));
        
        FEMInterpolation *interpolation = [[FEMInterpolation alloc] init];
        double point[3], sum, uvw[3];
        
        for (int t=0; t<mesh.numberOfBoundaryElements; t++) {
            
            if (elements[mesh.numberOfBulkElements+t].BoundaryInfo->Constraint != bcParams.tag) continue;
            
            left = elements[mesh.numberOfBulkElements+t].BoundaryInfo->Left;
            right = elements[mesh.numberOfBulkElements+t].BoundaryInfo->Right;
            
            // Here we really need both sides to be able to continue
            if (left == NULL || right == NULL) {
                parentMissing++;
                continue;
            }
            
            posSides = 0;
            actSides = 0;
            if (left != NULL) {
                posSides++;
                // TODO: add support for parallel run
            }
            if (right != NULL) {
                posSides++;
                // TODO: add support for parallel run
            }
            if (noHalo == YES && actSides == 0) continue;
            
            if (localConstraints == YES) {
                coeff = 1.0;
            } else {
                coeff = (1.0 * actSides) / (1.0 * posSides);
            }
            
            // Consistently choose the face with the old edges
            allLeft = YES;
            for (int i=0; i<left->Type.NumberOfNodes; i++) {
                if (left->NodeIndexes[i] >= noOrigNodes) {
                    allLeft = NO;
                    break;
                }
            }
            allRight = YES;
            for (int i=0; i<right->Type.NumberOfNodes; i++) {
                if (right->NodeIndexes[i] >= noOrigNodes) {
                    allRight = NO;
                    break;
                }
            }
            if (allLeft == YES) {
                // No ops!
            } else if (allRight == YES) {
                swap = left;
                left = right;
                right = swap;
            } else {
                // We already complained once
                continue;
            }

            oldFace = [self FEMMeshUtils_findFaceParent:left element:&elements[mesh.numberOfBulkElements+t] model:model];
            nn = elements[mesh.numberOfBulkElements+t].sizeNodeIndexes;
            memcpy(indexes, elements[mesh.numberOfBulkElements+t].NodeIndexes, nn*sizeof(int));
            for (int i=0; i<nn; i++) {
                elements[mesh.numberOfBulkElements+t].NodeIndexes[i] = nodePerm[indexes[i]] + noOrigNodes;
            }
            newFace = [self FEMMeshUtils_findFaceParent:right element:&elements[mesh.numberOfBulkElements+t] model:model];
            memcpy(elements[mesh.numberOfBulkElements+t].NodeIndexes, indexes, nn*sizeof(int));

            parentFound++;
            oldMap = [elementDescription getEdgeMap:oldFace->Type.ElementCode/100];
            newMap = [elementDescription getEdgeMap:newFace->Type.ElementCode/100];
            
            IP = GaussQuadrature(oldFace, NULL, NULL);
            for (int it=0; it<IP->n; it++) {
                u = IP->u[it];
                v = IP->v[it];
                w = IP->w[it];
                
                nn = oldFace->Type.NumberOfNodes;
                for (int i=0; i<nn; i++) {
                    elementNodes->x[i] = nodes->x[oldFace->NodeIndexes[i]];
                    elementNodes->y[i] = nodes->y[oldFace->NodeIndexes[i]];
                    elementNodes->z[i] = nodes->z[oldFace->NodeIndexes[i]];
                }
                stat = [integration setBasisForElement:oldFace elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                stat = [integration setBasisFirstDerivativeForElement:oldFace elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                stat = [integration setMetricDeterminantForElement:oldFace elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
                
                [elementDescription getEdgeBasisElement:oldFace wBasis:wBasis rotWBasis:rotWBasis basis:basis dBasisdx:basisFirstDerivative];
                
                point[0] = cblas_ddot(nn, integration.basis, 1, elementNodes->x, 1);
                point[1] = cblas_ddot(nn, integration.basis, 1, elementNodes->y, 1);
                point[2] = cblas_ddot(nn, integration.basis, 1, elementNodes->z, 1);
                
                nn = newFace->Type.NumberOfNodes;
                for (int i=0; i<nn; i++) {
                    elementNodes->x[i] = nodes->x[newFace->NodeIndexes[i]];
                    elementNodes->y[i] = nodes->y[newFace->NodeIndexes[i]];
                    elementNodes->z[i] = nodes->z[newFace->NodeIndexes[i]];
                }
                
                found = [interpolation isPointInElement:newFace elementNodes:elementNodes point:point localCoordinates:uvw globalEpsilon:NULL localEpsilon:NULL numericEpsilon:NULL globalDistance:NULL localDistance:NULL model:model elementDescription:elementDescription elementMaps:elementMaps edgeBasis:NULL];
                u = uvw[0]; v = uvw[1]; w = uvw[2];
                stat = [integration setBasisForElement:newFace elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                stat = [integration setBasisFirstDerivativeForElement:newFace elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w withBubbles:NO basisDegree:NULL];
                stat = [integration setMetricDeterminantForElement:newFace elementNodes:elementNodes inMesh:mesh firstEvaluationPoint:u secondEvaluationPoint:v thirdEvaluationPoint:w];
                detJ = integration.metricDeterminant;
                
                [elementDescription getEdgeBasisElement:newFace wBasis:wBasis2 rotWBasis:rotWBasis basis:basis dBasisdx:basisFirstDerivative];
                
                weight = detJ * IP->s[it] * coeff;
                
                // Go through combinations of edges and find the edges for which the
                // indexes are th same
                for (int i=0; i<oldFace->Type.NumberOfEdges; i++) {
                    e1 = oldFace->EdgeIndexes[i];
                    
                    if (eqInd[e1] < 0) {
                        eqInd[e1] = indp;
                        invPerm[indp] = n + e1;
                        indp++;
                    }
                    
                    if (setDiagEdges == YES) {
                        i1 = oldFace->NodeIndexes[oldMap[i][0]];
                        i1 = noOrigNodes + nodePerm[i1];
                        i2 = oldFace->NodeIndexes[oldMap[i][1]];
                        i2 = noOrigNodes + nodePerm[i2];
                        
                        for (j=0; j<newFace->Type.NumberOfEdges; j++) {
                            j1 = newFace->NodeIndexes[newMap[j][0]];
                            j2 = newFace->NodeIndexes[newMap[j][1]];
                            if ((i1 == j1 && i2 == j2) || (i1 == j2 && i2 == j1)) break;
                        }
                        sum = 0.0;
                        for (int k=0; k<3; k++) {
                            sum = sum + wBasis[i][k] * wBasis[i][k];
                        }
                        val = weight * sum;
                        if (fabs(val) >= 10.0*AEPS) [listmatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:eqInd[e1] andIndex:n+e1 value:val setValue:NULL];
                        
                        e2 = newFace->EdgeIndexes[j];
                        sum = 0.0;
                        for (int k=0; k<3; k++) {
                            sum = sum + wBasis[i][k] * wBasis2[j][k];
                        }
                        val = weight * sum;
                        if (fabs(val) >= 10.0*AEPS) [listmatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:eqInd[e1] andIndex:n+e2 value:-val setValue:NULL];
                    } else {
                        for (int j=0; j<newFace->Type.NumberOfEdges; j++) {
                            e2 = newFace->EdgeIndexes[j];
                            e12 = oldFace->EdgeIndexes[j];
                            
                            sum = 0.0;
                            for (int k=0; k<3; k++) {
                                sum = sum + wBasis[i][k] * wBasis[j][k];
                            }
                            val = weight * sum;
                            if (fabs(val) >= 10.0*AEPS) [listmatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:eqInd[e1] andIndex:n+e12 value:val setValue:NULL];
                            
                            sum = 0.0;
                            for (int k=0; k<3; k++) {
                                sum = sum + wBasis[i][k] * wBasis2[j][k];
                            }
                            val = weight * sum;
                            if (fabs(val) >= 10.0*AEPS) [listmatrix addToMatrixElement:matrixContainers->ListMatrix atIndex:eqInd[e1] andIndex:n+e2 value:-val setValue:NULL];
                        }
                    }
                }
            }
        }
        
        free_ivector(eqInd, 0, mesh.numberOfEdges-1);
        if (doNodes == NO && parentMissing > 0) {
            //TODO: add support for parallel run
        }
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: created projector for %d discontinuous edges.\n", indp-noDiscontNodes);
        [elementMaps deallocation];
    }
    
    // Convert from list matrix to CRS matrix format
    [listmatrix convertToCRSMatrix:projector];
    
    if (projector.numberOfRows > 0) {
        FEMMatrixCRS *crsMatrix = [[FEMMatrixCRS alloc] init];
        BOOL sortValues = YES;
        [crsMatrix sortMatrix:projector alsoValues:&sortValues];
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_weightedProjectorDiscontinousMesh: number of entries in projector matrix: %d.\n", matrixContainers->sizeCols);
    } else {
        [projector deallocation];
        projector = nil;
    }
    
jump:
    
    free_dvector(elementNodes->x, 0, mesh.maxElementDofs-1);
    free_dvector(elementNodes->y, 0, mesh.maxElementDofs-1);
    free_dvector(elementNodes->z, 0, mesh.maxElementDofs-1);
    free(elementNodes);
    
    free_ivector(indexes, 0, mesh.maxElementDofs-1);
    free_ivector(disContIndexes, 0, mesh.maxElementDofs-1);
    
    free_dmatrix(wBasis, 0, mesh.maxElementDofs-1, 0, 2);
    free_dmatrix(wBasis2, 0, mesh.maxElementDofs-1, 0, 2);
    free_dmatrix(rotWBasis, 0, mesh.maxElementDofs-1, 0, 2);
    
    [integration deallocation:mesh];
    
    return projector;
}

/*************************************************************************
 
    Create a nodal projector related to discontinuous interface

    Method corresponds to Elmer from git on October 27 2015

*************************************************************************/
-(FEMMatrix * __nonnull)FEMMeshUtils_nodalProjectorDiscontinuousMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model boundary:(int)bc listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    BOOL found;
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_nodalProjectorDiscontinuousMesh: creating nodal projector for discontinous boudnary.\n");
    
    if (mesh.discontinuousMesh == YES) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_nodalProjectorDiscontinuousMesh: discontinous mesh not created.\n");
        return nil;
    }
    
    int j = 0;
    for (FEMBoundaryCondition *boundary in model.boundaryConditions) {
        if ([listUtilities listGetLogical:model inArray:boundary.valuesList forVariable:@"discontinuous boundary" info:&found] == YES) j++;
    }
    // This is a temporal limitation
    if (j > 1) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_nodalProjectorDiscontinuousMesh: only one BC (not %d) for discontinuous boundary.\n", j);
    }
    
    int *nodePerm = mesh.getDiscontinousPerm;
    int n = mesh.numberOfNodes;
    int m = 0;
    for (int i=0; i<n; i++) {
        if (nodePerm[i] >= 0) m++;
    }
    
    FEMMatrix *projector = [[FEMMatrix alloc] init];
    projector.projectorType = PROJECTOR_TYPE_NODAL;
    projector.projectorBC = bc;
    
    matrixArraysContainer *matrixContainers = projector.getContainers;
    matrixContainers->Cols = intvec(0, m-1);
    matrixContainers->sizeCols = m;
    matrixContainers->Values = doublevec(0, m-1);
    matrixContainers->sizeValues = m;
    matrixContainers->Rows = intvec(0, (m+1)-1);
    matrixContainers->sizeRows = m + 1;
    matrixContainers->InvPerm = intvec(0, m-1);
    matrixContainers->sizeInvPerm = m;
    projector.numberOfRows = m;
    
    memset(matrixContainers->Values, 0.0, m*sizeof(double) );
    for (int i=0; i<m+1; i++) {
        matrixContainers->Rows[i] = i;
    }
    for (int i=0; i<n; i++) {
        j = nodePerm[i];
        if (j < 0) continue;
        matrixContainers->Cols[j] = n + j;
        matrixContainers->InvPerm[j] = i;
    }
    return projector;
}

/**************************************************************************************************************
 
    Create master and slave mesh for the interface in order to at a later stage create projector
    matrix to implement periodicity or mortar elements. The idea is to use a reduced set of elements 
    and thereby speed up the mapping process. Also this gives more flexibility in transformation
    operations since the nodes may be ereased after use. Interface meshes consist of boundary elements only.
 
    Method corresponds to Elmer from git on October 27 2015

**************************************************************************************************************/
-(BOOL)FEMMeshUtils_createInterfaceMeshesModel:(FEMModel * __nonnull)model mesh:(FEMMesh * __nonnull)mesh masterBoundary:(int)mbd targetBoundary:(int)trgt bMesh1:(FEMMesh * __nonnull)bMesh1 bMesh2:(FEMMesh * __nonnull)bMesh2 {
    
    int k1, k2, n, constraint, en, haloCount, ind;
    int *pPerm;
    FEMMesh *pMesh;
    Element_t *bMesh1Elements = NULL, *bMesh2Elements = NULL, *elements = NULL, *face = NULL, *left = NULL, *parent = NULL,
              *pMeshElements = NULL, *right = NULL;
    bool *activeNodes = NULL;
    BOOL any, checkForHalo, narrowHalo, noHalo, onTheFlyBC, success, targetActive, thisActive;
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_createInterfaceMeshesModel: making a list of elements at interface.\n");
    
    // BCs index should start from 0
    if (mbd < 0 || trgt < 0) {
        fatal("FEMMeshUtils:FEMMeshUtils_createInterfaceMeshesModel", "Invalid target boundaries.");
    }
    
    // If the target is larger than number of given boundaries given then
    // it has probably been created on-the-fly from a discontinuous boundary
    onTheFlyBC = ((trgt+1) > model.numberOfBoundaryConditions) ? YES : NO;
    
    // If parallel we may have some excess halo elements.
    // To eliminate them, mark the nodes that are associated to elements truly owned
    narrowHalo = NO;
    noHalo = NO;
    
    // TODO: add support for parallel run
    
    // This is just temporily set to false always until the logic has been tested
    checkForHalo = (narrowHalo == YES || noHalo == YES) ? YES : NO;
    
    elements = mesh.getElements;
    
    if (checkForHalo == YES) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_createInterfaceMeshesModel: checking for halo elements.\n");
        activeNodes = boolvec(0, mesh.numberOfNodes-1);
        haloCount = 0;
        memset(activeNodes, 0, mesh.numberOfNodes*sizeof(bool));
        for (int i=0; i<mesh.numberOfBoundaryElements; i++) {
            if (elements[mesh.numberOfBulkElements+i].Type.ElementCode <= 200) continue;
            
            left = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Left;
            if (left != NULL) {
                // TODO: add support for parallel run
            }
            right = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Right;
            if (right != NULL) {
                // TODO: add support for parallel run
            }
        }
        // No halo element found on the boundary so no need to check them later
        if (haloCount == 0) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_createInterfaceMeshesModel: found no halo elements to eliminate.\n");
            free_bvector(activeNodes, 0, mesh.numberOfNodes-1);
            checkForHalo = NO;
        } else {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_createInterfaceMeshesModel: number of halo to eliminate: %d.\n", haloCount);
        }
    }
    
    // Search elements in this boundary and its periodic
    // conterpart
    int n1 = 0;
    int n2 = 0;
    haloCount = 0;
    FEMBoundaryCondition *bcParamsMaster = (model.boundaryConditions)[mbd];
    FEMBoundaryCondition *bcParamsTarget = (model.boundaryConditions)[trgt];
    for (int i=0; i<mesh.numberOfBoundaryElements; i++) {
        if (elements[mesh.numberOfBulkElements+i].Type.ElementCode <= 200) continue;
        
        constraint = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Constraint;
        if (bcParamsMaster.tag == constraint) {
            if (checkForHalo == YES) {
                if (narrowHalo == YES) {
                    any = NO;
                    for (int j=0; j<elements[mesh.numberOfBulkElements+i].Type.NumberOfNodes; j++) {
                        if (activeNodes[elements[mesh.numberOfBulkElements+i].NodeIndexes[j]] == true) {
                            any = YES;
                            break;
                        }
                    }
                    if (any == YES) {
                        n1++;
                    } else {
                        haloCount++;
                    }
                } else if (noHalo == YES) {
                    thisActive = noHalo;
                    left = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Left;
                    if (left != NULL) {
                        // TODO: add support for parallel run
                    }
                    right = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Right;
                    if (right != NULL) {
                        // TODO: add support for parallel run
                    }
                    if (thisActive == YES) {
                        n1++;
                    } else {
                        haloCount++;
                    }
                }
            } else {
                n1++;
            }
        }
        
        if (onTheFlyBC == YES) {
            if ((trgt+1) == constraint) n2++;
        } else {
            if (bcParamsTarget.tag == constraint) n2++;
        }
    }
    
    if (checkForHalo == YES) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_createInterfaceMeshesModel: number of halo elements eliminated: %d.\n", haloCount);
    }
    
    if (n1 <= 0 || n2 <= 0) {
        // This is too conservative in parallel
        return success = NO;
    }
    
    // Intialize mesh structures for boundaries, this is
    // for getting the mesh projector
    bMesh1.parent = mesh;
    bMesh2.parent = mesh;
    
    bMesh1Elements = (Element_t*) malloc( sizeof(Element_t) * n1 );
    bMesh2Elements = (Element_t*) malloc( sizeof(Element_t) * n2 );
    
    int *perm1 = intvec(0, mesh.numberOfNodes-1);
    int *perm2 = intvec(0, mesh.numberOfNodes-1);
    
    // Fill in the mesh element structures with the boundary elements
    n1 = 0;
    n2 = 0;
    memset(perm1, -1, mesh.numberOfNodes*sizeof(int));
    memset(perm2, -1, mesh.numberOfNodes*sizeof(int));
    bMesh1.maxElementNodes = 0;
    bMesh2.maxElementNodes = 0;
    
    for (int i=0; i<mesh.numberOfBoundaryElements; i++) {
        
        constraint = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Constraint;
        thisActive = (bcParamsMaster.tag == constraint) ? YES : NO;
        if (thisActive == YES && checkForHalo == YES) {
            if (narrowHalo == YES) {
                any = NO;
                for (int j=0; j<elements[mesh.numberOfBulkElements+i].Type.NumberOfNodes; j++) {
                    if (activeNodes[elements[mesh.numberOfBulkElements+i].NodeIndexes[j]] == false) {
                        any = YES;
                        break;
                    }
                }
                if (any == YES) thisActive = NO;
            } else if (noHalo == YES) {
                thisActive = NO;
                left = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Left;
                if (left != NULL) {
                    // TODO: add support for parallel run
                }
                right = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Right;
                if (right != NULL) {
                    // TODO: add support for parallel run
                }
            }
        }
        
        if (onTheFlyBC == YES) {
            targetActive = ((trgt+1) == constraint) ? YES : NO;
        } else {
            targetActive = (bcParamsTarget.tag == constraint) ? YES : NO;
        }
        
        if (!(thisActive || targetActive)) continue;
        
        //Set the pointer accordingly so we need to code the complex stuff only once
        if (thisActive == YES) {
            ind = n1;
            pMesh = bMesh1;
            pMeshElements = bMesh1Elements;
            pPerm = perm1;
            n1++;
        } else {
            ind = n2;
            pMesh = bMesh2;
            pMeshElements = bMesh2Elements;
            pPerm = perm2;
            n2++;
        }
        
        n = elements[mesh.numberOfBulkElements+i].Type.NumberOfNodes;
        pMesh.maxElementNodes = max(pMesh.maxElementNodes, n);
        pMeshElements[ind] = elements[mesh.numberOfBulkElements+i];
        
        pMeshElements[ind].NodeIndexes = intvec(0, n-1);
        pMeshElements[ind].sizeNodeIndexes = n;
        
        if (mesh.numberOfFaces == 0 || mesh.numberOfEdges == 0) {
            memcpy(pMeshElements[ind].NodeIndexes, elements[mesh.numberOfBulkElements+i].NodeIndexes, n*sizeof(int));
        } else {
            // If we have edge dofs, we want the face element be associated with the
            // face list since that only has properlly defined edge indexes
            parent = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Left;
            if (parent == NULL) {
                parent = elements[mesh.numberOfBulkElements+i].BoundaryInfo->Right;
            }
            face = [self FEMMeshUtils_findFaceParent:parent element:&elements[mesh.numberOfBulkElements+i] model:model];
            memcpy(pMeshElements[ind].NodeIndexes, face->NodeIndexes, n*sizeof(int));
            
            // Set the element index to be face index as it may be needed
            // for the edge elements
            pMeshElements[ind].ElementIndex = face->ElementIndex;
            
            if (face->Pdefs != NULL) {
                pMeshElements[ind].Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
                memcpy(pMeshElements[ind].Pdefs, face->Pdefs, sizeof(PElementDefs_t));
            }
            
            en = face->Type.NumberOfEdges;
            pMeshElements[ind].EdgeIndexes = intvec(0, en-1);
            pMeshElements[ind].sizeEdgeIndexes = en;
            memcpy(pMeshElements[ind].EdgeIndexes, face->EdgeIndexes, en*sizeof(int));
        }
        
        for (int j=0; j<n; j++) {
            pPerm[elements[mesh.numberOfBulkElements+i].NodeIndexes[j]] = 0;
        }
    }
    
    // Fill in the mesh node structures with the boundary nodes
    bMesh1.numberOfBulkElements = n1;
    bMesh2.numberOfBulkElements = n2;
    
    bMesh1.numberOfNodes = 0;
    for (int i=0; i<mesh.numberOfNodes; i++) {
        if (perm1[i] >= 0) bMesh1.numberOfNodes++;
    }
    bMesh2.numberOfNodes = 0;
    for (int i=0; i<mesh.numberOfNodes; i++) {
        if (perm2[i] >= 0) bMesh2.numberOfNodes++;
    }
    
    // As there were some active boundary elements, this consition should
    // really never be possible
    if (bMesh1.numberOfNodes == 0 || bMesh2.numberOfNodes == 0) {
        fprintf(stderr, "FEMMeshUtils:FEMMeshUtils_createInterfaceMeshesModel: no active nodes on periodic boundary.\n");
        fatal("FEMMeshUtils:FEMMeshUtils_createInterfaceMeshesModel");
    }
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_createInterfaceMeshesModel: number of periodic nodes: %d, %d.\n", bMesh1.numberOfNodes, bMesh2.numberOfNodes);
    
    Nodes_t *meshNodes = mesh.getNodes;
    
    Nodes_t *mesh1Nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    mesh1Nodes->x = doublevec(0, bMesh1.numberOfNodes-1);
    mesh1Nodes->y = doublevec(0, bMesh1.numberOfNodes-1);
    mesh1Nodes->z = doublevec(0, bMesh1.numberOfNodes-1);
    mesh1Nodes->numberOfNodes =  bMesh1.numberOfNodes;
    
    Nodes_t *mesh2Nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    mesh2Nodes->x = doublevec(0, bMesh2.numberOfNodes-1);
    mesh2Nodes->y = doublevec(0, bMesh2.numberOfNodes-1);
    mesh2Nodes->z = doublevec(0, bMesh2.numberOfNodes-1);
    mesh2Nodes->numberOfNodes = bMesh2.numberOfNodes;
    
    int *invPermMesh1 = intvec(0, bMesh1.numberOfNodes-1);
    int *invPermMesh2 = intvec(0, bMesh2.numberOfNodes-1);
    memset(invPermMesh1, -1, bMesh1.numberOfNodes*sizeof(int));
    memset(invPermMesh2, -1, bMesh2.numberOfNodes*sizeof(int));
    
    // Now create the master and target meshes that only include the active elements
    k1 = 0; k2 = 0;
    for (int i=0; i<mesh.numberOfNodes; i++) {
        
        if (perm1[i] >= 0) {
            perm1[i] = k1;
            invPermMesh1[k1] = i;
            
            mesh1Nodes->x[k1] = meshNodes->x[i];
            mesh1Nodes->y[k1] = meshNodes->y[i];
            mesh1Nodes->z[k1] = meshNodes->z[i];
            k1++;
        }
        
        if (perm2[i] >= 0) {
            perm2[i] = k2;
            invPermMesh2[k2] = i;
            
            mesh2Nodes->x[k2] = meshNodes->x[i];
            mesh2Nodes->y[k2] = meshNodes->y[i];
            mesh2Nodes->z[k2] = meshNodes->z[i];
            k2++;
        }
    }
    
    // Finally, renumber the element node pointers to use
    // only boundary nodes
    for (int i=0; i<n1; i++) {
        for (int j=0; j<bMesh1Elements[i].Type.NumberOfNodes; j++) {
            bMesh1Elements[i].NodeIndexes[j] = perm1[bMesh1Elements[i].NodeIndexes[j]];
        }
    }
    for (int i=0; i<n2; i++) {
        for (int j=0; j<bMesh2Elements[i].Type.NumberOfNodes; j++) {
            bMesh2Elements[i].NodeIndexes[j] = perm2[bMesh2Elements[i].NodeIndexes[j]];
        }
    }
    
    [bMesh1 assignElements:bMesh1Elements];
    [bMesh1 assignNodes:mesh1Nodes];
    [bMesh1 assignInvPerm:invPermMesh1];
    
    [bMesh2 assignElements:bMesh2Elements];
    [bMesh2 assignNodes:mesh2Nodes];
    [bMesh2 assignInvPerm:invPermMesh2];
    
    free_ivector(perm1, 0, mesh.numberOfNodes-1);
    free_ivector(perm2, 0, mesh.numberOfNodes-1);
    
    if (checkForHalo == YES) free_bvector(activeNodes, 0, mesh.numberOfNodes-1);
    
    return success = YES;
}

/***************************************************************************
 
    Given a permutation map the (x,y,z) such that the projector can better 
    be applied. E.g. if boundary has constant x, take that as the last 
    coordinate.
 
    Method corresponds to Elmer from git on October 27 2015
 
***************************************************************************/
-(void)FEMMeshUtils_mapInterfaceCoordinateMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 BCParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    double *nodesX = NULL, *nodesY = NULL, *nodesZ = NULL;
    Nodes_t *bMeshNodes = NULL;
    listBuffer coordMap = { NULL, NULL, NULL, NULL, 0, 0, 0};
    BOOL all1, all2, all3, found;
    
    // Perform coordinate mapping
    found = [listUtilities listGetIntegerArray:model inArray:bcParams.valuesList forVariable:@"projector coordinate mapping" buffer:&coordMap];
    if (found == NO) return;
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_mapInterfaceCoordinateMesh1: performing coordinate mapping.\n");
    
    // Check if all values are different than 1
    all1 = YES;
    for (int i=0; i<3; i++) {
        if (coordMap.ivector[i] == 1) {
            all1 = NO;
            break;
        }
    }
    // Check if all values are different than 2
    all2 = YES;
    for (int i=0; i<3; i++) {
        if (coordMap.ivector[i] == 2) {
            all2 = NO;
            break;
        }
    }
    // Check if all values are different than 3
    all3 = YES;
    for (int i=0; i<3; i++) {
        if (coordMap.ivector[i] == 3) {
            all3 = NO;
            break;
        }
    }
    
    if (coordMap.m != 3 || (all1 == YES || all2 == YES || all3 == YES)) {
        fprintf(stderr, "FEMMeshUtils:FEMMeshUtils_mapInterfaceCoordinateMesh1: inconsistent coordinate mapping: \n");
        for (int i=0; i<coordMap.m; i++) {
            fprintf(stderr, "FEMMeshUtils:FEMMeshUtils_mapInterfaceCoordinateMesh1: %d\n", coordMap.ivector[i]);
        }
        fprintf(stderr, "FEMMeshUtils:FEMMeshUtils_mapInterfaceCoordinateMesh1: coordinate mapping should be a permutation of 1, 2 and 3.\n");
        fatal("FEMMeshUtils:FEMMeshUtils_mapInterfaceCoordinateMesh1");
    }
    
    for (int meshNo=1; meshNo<=2; meshNo++) {
        if (meshNo == 1) {
            bMeshNodes = bMesh1.getNodes;
        } else {
            bMeshNodes = bMesh2.getNodes;
        }
        
        if (coordMap.ivector[0] == 1) {
            nodesX = bMeshNodes->x;
        } else if (coordMap.ivector[0] == 2) {
             nodesX = bMeshNodes->y;
        } else {
             nodesX = bMeshNodes->z;
        }
        
        if (coordMap.ivector[1] == 1) {
            nodesY = bMeshNodes->x;
        } else if (coordMap.ivector[1] == 2) {
            nodesY = bMeshNodes->y;
        } else {
            nodesY = bMeshNodes->z;
        }
        
        if (coordMap.ivector[2] == 1) {
            nodesZ = bMeshNodes->x;
        } else if (coordMap.ivector[2] == 2) {
            nodesZ = bMeshNodes->y;
        } else {
            nodesZ = bMeshNodes->z;
        }
        
        bMeshNodes->x = nodesX;
        bMeshNodes->y = nodesY;
        bMeshNodes->z = nodesZ;
    }
    
    if (coordMap.ivector != NULL) {
        free_ivector(coordMap.ivector, 0, coordMap.m-1);
    }
}

/***********************************************************************************
 
    Simply fitting of cylinder into a point cloud. This is done in two phases.
    1) The axis of the cylinder is found by minimizing the \sum((n_i*t)^2)
      for each component of of t where n_i:s are the surface normals.
      This is fully generic and assumes no positions.
    2) The radius and center point of the cylinder are found by fitting a circle
      in the chosen plane to three representative points. Currently the fitting
      can only be done in x-y plane.
 
    Method corresponds to Elmer from git on October 27 2015
 
***********************************************************************************/
-(void)FEMMeshUtils_cylinderFitMesh:(FEMMesh * __nonnull)pMesh pParams:(FEMBoundaryCondition * __nonnull)pParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    int i, n, axisI, circleInd[3];
    double **a, axisNormal[3], coord[3], dist, d1, d2, circleCoord[3][3], ninj[3][3], maxDist, minDist, normals[3], sum, tangent1[3], tangent2[3];
    BOOL check;
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_cylinderFitMesh: trying to fit a cylinder to the surface patch.\n");
    
    FEMElementUtils *elementUtils = [[FEMElementUtils alloc] init];
    FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
    FEMLinearAlgebra *linearAlgebra = [[FEMLinearAlgebra alloc] init];
    
    Element_t *elements = pMesh.getElements;
    Nodes_t *meshNodes = pMesh.getNodes;
    
    memset(*ninj, 0.0, (3*3)*sizeof(double));
    
    Nodes_t *nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    nodes->x = doublevec(0, pMesh.maxElementNodes-1);
    nodes->y = doublevec(0, pMesh.maxElementNodes-1);
    nodes->z = doublevec(0, pMesh.maxElementNodes-1);
    
    a = doublematrix(0, 2, 0, 2);
    
    // If the initial mesh is in 2D, there is really no need to figure out the
    // direction of the rotational axis. It can only be aligned with the z-axis
    FEMMesh *mesh = (FEMMesh *)model.mesh;
    if (mesh.dimension == 2) {
        memset(axisNormal, 0.0, sizeof(axisNormal));
        axisNormal[2] = 1.0;
        goto jump;
    }
    
    // Compute the inner product of <N*N> for the elements
    check = NO;
    for (int t=0; t<pMesh.numberOfBulkElements; t++) {
        n = elements[t].Type.NumberOfNodes;
        for (i=0; i<n; i++) {
            nodes->x[i] = meshNodes->x[elements[t].NodeIndexes[i]];
            nodes->y[i] = meshNodes->y[elements[t].NodeIndexes[i]];
            nodes->z[i] = meshNodes->z[elements[t].NodeIndexes[i]];
        }
        
        [elementDescription normalVectorForBDElement:&elements[t] boundaryNodes:nodes mesh:pMesh paraU:NULL paraV:NULL check:&check normals:normals];
        for (i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                ninj[i][j] = ninj[i][j] + normals[i] * normals[j];
            }
        }
    }
        
    // Normalize by the number of bulk elements
    for (i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            ninj[i][j] = ninj[i][i] / pMesh.numberOfBulkElements;
        }
    }
    
    // The potential direction for the cylinder axis is the direction with least hits for the normal
    axisI = 0;
    for (i=1; i<3; i++) {
        if (ninj[i][i < ninj[axisI][axisI]]) axisI = i;
    }
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_cylinderFitMesh: axis coordinate set to be: %d.\n", axisI);
    
    // Keep the dominating direction fixed and iteratively solve the two other directories
    memset(axisNormal, 0.0, sizeof(axisNormal));
    axisNormal[axisI] = 1.0;
    
    // Basically we could solve from equation Ax=0 the tangent but only up to a constant.
    // Thus we enforce the axis direction to one by manipulation by the manipulation of
    // the matix equation, thereby can get a unique solution
    memcpy(*a, *ninj, (3*3)*sizeof(double));
    for (i=0; i<3; i++) {
        a[axisI][i] = 0.0;
    }
    a[axisI][axisI] = 1.0;
    [linearAlgebra invertMatrix:a ofSize:3];
    for (i=0; i<3; i++) {
        axisNormal[i] = a[i][axisI];
    }
    
    // Normalize the axis normal length to one
    sum = 0.0;
    for (i=0; i<3; i++) {
        sum = sum + pow(axisNormal[i], 2.0);
    }
    for (i=0; i<3; i++) {
        axisNormal[i] = axisNormal[i] /  sqrt(sum);
    }
    if (1.0 - fabs(axisNormal[2]) > 1.0e-5) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_cylinderFitMesh: the cylinder axis is not aligned with z-axis.\n");
    }
    
jump:
    
    [elementUtils tangentDirectionsForNormal:axisNormal tangent1:tangent1 tangent2:tangent2];
    
    // Finding three points with maximum distance in the tangent directions
    
    // First, find the single extremum point in the first tangent direction
    // Save the local coordinates in the N-T system of the cylinder
    minDist = HUGE_VAL;
    for (i=0; i<pMesh.numberOfNodes; i++) {
        coord[0] = meshNodes->x[i];
        coord[1] = meshNodes->y[i];
        coord[2] = meshNodes->z[i];
        
        d1 = 0.0;
        for (int j=0; j<3; j++) {
            d1 = d1 + tangent1[j] * coord[j];
        }
        if (d1 < minDist) {
            minDist = d1;
            circleInd[0] = i;
        }
    }
    
    i = circleInd[0];
    coord[0] = meshNodes->x[i];
    coord[1] = meshNodes->y[i];
    coord[2] = meshNodes->z[i];
    
    circleCoord[0][0] = 0.0;
    for (int i=0; i<3; i++) {
        circleCoord[0][0] = circleCoord[0][0] + tangent1[i] * coord[i];
    }
    
    circleCoord[0][1] = 0.0;
    for (int i=0; i<3; i++) {
        circleCoord[0][1] = circleCoord[0][1] + tangent2[i] * coord[i];
    }
    
    circleCoord[0][2] = 0.0;
    for (int i=0; i<3; i++) {
        circleCoord[0][2] = circleCoord[0][2] + axisNormal[i] * coord[i];
    }
    
    // Find two more points such that their minimum distance to the previous point(s)
    // is maximized. This takes some time but the further the nodes are apart the more
    // accurate it will be to fit the circle to the points. Also if there is just
    // a symmetric section of the cylinder it is important to find the points rigorously.
    for (int j=2; j<=3; j++) {
        // The maximum minimum distance of any node from the previously defined nodes
        maxDist = 0;
        for (i=0; i<pMesh.numberOfNodes; i++) {
            coord[0] = meshNodes->x[i];
            coord[1] = meshNodes->y[i];
            coord[2] = meshNodes->z[i];
            
            // Minimum distance from the previously defined nodes
            minDist = HUGE_VAL;
            for (int k=0; k<j-1; k++) {
                d1 = 0.0;
                for (int l=0; l<3; l++) {
                    d1 = d1 + tangent1[l] * coord[l];
                }
                d2 = 0.0;
                for (int l=0; l<3; l++) {
                    d2 = d2 + tangent2[l] * coord[l];
                }
                dist = pow((d1 - circleCoord[k][0]), 2.0) + pow((d2 - circleCoord[k][1]), 2.0);
                minDist = min(dist, minDist);
            }
            
            // If the minimum distance is greater than in any other node, choose this
            if (maxDist < minDist) {
                maxDist = minDist;
                circleInd[j-1] = i;
            }
        }
        
        // Ok, we have found the point, now set the circle coordinates
        i = circleInd[j-1];
        coord[0] = meshNodes->x[i];
        coord[1] = meshNodes->y[i];
        coord[2] = meshNodes->z[i];
        
        circleCoord[j-1][0] = 0.0;
        for (int i=0; i<3; i++) {
            circleCoord[j-1][0] = circleCoord[j-1][0] + tangent1[i] * coord[i];
        }

        circleCoord[j-1][1] = 0.0;
        for (int i=0; i<3; i++) {
            circleCoord[j-1][1] = circleCoord[j-1][1] + tangent2[i] * coord[i];
        }

        circleCoord[j-1][2] = 0.0;
        for (int i=0; i<3; i++) {
            circleCoord[j-1][2] = circleCoord[j-1][2] + axisNormal[i] * coord[i];
        }
    }
    
    // Given three nodes it is possible to analytically compute the center point and
    // radius of the cylinder from a 4x4 determinant equation. The matrices values
    // m1i are the determinants of the comatrices.
    for (int i=0; i<3; i++) {
        a[i][0] = circleCoord[i][0]; // x
        a[i][1] = circleCoord[i][1]; // y
        a[i][2] = 1.0;
    }
    double m11 = det3x3(a);
    
    for (int i=0; i<3; i++) {
        a[i][0] = pow(circleCoord[i][0], 2.0) + pow(circleCoord[i][1], 2.0); // x^2+y^2
        a[i][1] = circleCoord[i][1]; // y
        a[i][2] = 1.0;
    }
    double m12 = det3x3(a);
    
    for (int i=0; i<3; i++) {
        a[i][0] = pow(circleCoord[i][0], 2.0) + pow(circleCoord[i][1], 2.0); // x^2+y^2
        a[i][1] = circleCoord[i][0]; // x
        a[i][2] = 1.0;
    }
    double m13 = det3x3(a);
    
    for (int i=0; i<3; i++) {
        a[i][0] = pow(circleCoord[i][0], 2.0) + pow(circleCoord[i][1], 2.0); // x^2+y^2
        a[i][1] = circleCoord[i][0]; // x
        a[i][2] = circleCoord[i][1]; // y
    }
    //double m14 = det3x3(a);
    
    if (fabs(m11) < DBL_EPSILON) {
        fatal("FEMMeshUtils:FEMMeshUtils_cylinderFitMesh", "Points can not be a circle.");
    }
    
    double x0 = 0.5 * m12 / m11;
    double y0 = -0.5 * m13 / m11;
    //double rad = sqrt(pow(x0, 2.0) + pow(y0, 2.0) + m14 / m11);
    
    for (int i=0; i<3; i++) {
        coord[i] = x0 * tangent1[i] + y0 * tangent2[i];
    }
    
    [listUtilities addConstRealInClassList:pParams theVariable:@"rotational projector center x" withValue:&coord[0] orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:pParams theVariable:@"rotational projector center y" withValue:&coord[1] orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:pParams theVariable:@"rotational projector center z" withValue:&coord[2] orUsingBlock:nil string:nil];
    
    [listUtilities addConstRealInClassList:pParams theVariable:@"rotational projector normal x" withValue:&axisNormal[0] orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:pParams theVariable:@"rotational projector normal y" withValue:&axisNormal[1] orUsingBlock:nil string:nil];
    [listUtilities addConstRealInClassList:pParams theVariable:@"rotational projector normal z" withValue:&axisNormal[2] orUsingBlock:nil string:nil];
    
    free_dmatrix(a, 0, 2, 0, 2);
    free_dvector(nodes->x, 0, pMesh.maxElementNodes-1);
    free_dvector(nodes->y, 0, pMesh.maxElementNodes-1);
    free_dvector(nodes->z, 0, pMesh.maxElementNodes-1);
    free(nodes);
}

/****************************************************************************************

    Given two interface meshes for nonconforming rotating boundaries, make a coordinate
    transformation to (phi,z) level where the interpolation accuracy is not limited by 
    the curvilinear coordinates. Also ensure that the master nodes are manipulated so 
    that they for sure hit the target nodes.
 
    Method corresponds to Elmer from git on October 27 2015
 
****************************************************************************************/
-(void)FEMMeshUtils_rotationInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams cylindrical:(BOOL)cylindrical radius:(double)radius fullCircle:(BOOL * __nonnull)fullCircle model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    int n;
    double alpha, dfii1, dfii2, degOffset=0.0, eps_rad, err1, err2, fii, fii0=0.0, fmin, fmax, nSymmetry, normals[3], rad, tangent1[3], tangent2[3], x[3],
           xtmp[3], x0[3], x1_min[3], x1_max[3], x1r_min[3], x1r_max[3], x2_min[3], x2_max[3], x2r_min[3], x2r_max[3];
    FEMMesh *pMesh;
    Element_t *pMeshElements = NULL;
    Nodes_t *pMeshNodes = NULL;
    BOOL found, gotCenter, gotNormal, hit0, hit90, hit180, hit270, moveAngle = NO, setDegOffset;
    
    FEMElementUtils *elementUtils = [[FEMElementUtils alloc] init];
    
    // We choose degrees as they are more intuitive
    double rad2Deg = 180.0 / M_PI;
    int maxElementNodes = bMesh2.maxElementNodes;
    double *angles = doublevec(0, maxElementNodes-1);
    
    *fullCircle = NO;
    
    // Cylindrical projector is fitted always and rotational only when requested
    if ([listUtilities listGetLogical:model inArray:bcParams.valuesList forVariable:@"rotational projector center fit" info:&found] == YES ||
        cylindrical == YES) {
        if ([listUtilities listCheckPresentVariable:@"rotational projector center x" inArray:bcParams.valuesList] == NO) {
            [self FEMMeshUtils_cylinderFitMesh:bMesh1 pParams:bcParams model:model listUtilities:listUtilities];
        }
    }
    
    x0[0] = [listUtilities listGetConstReal:model inArray:bcParams.valuesList forVariable:@"rotational projector center x" info:&gotCenter minValue:NULL maxValue:NULL];
    x0[1] = [listUtilities listGetConstReal:model inArray:bcParams.valuesList forVariable:@"rotational projector center y" info:&found minValue:NULL maxValue:NULL];
    gotCenter = (gotCenter == YES || found == YES) ? YES : NO;
    x0[2] = [listUtilities listGetConstReal:model inArray:bcParams.valuesList forVariable:@"rotational projector center z" info:&found minValue:NULL maxValue:NULL];
    gotCenter = (gotCenter == YES || found == YES) ? YES : NO;
    
    normals[0] = [listUtilities listGetConstReal:model inArray:bcParams.valuesList forVariable:@"rotational projector normal x" info:&gotNormal minValue:NULL maxValue:NULL];
    normals[1] = [listUtilities listGetConstReal:model inArray:bcParams.valuesList forVariable:@"rotational projector normal y" info:&found minValue:NULL maxValue:NULL];
    gotNormal = (gotNormal == YES || found == YES) ? YES : NO;
    normals[2] = [listUtilities listGetConstReal:model inArray:bcParams.valuesList forVariable:@"rotational projector normal z" info:&found minValue:NULL maxValue:NULL];
    gotNormal = (gotNormal == YES || found == YES) ? YES : NO;
    
    if (gotNormal == YES) [elementUtils tangentDirectionsForNormal:normals tangent1:tangent1 tangent2:tangent2];
    
    // Go through master (k=1) and target mesh (k=2)
    for (int k=1; k<=2; k++) {
        
        // Potentially the projector may be set to rotate by just adding an offset
        // to the angle. This may depends on time
        if (k == 1) {
            degOffset = [listUtilities listGetConstReal:model inArray:bcParams.valuesList forVariable:@"rotational projector angle offset" info:&setDegOffset minValue:NULL maxValue:NULL];
        } else {
            setDegOffset = NO;
        }
        
        if (k == 1) {
            pMesh = bMesh1;
        } else {
            pMesh = bMesh2;
        }
        
        pMeshElements = pMesh.getElements;
        pMeshNodes = pMesh.getNodes;
        // Check the initial boundary boxes
        vDSP_minvD(pMeshNodes->x, 1, &x2_min[0], pMesh.numberOfNodes);
        vDSP_minvD(pMeshNodes->y, 1, &x2_min[1], pMesh.numberOfNodes);
        vDSP_minvD(pMeshNodes->z, 1, &x2_min[2], pMesh.numberOfNodes);
        
        vDSP_maxvD(pMeshNodes->x, 1, &x2_max[0], pMesh.numberOfNodes);
        vDSP_maxvD(pMeshNodes->y, 1, &x2_max[1], pMesh.numberOfNodes);
        vDSP_maxvD(pMeshNodes->z, 1, &x2_max[2], pMesh.numberOfNodes);
        
        if (k == 1) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: initial extrema for this bondary (x, ,y ,z).\n");
        } else if (k == 2) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: initial extrema for target bondary (x, ,y ,z).\n");
        }
        for (int i=0; i<3; i++) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: coordinate: %d: %f, %f.\n", i+1, x2_min[i], x2_max[i]);
        }
        
        // Memorize the bounding box of the master mesh
        if (k == 1) {
            memcpy(x1_min, x2_min, sizeof(x2_min));
            memcpy(x1_max, x2_max, sizeof(x2_max));
        }
        
        // Do the actual coordinate transformation
        n = pMesh.numberOfNodes;
        for (int i=0; i<n; i++) {
            x[0] = pMeshNodes->x[i];
            x[1] = pMeshNodes->y[i];
            x[2] = pMeshNodes->z[i];
            
            // Substract the center of axis
            if (gotCenter == YES) {
                for (int j=0; j<3; j++) {
                    x[j] = x[j] - x0[j];
                }
            }
            
            if (gotNormal == YES) {
                memcpy(xtmp, x, sizeof(x));
                x[0] = cblas_ddot(3, tangent1, 1, xtmp, 1);
                x[1] = cblas_ddot(3, tangent2, 1, xtmp, 1);
                x[2] = cblas_ddot(3, normals, 1, xtmp, 1);
            }
            
            // Set the angle to be the first coordinate as it may sometimes be the
            // only nonzero coordinate. z-coordinate is always unchanged
            alpha = rad2Deg * atan2(x[1], x[0]);
            rad = sqrt(pow(x[0], 2.0) + pow(x[1], 2.0));
            
            // Set the offset and revert then the angle to range [-180,180]
            if (setDegOffset == YES) {
                alpha = fmod((alpha + degOffset), 360.0);
                if (alpha > 180.0) alpha = alpha - 360.0;
            }
            
            pMeshNodes->x[i] = alpha;
            pMeshNodes->y[i] = x[2];
            pMeshNodes->z[i] = rad;
        }
        
        // For cyclindrical projector, follow exactly the same logic salve and master
        if (cylindrical == YES && k == 2) {
            if (moveAngle == YES) {
                fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: moving the second discontinuity to same angle.\n");
                for (int i=0; i<pMesh.numberOfNodes; i++) {
                    if (pMeshNodes->x[i] < fii0) pMeshNodes->x[i] = pMeshNodes->x[i] + 360.0;
                }
            }
        } else {
            // Let's see if we have a full angle to operate or not.
            // If not, then make the interval continuous.
            // Here we check only four critical angles: (0, 90. 180, 270) degs
            hit0 = NO; hit90 = NO; hit180 = NO; hit270 = NO;
            moveAngle = NO; fii = 0.0; fii0 = 0.0;
            
            for (int i=0; i<pMesh.numberOfBulkElements; i++) {
                n = pMeshElements[i].Type.NumberOfNodes;
                for (int j=0; j<n; j++) {
                    angles[j] = pMeshNodes->x[pMeshElements->NodeIndexes[j]];
                }
                vDSP_minvD(angles, 1, &fmin, n);
                vDSP_maxvD(angles, 1, &fmax, n);
                
                if (fmax - fmin > 180.0) {
                    hit180 = YES;
                } else {
                    if (fmax >= 0.0 && fmin <= 0.0) hit0 = YES;
                    if (fmax >= 90.0 && fmin <= 90.0) hit90 = YES;
                    if(fmax >= -90.0 && fmin <= -90.0) hit270 = YES;
                }
            }
            *fullCircle = (hit0 == YES && hit90 == YES && hit180 == YES && hit270 == YES) ? YES : NO;
            
            // Elininate the problematic discontinuity in case we have no full circle
            // The discontinuity will be moved to some of angles (-90, 0, 90)
            if (*fullCircle ==YES) {
                fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: cylindrical interface seems to be a full circle.\n");
            } else if (hit180 == YES) {
                moveAngle = YES;
                if (hit0 == NO) {
                    fii = 0.0;
                } else if (hit270 == NO) {
                    fii = -90.0;
                } else if (hit90 == NO) {
                    fii = 90.0;
                }
                for (int i=0; i<pMesh.numberOfNodes; i++) {
                    if (pMeshNodes->x[i] < fii) pMeshNodes->x[i] = pMeshNodes->x[i] + 360.0;
                }
                fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: moving discontinuity of angle to: %f.\n", fii);
                fii0 = fii;
            }
        }
        
        // Check the transformed bounding boxes
        vDSP_minvD(pMeshNodes->x, 1, &x2r_min[0], pMesh.numberOfNodes);
        vDSP_minvD(pMeshNodes->y, 1, &x2r_min[1], pMesh.numberOfNodes);
        vDSP_minvD(pMeshNodes->z, 1, &x2r_min[2], pMesh.numberOfNodes);
        
        vDSP_maxvD(pMeshNodes->x, 1, &x2r_max[0], pMesh.numberOfNodes);
        vDSP_maxvD(pMeshNodes->y, 1, &x2r_max[1], pMesh.numberOfNodes);
        vDSP_maxvD(pMeshNodes->z, 1, &x2r_max[2], pMesh.numberOfNodes);
        
        if (k == 1) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: transformed extrema for this boundary (phi, z, r).\n");
        } else if (k == 2) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: transformed extrema for target boundary (phi, z, r).\n");
        }
        for (int i=0; i<3; i++) {
             fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: coordinate: %d: %f, %f.\n", i+1, x2r_min[i], x2r_max[i]);
        }
        
        if (x2r_min[2] < DBL_EPSILON) {
            fatal("FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1", "Radius cannot be almost zero.");
        }
        
        // Memorize the bounding box for the 1st mesh
        if (k == 1) {
            memcpy(x1r_min, x2r_min, sizeof(x2r_min));
            memcpy(x1r_max, x2r_max, sizeof(x2r_max));
        }
    }
    
    eps_rad = 1.0e-3;
    
    // Choose radius to be max radius of thi boundary
    radius = x1r_max[2];
    
    err1 = (x1r_max[2] - x1r_min[2]) / radius;
    err2 = (x2r_max[2] - x2r_min[2]) / radius;
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: discrepancy from constant radius: %f.\n", err1);
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: discrepancy from constant radius: %f.\n", err2);
    if (err1 > eps_rad || err2 > eps_rad) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: discrepancy of radius is rather large.\n");
    }
    
    // Ok, so we have concluded that the interface has constant radius
    // therefore the constant radius may be removed from the mesh description.
    // Or perhaps we don't remove to allow more intelligent projector building
    // for contact mechanics.
    
    // Check whether the z-coordinate is constant or not.
    // Constant z-coordinate implies 1D system, otherwise 2D system.
    err1 = (x1r_max[1] - x1r_min[1]) / radius;
    err2 = (x2r_max[1] - x2r_min[1]) / radius;
    
    if (err1 < eps_rad || err2 < eps_rad) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: the effective interface meshes are 1D.\n");
        Nodes_t *bMesh1Nodes = bMesh1.getNodes;
        Nodes_t *bMesh2Nodes = bMesh2.getNodes;
        memset(bMesh1Nodes->y, 0.0, bMesh1.numberOfNodes*sizeof(double) );
        memset(bMesh2Nodes->y, 0.0, bMesh2.numberOfNodes*sizeof(double) );
    } else {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: the effective interface meshes are 2D.\n");
    }
    
    // Some pieces of the code can not ork with 1D meshes, this choice is ok for all steps
    bMesh1.dimension = 2;
    bMesh2.dimension = 2;
    
    // Cylindrical interfaces does not have symmetry as does the rotational
    if (cylindrical == YES || *fullCircle == YES) {
        free_dvector(angles, 0, maxElementNodes-1);
        return;
    };
    
    // If we are studying a symmeric segment, then analyze further the angle
    dfii1 = x1r_max[0] - x1r_min[0];
    dfii2 = x2r_max[0] - x2r_min[0];
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: this boundary dfii: %f.\n", dfii1);
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: this boundary dfii: %f.\n", dfii2);
    
    err1 = 2.0 * fabs(dfii1 - dfii2) / (dfii1 + dfii2);
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: discrepancy in dfii: %f.\n", err1);
    
    int i = [listUtilities listGetInteger:model inArray:bcParams.valuesList forVariable:@"rotational projector periods" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) {
        nSymmetry = 360.0 / dfii2;
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: suggested sections in target: %f.\n", nSymmetry);
        if (fabs(nSymmetry - round(nSymmetry)) > 0.01) {
            if (dfii1 < dfii2) {
                fprintf(stderr, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: you might try to switch master and target.\n");
            }
            fatal("FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1", "Check your settings, this can not be periodic.");
        }
        int value = round(nSymmetry);
        [listUtilities addIntegerInClassList:bcParams.valuesList theVariable:@"rotational projector periods" withValue:&value orUsingBlock:NULL];
    } else {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: using enforced number of periods: %d.\n", i);
        nSymmetry = 360.0 / dfii2;
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: suggested number of periods: %f.\n", nSymmetry);
    }
    
    free_dvector(angles, 0, maxElementNodes-1);
}

/*********************************************************************************
 
    Given two interface meshes for nonconforming radial boundaries, make
    a coordinate transformation to (r,z) level. This is always a symmetry 
    condition and can not be a contact condition.
 
    Method corresponds to Elmer from git on October 27 2015
 
*********************************************************************************/
-(void)FEMMeshUtils_radiusInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    double phi, phierr, r, x[3], x1_min[3], x1_max[3], x2_min[3], x2_max[3], z;
    FEMMesh *pMesh;
    Nodes_t *pMeshNodes = NULL;
    
    // We choose degrees a they are more intuitive
    double rad2Deg = 180.0 / M_PI;
    
    // Go through master (k=1) and target (k=2)
    for (int k=1; k<=2; k++) {
        
        if (k == 1) {
            pMesh = bMesh1;
        } else {
            pMesh = bMesh2;
        }
        pMeshNodes = pMesh.getNodes;
        
        for (int i=0; i<3; i++) {
            x2_min[i] = HUGE_VAL;
            x2_max[i] = -HUGE_VAL;
        }
        
        // Loop over all nodes
        for (int i=0; i<pMesh.numberOfNodes; i++) {
            x[0] = pMeshNodes->x[i];
            x[1] = pMeshNodes->y[i];
            x[2] = pMeshNodes->z[i];
            
            // Do the actual coordinate transformation
            r = sqrt(pow(x[0], 2.0) + pow(x[1], 2.0));
            phi = rad2Deg * atan2(x[1], x[0]);
            z = x[2];
            
            pMeshNodes->x[i] = r;
            pMeshNodes->y[i] = z;
            pMeshNodes->z[i] = 0.0;
            
            // This is just to cjeck a posteriori that the ranges are ok
            x2_min[0] = min(r, x2_min[0]);
            if (r > DBL_EPSILON) {
                x2_min[1] = min(phi, x2_min[1]);
            }
            x2_min[2] = min(z, x2_min[2]);
            
            x2_max[0] = max(r, x2_max[0]);
            if (r > DBL_EPSILON) {
                x2_max[1] = max(phi, x2_max[1]);
            }
            x2_max[2] = max(z, x2_max[2]);
        }
        
        // Memorize the bounding box of the master mesh
        if (k == 1) {
            memcpy(x1_min, x2_min, sizeof(x2_min));
            memcpy(x1_max, x2_max, sizeof(x2_max));
        }
        
        if (k == 1) {
            fprintf(stdout, "EMMeshUtils:FEMMeshUtils_radiusInterfaceMesh1: transformed extrema for this boundary (phi,r,z).\n");
        } else if (k == 2) {
            fprintf(stdout, "EMMeshUtils:FEMMeshUtils_radiusInterfaceMesh1: transformed extrema for target boundary (phi,r,z).\n");
        }
        for (int i=0; i<3; i++) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: coordinate: %d: %f, %f.\n", i+1, x2_min[i], x2_max[i]);
        }
        
        phierr = x2_max[1] - x2_min[1];
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: discrepancy from constant angle (degs): %f.\n", phierr);
    }
    
    // Error in radius
    // Choose radius to be max radius of either boundary
    double rad = max(x1_max[0], x2_max[0]);
    double err1 = fabs(x1_max[0] - x2_max[0]) / rad;
    double err2 = fabs(x1_min[0] - x2_min[0]) / rad;
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: discrepancy in maximum radius: %f.\n", err1);
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: discrepancy in minimum radius: %f.\n", err2);
    
    double eps_rad = 1.0e-3;
    if (err1 > eps_rad || err2 > eps_rad) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_rotationInterfaceMesh1: discrepancy of radius may be too large.\n");
    }
    
    // Some pieces of the code can not work with 1D meshes, this choice is ok for all steps
    bMesh1.dimension = 2;
    bMesh2.dimension = 2;
}

/***************************************************************
 
    Given two interface meshes flatten them to (x,y) plane.

    Method corresponds to Elmer from git on October 27 2015
 
***************************************************************/
-(void)FEMMeshUtils_flatInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    int minDiffI=0;
    double *coord = NULL, diff, maxDiff, minDiff, relDiff, relDiff1=0;
    FEMMesh *bMesh;
    Nodes_t *bMeshNodes = NULL;
    
    BOOL found;
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_flatInterfaceMesh1: flattening interface meshes to 2D.\n");
    
    int meshDim = model.dimension;
    int flatDim = [listUtilities listGetInteger:model inArray:bcParams.valuesList forVariable:@"flat projector coordinate" info:&found minValue:NULL maxValue:NULL];
    BOOL reduceDim = [listUtilities listGetInteger:model inArray:bcParams.valuesList forVariable:@"flat projector reduce dimension" info:&found minValue:NULL maxValue:NULL];
    
    if (found == NO) {
        double minVal, maxVal;
        for (int j=1; j<=2; j++) {
            
            if (j == 1) {
                bMesh = bMesh1;
            } else {
                bMesh = bMesh2;
            }
            bMeshNodes =  bMesh.getNodes;
            
            maxDiff = 0.0;
            minDiff = HUGE_VAL;
            
            for (int i=1; i<=meshDim; i++) {
                if (i == 1) {
                    coord = bMeshNodes->x;
                } else if (i == 2) {
                    coord = bMeshNodes->y;
                } else {
                    coord = bMeshNodes->z;
                }
                vDSP_minvD(coord, 1, &minVal, bMesh.numberOfNodes);
                vDSP_maxvD(coord, 1, &maxVal, bMesh.numberOfNodes);
                diff = maxVal - minVal;
                maxDiff = max(diff, maxDiff);
                if (diff < minDiff) {
                    minDiff = diff;
                    minDiffI = i;
                }
            }
            
            relDiff = minDiff / maxDiff;
            if (j == 1) {
                flatDim = minDiffI;
                relDiff1 = relDiff;
            } else if (j == 2) {
                if (relDiff < relDiff1) flatDim = minDiffI;
            }
        }
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_flatInterfaceMesh1: flat projector coordinate set to: %d.\n", flatDim);
        [listUtilities addIntegerInClassList:bcParams theVariable:@"flat projector coordinate" withValue:&flatDim orUsingBlock:nil];
    }
    
    for (int j=1; j<=2; j++) {
        
        if (j == 1) {
            bMesh = bMesh1;
        } else {
            bMesh = bMesh2;
        }
        bMeshNodes =  bMesh.getNodes;

        // Set the 3rd component to be the distance in the flat interface
        if (flatDim == 3) {
            continue;
        } else if (flatDim == 2) {
            coord = bMeshNodes->y;
            coord = bMeshNodes->y = bMeshNodes->z;
            coord = bMeshNodes->z = coord;
            if (meshDim == 2) memset( bMeshNodes->y, 0.0, bMesh.numberOfNodes*sizeof(double) );
        } else if (flatDim == 1) {
            coord = bMeshNodes->x;
            bMeshNodes->x = bMeshNodes->y;
            bMeshNodes->y =  bMeshNodes->z;
            bMeshNodes->z = coord;
            if (meshDim == 2) memset( bMeshNodes->y, 0.0, bMesh.numberOfNodes*sizeof(double) );
        }
        if (reduceDim == YES) memset( bMeshNodes->z, 0.0, bMesh.numberOfNodes*sizeof(double) );
        
        // Some pieces of the code can not work with 1D meshes, this choice is ok for ll steps
        bMesh.dimension = 2;
    }
}

/*******************************************************************
 
    Given two interface meshes flatten them into the plane that
    best fits either of the meshes.
 
    Method corresponds to Elmer from git on October 27 2015
 
*******************************************************************/
-(void)FEMMeshUtils_planeInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    int n;
    double coord[3], detJ, length, normals[3], normals0[3], normalSum[3], planeLess=0, planeLess1=0, **planeNormal, **planeNormal1, refSum, sum,
           tangent[3], tangent2[3];
    FEMMesh *bMesh;
    FEMNumericIntegration *integration;
    Element_t *bMeshElements = NULL;
    Nodes_t *bMeshNodes = NULL;
    listBuffer pNormals = { NULL, NULL, NULL, NULL, 0, 0, 0};
    GaussIntegrationPoints *IP = NULL;
    BOOL found, normal0Set, stat;
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1: flattening interface meshes to a place.\n");
    
    FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
    FEMElementUtils *elementUtils = [[FEMElementUtils alloc] init];
    
    int meshDim = model.dimension;
    found = [listUtilities listGetConstRealArray:model inArray:bcParams.valuesList forVariable:@"plane projector normal" buffer:&pNormals];
    
    // If the projector normal is not given, determine it first
    if (found == NO) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1: could not find place projector normal, so determining it.\n");
        Nodes_t *elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
        elementNodes->x = doublevec(0, MAX_ELEMENT_NODES-1);
        elementNodes->y = doublevec(0, MAX_ELEMENT_NODES-1);
        elementNodes->z = doublevec(0, MAX_ELEMENT_NODES-1);
        memset(elementNodes->x, 0.0, MAX_ELEMENT_NODES*sizeof(double) );
        memset(elementNodes->y, 0.0, MAX_ELEMENT_NODES*sizeof(double) );
        memset(elementNodes->z, 0.0, MAX_ELEMENT_NODES*sizeof(double) );
        
        planeNormal = doublematrix(0, 2, 0, 0);
        planeNormal1 = doublematrix(0, 2, 0, 0);
        
        BOOL check = NO;
        // Fit a plane to both datasets
        for (int j=1; j<=2; j++) {
            if (j == 1) {
                bMesh = bMesh1;
            } else {
                bMesh = bMesh2;
            }
            bMeshElements = bMesh.getElements;
            bMeshNodes = bMesh.getNodes;
            
            integration = [[FEMNumericIntegration alloc] init];
            if ([integration allocation:bMesh] == NO) fatal("FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1", "Allocation error in FEMNumericIntegration.");
            
            memset(normalSum, 0.0, sizeof(normalSum) );
            refSum = 0.0;
            normal0Set = NO;
            
            // We use the dot2Min and normal2 temporarily also for first mesh, with k=1
            for (int i=0; i<bMesh.numberOfBulkElements; i++) {
                n = bMeshElements[i].Type.NumberOfNodes;
                IP = GaussQuadrature(&bMeshElements[i], NULL, NULL);
                
                for (int k=0; k<n; k++) {
                    elementNodes->x[k] = bMeshNodes->x[bMeshElements[i].NodeIndexes[k]];
                    elementNodes->y[k] = bMeshNodes->y[bMeshElements[i].NodeIndexes[k]];
                    elementNodes->z[k] = bMeshNodes->z[bMeshElements[i].NodeIndexes[k]];
                }
                for (int nip=0; nip<IP->n; nip++) {
                    stat = [integration setBasisForElement:&bMeshElements[i] elementNodes:elementNodes inMesh:bMesh firstEvaluationPoint:IP->u[nip] secondEvaluationPoint:IP->v[nip] thirdEvaluationPoint:IP->w[nip] withBubbles:NO basisDegree:NULL];
                    stat = [integration setMetricDeterminantForElement:&bMeshElements[i] elementNodes:elementNodes inMesh:bMesh firstEvaluationPoint:IP->u[nip] secondEvaluationPoint:IP->v[nip] thirdEvaluationPoint:IP->w[nip]];
                    detJ = integration.metricDeterminant;
                    
                    [elementDescription normalVectorForBDElement:&bMeshElements[i] boundaryNodes:elementNodes mesh:bMesh paraU:&IP->u[nip] paraV:&IP->v[nip] check:&check normals:normals];
                    if (normal0Set == NO) {
                        memcpy(normals0, normals, sizeof(normals));
                        normal0Set = YES;
                    }
                    if (cblas_ddot(3, normals, 1, normals0, 1) < 0.0) {
                        for (int k=0; k<3; k++) {
                            normals[k] = -normals[k];
                        }
                    }
                    for (int k=0; k<3; k++) {
                        normalSum[k] = normalSum[k] + IP->s[nip] * detJ * normals[k];
                    }
                    refSum = refSum + IP->s[nip] * detJ;
                }
            }
            
            // Normalize the normal to unity length
            sum = 0.0;
            for (int k=0; k<3; k++) {
                sum = sum + pow(normalSum[k], 2.0);
            }
            length = sqrt(sum);
            for (int k=0; k<3; k++) {
                planeNormal[k][0] = normalSum[k] / length;
            }
            
            // Planeless is one if all the normals have the same direction
            planeLess = length / refSum;
            
            // Save the key parameters of the fist mesh
            if (j == 1) {
                memcpy(*planeNormal1, *planeNormal, (3*1)*sizeof(double));
                planeLess1 = planeLess;
            }
            
            [integration deallocation:bMesh];
        }
        
        // Choose the mesh for which is close to a plane
        if (planeLess1 > planeLess) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1: selecting slave normal.\n");
            memcpy(*planeNormal, *planeNormal1, (3*1)*sizeof(double));
        } else {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1: selecting master normal.\n");
            for (int i=0; i<3; i++) {
               planeNormal[i][0] = -planeNormal[i][0];
            }
        }
        for (int i=0; i<3; i++) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1: plane normal selected: %f.\n", planeNormal[i][0]);
        }
        [listUtilities addConstRealArrayInClassList:bcParams theVariable:@"plane projector normal" withValues:planeNormal size1:3 size2:1 orUsingBlock:nil string:nil];
        
        free_dmatrix(planeNormal, 0, 2, 0, 0);
        free_dmatrix(planeNormal1, 0, 2, 0, 0);
        
        free_dvector(elementNodes->x, 0, MAX_ELEMENT_NODES-1);
        free_dvector(elementNodes->y, 0, MAX_ELEMENT_NODES-1);
        free_dvector(elementNodes->z, 0, MAX_ELEMENT_NODES-1);
        free(elementNodes);
        
        found = [listUtilities listGetConstRealArray:model inArray:bcParams.valuesList forVariable:@"plane projector normal" buffer:&pNormals];
    }
    
    for (int i=0; i<3; i++) {
        normals[i] = pNormals.matrix[i][0];
    }
    [elementUtils tangentDirectionsForNormal:normals tangent1:tangent tangent2:tangent2];
    
    double minVal, maxVal;
    for (int j=1; j<=2; j++) {
        if (j == 1) {
            bMesh = bMesh1;
        } else {
            bMesh = bMesh2;
        }
        bMeshElements = bMesh.getElements;
        bMeshNodes = bMesh.getNodes;
        
        for (int i=0; i<bMesh.numberOfNodes; i++) {
            coord[0] = bMeshNodes->x[i];
            coord[1] = bMeshNodes->y[i];
            coord[2] = bMeshNodes->z[i];
            
            bMeshNodes->x[i] = cblas_ddot(3, coord, 1, tangent, 1);
            if (meshDim == 3) {
                bMeshNodes->y[i] = cblas_ddot(3, coord, 1, tangent2, 1);
            } else {
                bMeshNodes->y[i] = 0.0;
            }
            bMeshNodes->z[i] = cblas_ddot(3, coord, 1, normals, 1);
        }
        
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1: range for mesh: %d\n", j);
        vDSP_minvD(bMeshNodes->x, 1, &minVal, bMesh.numberOfNodes);
        vDSP_maxvD(bMeshNodes->x, 1, &maxVal, bMesh.numberOfNodes);
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1: x: %f %f\n", minVal, maxVal);
        vDSP_minvD(bMeshNodes->y, 1, &minVal, bMesh.numberOfNodes);
        vDSP_maxvD(bMeshNodes->y, 1, &maxVal, bMesh.numberOfNodes);
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1: y: %f %f\n", minVal, maxVal);
        vDSP_minvD(bMeshNodes->z, 1, &minVal, bMesh.numberOfNodes);
        vDSP_maxvD(bMeshNodes->z, 1, &maxVal, bMesh.numberOfNodes);
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_planeInterfaceMesh1: z: %f %f\n", minVal, maxVal);
        
        bMesh.dimension = 2;
    }
    
    if (pNormals.matrix != NULL) {
        free_dmatrix(pNormals.matrix, 0, pNormals.m-1, 0, pNormals.n-1);
    }
}

/*****************************************************************************************
 
    Given two interface meshes, check the angle between them using the normal
    vectors of the first element. Also check that all other elements are
    aligned with the first one. Only then is it possible to determine the angle.
 
    Method corresponds to Elmer from git on October 27 2015
 
*****************************************************************************************/
-(void)FEMMeshUtils_checkInterfaceAngleMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 angles:(double * __nonnull)angles gotAngles:(BOOL * __nonnull)gotAngles {
    
    int n;
    double alpha=0.0, dot1Min=0.0, dot2Min=0.0, normals[3], normal1[3], normal2[3];
    FEMMesh *pMesh;
    Element_t *pMeshElements = NULL;
    Nodes_t *pMeshNodes = NULL;
    
    FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
    
    // Currently check of the normal direction is not enforced since at this stage
    // the model nodes holding structure may not exist!
    // This means that there may be a 180 error in the directions.
    // Therefore an angle smaller than 180 is always chosen
    Nodes_t *elementNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    elementNodes->x = doublevec(0, max(bMesh1.maxElementNodes, bMesh2.maxElementNodes)-1);
    elementNodes->y = doublevec(0, max(bMesh1.maxElementNodes, bMesh2.maxElementNodes)-1);
    elementNodes->z = doublevec(0, max(bMesh1.maxElementNodes, bMesh2.maxElementNodes)-1);
    
    BOOL check = NO;
    for (int k=1; k<=2; k++) {
        if (k == 1) {
            pMesh = bMesh1;
        } else {
            pMesh = bMesh2;
        }
        pMeshElements = pMesh.getElements;
        pMeshNodes = pMesh.getNodes;
        
        // We use the dot2Min and normal2 temporarily also for first mesh, with k=1
        for (int i=0; i<pMesh.numberOfBoundaryElements; i++) {
            n = pMeshElements[i].Type.NumberOfNodes;
            for (int j=0; j<n; j++) {
                elementNodes->x[j] = pMeshNodes->x[pMeshElements[i].NodeIndexes[j]];
                elementNodes->y[j] = pMeshNodes->y[pMeshElements[i].NodeIndexes[j]];
                elementNodes->z[j] = pMeshNodes->z[pMeshElements[i].NodeIndexes[j]];
            }
            [elementDescription normalVectorForBDElement:&pMeshElements[i] boundaryNodes:elementNodes mesh:pMesh paraU:NULL paraV:NULL check:&check normals:normals];
            
            if (i == 0) {
                memcpy(normal2, normals, sizeof(normals));
                dot2Min = 1.0;
            } else {
                dot2Min = min(dot2Min, cblas_ddot(3, normals, 1, normal2, 1));
            }
        }
        
        if (k == 1) {
            memcpy(normal1, normal2, sizeof(normal2));
            dot1Min = dot2Min;
        }
    }
    
    BOOL constantNormals = ((1.0 - dot1Min < 1.0e-6) && (1.0 - dot2Min < 1.0e-6)) ? YES : NO;
    if (constantNormals == YES) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_checkInterfaceAngleMesh1: master normal: \n");
        for (int i=0; i<3; i++) {
            fprintf(stdout, "%f\n", normal1[i]);
        }
        
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_checkInterfaceAngleMesh1: initial target normal: \n");
        for (int i=0; i<3; i++) {
            fprintf(stdout, "%f\n", normal2[i]);
        }
        
        // The full angle between the two normals
        alpha = acos(cblas_ddot(3, normal1, 1, normal2, 1)) * 180.0 / M_PI;
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_checkInterfaceAngleMesh1: suggested angle between two normals in degs (+/- 180): %f.\n", alpha);
    } else {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_checkInterfaceAngleMesh1: could not suggest angle.\n");
    }
    
    *gotAngles = NO;
    if (constantNormals == NO) {
        fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_checkInterfaceAngleMesh1: normals are not constant, can not test for rotation.\n");
    } else if (alpha > DBL_EPSILON) {
        for (int i=0; i<3; i++) {
            if (fabs(normal1[i] - normal2[i]) < DBL_EPSILON) {
                *gotAngles = YES;
                fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_checkInterfaceAngleMesh1: rotation around axis %i in degs: %f\n", i+1, alpha);
                angles[i] = alpha;
                break;
            }
        }
        if (*gotAngles == NO) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_checkInterfaceAngleMesh1: could not define axis, improve algorithm.\n");
        }
    }
    
    free_dvector(elementNodes->x, 0,  max(bMesh1.maxElementNodes, bMesh2.maxElementNodes)-1);
    free_dvector(elementNodes->y, 0,  max(bMesh1.maxElementNodes, bMesh2.maxElementNodes)-1);
    free_dvector(elementNodes->z, 0,  max(bMesh1.maxElementNodes, bMesh2.maxElementNodes)-1);
    free(elementNodes);
}

/******************************************************************
 
    Given two meshes that should occupy the same domain in space,
    use rotation, scaling and translation to achive this goal.
 
    Method corresponds to Elmer from git on October 27 2015
 
******************************************************************/
-(void)FEMMeshUtils_overlayInterfaceMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 bcParams:(FEMBoundaryCondition * __nonnull)bcParams model:(FEMModel * __nonnull)model listUtilities:(FEMListUtilities * __nonnull)listUtilities {
    
    double alpha, angles[3], scl[3], x[4], x1_min[3], x1_max[3], x2_min[3], x2_max[3], x2r_min[3], x2r_max[3], y[4];
    Nodes_t *bMesh1Nodes = NULL, *bMesh2Nodes = NULL;
    listBuffer pArray = { NULL, NULL, NULL, NULL, 0, 0, 0};
    BOOL found;
    
    bMesh1Nodes = bMesh1.getNodes;
    bMesh2Nodes = bMesh2.getNodes;
    
    // First check the bounding boxes
    vDSP_minvD(bMesh1Nodes->x, 1, &x1_min[0], bMesh1.numberOfNodes);
    vDSP_minvD(bMesh1Nodes->y, 1, &x1_min[1], bMesh1.numberOfNodes);
    vDSP_minvD(bMesh1Nodes->z, 1, &x1_min[2], bMesh1.numberOfNodes);
    
    vDSP_maxvD(bMesh1Nodes->x, 1, &x1_max[0], bMesh1.numberOfNodes);
    vDSP_maxvD(bMesh1Nodes->y, 1, &x1_max[1], bMesh1.numberOfNodes);
    vDSP_maxvD(bMesh1Nodes->z, 1, &x1_max[2], bMesh1.numberOfNodes);
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: minimum values for this periodic BC: %e %e %e\n", x1_min[0], x1_min[1], x1_min[2]);
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: maximum values for this periodic BC: %e %e %e\n", x1_max[0], x1_max[1], x1_max[2]);
    
    vDSP_minvD(bMesh2Nodes->x, 1, &x2_min[0], bMesh2.numberOfNodes);
    vDSP_minvD(bMesh2Nodes->y, 1, &x2_min[1], bMesh2.numberOfNodes);
    vDSP_minvD(bMesh2Nodes->z, 1, &x2_min[2], bMesh2.numberOfNodes);

    vDSP_maxvD(bMesh2Nodes->x, 1, &x2_max[0], bMesh2.numberOfNodes);
    vDSP_maxvD(bMesh2Nodes->y, 1, &x2_max[1], bMesh2.numberOfNodes);
    vDSP_maxvD(bMesh2Nodes->z, 1, &x2_max[2], bMesh2.numberOfNodes);
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: minimum values for target periodic BC: %e %e %e\n", x2_min[0], x2_min[1], x2_min[2]);
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: maximum values for target periodic BC: %e %e %e\n", x2_max[0], x2_max[1], x2_max[2]);
    
    double **trfMatrix = doublematrix(0, 3, 0, 3);
    
    // If whole transformation matrix given, it will be used directly
    found = [listUtilities listGetConstRealArray:model inArray:bcParams.valuesList forVariable:@"periodic bc matrix" buffer:&pArray];
    if (found == YES) {
        for (int i=0; i<pArray.m; i++) {
            for (int j=0; j<pArray.n; j++) {
                trfMatrix[i][j] = pArray.matrix[j][i];
            }
        }
    } else {
        // Otherwise check for rotation, scaling and translation
        
        // Initialize the mapping matrices
        double **identity = doublematrix(0, 3, 0, 3);
        double **trsMatrix = doublematrix(0, 3, 0, 3);
        double **rotMatrix = doublematrix(0, 3, 0, 3);
        double **sclMatrix = doublematrix(0, 3, 0, 3);
        
        memset( *identity, 0.0, (4*4)*sizeof(double) );
        for (int i=0; i<4; i++) {
            identity[i][i] = 1.0;
        }
        memcpy(*trsMatrix, *identity, (4*4)*sizeof(double));
        memcpy(*rotMatrix, *identity, (4*4)*sizeof(double));
        memcpy(*sclMatrix, *identity, (4*4)*sizeof(double));
        
        // Rotations:
        // These are called first since they are not accounted for in
        // the automatic scaling and translation
        memset(angles, 0.0, sizeof(angles) );
        BOOL gotRotate = [listUtilities listGetConstRealArray:model inArray:bcParams.valuesList forVariable:@"periodic bc rotate" buffer:&pArray];
        if (gotRotate == YES) {
            for (int i=0; i<3; i++) {
                angles[i] = pArray.matrix[i][0];
            }
        } else {
            if ([listUtilities listGetLogical:model inArray:bcParams.valuesList forVariable:@"periodic bc rotate automatic" info:&found] == YES) {
                [self FEMMeshUtils_checkInterfaceAngleMesh1:bMesh1 mesh2:bMesh2 angles:angles gotAngles:&gotRotate];
            }
        }
        
        if (gotRotate == YES) {
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: rotating target with: %f %f %f\n", angles[0], angles[1], angles[2]);
            double **C = doublematrix(0, 3, 0, 3);
            for (int i=0; i<3; i++) {
                alpha = angles[i] * M_PI / 180.0;
                if (fabs(alpha) < DBL_MIN) continue;
                memcpy(*trfMatrix, *identity, (4*4)*sizeof(double));
                switch (i) {
                    case 0:
                        trfMatrix[1][1] = cos(alpha);
                        trfMatrix[1][2] = -sin(alpha);
                        trfMatrix[2][1] = sin(alpha);
                        trfMatrix[2][2] = cos(alpha);
                        break;
                    case 1:
                        trfMatrix[0][0] = cos(alpha);
                        trfMatrix[0][2] = -sin(alpha);
                        trfMatrix[2][0] = sin(alpha);
                        trfMatrix[2][2] = cos(alpha);
                        break;
                    case 2:
                        trfMatrix[0][0] = cos(alpha);
                        trfMatrix[0][1] = -sin(alpha);
                        trfMatrix[1][0] = sin(alpha);
                        trfMatrix[1][1] = cos(alpha);
                        break;
                }
                memset( *C, 0.0, (4*4)*sizeof(double) );
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, 4, 4, 1.0, *rotMatrix, 4, *trfMatrix, 4, 0.0, *C, 4);
                memcpy(*rotMatrix, *C, (4*4)*sizeof(double));
            }
            for (int i=0; i<bMesh2.numberOfNodes; i++) {
                x[0] = bMesh2Nodes->x[i];
                x[1] = bMesh2Nodes->y[i];
                x[2] = bMesh2Nodes->z[i];
                x[3] = 1.0;
                memset(y, 0.0, sizeof(y) );
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 4, 4, 1.0, *rotMatrix, 4, x, 1, 0.0, y, 1);
                bMesh2Nodes->x[i] = y[0];
                bMesh2Nodes->y[i] = y[1];
                bMesh2Nodes->z[i] = y[2];
            }
            vDSP_minvD(bMesh2Nodes->x, 1, &x2r_min[0], bMesh2.numberOfNodes);
            vDSP_minvD(bMesh2Nodes->y, 1, &x2r_min[1], bMesh2.numberOfNodes);
            vDSP_minvD(bMesh2Nodes->z, 1, &x2r_min[2], bMesh2.numberOfNodes);

            vDSP_maxvD(bMesh2Nodes->x, 1, &x2r_max[0], bMesh2.numberOfNodes);
            vDSP_maxvD(bMesh2Nodes->y, 1, &x2r_max[1], bMesh2.numberOfNodes);
            vDSP_maxvD(bMesh2Nodes->z, 1, &x2r_max[2], bMesh2.numberOfNodes);
            
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: minimum values for rotated target: %e %e %e\n", x2r_min[0], x2r_min[1], x2r_min[2]);
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: maximum values for rotated target: %e %e %e\n", x2r_max[0], x2r_max[1], x2r_max[2]);
            free_dmatrix(C, 0, 3, 0, 3);
        } else {
            memcpy(x2r_min, x2_min, sizeof(x2_min));
            memcpy(x2r_max, x2_max, sizeof(x2_max));
        }
        
        // Scaling:
        // This is either given or enforced by requiring bounding boxes to be of the same size
        found = [listUtilities listGetConstRealArray:model inArray:bcParams.valuesList forVariable:@"periodic bc scale" buffer:&pArray];
        if (found == YES) {
            for (int i=0; i<pArray.m; i++) {
                sclMatrix[i][i] = pArray.matrix[i][0];
            }
        } else {
            // Define scaling from the bounding boxes.
            // This assumes isotropic scaling since component-wise scaling
            // was prone to errors
            double s1 = 0.0;
            double s2 = 0.0;
            for (int i=0; i<3; i++) {
                s1 = s1 + pow(x1_max[i] - x1_min[i], 2.0);
            }
            for (int i=0; i<3; i++) {
                s2 = s2 + pow(x2r_max[i] - x2r_min[i], 2.0);
            }
            if (s2 > DBL_EPSILON) {
                for (int i=0; i<3; i++) {
                    scl[i] = sqrt(s1 / s2);
                }
            } else {
                for (int i=0; i<3; i++) {
                    scl[i] = 1.0;
                }
            }
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: scaling with: %e %e %e\n", scl[0], scl[1], scl[2]);
            for (int i=0; i<3; i++) {
                sclMatrix[i][i] = scl[i];
            }
        }
        
        // Translations:
        // And finally define translations
        found = [listUtilities listGetConstRealArray:model inArray:bcParams.valuesList forVariable:@"periodic bc translate" buffer:&pArray];
        if (found == YES) {
            for (int i=0; i<pArray.m; i++) {
                trsMatrix[3][i] = pArray.matrix[i][0];
            }
        } else {
            // Define translations so that the lower left corner is the same
            for (int i=0; i<3; i++) {
                trsMatrix[3][i] = x1_min[i] - sclMatrix[i][i] * x2r_min[i];
            }
            fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: translation: %e %e %e\n", trsMatrix[3][0], trsMatrix[3][1], trsMatrix[3][2]);
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 4, 4, 4, 1.0, *sclMatrix, 4, *trsMatrix, 4, 0.0, *trfMatrix, 4);
        
        free_dmatrix(identity, 0, 3, 0, 3);
        free_dmatrix(trsMatrix, 0, 3, 0, 3);
        free_dmatrix(rotMatrix, 0, 3, 0, 3);
        free_dmatrix(sclMatrix, 0, 3, 0, 3);
    }
    
    //Now transform the coordinates
    for (int i=0; i<bMesh2.numberOfNodes; i++) {
        x[0] = bMesh2Nodes->x[i];
        x[1] = bMesh2Nodes->y[i];
        x[2] = bMesh2Nodes->z[i];
        x[3] = 1.0;
        memset(y, 0.0, sizeof(y) );
        // The following does MATMUL(x,trfMatrix)
        cblas_dgemv(CblasRowMajor, CblasTrans, 4, 4, 1.0, *trfMatrix, 4, x, 1, 0.0, y, 1);
        bMesh2Nodes->x[i] = y[0] / y[3];
        bMesh2Nodes->y[i] = y[1] / y[3];
        bMesh2Nodes->z[i] = y[2] / y[3];
    }
    
    vDSP_minvD(bMesh2Nodes->x, 1, &x2r_min[0], bMesh2.numberOfNodes);
    vDSP_minvD(bMesh2Nodes->y, 1, &x2r_min[1], bMesh2.numberOfNodes);
    vDSP_minvD(bMesh2Nodes->z, 1, &x2r_min[2], bMesh2.numberOfNodes);

    vDSP_maxvD(bMesh2Nodes->x, 1, &x2r_max[0], bMesh2.numberOfNodes);
    vDSP_maxvD(bMesh2Nodes->y, 1, &x2r_max[1], bMesh2.numberOfNodes);
    vDSP_maxvD(bMesh2Nodes->z, 1, &x2r_max[2], bMesh2.numberOfNodes);
    
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: minimum values for transformed target: %e %e %e\n", x2r_min[0], x2r_min[1], x2r_min[2]);
    fprintf(stdout, "FEMMeshUtils:FEMMeshUtils_overlayInterfaceMesh1: maximum values for transformed target: %e %e %e\n", x2r_max[0], x2r_max[1], x2r_max[2]);

    free_dmatrix(trfMatrix, 0, 3, 0, 3);
    
    if (pArray.matrix != NULL) {
        free_dmatrix(pArray.matrix, 0, pArray.m-1, 0, pArray.n-1);
    }
}

/**************************************************************************************
 
    Create a projector for mapping between interfaces using the Galerkin method
    A temporal mesh structure with a node for each Gaussian integration point is
    created. Then this projector matrix is transferred to a projector on the nodal
    coordinates.
 
    Method corresponds to Elmer from git on October 27 2015
 
**************************************************************************************/
-(FEMMatrix * __nonnull)FEMMeshUtils_nodalProjectorMesh2:(FEMMesh * __nonnull)bMesh2 mesh1:(FEMMesh * __nonnull)bMesh1 useQuadrantTree:(BOOL)useQuadrantTree repeating:(BOOL)repeating antiRepeating:(BOOL)antiRepeating model:(FEMModel * __nonnull)model {
    
    FEMMatrix *projector;
    
    int *invPerm1 = bMesh1.getInvPerm;
    int *invPerm2 = bMesh2.getInvPerm;
    
    // Set the nodes of mesh1 to be the interval defined by mesh2
    BOOL *mirrorNode = NULL;
    if (repeating == YES) {
        if (antiRepeating == YES) {
            mirrorNode = (BOOL*)malloc(sizeof(BOOL) * bMesh1.numberOfNodes );
            memset(mirrorNode, 0, bMesh1.numberOfNodes*sizeof(BOOL) );
        }
        int size = bMesh1.numberOfNodes;
        [self preRotationalProjectorMesh1:bMesh1 mesh2:bMesh2 mirrorNode:mirrorNode sizeMirrorNode:&size];
    }
    
    // Create the projector usig the nodal points
    // This corresponds to numerical integration of the collocation method
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    projector = [utilities meshProjectorMesh1:bMesh2 mesh2:bMesh1 model:model useQuadrantTree:&useQuadrantTree transpose:NULL];
    projector.projectorType = PROJECTOR_TYPE_NODAL;
    
    matrixArraysContainer *projectorContainers = projector.getContainers;
    int *cols = projectorContainers->Cols;
    int *rows = projectorContainers->Rows;
    
    // One needs to change the sign of the projector for the mirror nodes
    if (antiRepeating == YES) {
        int size =  bMesh1.numberOfNodes;
        [self postRotationalProjector:projector mirrorNode:mirrorNode sizeMirrorNode:&size];
        if (mirrorNode != NULL) free(mirrorNode);
    }
    
    // Now return from the indexes of the interfaces mesh system to the
    // original mesh system
    projectorContainers->InvPerm = intvec(0, bMesh1.numberOfNodes-1);
    projectorContainers->sizeInvPerm = bMesh1.numberOfNodes;
    memcpy(projectorContainers->InvPerm, invPerm1, bMesh1.numberOfNodes*sizeof(int));
    
    int k;
    for (int i=0; i<projector.numberOfRows; i++) {
        for (int j=rows[i]; j<=rows[i+1]-1; j++) {
            k = cols[j];
            if (k >= 0) cols[j] = invPerm2[k];
        }
    }
    
    return projector;
}

/*************************************************************
 
    Set projector abs(rowsum) to unity
 
    Method corresponds to Elmer from git on October 27 2015

*************************************************************/
 -(void)FEMMeshUtils_setProjectorRowSum:(FEMMatrix * __nonnull)projector {
     
     double rowsum;
    
     matrixArraysContainer *projectorContainers = projector.getContainers;
     
     for (int i=0; i<projector.numberOfRows; i++) {
         rowsum = 0.0;
         for (int j=projectorContainers->Rows[i]; j<=projectorContainers->Rows[i+1]-1; j++) {
             rowsum = rowsum + fabs(projectorContainers->Values[j]);
         }
         for (int j=projectorContainers->Rows[i]; j<=projectorContainers->Rows[i+1]-1; j++) {
             projectorContainers->Values[j] = projectorContainers->Values[j] / rowsum;
         }
     }
}

#pragma mark Public methods

- (id)init
{
    self = [super init];
    if (self) {
        //Initialize here
    }
    
    return self;
}

/************************************************************************

    Find 2D mesh edges.
 
************************************************************************/
-(void)findEdges2DInMesh:(FEMMesh * __nonnull)mesh {
    
    int i, j, k, n=0, numbOfEdges;
    int node1, node2, edge=-1, swap, degree;
    HashTable_t *hashTable = NULL;
    HashEntry_t *hashPtr = NULL, *hashPtr1 = NULL;
    Element_t *elements = NULL, *element = NULL, *edges = NULL;
    ElementType_t *elmType = NULL;
    FEMElementDescription *elementDescription;
    BOOL found;
    
    elementDescription = [FEMElementDescription sharedElementDescription];
    
    elements = mesh.getElements;
    edges = mesh.getEdges;
    edges = (Element_t*)malloc(sizeof(Element_t) * (4*mesh.numberOfBulkElements) );
    
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        
        switch (element->Type.ElementCode / 100) {
            case 3:
                n = 3;
                break;
            case 4:
                n = 4;
                break;
        }
        
        if (element->EdgeIndexes == NULL) element->EdgeIndexes = intvec(0, n-1);
        memset( element->EdgeIndexes, -1, n*sizeof(int) );
    }
    
    hashTable = (HashTable_t*)malloc(sizeof(HashTable_t) * (mesh.numberOfNodes) );
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashTable[i].head = NULL;
    }
    
    // Loop over elements
    numbOfEdges = 0;
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        
        switch (element->Type.ElementCode / 100) {
            case 3:
                n = 3;
                break;
            case 4:
                n = 4;
                break;
        }
        
        // Loop over edge of every element
        for (k=0; k<n; k++) {
            // We use min(node1, node2) as the hash table key
            node1 = element->NodeIndexes[k];
            if (k < n-1) {
                node2 = element->NodeIndexes[k+1];
            } else {
                node2 = element->NodeIndexes[0];
            }
            
            if (node2 < node1) {
                swap = node1;
                node1 = node2;
                node2 = swap;
            }
            
            // Look for the edge for the hash table
            hashPtr = hashTable[node1].head;
            found = NO;
            while (hashPtr != NULL) {
                if (hashPtr->node1 == node2) {
                    found = YES;
                    edge = hashPtr->edge;
                    break;
                }
                hashPtr = hashPtr->next;
            }
            
            // Existing edge, update structures
            if (found == YES) {
                element->EdgeIndexes[k] = edge;
                edges[edge].BoundaryInfo->Right = element;
            } else {
                
                //Edge not there yet, create it
                edge = numbOfEdges;
                degree = element->Type.BasisFunctionDegree;
                
                edges[edge].ElementIndex = edge+1;
                edges[edge].NodeIndexes = intvec(0, (degree+1)-1);
                edges[edge].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(edges[edge].BoundaryInfo);
                BOOL noStab = NO;
                elmType = [elementDescription getElementType:(201+degree) inMesh:mesh stabilization:&noStab];
                edges[edge].Type = *elmType;
                
                edges[edge].NodeIndexes[0] = element->NodeIndexes[k];
                if (k < n-1) {
                    edges[edge].NodeIndexes[1] = element->NodeIndexes[k+1];
                } else {
                    edges[edge].NodeIndexes[1] = element->NodeIndexes[0];
                }
                
                for (j=1; j<degree; j++) {
                    edges[edge].NodeIndexes[j+1] = element->NodeIndexes[k+n+j-1];
                }
                
                // Create P element definitions if necessary
                if (element->Pdefs != NULL) {
                    edges[edge].Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
                    edges[edge].Pdefs->p = 0;
                } else {
                    edges[edge].Pdefs = NULL;
                }
                
                edges[edge].NDOFs = 0;
                if (element->NDOFs != 0) edges[edge].NDOFs = edges[edge].Type.NumberOfNodes;
                edges[edge].BDOFs = 0;
                edges[edge].DGDOFs = 0;
                edges[edge].EdgeIndexes = NULL;
                edges[edge].FaceIndexes = NULL;
                
                element->EdgeIndexes[k] = edge;
                
                edges[edge].BoundaryInfo->Left = element;
                edges[edge].BoundaryInfo->Right = NULL;
                
                // Update the hash table
                hashPtr = (HashEntry_t *)malloc(sizeof(HashEntry_t));
                hashPtr->edge = edge;
                hashPtr->node1 = node2;
                hashPtr->next = hashTable[node1].head;
                hashTable[node1].head = hashPtr;
                
                numbOfEdges++;
            }
        }
    }
    
    mesh.numberOfEdges = numbOfEdges;
    
    // Delete the hash table
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashPtr = hashTable[i].head;
        while (hashPtr != NULL) {
            hashPtr1 = hashPtr->next;
            free(hashPtr);
            hashPtr = hashPtr1;
        }
    }
    free(hashTable);
}

/************************************************************************
 
 Find 3D mesh edges.
 
************************************************************************/
-(void)findEdges3DInMesh:(FEMMesh * __nonnull)mesh {
    
    int i, ii, jj, k, n=0, numbOfEdges;
    int n1, n2;
    int node1, node2, edge=-1, degree;
    int **edgeMap = NULL, **faceEdgeMap = NULL;
    int **tetraEdgeMap, **tetraFaceMap, **brickEdgeMap, **wedgeEdgeMap, **pyramidEdgeMap;
    int **tetraFaceEdgeMap, **brickFaceEdgeMap, **wedgeFaceEdgeMap, **pyramidFaceEdgeMap;
    HashTable_t *hashTable = NULL;
    HashEntry_t *hashPtr = NULL, *hashPtr1 = NULL;
    Element_t *elements = NULL, *element = NULL, *edges = NULL, *faces = NULL;
    ElementType_t *elmType = NULL;
    FEMPElementMaps *elementMaps;
    FEMElementDescription *elementDescription;
    BOOL found;

    tetraEdgeMap = intmatrix(0, 5, 0, 2);
    tetraFaceMap = intmatrix(0, 3, 0, 5);
    brickEdgeMap = intmatrix(0, 5, 0, 8);
    wedgeEdgeMap = intmatrix(0, 4, 0, 3);
    pyramidEdgeMap = intmatrix(0, 4, 0, 7);
    
    tetraFaceEdgeMap = intmatrix(0, 3, 0, 2);
    brickFaceEdgeMap = intmatrix(0, 7, 0, 3);
    wedgeFaceEdgeMap = intmatrix(0, 5, 0, 3);
    pyramidFaceEdgeMap = intmatrix(0, 4, 0, 3);
    
    tetraFaceMap[0][0] = 0; tetraFaceMap[0][1] = 1; tetraFaceMap[0][2] = 2; tetraFaceMap[0][3] = 4; tetraFaceMap[0][4] = 5; tetraFaceMap[0][5] = 6;
    tetraFaceMap[1][0] = 0; tetraFaceMap[1][1] = 1; tetraFaceMap[1][2] = 3; tetraFaceMap[1][3] = 4; tetraFaceMap[1][4] = 8; tetraFaceMap[1][5] = 7;
    tetraFaceMap[2][0] = 1; tetraFaceMap[2][1] = 2; tetraFaceMap[2][2] = 3; tetraFaceMap[2][3] = 5; tetraFaceMap[2][4] = 9; tetraFaceMap[2][5] = 8;
    tetraFaceMap[3][0] = 2; tetraFaceMap[3][1] = 0; tetraFaceMap[3][2] = 3; tetraFaceMap[3][3] = 6; tetraFaceMap[3][4] = 7; tetraFaceMap[3][5] = 9;
    
    tetraFaceEdgeMap[0][0] = 0; tetraFaceEdgeMap[0][1] = 1; tetraFaceEdgeMap[0][2] = 2;
    tetraFaceEdgeMap[1][0] = 0; tetraFaceEdgeMap[1][1] = 4; tetraFaceEdgeMap[1][2] = 3;
    tetraFaceEdgeMap[2][0] = 1; tetraFaceEdgeMap[2][1] = 5; tetraFaceEdgeMap[2][2] = 4;
    tetraFaceEdgeMap[3][0] = 2; tetraFaceEdgeMap[3][1] = 3; tetraFaceEdgeMap[3][2] = 5;
    
    tetraEdgeMap[0][0] = 0; tetraEdgeMap[0][1] = 1; tetraEdgeMap[0][2] = 4;
    tetraEdgeMap[1][0] = 1; tetraEdgeMap[1][1] = 2; tetraEdgeMap[1][2] = 5;
    tetraEdgeMap[2][0] = 2; tetraEdgeMap[2][1] = 0; tetraEdgeMap[2][2] = 6;
    tetraEdgeMap[3][0] = 0; tetraEdgeMap[3][1] = 3; tetraEdgeMap[3][2] = 7;
    tetraEdgeMap[4][0] = 1; tetraEdgeMap[4][1] = 3; tetraEdgeMap[4][2] = 8;
    tetraEdgeMap[5][0] = 2; tetraEdgeMap[5][1] = 3; tetraEdgeMap[5][2] = 9;
    
    pyramidEdgeMap[0][0] = 0; pyramidEdgeMap[0][1] = 1; pyramidEdgeMap[0][2] = 0;
    pyramidEdgeMap[1][0] = 1; pyramidEdgeMap[1][1] = 2; pyramidEdgeMap[1][2] = 0;
    pyramidEdgeMap[2][0] = 2; pyramidEdgeMap[2][1] = 3; pyramidEdgeMap[2][2] = 0;
    pyramidEdgeMap[3][0] = 3; pyramidEdgeMap[3][1] = 0; pyramidEdgeMap[3][2] = 0;
    pyramidEdgeMap[4][0] = 0; pyramidEdgeMap[4][1] = 4; pyramidEdgeMap[4][2] = 0;
    pyramidEdgeMap[5][0] = 1; pyramidEdgeMap[5][1] = 4; pyramidEdgeMap[5][2] = 0;
    pyramidEdgeMap[6][0] = 2; pyramidEdgeMap[6][1] = 4; pyramidEdgeMap[6][2] = 0;
    pyramidEdgeMap[7][0] = 3; pyramidEdgeMap[7][1] = 4; pyramidEdgeMap[7][2] = 0;
    
    pyramidFaceEdgeMap[0][0] = 0; pyramidFaceEdgeMap[0][1] = 1; pyramidFaceEdgeMap[0][2] = 2; pyramidFaceEdgeMap[0][3] = 3;
    pyramidFaceEdgeMap[1][0] = 0; pyramidFaceEdgeMap[1][1] = 5; pyramidFaceEdgeMap[1][2] = 4; pyramidFaceEdgeMap[1][3] = -1;
    pyramidFaceEdgeMap[2][0] = 1; pyramidFaceEdgeMap[2][1] = 6; pyramidFaceEdgeMap[2][2] = 5; pyramidFaceEdgeMap[2][3] = -1;
    pyramidFaceEdgeMap[3][0] = 2; pyramidFaceEdgeMap[3][1] = 7; pyramidFaceEdgeMap[3][2] = 6; pyramidFaceEdgeMap[3][3] = -1;
    pyramidFaceEdgeMap[4][0] = 3; pyramidFaceEdgeMap[4][1] = 4; pyramidFaceEdgeMap[4][2] = 7; pyramidFaceEdgeMap[4][3] = -1;
    
    wedgeEdgeMap[0][0] = 0; wedgeEdgeMap[0][1] = 1; wedgeEdgeMap[0][2] = 0;
    wedgeEdgeMap[1][0] = 1; wedgeEdgeMap[1][1] = 2; wedgeEdgeMap[1][2] = 0;
    wedgeEdgeMap[2][0] = 0; wedgeEdgeMap[2][1] = 2; wedgeEdgeMap[2][2] = 0;
    wedgeEdgeMap[3][0] = 3; wedgeEdgeMap[3][1] = 4; wedgeEdgeMap[3][2] = 0;
    wedgeEdgeMap[4][0] = 4; wedgeEdgeMap[4][1] = 5; wedgeEdgeMap[4][2] = 0;
    wedgeEdgeMap[5][0] = 5; wedgeEdgeMap[5][1] = 3; wedgeEdgeMap[5][2] = 0;
    wedgeEdgeMap[6][0] = 0; wedgeEdgeMap[6][1] = 3; wedgeEdgeMap[6][2] = 0;
    wedgeEdgeMap[7][0] = 1; wedgeEdgeMap[7][1] = 4; wedgeEdgeMap[7][2] = 0;
    wedgeEdgeMap[8][0] = 2; wedgeEdgeMap[8][1] = 5; wedgeEdgeMap[8][2] = 0;
    
    wedgeFaceEdgeMap[0][0] = 0; wedgeFaceEdgeMap[0][1] = 1; wedgeFaceEdgeMap[0][2] = 2; wedgeFaceEdgeMap[0][3] = -1;
    wedgeFaceEdgeMap[1][0] = 3; wedgeFaceEdgeMap[1][1] = 4; wedgeFaceEdgeMap[1][2] = 5; wedgeFaceEdgeMap[1][3] = -1;
    wedgeFaceEdgeMap[2][0] = 0; wedgeFaceEdgeMap[2][1] = 7; wedgeFaceEdgeMap[2][2] = 3; wedgeFaceEdgeMap[2][3] = 6;
    wedgeFaceEdgeMap[3][0] = 1; wedgeFaceEdgeMap[3][1] = 8; wedgeFaceEdgeMap[3][2] = 4; wedgeFaceEdgeMap[3][3] = 7;
    wedgeFaceEdgeMap[4][0] = 2; wedgeFaceEdgeMap[4][1] = 6; wedgeFaceEdgeMap[4][2] = 5; wedgeFaceEdgeMap[4][3] = 8;
    
    brickEdgeMap[0][0] = 0; brickEdgeMap[0][1] = 1; brickEdgeMap[0][2] = 8;
    brickEdgeMap[1][0] = 1; brickEdgeMap[1][1] = 2; brickEdgeMap[1][2] = 9;
    brickEdgeMap[2][0] = 3; brickEdgeMap[2][1] = 2; brickEdgeMap[2][2] = 10;
    brickEdgeMap[3][0] = 0; brickEdgeMap[3][1] = 3; brickEdgeMap[3][2] = 11;
    brickEdgeMap[4][0] = 4; brickEdgeMap[4][1] = 5; brickEdgeMap[4][2] = 12;
    brickEdgeMap[5][0] = 5; brickEdgeMap[5][1] = 6; brickEdgeMap[5][2] = 13;
    brickEdgeMap[6][0] = 7; brickEdgeMap[6][1] = 6; brickEdgeMap[6][2] = 14;
    brickEdgeMap[7][0] = 4; brickEdgeMap[7][1] = 7; brickEdgeMap[7][2] = 15;
    brickEdgeMap[8][0] = 0; brickEdgeMap[8][1] = 4; brickEdgeMap[8][2] = 16;
    brickEdgeMap[9][0] = 1; brickEdgeMap[9][1] = 5; brickEdgeMap[9][2] = 17;
    brickEdgeMap[10][0] = 2; brickEdgeMap[10][1] = 6; brickEdgeMap[10][2] = 18;
    brickEdgeMap[11][0] = 3; brickEdgeMap[11][1] = 7; brickEdgeMap[11][2] = 19;
    
    brickFaceEdgeMap[0][0] = 0; brickFaceEdgeMap[0][1] = 1; brickFaceEdgeMap[0][2] = 2; brickFaceEdgeMap[0][3] = 3;
    brickFaceEdgeMap[1][0] = 4; brickFaceEdgeMap[1][1] = 5; brickFaceEdgeMap[1][2] = 6; brickFaceEdgeMap[1][3] = 7;
    brickFaceEdgeMap[2][0] = 0; brickFaceEdgeMap[2][1] = 9; brickFaceEdgeMap[2][2] = 4; brickFaceEdgeMap[2][3] = 8;
    brickFaceEdgeMap[3][0] = 1; brickFaceEdgeMap[3][1] = 10; brickFaceEdgeMap[3][2] = 5; brickFaceEdgeMap[3][3] = 9;
    brickFaceEdgeMap[4][0] = 2; brickFaceEdgeMap[4][1] = 11; brickFaceEdgeMap[4][2] = 6; brickFaceEdgeMap[4][3] = 10;
    brickFaceEdgeMap[5][0] = 3; brickFaceEdgeMap[5][1] = 8; brickFaceEdgeMap[5][2] = 7; brickFaceEdgeMap[5][3] = 11;
    
    elementDescription = [FEMElementDescription sharedElementDescription];
    elementMaps = [[FEMPElementMaps alloc] init];
    
    elements = mesh.getElements;
    faces = mesh.getFaces;
    edges = (Element_t*)malloc(sizeof(Element_t) * (12*mesh.numberOfBulkElements) );
    
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        if (element->EdgeIndexes == NULL) element->EdgeIndexes = intvec(0, 11);
        memset( element->EdgeIndexes, -1, 12*sizeof(int) );
    }
    
    hashTable = (HashTable_t*)malloc(sizeof(HashTable_t) * (mesh.numberOfNodes) );
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashTable[i].head = NULL;
    }

    // Loop over elements
    numbOfEdges = 0;
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        
        // For p elements, mappings are different
        if (element->Pdefs != NULL) {
            [elementMaps getEdgeMapForElement:element edgeMap:edgeMap];
            [elementMaps getFaceEdgeMapForElement:element faceEdgeMap:faceEdgeMap];
            n = element->Type.NumberOfEdges;
        } else {
            switch (element->Type.ElementCode / 100) {
                case 5:
                    n = 6;
                    edgeMap = tetraEdgeMap;
                    faceEdgeMap = tetraFaceEdgeMap;
                    break;
                case 6:
                    n = 8;
                    edgeMap = pyramidEdgeMap;
                    faceEdgeMap = pyramidFaceEdgeMap;
                    break;
                case 7:
                    n = 9;
                    edgeMap = wedgeEdgeMap;
                    faceEdgeMap = wedgeFaceEdgeMap;
                    break;
                case 8:
                    n = 12;
                    edgeMap = brickEdgeMap;
                    faceEdgeMap = brickFaceEdgeMap;
                    break;
                default:
                    fprintf(stderr, "FEMMeshUtils:findEdges3DInMesh: element type %d not implemented.\n", element->Type.ElementCode);
                    fatal("FEMMeshUtils:findEdges3DInMesh");
                    break;
            }
        }
        
        // Loop over every edge of every element
        for (k=0; k<n; k++) {
         
            // Use min(node1,node2) as key to hash table
            n1 = element->NodeIndexes[edgeMap[k][0]];
            n2 = element->NodeIndexes[edgeMap[k][1]];
            if (n1 < n2) {
                node1 = n1;
                node2 = n2;
            } else {
                node1 = n2;
                node2 = n1;
            }
            
            // Look the edge from the hash table
            hashPtr = hashTable[node1].head;
            found = NO;
            while (hashPtr != NULL) {
                if (hashPtr->node1 == node2) {
                    found = YES;
                    edge = hashPtr->edge;
                    break;
                }
                hashPtr = hashPtr->next;
            }
            
            // Existing edge, update structures
            if (found) {
                element->EdgeIndexes[k] = edge;
                
                // Mark as an edge of pyramid square face
                if ([elementMaps isPPyramid:element] == YES && k < 4) {
                    edges[edge].Pdefs->PyramidQuadEdge = true;
                }
                
                if (faces != NULL) {
                    for (ii=0; ii<element->Type.NumberOfFaces; ii++) {
                        if (faces[element->FaceIndexes[ii]].EdgeIndexes == NULL) {
                            faces[element->FaceIndexes[ii]].EdgeIndexes = intvec(0, faces[element->FaceIndexes[ii]].Type.NumberOfEdges-1);
                            memset( faces[element->FaceIndexes[ii]].EdgeIndexes, -1, faces[element->FaceIndexes[ii]].Type.NumberOfEdges*sizeof(int) );
                        }
                        for (jj=0; jj<faces[element->FaceIndexes[ii]].Type.NumberOfEdges; jj++) {
                            if (faceEdgeMap[ii][jj] == k) {
                                faces[element->FaceIndexes[ii]].EdgeIndexes[jj] = edge;
                                if (edges[edge].BoundaryInfo->Left == NULL) {
                                    edges[edge].BoundaryInfo->Left = &faces[element->FaceIndexes[ii]];
                                } else {
                                    edges[edge].BoundaryInfo->Right = &faces[element->FaceIndexes[ii]];
                                }
                            }
                        }
                    }
                }
            } else {
                
                //Edge not yet there, create it
                edge = numbOfEdges;
                edges[edge].ElementIndex = edge+1;
                degree = element->Type.BasisFunctionDegree;
                
                // Edge is always a line segment with deg+1 nodes
                BOOL noStab = NO;
                elmType = [elementDescription getElementType:(201+degree) inMesh:mesh stabilization:&noStab];
                edges[edge].Type = *elmType;
                
                edges[edge].NDOFs = 0;
                if (element->NDOFs != 0) edges[edge].NDOFs = edges[edge].Type.NumberOfNodes;
                edges[edge].BDOFs = 0;
                edges[edge].DGDOFs = 0;
                edges[edge].EdgeIndexes = NULL;
                edges[edge].FaceIndexes = NULL;
                
                edges[edge].NodeIndexes = intvec(0, (degree+1)-1);
                for (n2=0; n2<degree+1; n2++) {
                    edges[edge].NodeIndexes[n2] = element->NodeIndexes[edgeMap[k][n2]];
                }
                
                element->EdgeIndexes[k] = edge;
                edges[edge].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(edges[edge].BoundaryInfo);
                edges[edge].BoundaryInfo->Left = NULL;
                edges[edge].BoundaryInfo->Right = NULL;
                
                // Allocate p element definitions
                if (element->Pdefs != NULL) {
                    edges[edge].Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
                    edges[edge].Pdefs->p = 0;
                    edges[edge].Pdefs->PyramidQuadEdge = false;
                    // Here mark edge as edge of pyramid if needed (or set as not)
                    if ([elementMaps isPPyramid:element] == YES && k < 4) {
                        edges[edge].Pdefs->PyramidQuadEdge = true;
                    }
                } else {
                    edges[edge].Pdefs = NULL;
                }
                
                if (faces != NULL) {
                    for (ii=0; ii<element->Type.NumberOfFaces; ii++) {
                        if (faces[element->FaceIndexes[ii]].EdgeIndexes == NULL) {
                            faces[element->FaceIndexes[ii]].EdgeIndexes = intvec(0, faces[element->FaceIndexes[ii]].Type.NumberOfEdges-1);
                            memset(faces[element->FaceIndexes[ii]].EdgeIndexes, -1, faces[element->FaceIndexes[ii]].Type.NumberOfEdges*sizeof(int) );
                        }
                        for (jj=0; jj<faces[element->FaceIndexes[ii]].Type.NumberOfEdges; jj++) {
                            if (faceEdgeMap[ii][jj] == k) {
                                faces[element->FaceIndexes[ii]].EdgeIndexes[jj] = edge;
                                if (edges[edge].BoundaryInfo->Left == NULL) {
                                    edges[edge].BoundaryInfo->Left = &faces[element->FaceIndexes[ii]];
                                } else {
                                    edges[edge].BoundaryInfo->Right = &faces[element->FaceIndexes[ii]];
                                }
                            }
                        }
                    }
                }
                
                // Update the hash table
                hashPtr = (HashEntry_t *)malloc(sizeof(HashEntry_t));
                hashPtr->edge = edge;
                hashPtr->node1 = node2;
                hashPtr->next = hashTable[node1].head;
                hashTable[node1].head = hashPtr;

                
                numbOfEdges++;
            }
        }
    }
    
    [mesh assignEdges:edges];
    mesh.numberOfEdges = numbOfEdges;
    
    // Delete the hash table
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashPtr = hashTable[i].head;
        while (hashPtr != NULL) {
            hashPtr1 = hashPtr->next;
            free(hashPtr);
            hashPtr = hashPtr1;
        }
    }
    free(hashTable);
    
    if (faces != NULL) [self FEMMeshUtils_fixFaceEdges:mesh];
    
    free_imatrix(tetraEdgeMap, 0, 5, 0, 2);
    free_imatrix(tetraFaceMap, 0, 3, 0, 5);
    free_imatrix(brickEdgeMap, 0, 5, 0, 8);
    free_imatrix(wedgeEdgeMap, 0, 4, 0, 3);
    free_imatrix(pyramidEdgeMap, 0, 4, 0, 7);
    
    free_imatrix(tetraFaceEdgeMap, 0, 3, 0, 2);
    free_imatrix(brickFaceEdgeMap, 0, 7, 0, 3);
    free_imatrix(wedgeFaceEdgeMap, 0, 5, 0, 3);
    free_imatrix(pyramidFaceEdgeMap, 0, 4, 0, 3);

    [elementMaps deallocation];
}

/************************************************************************
 
 Find 3D mesh faces.
 
 ************************************************************************/
-(void)findFaces3DInMesh:(FEMMesh * __nonnull)mesh {
    
    int i, j, k, n, numbOfFaces;
    int n1 = 0, n2;
    int node1, node2, node3, face=-1, degree;
    int **faceMap = NULL;
    int **tetraFaceMap, **brickFaceMap, **wedgeFaceNap, **pyramidFaceMap;
    int *nf;
    HashTable_t *hashTable = NULL;
    HashEntry_t *hashPtr = NULL, *hashPtr1 = NULL;
    Element_t *elements = NULL, *element = NULL, *faces = NULL;
    ElementType_t *elmType = NULL;
    FEMPElementMaps *elementMaps;
    FEMElementDescription *elementDescription;
    BOOL found;
    
    elementMaps = [[FEMPElementMaps alloc] init];
    elementDescription = [FEMElementDescription sharedElementDescription];
    
    tetraFaceMap = intmatrix(0, 3, 0, 5);
    brickFaceMap = intmatrix(0, 5, 0, 8);
    wedgeFaceNap = intmatrix(0, 4, 0, 3);
    pyramidFaceMap = intmatrix(0, 4, 0, 7);
    nf = intvec(0, 3);
    
    tetraFaceMap[0][0] = 0; tetraFaceMap[0][1] = 1; tetraFaceMap[0][2] = 2; tetraFaceMap[0][3] = 4; tetraFaceMap[0][4] = 5; tetraFaceMap[0][5] = 6;
    tetraFaceMap[1][0] = 0; tetraFaceMap[1][1] = 1; tetraFaceMap[1][2] = 3; tetraFaceMap[1][3] = 4; tetraFaceMap[1][4] = 8; tetraFaceMap[1][5] = 7;
    tetraFaceMap[2][0] = 1; tetraFaceMap[2][1] = 2; tetraFaceMap[2][2] = 3; tetraFaceMap[2][3] = 5; tetraFaceMap[2][4] = 9; tetraFaceMap[2][5] = 8;
    tetraFaceMap[3][0] = 2; tetraFaceMap[3][1] = 0; tetraFaceMap[3][2] = 3; tetraFaceMap[3][3] = 6; tetraFaceMap[3][4] = 7; tetraFaceMap[3][5] = 9;
    
    wedgeFaceNap[0][0] = 0; wedgeFaceNap[0][1] = 1; wedgeFaceNap[0][2] = 2; wedgeFaceNap[0][3] = -1;
    wedgeFaceNap[1][0] = 3; wedgeFaceNap[1][1] = 4; wedgeFaceNap[1][2] = 5; wedgeFaceNap[1][3] = -1;
    wedgeFaceNap[2][0] = 0; wedgeFaceNap[2][1] = 1; wedgeFaceNap[2][2] = 4; wedgeFaceNap[2][3] = 3;
    wedgeFaceNap[3][0] = 2; wedgeFaceNap[3][1] = 1; wedgeFaceNap[3][2] = 4; wedgeFaceNap[3][3] = 5;
    wedgeFaceNap[4][0] = 2; wedgeFaceNap[4][1] = 0; wedgeFaceNap[4][2] = 3; wedgeFaceNap[4][3] = 5;
    
    pyramidFaceMap[0][0] = 0; pyramidFaceMap[0][1] = 1; pyramidFaceMap[0][2] = 2; pyramidFaceMap[0][3] = 3; pyramidFaceMap[0][4] = 5;
    pyramidFaceMap[0][5] = 6; pyramidFaceMap[0][6] = 7; pyramidFaceMap[0][7] = 8;
    pyramidFaceMap[1][0] = 0; pyramidFaceMap[1][1] = 1; pyramidFaceMap[1][2] = 4; pyramidFaceMap[1][3] = 5; pyramidFaceMap[1][4] = 10;
    pyramidFaceMap[1][5] = 9; pyramidFaceMap[1][6] = -1; pyramidFaceMap[1][7] = -1;
    pyramidFaceMap[2][0] = 1; pyramidFaceMap[2][1] = 2; pyramidFaceMap[2][2] = 4; pyramidFaceMap[2][3] = 6; pyramidFaceMap[2][4] = 11;
    pyramidFaceMap[2][5] = 10; pyramidFaceMap[2][6] = -1; pyramidFaceMap[2][7] = -1;
    pyramidFaceMap[3][0] = 2; pyramidFaceMap[3][1] = 3; pyramidFaceMap[3][2] = 4; pyramidFaceMap[3][3] = 7; pyramidFaceMap[3][4] = 12;
    pyramidFaceMap[3][5] = 11; pyramidFaceMap[3][6] = -1; pyramidFaceMap[3][7] = -1;
    pyramidFaceMap[4][0] = 3; pyramidFaceMap[4][1] = 0; pyramidFaceMap[4][2] = 4; pyramidFaceMap[4][3] = 8; pyramidFaceMap[4][4] = 9;
    pyramidFaceMap[4][5] = 12; pyramidFaceMap[4][6] = -1; pyramidFaceMap[4][7] = -1;
    
    brickFaceMap[0][0] = 0; brickFaceMap[0][1] = 1; brickFaceMap[0][2] = 2; brickFaceMap[0][3] = 3; brickFaceMap[0][4] = 8; brickFaceMap[0][5] = 9;
    brickFaceMap[0][6] = 10; brickFaceMap[0][7] = 11; brickFaceMap[0][8] = 24;
    brickFaceMap[1][0] = 4; brickFaceMap[1][1] = 5; brickFaceMap[1][2] = 6; brickFaceMap[1][3] = 7; brickFaceMap[1][4] = 16; brickFaceMap[1][5] = 17;
    brickFaceMap[1][6] = 18; brickFaceMap[1][7] = 19; brickFaceMap[1][8] = 25;
    brickFaceMap[2][0] = 0; brickFaceMap[2][1] = 1; brickFaceMap[2][2] = 5; brickFaceMap[2][3] = 4; brickFaceMap[2][4] = 8; brickFaceMap[2][5] = 13;
    brickFaceMap[2][6] = 16; brickFaceMap[2][7] = 12; brickFaceMap[2][8] = 20;
    brickFaceMap[3][0] = 1; brickFaceMap[3][1] = 2; brickFaceMap[3][2] = 6; brickFaceMap[3][3] = 5; brickFaceMap[3][4] = 9; brickFaceMap[3][5] = 14;
    brickFaceMap[3][6] = 16; brickFaceMap[3][7] = 13; brickFaceMap[3][8] = 21;
    brickFaceMap[4][0] = 2; brickFaceMap[4][1] = 3; brickFaceMap[4][2] = 7; brickFaceMap[4][3] = 6; brickFaceMap[4][4] = 10; brickFaceMap[4][5] = 15;
    brickFaceMap[4][6] = 18; brickFaceMap[4][7] = 14; brickFaceMap[4][8] = 22;
    brickFaceMap[5][0] = 3; brickFaceMap[5][1] = 0; brickFaceMap[5][2] = 4; brickFaceMap[5][3] = 7; brickFaceMap[5][4] = 11; brickFaceMap[5][5] = 12;
    brickFaceMap[5][6] = 19; brickFaceMap[5][7] = 15; brickFaceMap[5][8] = 23;
    
    elements = mesh.getElements;
    faces = (Element_t*)malloc(sizeof(Element_t) * (6*mesh.numberOfBulkElements) );
    
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        if (element->FaceIndexes == NULL) element->FaceIndexes = intvec(0, 5);
        memset( element->FaceIndexes, -1, 6*sizeof(int) );
    }
    
    hashTable = (HashTable_t*)malloc(sizeof(HashTable_t) * (mesh.numberOfNodes) );
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashTable[i].head = NULL;
    }
    
    // Loop over elements
    numbOfFaces = 0;
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        element = &elements[i];
        
        // For p elements, mappings are different
        if (element->Pdefs != NULL) {
            [elementMaps getFaceMapForElement:element faceMap:faceMap];
            n = element->Type.NumberOfFaces;
        } else {
            switch (element->Type.ElementCode / 100) {
                case 5:
                    n = 4;
                    faceMap = tetraFaceMap;
                    break;
                case 6:
                    n = 5;
                    faceMap = pyramidFaceMap;
                    break;
                case 7:
                    n = 5;
                    faceMap = wedgeFaceNap;
                    break;
                case 8:
                    n = 6;
                    faceMap = brickFaceMap;
                    break;
                default:
                    continue;
                    break;
            }
        }
        
        // Loop over every face of every element
        for (k=0; k<n; k++) {
            
            // We use min(node1,node2,node3) as the hash table key
            switch (element->Type.ElementCode / 100) {
                case 5:
                    
                    // Tetras:
                    // =======
                    for (j=0; j<3; j++) {
                        nf[j] = element->NodeIndexes[faceMap[k][j]];
                    }
                    vDSP_vsort((float *)nf, 3, 1);
                    break;
                    
                case 6:
                    
                    // Pyramids:
                    // ========
                    if (k == 0) {
                        for (j=0; j<4; j++) {
                            nf[j] = element->NodeIndexes[faceMap[k][j]];
                        }
                        vDSP_vsort((float *)nf, 4, 1);
                    } else {
                        for (j=0; j<3; j++) {
                            nf[j] = element->NodeIndexes[faceMap[k][j]];
                        }
                        vDSP_vsort((float *)nf, 3, 1);
                    }
                    break;
                    
                case 7:
                    
                    // Wedges:
                    // ======
                    if (k <= 1) {
                        for (j=0; j<3; j++) {
                            nf[j] = element->NodeIndexes[faceMap[k][j]];
                        }
                        vDSP_vsort((float *)nf, 3, 1);
                    } else {
                        for (j=0; j<4; j++) {
                            nf[j] = element->NodeIndexes[faceMap[k][j]];
                        }
                        vDSP_vsort((float *)nf, 4, 1);
                    }
                    break;
                    
                case 8:
                    
                    // Bricks:
                    // ======
                    for (j=0; j<4; j++) {
                        nf[j] = element->NodeIndexes[faceMap[k][j]];
                    }
                    vDSP_vsort((float *)nf, 4, 1);
                    break;
                    
                default:
                    fprintf(stderr, "FEMMeshUtils:findFaces3DForMesh: element type %d not implemented.\n", element->Type.ElementCode);
                    fatal("FEMMeshUtils:findFaces3DForMesh");
                    break;
            }
            
            node1 = nf[0];
            node2 = nf[1];
            node3 = nf[2];
            
            // Look the face from the hash table
            hashPtr = hashTable[node1].head;
            found = NO;
            while (hashPtr != NULL) {
                if (hashPtr->node1 == node2 && hashPtr->node2 == node3) {
                    found = YES;
                    face = hashPtr->face;
                    break;
                }
                hashPtr = hashPtr->next;
            }
            
            // Existing face, update structures
            if (found == YES) {
                element->FaceIndexes[k] = face;
                faces[face].BoundaryInfo->Right = element;
            } else {
                
                // Face not yet there, create it
                face = numbOfFaces;
                faces[face].ElementIndex = face+1;
                
                degree = element->Type.BasisFunctionDegree;
                switch (element->Type.ElementCode / 100) {
                    case 5:
                        // Tetras:
                        // ======
                        switch (degree) {
                            case 1:
                                n1 = 3;
                                break;
                            case 2:
                                n1 = 6;
                                break;
                            case 3:
                                n1 = 10;
                                break;
                        }
                        n1 = 3; // TODO: The switch statement above is then useless?
                        BOOL noStab = NO;
                        elmType = [elementDescription getElementType:(300+n1) inMesh:mesh stabilization:&noStab];
                        faces[face].Type = *elmType;
                        break;
                        
                    case 6:
                        // Pyramid (only 605 supported)
                        // ============================
                        if (k == 0) {
                            n1 = 4 ;
                            BOOL noStab = NO;
                            elmType = [elementDescription getElementType:(400+n1) inMesh:mesh stabilization:&noStab];
                            faces[face].Type = *elmType;
                        } else {
                            n1 = 3;
                            BOOL noStab = NO;
                            elmType = [elementDescription getElementType:(300+n1) inMesh:mesh stabilization:&noStab];
                            faces[face].Type = *elmType;;
                        }
                        break;
                        
                    case 7:
                        // Wedges (only 706 supported)
                        // ===========================
                        if (k <= 1) {
                            n1 = 3;
                            BOOL noStab = NO;
                            elmType = [elementDescription getElementType:303 inMesh:mesh stabilization:&noStab];
                            face[faces].Type = *elmType;
                        } else {
                            n1 = 4;
                            BOOL noStab = NO;
                            elmType = [elementDescription getElementType:404 inMesh:mesh stabilization:&noStab];
                            faces[face].Type = *elmType;
                        }
                        break;
                        
                    case 8:
                        // Bricks
                        // ======
                        switch (element->Type.NumberOfNodes) {
                            case 8:
                                n1 = 4;
                                break;
                            case 20:
                                n1 = 8;
                                break;
                            case 27:
                                n1 = 9;
                                break;
                        }
                        noStab = NO;
                        elmType = [elementDescription getElementType:(400+n1) inMesh:mesh stabilization:&noStab];
                        faces[face].Type = *elmType;
                        break;
                        
                    default:
                        fprintf(stderr, "FEMMeshUtils:findFaces3DForMesh: element type %d not implemented.\n", element->Type.ElementCode);
                        fatal("FEMMeshUtils:findFaces3DForMesh");
                        break;
                }
                
                // Create P element definitions if necessary
                if (element->Pdefs != NULL) {
                    faces[face].Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
                    faces[face].Pdefs->p = 0;
                } else {
                    faces[face].Pdefs = NULL;
                }
                
                faces[face].NDOFs = 0;
                if (element->NDOFs != 0) faces[face].NDOFs = faces[face].Type.NumberOfNodes;
                faces[face].BDOFs = 0;
                faces[face].DGDOFs = 0;
                faces[face].EdgeIndexes = NULL;
                faces[face].FaceIndexes = NULL;
                
                faces[face].NodeIndexes = intvec(0, n1-1);
                for (n2=0; n2<n1; n2++) {
                    faces[face].NodeIndexes[n2] = element->NodeIndexes[faceMap[k][n2]];
                }
                
                element->FaceIndexes[k] = face;
                
                faces[face].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(faces[face].BoundaryInfo);
                faces[face].BoundaryInfo->Left = element;
                faces[face].BoundaryInfo->Right = NULL;
                
                // Update the hash table
                hashPtr = (HashEntry_t *)malloc(sizeof(HashEntry_t));
                hashPtr->face = face;
                hashPtr->node1 = node2;
                hashPtr->node2 = node3;
                hashPtr->next = hashTable[node1].head;
                hashTable[node1].head = hashPtr;
                
                numbOfFaces++;
            }
        }
    }
    
    [mesh assignFaces:faces];
    mesh.numberOfFaces = numbOfFaces;
    
    // Delete the hash table
    for (i=0; i<mesh.numberOfNodes; i++) {
        hashPtr = hashTable[i].head;
        while (hashPtr != NULL) {
            hashPtr1 = hashPtr->next;
            free(hashPtr);
            hashPtr = hashPtr1;
        }
    }
    free(hashTable);
    
    free_imatrix(tetraFaceMap, 0, 3, 0, 5);
    free_imatrix(brickFaceMap, 0, 5, 0, 8);
    free_imatrix(wedgeFaceNap, 0, 4, 0, 3);
    free_imatrix(pyramidFaceMap, 0, 4, 0, 7);
    free_ivector(nf, 0, 3);
    
    [elementMaps deallocation];
}

/************************************************************************
 
    Generate element edge (face in 3D) tables for given mesh.
    Only for triangles and tetras. If mesh already has edges, do nothing

************************************************************************/
-(void)findEdgesForMesh:(FEMMesh * __nonnull)mesh findEdges:(BOOL * __nullable)present {
    
    Element_t *edges = NULL, *faces = NULL;
    BOOL findEdges3D;
    
    if (present != NULL) {
        findEdges3D = *present;
    } else {
        findEdges3D = YES;
    }
    
    edges = mesh.getEdges;
    faces = mesh.getFaces;
    
    switch (mesh.dimension) {
        case 2:
            if (edges == NULL) [self findEdges2DInMesh:mesh];
            break;
        case 3:
            if (faces == NULL) [self findFaces3DInMesh:mesh];
            if (findEdges3D == YES) {
                if (edges == NULL) [self findEdges3DInMesh:mesh];
            }
            break;
    }
    
    [self FEMMeshUtils_assignConstraints:mesh];
}

/*************************************************************************
    
    Assign local number of edge to given boundary element. Also copies all
    p element attributes from element edge to boundary edge

*************************************************************************/
-(void)assignLocalNumberToEdgeElement:(Element_t * __nonnull)edge fromElement:(Element_t * __nonnull)element inMesh:(FEMMesh * __nonnull)mesh {
    
    int i, j, n, edgeNumber, numbEdges, bMap[4];
    Element_t *entity = NULL;
    FEMPElementMaps *pMaps;
    
    pMaps = [[FEMPElementMaps alloc] init];
    
    numbEdges = 0;
    switch (element->Type.dimension) {
        case 2:
            numbEdges = element->Type.NumberOfEdges;
            break;
        case 3:
            numbEdges = element->Type.NumberOfFaces;
            break;
        default:
            fprintf(stdout, "FEMMeshUtils:assignLocalNumberToEdgeElement: unsupported dimension.");
            return;
            break;
    }
    
    // For each edge or face in element, try to find local number
    for (edgeNumber=0; edgeNumber<numbEdges; edgeNumber++) {
        // If edges have been created, stop search. Actually, this should not happen
        if (element->EdgeIndexes == NULL) {
            fprintf(stdout, "FEMMeshUtils:assignLocalNumberToEdgeElement: edges have not been creates. Returning now.");
            return;
        }
        
        entity = [self FEMMeshUtils_getEntityForElement:element edge:edgeNumber inMesh:mesh];
        
        // Edge element not found. This should not be possible, unless there is
        // an error in the mesh read in process
        if (entity == NULL) {
            fprintf(stdout, "FEMMeshUtils:assignLocalNumberToEdgeElement: edge element not found.");
            return;
        }
        
        n = 0;
        // For each element node
        for (i=0; i<entity->Type.NumberOfNodes; i++) {
            // For each node in edge element
            for (j=0; j<edge->Type.NumberOfNodes; j++) {
                // If entity and edge element node match, increment counter
                if (entity->NodeIndexes[i] == edge->NodeIndexes[j]) n++;
            }
        }
        
        // If all nodes are on boundary, edge was found
        if (n == edge->Type.NumberOfNodes) {
            edge->Pdefs->LocalNumber = edgeNumber;
            
            // Change ordering of global nodes to match that of element
            [pMaps getBoundaryMapForElement:element localNumber:edgeNumber resultMap:bMap];
            for (i=0; i<n; i++) {
                edge->NodeIndexes[i] = element->NodeIndexes[bMap[i]];
            }
            
            // Copy attributes of edge element to boundary element
            // Misc attributes
            edge->Pdefs->isEdge = entity->Pdefs->isEdge;
            
            // Gauss points
            edge->Pdefs->GaussPoints = entity->Pdefs->GaussPoints;
            
            // Element p (and boundary bubble dofs)
            edge->BDOFs = entity->BDOFs;
            edge->Pdefs->p = entity->Pdefs->p;
            
            // If this boundary has edges copy edge indexes
            if (entity->EdgeIndexes != NULL) {
                // Allocate element edges to element
                n = entity->Type.NumberOfEdges;
                [pMaps getFaceEdgeMapForElement:element index:edgeNumber resultMap:bMap];
             
                if (edge->EdgeIndexes != NULL) {
                    free_ivector(edge->EdgeIndexes, 0, edge->sizeEdgeIndexes-1);
                }
                edge->EdgeIndexes = intvec(0, n-1);
                edge->sizeEdgeIndexes = n;
                // Copy edges from edge to boundary edge
                for (i=0; i<n; i++) {
                    edge->EdgeIndexes[i] = element->EdgeIndexes[bMap[i]];
                }
            }
            
            // Edge fields copied and local edge found so return
            return;
        }
    }
    
    // If we are here, local number not found
    fprintf(stdout, "FEMMeshUtils:assignLocalNumberToEdgeElement: unable to find local edge.");
}

/************************************************************************************************************
 
    Create a projector between Master and Target boundaries. The projector may be a nodal projector x=Px or 
    a weigted Galerking projector such that Qx=Px. In the first case the projector will be P and in the 
    second case [Q-P].
 
    Method corresponds to Elmer from git on October 27 2015

************************************************************************************************************/
-(FEMMatrix * __nullable)periodicProjectorInModel:(FEMModel * __nonnull)model forMesh:(FEMMesh * __nonnull)mesh masterBoundary:(int)mbd targetBoundary:(int)trgt galerking:(BOOL * __nullable)galerkin {
    
    int dim;
    FEMMatrix *projector;
    FEMMesh *bMesh1, *bMesh2;
    BOOL antiradial, antirotational, antisliding, cylindrical, doEdges, doNodes, flat, found, intGalerkin, levelProj, plane, radial, rotational,
         sliding, success, useExtProjector;
    
    // BCs index should start from 0
    if (mbd < 0) return nil;
    
    FEMCore *core = [FEMCore sharedCore];
    FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
    FEMSolution *solution = (FEMSolution *)model.solution;
    
    dim = model.dimension;
    
    [listUtilities resetTimer:@"periodicProjector" model:model];
    
    FEMBoundaryCondition *boundaryConditionAtId = (model.boundaryConditions)[mbd];
    FEMMesh *pMesh = mesh;
    
    // Whether to choose nodal or Galerkin projector is determined by an optional
    // flag. The default is the nodal projector
    if (galerkin != NULL) {
        intGalerkin = *galerkin;
    } else {
        intGalerkin = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"galerkin projector" info:&found];
    }
    
    // If the boundary is discontinous, then we have the luxury of creating the projector
    // very cheeply using the permutation vector. This does not need the target as the
    // boundary is self-contained.
    if ([listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"discontinuous boundary" info:&found] == YES && mesh.isDiscontinuousMesh == YES) {
        if (intGalerkin == YES) {
            projector = [self FEMMeshUtils_weightedProjectorDiscontinousMesh:pMesh model:model boundary:mbd listUtilities:listUtilities];
        } else {
            projector = [self FEMMeshUtils_nodalProjectorDiscontinuousMesh:pMesh model:model boundary:mbd listUtilities:listUtilities];
        }
        
        if (projector == nil) return nil;
        goto jump;
    }
    
    // BCs index should start from 0
    if (trgt < 0) return nil;
    
    // Create the mesh projector, and if needed also eliminate the ghost nodes.
    // There are two choices of projector: a nodal P in x = PX, and a Galerkin projector
    // [Q-P] in Qx = Px. The projector is assumed to be either a rotational projector with no
    // translation and roration, or a generic one with possible cooridinate mapping
    fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: -----------------------------------------------------------\n");
    fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: creating projector between BCs %d and %d.\n", mbd+1, trgt+1);
    fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: -----------------------------------------------------------\n");
    
    // Create temporal mesh structures that are utilized when making the
    // projector between the two boundaries
    bMesh1 = [[FEMMesh alloc] init];
    bMesh2 = [[FEMMesh alloc] init];
    
    success = [self FEMMeshUtils_createInterfaceMeshesModel:model mesh:mesh masterBoundary:mbd targetBoundary:trgt bMesh1:bMesh1 bMesh2:bMesh2];
    if (success == NO) {
        [bMesh1 deallocation];
        [bMesh2 deallocation];
        return nil;
    }
    
    // Do we have external method to take care of the projection matrix creation
    useExtProjector = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"external projector" info:&found];
    
    // If requested, map the interface coordinate from (x,y,z) to any permutation of these
    [self FEMMeshUtils_mapInterfaceCoordinateMesh1:bMesh1 mesh2:bMesh2 BCParams:boundaryConditionAtId model:model listUtilities:listUtilities];
    
    Element_t *elements = mesh.getElements;
    Nodes_t *bMesh1Nodes = bMesh1.getNodes;
    Nodes_t *bMesh2Nodes = bMesh2.getNodes;
    
    // Check whether to use (anti) rotational projector.
    // We don't really know on which side the projector was called so
    // let's check both sides
    rotational = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"rotational projector" info:&found];
    antirotational = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"anti rotational projector" info:&found];
    if (antirotational == YES) rotational = YES;
    
    cylindrical = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"cylindrical projector" info:&found];
    
    radial = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"radial projector" info:&found];
    antiradial = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"anti radial projector" info:&found];
    if (antiradial == YES) radial = YES;
    
    sliding = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"sliding projector" info:&found];
    antisliding = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"anti sliding projector" info:&found];
    if (antisliding == YES) sliding = YES;
    
    flat = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"flat projector" info:&found];
    plane = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"plane projector" info:&found];
    
    if (radial == YES) fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: enforcig radial projector.\n");
    if (sliding == YES) fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: enforcing sliding projector.\n");
    if (cylindrical == YES) fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: enforcing cyclindrical projector.\n");
    if (rotational == YES) fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: enforcing rotational projector.\n");
    if (flat == YES) fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: enforcing flat projector.\n");
    if (plane == YES) fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: enforcing plane projector.\n");
    
    double nodeScale = [listUtilities listGetConstReal:model inArray:boundaryConditionAtId.valuesList forVariable:@"mortar bc scaling" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) {
        if (antiradial == YES) {
            nodeScale = -1.0;
        } else {
            nodeScale = 1.0;
        }
    }
    
    BOOL nodalJump = [listUtilities listCheckPrefix:@"mortar bc coefficient" inArray:boundaryConditionAtId.valuesList];
    if (nodalJump == NO) {
        nodalJump = [listUtilities listCheckPrefix:@"mortar bc resistivity" inArray:boundaryConditionAtId.valuesList];
    }
    
    // There are tailored projectors for simplified interfaces
    
    // Stride projector is obsolite and has been eliminated
    if ([listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"stride projector" info:&found] == YES) {
        [listUtilities addLogicalInClassList:boundaryConditionAtId theVariable:@"level projector" withValue:YES];
        [listUtilities addLogicalInClassList:boundaryConditionAtId theVariable:@"level projector strong" withValue:YES];
        fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: enforcing level projector instead of old stride projector.\n");
    }
    
    levelProj = [listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"level projector" info:&found];
    if (rotational == YES || cylindrical == YES || radial == YES || flat == YES || plane == YES ) {
        if (found == NO) {
            fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: enforcing level projector to YES with dimensional reduction.\n");
            levelProj = YES;
        } else if (levelProj == NO) {
            // If we have dimensionally reduced projector but don't use level projector to
            // integrate over it, then ensure that the third coordinate is set to zero
            memset(bMesh1Nodes->z, 0.0, bMesh1.numberOfNodes*sizeof(double) );
            memset(bMesh2Nodes->z, 0.0, bMesh2.numberOfNodes*sizeof(double) );
        }
    }
    
    if (levelProj == YES) {
        if ((solution.solutionInfo)[@"projector skip nodes"] != nil) {
            if ([(solution.solutionInfo)[@"projector skip nodes"] boolValue] == YES) doNodes = NO;
        } else {
            if ([listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"projector skip nodes" info:&found] == YES) {
                doNodes = NO;
            } else {
                doNodes = (mesh.numberOfNodes > 0) ? YES : NO;
            }
        }
        
        if ((solution.solutionInfo)[@"projector skip edges"] != nil) {
            if ([(solution.solutionInfo)[@"projector skip edges"] boolValue] == YES) doEdges = NO;
        } else {
            if ([listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"projector skip edges" info:&found] == YES) {
                doEdges = NO;
            } else {
                // We are conservative here since there may be edges in 2D which
                // still can not be used for creating the projector
                doEdges = (mesh.numberOfNodes > 0 && mesh.dimension == 3 && dim == 3) ? YES :NO;
                
                // Ensure that there is no p-elements that made us think that we have edges.
                // Here we assume that if there is any p-elements, then also the first element is such
                if (doEdges == YES) {
                    if ([core isPElement:&elements[0]] == YES) {
                        doEdges = NO;
                        fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: edge projector will not be created for p-element mesh.\n");
                    }
                }
            }
        }
    }
    
    // If the interface is rotational, move to (phi, z) plane and alter the phi coordinates
    // so that the meshes coincide.
    // Otherwise make the two meshes to coinside using rotation, translation and scaling
    double radius = 1.0;
    BOOL fullCircle = NO;
    if (rotational == YES || cylindrical == YES) {
        [self FEMMeshUtils_rotationInterfaceMesh1:bMesh1 mesh2:bMesh2 bcParams:boundaryConditionAtId cylindrical:cylindrical radius:radius fullCircle:&fullCircle model:model listUtilities:listUtilities];
    } else if (radial == YES) {
        [self FEMMeshUtils_radiusInterfaceMesh1:bMesh1 mesh2:bMesh2 bcParams:boundaryConditionAtId model:model listUtilities:listUtilities];
    } else if (flat == YES) {
        [self FEMMeshUtils_flatInterfaceMesh1:bMesh1 mesh2:bMesh2 bcParams:boundaryConditionAtId model:model listUtilities:listUtilities];
    } else if (plane == YES) {
        [self FEMMeshUtils_planeInterfaceMesh1:bMesh1 mesh2:bMesh2 bcParams:boundaryConditionAtId model:model listUtilities:listUtilities];
    } else if (sliding == NO) {
        [self FEMMeshUtils_overlayInterfaceMesh1:bMesh1 mesh2:bMesh2 bcParams:boundaryConditionAtId model:model listUtilities:listUtilities];
    }
    
    BOOL repeating = ((rotational == YES  && fullCircle == NO) || sliding == YES) ? YES : NO;
    BOOL antiRepeating = ((antirotational == YES && fullCircle == NO) || antisliding == YES) ? YES : NO;
    
    if (useExtProjector == YES) {
        // TODO: Implement this case
    } else if (levelProj == YES) {
        // TODO: Implement this case
    } else {
        if (fullCircle == YES) {
            fatal("FEMMeshUtils:periodicProjectorInModel", "A full circle can not be dealt with the generic projector.");
        }
        BOOL useQuandrantTree = [listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"use quadrant tree" info:&found];
        if (found == NO) useQuandrantTree = YES;
        if (intGalerkin == YES) {
            int *invPerm1 = bMesh1.getInvPerm;
            int *invPerm2 = bMesh2.getInvPerm;
            FEMInterpolateMeshToMesh *interpolationMeshToMesh = [[FEMInterpolateMeshToMesh alloc] init];
            projector = [interpolationMeshToMesh weightedProjectorMesh2:bMesh2 mesh1:bMesh1 inversePermutation2:invPerm2 sizeInversePermutation2:bMesh2.numberOfNodes inversePermutation1:invPerm1 sizeInversePermutation1:bMesh1.numberOfNodes useQuadrantTree:useQuandrantTree repeating:repeating antiRepeating:antiRepeating periodicScale:nodeScale nodalJump:nodalJump model:model];
        } else {
            projector = [self FEMMeshUtils_nodalProjectorMesh2:bMesh2 mesh1:bMesh1 useQuadrantTree:useQuandrantTree repeating:repeating antiRepeating:antiRepeating model:model];
        }
    }
    
    // Deallocate meshes
    [bMesh1 deallocation];
    [bMesh2 deallocation];
    
jump:
    projector.projectorBC = mbd;
    matrixArraysContainer *projectorContainers = projector.getContainers;
    
    if ([listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"projector set rowsum" info:&found] == YES) {
        [self FEMMeshUtils_setProjectorRowSum:projector];
    }
    
    double coeff = [listUtilities listGetConstReal:model inArray:boundaryConditionAtId.valuesList forVariable:@"projector multiplier" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"projector multiplier" info:&found minValue:NULL maxValue:NULL];
    if (found == YES) {
        for (int i=0; i<projectorContainers->sizeValues; i++) {
            projectorContainers->Values[i] = coeff * projectorContainers->Values[i];
        }
    }
    
    if ([listUtilities listGetLogical:model inArray:boundaryConditionAtId.valuesList forVariable:@"save projector" info:&found] == YES) {
        [self saveProjector:projector saveRowSum:YES prefix:[@"p" stringByAppendingString:[NSString stringWithFormat:@"%d",mbd+1]] invPerm:NULL];
        // Dual projector if it exists
        if (projector.ematrix != nil) {
            [self saveProjector:projector.ematrix saveRowSum:YES prefix:[@"pd" stringByAppendingString:[NSString stringWithFormat:@"%d",mbd+1]] invPerm:projectorContainers->InvPerm];
        }
        // Biorthogonal projector if it exists
        if (projector.child != nil) {
            [self saveProjector:projector.child saveRowSum:YES prefix:[@"pb" stringByAppendingString:[NSString stringWithFormat:@"%d",mbd+1]] invPerm:projectorContainers->InvPerm];
        }
        fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: save projector and stop.\n");
    }
    
    BOOL deleteTimer = YES;
    [listUtilities checkTimer:@"periodicProjector" deleteTimer:&deleteTimer resetTimer:NULL model:model];
    fprintf(stdout, "FEMMeshUtils:periodicProjectorInModel: projector created, now exiting...\n");
    
    return projector;
}

/********************************************************************************************
 
    Split a mesh equally to smaller pieces by performing a uniform split. Also known as
    mesh multiplication. A 2D element splits into 4 elemenrs of same form and a #d element
    into 8 elements.
    Currently works only for linear elements
 
********************************************************************************************/
-(FEMMesh * __nullable)splitMeshEqual:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model nodal:(double * __nullable)h sizeNodal:(int * __nullable)sizeNodal {
    
    int i, j=0, k, l, n, newElCnt, nodeCnt, edgeCnt, node, parentID, diag, minVal, maxVal;
    int faceNumber, edge1, edge2, edge3, edge4, node12, node23, node34, node41, node31;
    int n1, n2, n3=0, *eoldNodes, *faceNodes, *edgeNodes;        // Only linear so far
    int **child;
    double *nodals, *xh = NULL;
    double dxyz[3][3], dist[3], r, s, t, h1, h2, sum;
    BOOL found;
    Element_t *elements = NULL, *newElements = NULL, *faces = NULL, *edges = NULL, *eptr=NULL, *eParent = NULL;
    Nodes_t *nodes, *newNodes;
    FEMMesh *newMesh;
    FEMBoundaryCondition *boundaryConditionAtId;
    FEMElementDescription *elementDescription;
    FEMListUtilities *listUtil;
    
    if (mesh == nil) return nil;
    
    newMesh = [[FEMMesh alloc] init];
    
    elementDescription = [FEMElementDescription sharedElementDescription];
    
    [self findEdgesForMesh:mesh findEdges:NULL];
    
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: ********** Old mesh **********\n");
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: nodes: %d.\n", mesh.numberOfNodes);
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: bulk elements: %d.\n", mesh.numberOfBulkElements);
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: boundary elements: %d.\n", mesh.numberOfBoundaryElements);
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: edges: %d.\n", mesh.numberOfEdges);
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: faces: %d.\n", mesh.numberOfFaces);
    
    // Update nodal coordinates
    nodeCnt = mesh.numberOfNodes + mesh.numberOfEdges;
    
    elements = mesh.getElements;
    faces = mesh.getFaces;
    edges = mesh.getEdges;
    nodes = mesh.getNodes;
    
    // For bricks, count faces
    for (i=0; i<mesh.numberOfFaces; i++) {
        if (faces[i].Type.NumberOfNodes == 4) nodeCnt++;
    }
    
    // For quads and bricks, count centerpoints
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        switch (elements[i].Type.ElementCode / 100) {
            case 4:
            case 8:
                nodeCnt++;
                break;
        }
    }
    
    // New mesh node coordinate arrays
    newNodes = newMesh.getNodes;
    newNodes->x = doublevec(0, nodeCnt-1);
    newNodes->y = doublevec(0, nodeCnt-1);
    newNodes->z = doublevec(0, nodeCnt-1);
    
    // New mesh includes old mesh nodes
    for (i=0; i<mesh.numberOfNodes; i++) {
        newNodes->x[i] = nodes->x[i];
        newNodes->y[i] = nodes->y[i];
        newNodes->z[i] = nodes->z[i];
    }
    
    if (h != NULL) {
        xh = doublevec(0, nodeCnt-1);
        for (i=0; i<*sizeNodal; i++) {
            xh[i] = h[i];
        }
    }
    
    // Add edge centers
    for (i=0; i<mesh.numberOfEdges; i++) {
        k = edges[i].Type.NumberOfNodes;
        j = i + mesh.numberOfNodes;
        if (h != NULL) {
            h1 = h[edges[i].NodeIndexes[0]];
            h2 = h[edges[i].NodeIndexes[1]];
            r = 1.0 / (1.0 + h1/h2);
            newNodes->x[j] = r*nodes->x[edges[i].NodeIndexes[0]] + (1.0-r)*nodes->x[edges[i].NodeIndexes[1]];
            newNodes->y[j] = r*nodes->y[edges[i].NodeIndexes[0]] + (1.0-r)*nodes->y[edges[i].NodeIndexes[1]];
            newNodes->z[j] = r*nodes->z[edges[i].NodeIndexes[0]] + (1.0-r)*nodes->z[edges[i].NodeIndexes[1]];
            xh[j] = r*h1+(1.0-r)*h2;
        } else {
            sum = 0.0;
            for (l=0; l<edges[i].Type.NumberOfNodes; l++) {
                sum = sum + nodes->x[edges[i].NodeIndexes[l]];
            }
            newNodes->x[j] = sum / k;
            
            sum = 0.0;
            for (l=0; l<edges[i].Type.NumberOfNodes; l++) {
                sum = sum + nodes->y[edges[i].NodeIndexes[l]];
            }
            newNodes->y[j] = sum / k;
            
            sum = 0.0;
            for (l=0; l<edges[i].Type.NumberOfNodes; l++) {
                sum = sum + nodes->z[edges[i].NodeIndexes[l]];
            }
            newNodes->z[j] = sum / k;
        }
    }
    
    nodals = doublevec(0, mesh.maxElementNodes-1);
    memset( nodals, 0.0, mesh.maxElementNodes*sizeof(double) );
    
    // Add face centers for bricks
    for (i=0; i<mesh.numberOfFaces; i++) {
        k = faces[i].Type.NumberOfNodes;
        if (k == 4) {
            j = i + mesh.numberOfNodes + mesh.numberOfEdges;
            if (h != NULL) {
                n = mesh.numberOfNodes;
                h1 = xh[n+faces[i].EdgeIndexes[1]];
                h2 = xh[n+faces[i].EdgeIndexes[3]];
                r = 2.0/(1.0 + h1/h2) - 1.0 ;
                h1 = xh[n+faces[i].EdgeIndexes[2]];
                h2 = xh[n+faces[i].EdgeIndexes[0]];
                s = 2.0/(1.0 + h1/h2) - 1.0;
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    nodals[l] = nodes->x[faces[i].NodeIndexes[l]];
                }
                newNodes->x[j] = [elementDescription interpolate2DInElement:&faces[i] nodalValues:nodals evaluatedAt:r andAt:s];
                
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    nodals[l] = nodes->y[faces[i].NodeIndexes[l]];
                }
                newNodes->y[j] = [elementDescription interpolate2DInElement:&faces[i] nodalValues:nodals evaluatedAt:r andAt:s];
                
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    nodals[l] = nodes->z[faces[i].NodeIndexes[l]];
                }
                newNodes->z[j] = [elementDescription interpolate2DInElement:&faces[i] nodalValues:nodals evaluatedAt:r andAt:s];
                
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    nodals[l] = h[faces[i].NodeIndexes[l]];
                }
                xh[j] = [elementDescription interpolate2DInElement:&faces[i] nodalValues:nodals evaluatedAt:r andAt:s];
            } else {
                sum = 0.0;
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    sum = sum + nodes->x[faces[i].NodeIndexes[l]];
                }
                newNodes->x[j] = sum / k;
                
                sum = 0.0;
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    sum = sum + nodes->y[faces[i].NodeIndexes[l]];
                }
                newNodes->y[j] = sum / k;
                
                sum = 0.0;
                for (l=0; l<faces[i].Type.NumberOfNodes; l++) {
                    sum = sum + nodes->z[faces[i].NodeIndexes[l]];
                }
                newNodes->z[j] = sum / k;
            }
        }
    }
    
    memset( nodals, 0.0, mesh.maxElementNodes*sizeof(double) );
    
    // Add centerpoint for quads & bricks
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        k = elements[i].Type.NumberOfNodes;
        switch (elements[i].Type.ElementCode / 100) {
            case 4:
                j++;
                if (h != NULL) {
                    n = mesh.numberOfNodes;
                    h1 = xh[n+elements[i].EdgeIndexes[1]];
                    h2 = xh[n+elements[i].EdgeIndexes[3]];
                    r = 2.0/(1.0 + h1/h2) - 1.0;
                    h1 = xh[n+elements[i].EdgeIndexes[2]];
                    h2 = xh[n+elements[i].EdgeIndexes[0]];
                    s = 2.0/(1.0 + h1/h2) - 1.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->x[elements[i].NodeIndexes[l]];
                    }
                    newNodes->x[j] = [elementDescription interpolate2DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s];
                    
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->y[elements[i].NodeIndexes[l]];
                    }
                    newNodes->y[j] = [elementDescription interpolate2DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s];
                    
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->z[elements[i].NodeIndexes[l]];
                    }
                    newNodes->z[j] = [elementDescription interpolate2DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s];
                } else {
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->x[elements[i].NodeIndexes[l]];
                    }
                    newNodes->x[j] = sum / k;
                    
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->y[elements[i].NodeIndexes[l]];
                    }
                    newNodes->y[j] = sum / k;
                    
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->z[elements[i].NodeIndexes[l]];
                    }
                    newNodes->z[j] = sum / k;
                }
                break;
                
            case 8:
                j++;
                if (h != NULL) {
                    n = mesh.numberOfNodes + mesh.numberOfEdges;
                    h1 = xh[n+elements[i].FaceIndexes[3]];
                    h2 = xh[n+elements[i].FaceIndexes[5]];
                    r = 2.0/(1.0 + h1/h2) - 1.0;
                    
                    h1 = xh[n+elements[i].FaceIndexes[4]];
                    h2 = xh[n+elements[i].FaceIndexes[2]];
                    s = 2.0/(1.0 + h1/h2) - 1.0;
                    
                    h1 = xh[n+elements[i].FaceIndexes[1]];
                    h2 = xh[n+elements[i].FaceIndexes[0]];
                    t = 2.0/(1.0 + h1/h2) - 1.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->x[elements[i].NodeIndexes[l]];
                    }
                    newNodes->x[j] = [elementDescription interpolate3DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s andAt:t];
                    
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->y[elements[i].NodeIndexes[l]];
                    }
                    newNodes->y[j] = [elementDescription interpolate3DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s andAt:t];
                    
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        nodals[l] = nodes->z[elements[i].NodeIndexes[l]];
                    }
                    newNodes->z[j] = [elementDescription interpolate3DInElement:&elements[i] nodalValues:nodals evaluatedAt:r andAt:s andAt:t];
                } else {
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->x[elements[i].NodeIndexes[l]];
                    }
                    newNodes->x[j] = sum / k;
                    
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->y[elements[i].NodeIndexes[l]];
                    }
                    newNodes->y[j] = sum / k;
                    
                    sum = 0.0;
                    for (l=0; l<elements[i].Type.NumberOfNodes; l++) {
                        sum = sum + nodes->z[elements[i].NodeIndexes[l]];
                    }
                    newNodes->z[j] = sum / k;
                }
                break;
        }
    }
    
    if (xh != NULL) {
        free_dvector(xh, 0, nodeCnt-1);
    }
    free_dvector(nodals, 0, mesh.maxElementNodes-1);
    
    // Update new mesh node count
    newMesh.numberOfEdges = 0;
    newMesh.numberOfFaces = 0;
    newMesh.maxBdofs = mesh.maxBdofs;
    newMesh.maxEdgeDofs = mesh.maxEdgeDofs;
    newMesh.maxFaceDofs = mesh.maxFaceDofs;
    newMesh.maxElementDofs = mesh.maxElementDofs;
    newMesh.dimension = mesh.dimension;
    
    newMesh.numberOfNodes = nodeCnt;
    newNodes->numberOfNodes = nodeCnt;
    
    // Update bulk elements
    
    // First count new elements
    newElCnt = 0;
    for (i=0; i<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; i++) {
        
        switch (elements[i].Type.ElementCode / 100) {
            // Each element will be divided into 2^dim new elements
            case 2:
                newElCnt = newElCnt + 2;  // Lines
                break;
            case 3:
                newElCnt = newElCnt + 4;  // Trias
                break;
            case 4:
                newElCnt = newElCnt + 4;  // Quads
                break;
            case 5:
                newElCnt = newElCnt + 8;  // Tetras
                break;
            case 8:
                newElCnt = newElCnt + 8;  // Hexas
                break;
        }
    }
    newElements = newMesh.getElements;
    newElements = (Element_t*) malloc( sizeof(Element_t) * newElCnt);
    initElements(newElements, newElCnt);
    for (i=0; i<newElCnt; i++) {
        newElements[i].EdgeIndexes = NULL;
        newElements[i].FaceIndexes = NULL;
    }
    
    child = intmatrix(0, mesh.numberOfBulkElements-1, 0, 7);
    newElCnt = 0;
    nodeCnt = mesh.numberOfNodes;
    edgeCnt = mesh.numberOfEdges;
    
    // Index to old quad/hexa centerpoint node in the new nodal arrays
    node = nodeCnt + mesh.numberOfEdges + mesh.numberOfFaces;

    // Now update all new mesh elements
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        
        switch (elements[i].Type.ElementCode) {
            case 303:
                // Split triangle to four triangles from
                // edge centerpoints
                
                // 1st new element
                child[i][0] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                
                // 2nd new element
                child[i][1] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElCnt++;
                
                // 3rd new element
                child[i][2] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                
                // 4th new element
                child[i][3] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                break;
                
            case 404:
                // Index to old quad centerpoint node in the
                // new mesh nodal arrays
                // node = node not node = node + 1 as it's the case in Elmer
                
                // Spit quad to four new quads from edge
                // centerpoints and centerpoint of the element
                
                // 1st new element
                child[i][0] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = node;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElCnt++;
                
                // 2nd new element
                child[i][1] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = node;
                newElCnt++;
                
                // 3rd new element
                child[i][2] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = node;
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                
                // 4th new element
                child[i][3] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = node;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].NodeIndexes[3];
                newElCnt++;
                break;
                
            case 504:
                // Split tetra to 8 new elements from
                // corners and edge centerpoints
                
                // 1st new element
                child[i][0] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElCnt++;
                
                // 2nd new element
                child[i][1] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[4] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElCnt++;
                
                // 3rd new element
                child[i][2] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[5] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElCnt++;
                
                // 4th new element
                child[i][3] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[3];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[5] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[4] + nodeCnt;
                newElCnt++;
                
                // Then the annoying part; we will have to split the remaining octahedron into four elements. This can be done in three ways
                // of which only one preserves the minimum angle condition (Delaunay splitting)
                dxyz[0][0] = newNodes->x[elements[i].EdgeIndexes[3] + nodeCnt] - newNodes->x[elements[i].EdgeIndexes[1] + nodeCnt];
                dxyz[1][0] = newNodes->y[elements[i].EdgeIndexes[3] + nodeCnt] - newNodes->y[elements[i].EdgeIndexes[1] + nodeCnt];
                dxyz[2][0] = newNodes->z[elements[i].EdgeIndexes[3] + nodeCnt] - newNodes->z[elements[i].EdgeIndexes[1] + nodeCnt];
                
                dxyz[0][1] = newNodes->x[elements[i].EdgeIndexes[4] + nodeCnt] - newNodes->x[elements[i].EdgeIndexes[2] + nodeCnt];
                dxyz[1][1] = newNodes->y[elements[i].EdgeIndexes[4] + nodeCnt] - newNodes->y[elements[i].EdgeIndexes[2] + nodeCnt];
                dxyz[2][1] = newNodes->z[elements[i].EdgeIndexes[4] + nodeCnt] - newNodes->z[elements[i].EdgeIndexes[2] + nodeCnt];
                
                dxyz[0][2] = newNodes->x[elements[i].EdgeIndexes[5] + nodeCnt] - newNodes->x[elements[i].EdgeIndexes[0] + nodeCnt];
                dxyz[1][2] = newNodes->y[elements[i].EdgeIndexes[5] + nodeCnt] - newNodes->y[elements[i].EdgeIndexes[0] + nodeCnt];
                dxyz[2][2] = newNodes->z[elements[i].EdgeIndexes[5] + nodeCnt] - newNodes->z[elements[i].EdgeIndexes[0] + nodeCnt];
                
                dist[0] = sqrt( pow(dxyz[0][0], 2.0) + pow(dxyz[1][0], 2.0) + pow(dxyz[2][0], 2.0) );
                dist[1] = sqrt( pow(dxyz[0][1], 2.0) + pow(dxyz[1][1], 2.0) + pow(dxyz[2][1], 2.0) );
                dist[2] = sqrt( pow(dxyz[0][2], 2.0) + pow(dxyz[1][2], 2.0) + pow(dxyz[2][2], 2.0) );
                
                diag = 1; // The default diaginal for splitting is between edges 2-4
                if (dist[1] < dist[0] && dist[1] < dist[2]) diag = 2; // Edges 3-5
                if (dist[2] < dist[0] && dist[2] < dist[1]) diag = 3; // Edges 1-6
                
                switch (diag) {
                    case 1:
                        // 5th new element
                        child[i][4] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElCnt++;
                        
                        // 6th new element
                        child[i][5] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElCnt++;
                        
                        // 7th new element
                        child[i][6] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElCnt++;
                        
                        // 8th new element
                        child[i][7] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElCnt++;
                        break;
                        
                    case 2:
                        // 5th new element
                        child[i][4] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElCnt++;
                        
                        // 6th new element
                        child[i][5] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElCnt++;
                        
                        // 7th new element
                        child[i][6] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElCnt++;
                        
                        // 8th new element
                        child[i][7] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElCnt++;
                        break;
                        
                    case 3:
                        // 5th new element
                        child[i][4] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElCnt++;
                        
                        // 6th new element
                        child[i][5] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElCnt++;
                        
                        // 7th new element
                        child[i][6] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[4] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElCnt++;
                        
                        // 8th new element
                        child[i][7] = newElCnt;
                        newElements[newElCnt] = elements[i];
                        newElements[newElCnt].ElementIndex = newElCnt+1;
                        newElements[newElCnt].NodeIndexes = intvec(0, 3);
                        newElements[newElCnt].sizeNodeIndexes = 4;
                        newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[2] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[3] + nodeCnt;
                        newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[5] + nodeCnt;
                        newElCnt++;
                        break;
                }
                
            case 808:
                // Index to old quad centerpoints node in the new mesh nodal arrays
                // node = node not node = node + 1 as it's the case in Elmer
                
                // Split brick to 8 new bricks from edge centerpoints
                // and centerpoint of the element
                
                // 1st new element
                child[i][0] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].FaceIndexes[0] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].EdgeIndexes[8] + nodeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].FaceIndexes[2] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[6] = node;
                newElements[newElCnt].NodeIndexes[7] = elements[i].FaceIndexes[5] + nodeCnt + edgeCnt;
                newElCnt++;
                
                // 2nd new element
                child[i][1] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[0] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].FaceIndexes[0] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].FaceIndexes[2] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].EdgeIndexes[9] + nodeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].FaceIndexes[3] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[7] = node;
                newElCnt++;
                
                // 3rd new element
                child[i][2] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[3] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].FaceIndexes[0] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].NodeIndexes[3];
                newElements[newElCnt].NodeIndexes[4] = elements[i].FaceIndexes[5] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[5] = node;
                newElements[newElCnt].NodeIndexes[6] = elements[i].FaceIndexes[4] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].EdgeIndexes[11] + nodeCnt;
                newElCnt++;
                
                // 4th new element
                child[i][3] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].FaceIndexes[0] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[1] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[2] + nodeCnt;
                newElements[newElCnt].NodeIndexes[4] = node;
                newElements[newElCnt].NodeIndexes[5] = elements[i].FaceIndexes[3] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].EdgeIndexes[10] + nodeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].FaceIndexes[4] + nodeCnt + edgeCnt;
                newElCnt++;
                
                // 5th new element
                child[i][4] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].EdgeIndexes[8] + nodeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].FaceIndexes[2] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[2] = node;
                newElements[newElCnt].NodeIndexes[3] = elements[i].FaceIndexes[5] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].NodeIndexes[4];
                newElements[newElCnt].NodeIndexes[5] = elements[i].EdgeIndexes[4] + nodeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].FaceIndexes[1] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].EdgeIndexes[7] + nodeCnt;
                newElCnt++;
                
                // 6th new element
                child[i][5] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].FaceIndexes[2] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[1] = elements[i].EdgeIndexes[9] + nodeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].FaceIndexes[3] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[3] = node;
                newElements[newElCnt].NodeIndexes[4] = elements[i].EdgeIndexes[4] + nodeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].NodeIndexes[5];
                newElements[newElCnt].NodeIndexes[6] = elements[i].EdgeIndexes[5] + nodeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].FaceIndexes[1] + nodeCnt + edgeCnt;
                newElCnt++;
                
                // 7th new element
                child[i][6] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = elements[i].FaceIndexes[5] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[1] = node;
                newElements[newElCnt].NodeIndexes[2] = elements[i].FaceIndexes[4] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].EdgeIndexes[11] + nodeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].EdgeIndexes[7] + nodeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].FaceIndexes[1] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].EdgeIndexes[6] + nodeCnt;
                newElements[newElCnt].NodeIndexes[7] = elements[i].NodeIndexes[7];
                newElCnt++;
                
                // 8th new element
                child[i][7] = newElCnt;
                newElements[newElCnt] = elements[i];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 7);
                newElements[newElCnt].sizeNodeIndexes = 8;
                newElements[newElCnt].NodeIndexes[0] = node;
                newElements[newElCnt].NodeIndexes[1] = elements[i].FaceIndexes[3] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[2] = elements[i].EdgeIndexes[10] + nodeCnt;
                newElements[newElCnt].NodeIndexes[3] = elements[i].FaceIndexes[4] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[4] = elements[i].FaceIndexes[1] + nodeCnt + edgeCnt;
                newElements[newElCnt].NodeIndexes[5] = elements[i].EdgeIndexes[5] + nodeCnt;
                newElements[newElCnt].NodeIndexes[6] = elements[i].NodeIndexes[6];
                newElements[newElCnt].NodeIndexes[7] = elements[i].EdgeIndexes[6] + nodeCnt;
                newElCnt++;
                break;
                
            default:
                fprintf(stderr, "FEMMeshUtils:splitMeshEqual: element type %d not supported by the multigrid solver.\n", elements[i].Type.ElementCode);
                fatal("FEMMeshUtils:splitMeshEqual");
        }
    }
    
    
    // Update new mesh element counts
    newMesh.numberOfBulkElements = newElCnt;
    
    // Update boundary elements
    // Note:  internal boundaries not taken care of....!!!
    eoldNodes = intvec(0, 3);
    faceNodes = intvec(0, 3);
    edgeNodes = intvec(0, 1);
    for (i=0; i<mesh.numberOfBoundaryElements; i++) {
                
        // Get parent of the boundary element
        eParent = elements[i+mesh.numberOfBulkElements].BoundaryInfo->Left;
        if (eParent == NULL) eParent = elements[i+mesh.numberOfBulkElements].BoundaryInfo->Right;
        if (eParent == NULL) continue;
        
        parentID = eParent->ElementIndex-1;
        
        switch (elements[i+mesh.numberOfBulkElements].Type.ElementCode / 100) {
            case 2:
                // Line segments
                
                // Which edge of the parent element are we?
                for (edge1=0; edge1<eParent->sizeEdgeIndexes; edge1++) {
                    if ((elements[i+mesh.numberOfBulkElements].NodeIndexes[0] == edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0]
                         && elements[i+mesh.numberOfBulkElements].NodeIndexes[1] == edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]) ||
                        (elements[i+mesh.numberOfBulkElements].NodeIndexes[1] == edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0] &&
                         elements[i+mesh.numberOfBulkElements].NodeIndexes[0] == edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1])) break;
                }
                
                // Index of the old edge centerpoint in the
                // new mesh nodal arrays
                node = eParent->EdgeIndexes[edge1] + mesh.numberOfNodes;
                
                // 1st new element
                newElements[newElCnt] = elements[i+mesh.numberOfBulkElements];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 1);
                newElements[newElCnt].sizeNodeIndexes = 2;
                newElements[newElCnt].NodeIndexes[0] = elements[i+mesh.numberOfBulkElements].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = node;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfBulkElements].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the
                // children of the old mesh parent element
                for (j=0; j<4; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    found = NO;
                    for (k=0; k<n-1; k++) {
                        if ((newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[k] &&
                            newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[k+1]) ||
                            (newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[k] &&
                             newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[k+1])) {
                                found = YES;
                                break;
                        }
                    }
                    if ((newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[n-1] &&
                         newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[0]) ||
                        (newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[n-1] &&
                         newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[0])) {
                            found = YES;
                    }
                    if (found == YES) break;
                }
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 2nd new element
                newElements[newElCnt] = elements[i+mesh.numberOfBulkElements];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 1);
                newElements[newElCnt].sizeNodeIndexes = 2;
                newElements[newElCnt].NodeIndexes[0] = node;
                newElements[newElCnt].NodeIndexes[1] = elements[i+mesh.numberOfBulkElements].NodeIndexes[1];
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfBulkElements].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the
                // children of the old mesh parent element
                for (j=0; j<4; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    found = NO;
                    for (k=0; k<n-1; k++) {
                        if ((newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[k] &&
                             newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[k+1]) ||
                            (newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[k] &&
                             newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[k+1])) {
                                found = YES;
                                break;
                        }
                    }
                    if ((newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[n-1] &&
                         newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[0]) ||
                        (newElements[newElCnt].NodeIndexes[1] == newElements[child[parentID][j]].NodeIndexes[n-1] &&
                         newElements[newElCnt].NodeIndexes[0] == newElements[child[parentID][j]].NodeIndexes[0])) {
                            found = YES;
                    }
                    if (found == YES) break;
                }
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                break;
                
            case 3:
                // Trias
                
                // On which face of the parent element are we?
                for (j=0; j<3; j++) {
                    eoldNodes[j] = elements[i+mesh.numberOfBulkElements].NodeIndexes[j];
                }
                vDSP_vsort((float *)eoldNodes, 3, 1);
                
                for (faceNumber=0; faceNumber<eParent->sizeFaceIndexes; faceNumber++) {
                    for (j=0; j<3; j++) {
                        faceNodes[j] = faces[eParent->FaceIndexes[faceNumber]].NodeIndexes[j];
                    }
                    vDSP_vsort((float *)faceNodes, 3, 1);
                    
                    if (eoldNodes[0] == faceNodes[0] && eoldNodes[1] == faceNodes[1] && eoldNodes[2] == faceNodes[2]) break;
                }
                
                // Then what are the edges on this face?
                
                // First edge
                eoldNodes[0] = min(elements[i+mesh.numberOfBulkElements].NodeIndexes[0], elements[i+mesh.numberOfBulkElements].NodeIndexes[1]);
                eoldNodes[1] = max(elements[i+mesh.numberOfBulkElements].NodeIndexes[0], elements[i+mesh.numberOfBulkElements].NodeIndexes[1]);
                for (edge1=0; edge1<eParent->sizeEdgeIndexes; edge1++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Second edge
                eoldNodes[0] = min(elements[i+mesh.numberOfBulkElements].NodeIndexes[1], elements[i+mesh.numberOfBulkElements].NodeIndexes[2]);
                eoldNodes[1] = max(elements[i+mesh.numberOfBulkElements].NodeIndexes[1], elements[i+mesh.numberOfBulkElements].NodeIndexes[2]);
                for (edge2=0; edge2<eParent->sizeEdgeIndexes; edge2++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge2]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge2]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge2]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge2]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Third edge
                eoldNodes[0] = min(elements[i+mesh.numberOfBulkElements].NodeIndexes[2], elements[i+mesh.numberOfBulkElements].NodeIndexes[0]);
                eoldNodes[1] = max(elements[i+mesh.numberOfBulkElements].NodeIndexes[2], elements[i+mesh.numberOfBulkElements].NodeIndexes[0]);
                for (edge3=0; edge3<eParent->sizeEdgeIndexes; edge3++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge3]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge3]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge3]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge3]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Index of the old face and edge centerpoints
                // in the new mesh nodal arrays
                node12 = eParent->EdgeIndexes[edge1] + mesh.numberOfNodes;
                node23 = eParent->EdgeIndexes[edge2] + mesh.numberOfNodes;
                node31 = eParent->EdgeIndexes[edge3] + mesh.numberOfNodes;
                
                // 1st new element
                newElements[newElCnt] = elements[i+mesh.numberOfBulkElements];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = elements[i+mesh.numberOfBulkElements].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = node12;
                newElements[newElCnt].NodeIndexes[2] = node31;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfBulkElements].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0; // Count matches
                    for (n1=0; n1<3; n1++) {
                        for (n2=0; n2<4; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) fatal("FEMMeshUtils:splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 2nd new element
                newElements[newElCnt] = elements[i+mesh.numberOfBulkElements];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = node12;
                newElements[newElCnt].NodeIndexes[1] = elements[i+mesh.numberOfBulkElements].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = node23;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfBulkElements].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<3; n1++) {
                        for (n2=0; n2<4; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) fatal("FEMMeshUtils:splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 3rd new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = node12;
                newElements[newElCnt].NodeIndexes[1] = node23;
                newElements[newElCnt].NodeIndexes[2] = node31;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among the children
                // of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<3; n1++) {
                        for (n2=0; n2<4; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) fatal("FEMMeshUtils:splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 4th new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 2);
                newElements[newElCnt].sizeNodeIndexes = 3;
                newElements[newElCnt].NodeIndexes[0] = node31;
                newElements[newElCnt].NodeIndexes[1] = node23;
                newElements[newElCnt].NodeIndexes[2] = elements[i+mesh.numberOfNodes].NodeIndexes[2];
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mesh parent element among
                // the children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<3; n1++) {
                        for (n2=0; n2<4; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) fatal("FEMMeshUtils:splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                break;
                
            case 4:
                // Quads
                
                // On which face of the parent element are we?
                for (j=0; j<4; j++) {
                    eoldNodes[j] = elements[i+mesh.numberOfNodes].NodeIndexes[j];
                }
                vDSP_vsort((float *)eoldNodes, 4, 1);
                
                for (faceNumber=0; faceNumber<eParent->sizeFaceIndexes; faceNumber++) {
                    for (j=0; j<4; j++) {
                        faceNodes[j] = faces[eParent->FaceIndexes[faceNumber]].NodeIndexes[j];
                    }
                    vDSP_vsort((float *)faceNodes, 4, 1);
                    
                    if (eoldNodes[0] == faceNodes[0] && eoldNodes[1] == faceNodes[1] && eoldNodes[2] == faceNodes[2]) break;
                }
                
                // Then what are the edges on this face?
                
                // First edge
                eoldNodes[0] = min(elements[i+mesh.numberOfNodes].NodeIndexes[0], elements[i+mesh.numberOfNodes].NodeIndexes[1]);
                eoldNodes[1] = max(elements[i+mesh.numberOfNodes].NodeIndexes[0], elements[i+mesh.numberOfNodes].NodeIndexes[1]);
                for (edge1=0; edge1<eParent->sizeEdgeIndexes; edge1++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge1]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge1]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Second edge
                eoldNodes[0] = min(elements[i+mesh.numberOfNodes].NodeIndexes[1], elements[i+mesh.numberOfNodes].NodeIndexes[2]);
                eoldNodes[1] = max(elements[i+mesh.numberOfNodes].NodeIndexes[1], elements[i+mesh.numberOfNodes].NodeIndexes[2]);
                for (edge2=0; edge2<eParent->sizeEdgeIndexes; edge2++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge2]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge2]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge2]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge2]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Third edge
                eoldNodes[0] = min(elements[i+mesh.numberOfNodes].NodeIndexes[2], elements[i+mesh.numberOfNodes].NodeIndexes[3]);
                eoldNodes[1] = max(elements[i+mesh.numberOfNodes].NodeIndexes[2], elements[i+mesh.numberOfNodes].NodeIndexes[3]);
                for (edge3=0; edge3<eParent->sizeEdgeIndexes; edge3++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge3]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge3]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge3]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge3]].NodeIndexes[1]);
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Fourth edge
                eoldNodes[0] = min(elements[i+mesh.numberOfNodes].NodeIndexes[3], elements[i+mesh.numberOfNodes].NodeIndexes[0]);
                eoldNodes[1] = max(elements[i+mesh.numberOfNodes].NodeIndexes[3], elements[i+mesh.numberOfNodes].NodeIndexes[0]);
                for (edge4=0; edge4<eParent->sizeEdgeIndexes; edge4++) {
                    edgeNodes[0] = min(edges[eParent->EdgeIndexes[edge4]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge4]].NodeIndexes[1]);
                    edgeNodes[1] = max(edges[eParent->EdgeIndexes[edge4]].NodeIndexes[0], edges[eParent->EdgeIndexes[edge4]].NodeIndexes[1]);
                    
                    if (eoldNodes[0] == edgeNodes[0] && eoldNodes[1] == edgeNodes[1]) break;
                }
                
                // Index of the old face and edge centerpoints
                // in the new mesh nodal arrays
                node = eParent->FaceIndexes[faceNumber] + // Face mid-point
                        mesh.numberOfNodes + mesh.numberOfEdges;
                node12 = eParent->EdgeIndexes[edge1] + mesh.numberOfNodes;
                node23 = eParent->EdgeIndexes[edge2] + mesh.numberOfNodes;
                node34 = eParent->EdgeIndexes[edge3] + mesh.numberOfNodes;
                node41 = eParent->EdgeIndexes[edge4] + mesh.numberOfNodes;
                
                // 1st new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = elements[i+mesh.numberOfNodes].NodeIndexes[0];
                newElements[newElCnt].NodeIndexes[1] = node12;
                newElements[newElCnt].NodeIndexes[2] = node;
                newElements[newElCnt].NodeIndexes[3] = node41;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mech parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<4; n1++) {
                        for (n2=0; n2<8; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) fatal("FEMMeshUtils:splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 2nd new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = node12;
                newElements[newElCnt].NodeIndexes[1] = elements[i+mesh.numberOfNodes].NodeIndexes[1];
                newElements[newElCnt].NodeIndexes[2] = node23;
                newElements[newElCnt].NodeIndexes[3] = node;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mech parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<4; n1++) {
                        for (n2=0; n2<8; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) fatal("FEMMeshUtils:splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 3rd new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = node41;
                newElements[newElCnt].NodeIndexes[1] = node;
                newElements[newElCnt].NodeIndexes[2] = node34;
                newElements[newElCnt].NodeIndexes[3] = elements[i+mesh.numberOfNodes].NodeIndexes[3];
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mech parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0;  // Count matches
                    for (n1=0; n1<4; n1++) {
                        for (n2=0; n2<8; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2)  break;
                }
                if (n3 < 3) fatal("FEMMeshUtils:splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                // 4th new element
                newElements[newElCnt] = elements[i+mesh.numberOfNodes];
                newElements[newElCnt].ElementIndex = newElCnt+1;
                newElements[newElCnt].NodeIndexes = intvec(0, 3);
                newElements[newElCnt].sizeNodeIndexes = 4;
                newElements[newElCnt].NodeIndexes[0] = node;
                newElements[newElCnt].NodeIndexes[1] = node23;
                newElements[newElCnt].NodeIndexes[2] = elements[i+mesh.numberOfNodes].NodeIndexes[2];
                newElements[newElCnt].NodeIndexes[3] = node34;
                newElements[newElCnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
                initBoundaryInfo(newElements[newElCnt].BoundaryInfo);
                newElements[newElCnt].BoundaryInfo = elements[i+mesh.numberOfNodes].BoundaryInfo;
                newElements[newElCnt].BoundaryInfo->Left = NULL;
                newElements[newElCnt].BoundaryInfo->Right = NULL;
                
                // Search the new mech parent element among the
                // children of the old mesh parent element
                for (j=0; j<8; j++) {
                    eptr = &newElements[child[parentID][j]];
                    n = newElements[child[parentID][j]].Type.NumberOfNodes;
                    n3 = 0; // Count matches
                    for (n1=0; n1<4; n1++) {
                        for (n2=0; n2<8; n2++) {
                            if (newElements[newElCnt].NodeIndexes[n1] == newElements[child[parentID][j]].NodeIndexes[n2]) n3++;
                        }
                    }
                    if (n3 > 2) break;
                }
                if (n3 < 3) fatal("FEMMeshUtils:splitMeshEqual", "Parent element not found.");
                newElements[newElCnt].BoundaryInfo->Left = eptr;
                newElCnt++;
                
                break;
        }
        
    }
    free_ivector(eoldNodes, 0, 3);
    free_ivector(faceNodes, 0, 3);
    free_ivector(edgeNodes, 0, 1);
    
    // Update new mesh boundary element counts
    newMesh.numberOfBoundaryElements = newElCnt - newMesh.numberOfBulkElements;
    newMesh.maxElementDofs = mesh.maxElementDofs;
    newMesh.maxElementNodes = mesh.maxElementNodes;
    
    for (i=0; i<newMesh.numberOfBulkElements+newMesh.numberOfBoundaryElements; i++) {
        newElements[i].EdgeIndexes = NULL;
        newElements[i].FaceIndexes = NULL;
    }
    
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: ********** New mesh **********\n");
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: nodes: %d.\n", newMesh.numberOfNodes);
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: bulk elements: %d.\n", newMesh.numberOfBulkElements);
    fprintf(stdout, "FEMMeshUtils:splitMeshEqual: boundary elements: %d.\n", newMesh.numberOfBoundaryElements);
    
    // TODO: add support for parallel run
    
    // If periodic BC given, compute boundary mesh projector
    listUtil = [FEMListUtilities sharedListUtilities];
    for (i=0; i<model.numberOfBoundaryConditions; i++) {
        boundaryConditionAtId = (model.boundaryConditions)[i];
        minVal = 1;
        maxVal = model.numberOfBoundaryConditions;
        k = [listUtil listGetInteger:model inArray:boundaryConditionAtId.valuesList forVariable:@"periodic bc" info:&found minValue:&minVal maxValue:&maxVal];
        boundaryConditionAtId.pMatrix = [self periodicProjectorInModel:model forMesh:mesh masterBoundary:i targetBoundary:k-1 galerking:NULL];
    }
    
    free_imatrix(child, 0, mesh.numberOfBulkElements-1, 0, 7);
    
    [mesh deallocateMeshEdgeTables];
    [mesh deallocateMeshFaceTables];
    
    return newMesh;
}

-(void)setStabilizationParametersInMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model {
    
    int i, j, n;
    Element_t *elements = NULL;
    Nodes_t *nodes = NULL, *meshNodes = NULL;
    FEMElementDescription *elementDescription;
    
    for (FEMSolution *solution in model.solutions) {
        if (mesh == solution.mesh) {
            mesh.stabilize = (mesh.stabilize == YES || [(solution.solutionInfo)[@"stabilize"] boolValue] == YES) ? YES : NO;
            mesh.stabilize = (mesh.stabilize == YES || [(solution.solutionInfo)[@"stabilization method"] isEqualToString:@"vms"] == YES) ? YES : NO;
            mesh.stabilize = (mesh.stabilize == YES || [(solution.solutionInfo)[@"stabilization method"] isEqualToString:@"stabilized"] == YES) ? YES : NO;
        }
    }
    
    nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    nodes->x = doublevec(0, mesh.maxElementNodes-1);
    nodes->y = doublevec(0, mesh.maxElementNodes-1);
    nodes->z = doublevec(0, mesh.maxElementNodes-1);
    
    elements = mesh.getElements;
    meshNodes = mesh.getNodes;
    
    elementDescription = [FEMElementDescription sharedElementDescription];
    
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        n = elements[i].Type.NumberOfNodes;
        for (j=0; j<n; j++) {
            nodes->x[j] = meshNodes->x[elements[i].NodeIndexes[j]];
            nodes->y[j] = meshNodes->y[elements[i].NodeIndexes[j]];
            nodes->z[j] = meshNodes->z[elements[i].NodeIndexes[j]];
        }
        if (mesh.stabilize == YES) {
            [elementDescription computeStabilizationParameterInElement:&elements[i] nodes:nodes mesh:mesh numberOfNodes:n mk:&elements[i].StabilizationMK hk:&elements[i].hK];
        } else {
            [elementDescription elementDiameter:&elements[i] nodes:nodes];
        }
    }
    
    free_dvector(nodes->x, 0, mesh.maxElementNodes-1);
    free_dvector(nodes->y, 0, mesh.maxElementNodes-1);
    free_dvector(nodes->z, 0, mesh.maxElementNodes-1);
}

-(void)setCurrentMesh:(FEMMesh * __nonnull)mesh inModel:(FEMModel * __nonnull)model {
    
    Element_t *meshElements = NULL;
    Nodes_t *meshNodes = NULL;
    
    model.variables = mesh.variables;
    model.mesh = mesh;
    
    meshElements = mesh.getElements;
    meshNodes = mesh.getNodes;
    
    [model setNodes:meshNodes];
    model.numberOfNodes = mesh.numberOfNodes;
    
    [model SetElements:meshElements];
    model.maxElementNodes = mesh.maxElementNodes;
    model.numberOfBulkElements = mesh.numberOfBulkElements;
    model.numberOfBoundaryElements = mesh.numberOfBoundaryElements;
}

-(void)updateMesh:(FEMMesh * __nonnull)mesh inSolution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model {
    
    int i, j, k, n, dofs, n1, n2;
    int *permutation;
    double *work;
    BOOL onlySearch, found, optimizeBandwidth;
    NSString *str;
    FEMUtilities *utilities;
    FEMElementUtils *elementUtils;
    FEMMatrix *matrix;
    variableArraysContainer *saveVarContainers = NULL, *varContainers = NULL;
    matrixArraysContainer *matContainers = NULL, *solutionMatContainers = NULL;
    
    saveVarContainers = solution.variable.getContainers;
    dofs = solution.variable.dofs;
    
    solution.mesh = mesh;
    [self setCurrentMesh:mesh inModel:model];
    
    // Create matrix and variable structures for
    // current equation on the mesh
    utilities = [[FEMUtilities alloc] init];
    onlySearch = NO;
    solution.variable = [utilities getVariableFrom:mesh.variables model:model name:solution.variable.name onlySearch:&onlySearch maskName:NULL info:&found];
    
    varContainers = solution.variable.getContainers;
    permutation = intvec(0, varContainers->sizePerm-1);
    
    if ((solution.solutionInfo)[@"optimize bandwidth"] != nil) {
        optimizeBandwidth = [(solution.solutionInfo)[@"optimize bandwidth"] boolValue];
    } else optimizeBandwidth = YES;
    
    elementUtils = [[FEMElementUtils alloc] init];
    str = (solution.solutionInfo)[@"equation"];
    matrix = [elementUtils createMatrixInModel:model forSolution:solution mesh:mesh dofs:dofs permutation:permutation sizeOfPermutation:varContainers->sizePerm matrixFormat:MATRIX_CRS optimizeBandwidth:optimizeBandwidth equationName:str discontinuousGalerkinSolution:NULL globalBubbles:NULL nodalDofsOnly:NULL projectorDofs:NULL];
    
    if ((solution.solutionInfo)[@"linear system symmetric"] != nil) {
        matrix.symmetric = [(solution.solutionInfo)[@"linear system symmetric"] boolValue];
    }
    
    if ((solution.solutionInfo)[@"lumped mass matrix"] != nil) {
        matrix.lumped = [(solution.solutionInfo)[@"lumped mass matrix"] boolValue];
    }
    
    work = doublevec(0, varContainers->sizeValues-1);
    memcpy(work, varContainers->Values, varContainers->sizeValues*sizeof(double));
    for (k=1; k<=dofs; k++) {
        for (i=0; i<varContainers->sizePerm; i++) {
            if (permutation[i] >= 0) {
                varContainers->Values[dofs*(permutation[i]+1)-k] = work[dofs*(varContainers->Perm[i]+1)-k];
            }
        }
    }
    
    if (varContainers->PrevValues != NULL) {
        for (j=0; j<varContainers->size2PrevValues; j++) {
            for (i=0; i<varContainers->size1PrevValues; i++) {
                work[i] = varContainers->PrevValues[i][j];
            }
             for (k=1; k<=dofs; k++) {
                 for (i=0; i<varContainers->sizePerm; i++) {
                     if (permutation[i] >= 0) {
                         varContainers->PrevValues[dofs*(permutation[i]+1)-k][j] = work[dofs*(varContainers->Perm[i]+1)-k];
                     }
                 }
             }
        }
    }
    free_dvector(work, 0, varContainers->sizeValues-1);
    
    memcpy(varContainers->Perm, permutation, varContainers->sizePerm*sizeof(int));
    solution.variable.solution = solution;
    free_ivector(permutation, 0, varContainers->sizePerm-1);
    
    matContainers = matrix.getContainers;
    matContainers->RHS = doublevec(0, matrix.numberOfRows-1);
    matContainers->sizeRHS = matrix.numberOfRows;
    
    solutionMatContainers = solution.matrix.getContainers;
    
    if (saveVarContainers->EigenValues != NULL) {
        n = saveVarContainers->sizeEigenValues;
        
        if (n > 0) {
            solution.nOfEigenValues = n;
            varContainers->EigenValues = cdoublevec(0, n-1);
            varContainers->sizeEigenValues = n;
            varContainers->EigenVectors = cdoublematrix(0, n-1, 0, varContainers->sizeValues-1);
            varContainers->size1EigenVectors = n;
            varContainers->size2EigenVectors = varContainers->sizeValues;
            
            for (i=0; i<n; i++) {
                varContainers->EigenValues[i] = 0.0;
                for (j=0; j<varContainers->sizeValues; j++) {
                    varContainers->EigenVectors[i][j] = 0.0;
                }
            }
            matContainers->MassValues = doublevec(0, matContainers->sizeValues-1);
            matContainers->sizeMassValues = matContainers->sizeValues;
            memset( matContainers->MassValues, 0.0, matContainers->sizeMassValues*sizeof(double) );
        }
    } else if (solutionMatContainers->Force != NULL) {
        n1 = matrix.numberOfRows;
        n2 = solutionMatContainers->size2Force;
        matContainers->Force = doublematrix(0, n1-1, 0, n2-1);
        matContainers->size1force = n1;
        matContainers->size2Force = n2;
        memset( *matContainers->Force, 0.0, (n1*n2)*sizeof(double) );
    }
    
    solution.matrix = matrix;
    solution.mesh.changed = YES;
}

-(void)allocatePDefinitionsElement:(Element_t * __nonnull)element {
    
    element->Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
    if (element->Pdefs == NULL) fatal("FEMMeshUtils:allocatePDefinitionsElement", "Unable to allocate memory.");
    
    // Initialize fields
    element->Pdefs->p = 0;
    element->Pdefs->TetraType = 0;
    element->Pdefs->isEdge = false;
    element->Pdefs->PyramidQuadEdge = false;
    element->Pdefs->LocalNumber = 0;
    element->Pdefs->GaussPoints = 0;
}

-(void)setEdgeFaceDofsMesh:(FEMMesh * __nonnull)mesh edgeDofs:(int * __nullable)edgeDofs faceDofs:(int * __nullable)faceDofs {
    
    Element_t *elements = NULL, *edges = NULL, *faces = NULL;
    
    [self findEdgesForMesh:mesh findEdges:NULL];
    
    FEMPElementMaps *elementMaps = [[FEMPElementMaps alloc] init];
    
    elements = mesh.getElements;
    edges = mesh.getEdges;
    faces = mesh.getFaces;
    
    // Set edge and face polynomial degree and degrees of freedom for
    // all elements
    for (int i=0; i<mesh.numberOfBulkElements; i++) {
        // Iterate each edge of element
        for (int j=0; j<elements[i].Type.NumberOfEdges; j++) {
            // Set attributes of p element edges
            if (elements[i].Pdefs != NULL) {
                // Set edge polynomial degree and dofs
                edges[elements[i].EdgeIndexes[j]].Pdefs->p = max(elements[i].Pdefs->p, edges[elements[i].EdgeIndexes[j]].Pdefs->p);
                edges[elements[i].EdgeIndexes[j]].BDOFs = max(edges[elements[i].EdgeIndexes[j]].BDOFs, edges[elements[i].EdgeIndexes[j]].Pdefs->p-1);
                edges[elements[i].EdgeIndexes[j]].Pdefs->isEdge = true;
                // Get Gauss points for edge. If no dofs, two gauss points are
                // still needed for integration of linear equation
                edges[elements[i].EdgeIndexes[j]].Pdefs->GaussPoints = (int)pow((double)(edges[elements[i].EdgeIndexes[j]].BDOFs+2), (double)edges[elements[i].EdgeIndexes[j]].Type.dimension);
                
                if (edges[elements[i].EdgeIndexes[j]].BoundaryInfo->Left != NULL) {
                    [self assignLocalNumberToEdgeElement:&edges[elements[i].EdgeIndexes[j]] fromElement:edges[elements[i].EdgeIndexes[j]].BoundaryInfo->Left inMesh:mesh];
                } else {
                    [self assignLocalNumberToEdgeElement:&edges[elements[i].EdgeIndexes[j]] fromElement:edges[elements[i].EdgeIndexes[j]].BoundaryInfo->Right inMesh:mesh];
                }
            }
            // Other element types, which need edge dofs
            else if (edgeDofs != NULL) {
                edges[elements[i].EdgeIndexes[j]].BDOFs = max(edgeDofs[i], edges[elements[i].EdgeIndexes[j]].BDOFs);
            }
            
            // Get maximum dof for edges
            mesh.maxEdgeDofs = max(edges[elements[i].EdgeIndexes[j]].BDOFs, mesh.maxEdgeDofs);
        }
        
        // Iterate each face of element
        for (int j=0; j<elements[i].Type.NumberOfFaces; j++) {
            // Set attributes of p element faces
            if (elements[i].Pdefs != NULL) {
                // Set face polynomial degree and dofs
                faces[elements[i].FaceIndexes[j]].Pdefs->p = max(elements[i].Pdefs->p, faces[elements[i].FaceIndexes[j]].Pdefs->p);
                // Get number of face dofs
                faces[elements[i].FaceIndexes[j]].BDOFs = max(faces[elements[i].FaceIndexes[j]].BDOFs, [elementMaps getFaceDofsForElement:&elements[i] polyDegree:faces[elements[i].FaceIndexes[j]].Pdefs->p faceNumber:j]);
                faces[elements[i].FaceIndexes[j]].Pdefs->isEdge = true;
                faces[elements[i].FaceIndexes[j]].Pdefs->GaussPoints = [elementMaps getNumberOfGaussPointsForElement:&faces[elements[i].FaceIndexes[j]] inMesh:mesh];
                if (faces[elements[i].FaceIndexes[j]].BoundaryInfo->Left != NULL) {
                    [self assignLocalNumberToEdgeElement:&faces[elements[i].FaceIndexes[j]] fromElement:faces[elements[i].FaceIndexes[j]].BoundaryInfo->Left inMesh:mesh];
                } else {
                    [self assignLocalNumberToEdgeElement:&faces[elements[i].FaceIndexes[j]] fromElement:faces[elements[i].FaceIndexes[j]].BoundaryInfo->Right inMesh:mesh];
                }
            } else if (faceDofs != NULL) {
                faces[elements[i].FaceIndexes[j]].BDOFs = max(faceDofs[i], faces[elements[i].FaceIndexes[j]].BDOFs);
            }
            
            // Get maximum dof for faces
            mesh.maxFaceDofs = max(faces[elements[i].FaceIndexes[j]].BDOFs, mesh.maxFaceDofs);
        }
    }
    
    // Set local edges for boundary elements
    for (int i=mesh.numberOfBulkElements; i<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; i++) {
        // Here set local number and copy attributes to this boundary element for left parent
        if (elements[i].BoundaryInfo->Left != NULL) {
            // Local edges are only assigned for p elements
            if (elements[i].BoundaryInfo->Left->Pdefs != NULL) {
                [self allocatePDefinitionsElement:&elements[i]];
                elements[i].Pdefs->isEdge = true;
                [self assignLocalNumberToEdgeElement:&elements[i] fromElement:elements[i].BoundaryInfo->Left inMesh:mesh];
            }
        }
        
        // Here set local number and copy attributes to this boundary element for right parent
        if (elements[i].BoundaryInfo->Right != NULL) {
            // Local edges are only assigned for pe elements
            if (elements[i].BoundaryInfo->Right->Pdefs != NULL) {
                [self allocatePDefinitionsElement:&elements[i]];
                elements[i].Pdefs->isEdge = true;
                [self assignLocalNumberToEdgeElement:&elements[i] fromElement:elements[i].BoundaryInfo->Right inMesh:mesh];
            }
        }
    }
    [elementMaps deallocation];
}

-(void)setMaximumDofsMesh:(FEMMesh * __nonnull)mesh {
    
    Element_t *elements = mesh.getElements;
    FEMPElementMaps *elementMaps = [[FEMPElementMaps alloc] init];
    
    // Set Gauss points for each p element
    for (int i=0; i<mesh.numberOfBulkElements; i++) {
        if (elements[i].Pdefs != NULL) {
            elements[i].Pdefs->GaussPoints = [elementMaps getNumberOfGaussPointsForElement:&elements[i] inMesh:mesh];
        }
        
        // Set maximum element dofs here (because element size may have changed when
        // edges and faces have been set). This is the absolute worst case.
        // Element which has maxElementDofs may not even be present as a real element
        mesh.maxElementDofs = max(mesh.maxElementDofs,
                                  elements[i].Type.NumberOfNodes +
                                  elements[i].Type.NumberOfEdges*mesh.maxEdgeDofs +
                                  elements[i].Type.NumberOfFaces*mesh.maxFaceDofs +
                                  elements[i].BDOFs, elements[i].DGDOFs);
        mesh.maxBdofs = max(elements[i].BDOFs, mesh.maxBdofs);
    }
    
    for (int i=0; i<mesh.numberOfBulkElements; i++) {
        if (elements[i].BDOFs > 0) {
            elements[i].BubbleIndexes = intvec(0, elements[i].BDOFs-1);
            elements[i].sizeBubbleIndexes = elements[i].BDOFs;
            for (int j=0; j<elements[i].BDOFs; j++) {
                elements[i].BubbleIndexes[j] = mesh.maxBdofs*i+j;
            }
        }
    }
    [elementMaps deallocation];
}

/*************************************************************************
 
    Given a 2D mesh, extrude it to 3D. The 3rd coordinate will always be
    at the interval [0,1]. Therefore, the adaptation for different shapes
    must be done with StructuredMeshMapper or some similar utility.
    The top and bottom surface will be assigned Boundary Conditions tags
    with indexes one larger than the maximum used on by the 2D mesh.
 
    Method corresponds to Elmer from git on October 27 2015

*************************************************************************/
-(FEMMesh * __nonnull)extrudeMesh:(FEMMesh * __nonnull)mesh inLevels:(int)inLevels model:(FEMModel * __nonnull)model {
    
    int j, k, l, n, bcid, cnt, dg_n, extrudedCoord, ind[8], ln, max_body, nnodes;
    double currentCoord, w;
    double *activeCoord = NULL, *wTable;
    Element_t *elementsIn = NULL, *elementsOut = NULL;
    ElementType_t *elmType = NULL;
    Nodes_t *nodesIn = NULL, *nodesOut = NULL;
    FEMBoundaryCondition *boundaryConditionAtId;
    BOOL found, needEdges;
    
    FEMListUtilities * listUtilities = [FEMListUtilities sharedListUtilities];
    FEMMesh *meshOut = [[FEMMesh alloc] init];
    
    // TODO: Add support for parallel run
    
    // Generate volume nodal points
    n = mesh.numberOfNodes;
    nnodes = (inLevels + 2) * n;
    nodesIn = mesh.getNodes;
    
    nodesOut = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(nodesOut);
    nodesOut->x = doublevec(0, nnodes-1);
    nodesOut->y = doublevec(0, nnodes-1);
    nodesOut->z = doublevec(0, nnodes-1);
    nodesOut->numberOfNodes = nnodes;
    
    int gelements = mesh.numberOfBulkElements;
    // TODO: Add support for parallel run
    
    // Create the division for the 1D unit mesh
    wTable = doublevec(0, inLevels+1);
    [self FEMMeshUtils_unitSegmentDivisionTable:wTable levels:inLevels+1 model:model listUtilities:listUtilities];
    
    int minv = 1;
    int maxv = 3;
    extrudedCoord = [listUtilities listGetInteger:model inArray:model.simulation.valuesList forVariable:@"extruded coordinate index" info:&found minValue:&minv maxValue:&maxv];
    if (found == NO) extrudedCoord = 3;
    
    if (extrudedCoord == 1) {
        activeCoord = nodesOut->x;
    } else if (extrudedCoord == 2) {
        activeCoord = nodesOut->y;
    } else if (extrudedCoord == 3) {
        activeCoord = nodesOut->z;
    }
    
    BOOL preserveBaseline = [listUtilities listGetLogical:model inArray:model.simulation.valuesList forVariable:@"preserve baseline" info:&found];
    if (found == NO) preserveBaseline = NO;
    
    double minCoord = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"extruded minimum coordinate" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) minCoord = 0.0;
    
    double maxCoord = [listUtilities listGetConstReal:model inArray:model.simulation.valuesList forVariable:@"extruded maximum coordinate" info:&found minValue:NULL maxValue:NULL];
    if (found == NO) maxCoord = 1.0;
    
    cnt = 0;
    for (int i=0; i<=inLevels+1; i++) {
        w = wTable[i];
        currentCoord = w * maxCoord + (1.0 - w) * minCoord;
        for (j=0; j<mesh.numberOfNodes; j++) {
            nodesOut->x[cnt] = nodesIn->x[j];
            nodesOut->y[cnt] = nodesIn->y[j];
            nodesOut->z[cnt] = nodesIn->z[j];
            
            // Override the coordinate in the extruded direction by the value on the layer
            activeCoord[cnt] = currentCoord;
            
            // TODO: Add support for parallel run
            
            cnt++;
        }
    }
    meshOut.numberOfNodes = cnt;
    
    elementsIn = mesh.getElements;
    // Count 101 elements:
    // (this requires an extra layer)
    int cnt101 = 0;
    for (int i=mesh.numberOfBulkElements; i<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; i++) {
        if (elementsIn[i].Type.ElementCode == 101) cnt101++;
    }
    
    n = mesh.numberOfElements;
    if (preserveBaseline == YES) {
        elementsOut = (Element_t*) malloc( sizeof(Element_t) * (n*(inLevels+3) + mesh.numberOfBoundaryElements + cnt101) );
        initElements(elementsOut, (n*(inLevels+3) + mesh.numberOfBoundaryElements + cnt101));
    } else {
        elementsOut = (Element_t*) malloc( sizeof(Element_t) * (n*(inLevels+3) + cnt101) );
        initElements(elementsOut, (n*(inLevels+3) + cnt101));
    }
    
    // Generate volume bulk elements
    FEMElementDescription *elementDescription = [FEMElementDescription sharedElementDescription];
    meshOut.maxElementNodes = 0;
    cnt = 0, dg_n = 0;
    needEdges = NO;
    n = mesh.numberOfNodes;
    for (int i=0; i<=inLevels; i++) {
        for (j=0; j<mesh.numberOfBulkElements; j++) {
            
            memcpy(&elementsOut[cnt], &elementsIn[j], sizeof(Element_t));
            
            ln = 0;
            for (k=0; k<elementsIn[j].Type.NumberOfNodes; k++) {
                ind[ln] = elementsIn[j].NodeIndexes[k]+i*n;
                ln++;
            }
            for (k=0; k<elementsIn[j].Type.NumberOfNodes; k++) {
                ind[ln] = elementsIn[j].NodeIndexes[k]+(i+1)*n;
                ln++;
            }
            elementsOut[cnt].NDOFs = ln;
            meshOut.maxElementNodes = max(meshOut.maxElementNodes, ln);
            
            switch (ln) {
                case 6:
                    elmType = NULL;
                    elmType = [elementDescription getElementType:706 inMesh:mesh stabilization:NULL];
                    elementsOut[cnt].Type = *elmType;
                    break;
                case 8:
                    elmType = NULL;
                    elmType = [elementDescription getElementType:808 inMesh:mesh stabilization:NULL];
                    elementsOut[cnt].Type = *elmType;
                    break;
            }
            
            elementsOut[cnt].GElementIndex = elementsIn[j].GElementIndex + gelements * i;
            
            elementsOut[cnt].ElementIndex = cnt + 1;
            elementsOut[cnt].NodeIndexes = intvec(0, ln-1);
            elementsOut[cnt].sizeNodeIndexes = ln;
            memcpy(elementsOut[cnt].NodeIndexes, ind, ln*sizeof(int));
            elementsOut[cnt].DGIndexes = NULL;
            elementsOut[cnt].EdgeIndexes = NULL;
            elementsOut[cnt].FaceIndexes = NULL;
            elementsOut[cnt].BubbleIndexes = NULL;
            
            k = elementsOut[cnt].DGDOFs;
            if (k > 0) {
                elementsOut[cnt].DGDOFs = elementsOut[cnt].Type.NumberOfNodes;
                k = elementsOut[cnt].DGDOFs;
                elementsOut[cnt].DGIndexes = intvec(0, k-1);
                elementsOut[cnt].sizeDGIndexes = k;
                for (int l=0; l<k; l++) {
                    elementsOut[cnt].DGIndexes[l] = dg_n;
                    dg_n++;
                }
                needEdges = YES;
            }
            
            if (elementsIn[j].Pdefs != NULL) {
                needEdges = YES;
                elementsOut[cnt].Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
                memcpy(elementsOut[cnt].Pdefs, elementsIn[j].Pdefs, sizeof(PElementDefs_t));
            }
            cnt++;
        }
    }
    meshOut.numberOfBulkElements = cnt;
    
    int max_bid = 0;
    int max_baseline_bid = 0;
    
    if (preserveBaseline == YES) {
        for (int j=0; j<mesh.numberOfBoundaryElements; j++) {
            k = j + mesh.numberOfBulkElements;
            
            memcpy(&elementsOut[cnt], &elementsIn[k], sizeof(Element_t));
            
            elementsOut[cnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
            memcpy(elementsOut[cnt].BoundaryInfo, elementsIn[k].BoundaryInfo, sizeof(BoundaryInfo_t));
            
            max_bid = max(max_bid, elementsIn[k].BoundaryInfo->Constraint);
            
            if (elementsIn[k].BoundaryInfo->Left != NULL) {
                l = elementsIn[k].BoundaryInfo->Left->ElementIndex-1;
                elementsOut[cnt].BoundaryInfo->Left = &elementsOut[mesh.numberOfBulkElements*(inLevels+1) + (inLevels+2)*mesh.numberOfBoundaryElements+l];
            }
            if (elementsIn[k].BoundaryInfo->Right != NULL) {
                l = elementsIn[k].BoundaryInfo->Right->ElementIndex-1;
                elementsOut[cnt].BoundaryInfo->Right = &elementsOut[mesh.numberOfBulkElements*(inLevels+1) + (inLevels+2)*mesh.numberOfBoundaryElements+l];
            }
            
            if (elementsIn[k].Type.ElementCode >= 200) {
                elementsOut[cnt].NDOFs = 2;
                elementsOut[cnt].NodeIndexes = intvec(0, 1);
                elementsOut[cnt].sizeNodeIndexes = 2;
                ind[0] = elementsIn[k].NodeIndexes[0];
                ind[1] = elementsIn[k].NodeIndexes[1];
                memcpy(elementsOut[cnt].NodeIndexes, ind, 2*sizeof(int));
                elmType = NULL;
                elmType = [elementDescription getElementType:202 inMesh:mesh stabilization:NULL];
                elementsOut[cnt].Type = *elmType;
            } else {
                elementsOut[cnt].NDOFs = 1;
                l = elementsIn[k].sizeNodeIndexes;
                elementsOut[cnt].NodeIndexes = intvec(0, l-1);
                elementsOut[cnt].sizeNodeIndexes = l;
                memcpy(elementsOut[cnt].NodeIndexes, elementsIn[k].NodeIndexes, l*sizeof(int));
                elementsOut[cnt].Type = elementsIn[k].Type;
            }
            elementsOut[cnt].DGDOFs = 0;
            elementsOut[cnt].DGIndexes = NULL;
            elementsOut[cnt].ElementIndex = cnt + 1;
            elementsOut[cnt].Pdefs = NULL;
            elementsOut[cnt].EdgeIndexes = NULL;
            elementsOut[cnt].FaceIndexes = NULL;
            elementsOut[cnt].BubbleIndexes = NULL;
            
            cnt++;
        }
        // TODO: Add support for parallel run
        max_baseline_bid = max_bid;
    }
    
    // Add side boundaries with the bottom mesh boundary id's
    // (or shift ids if preserving the baseline boundary)
    for (int i=0; i<=inLevels; i++) {
        for (j=0; j<mesh.numberOfBoundaryElements; j++) {
            k = j + mesh.numberOfBulkElements;
            
            memcpy(&elementsOut[cnt], &elementsIn[k], sizeof(Element_t));
            elementsOut[cnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
            memcpy(elementsOut[cnt].BoundaryInfo, elementsIn[k].BoundaryInfo, sizeof(BoundaryInfo_t));
            
            elementsOut[cnt].BoundaryInfo->Constraint = elementsOut[cnt].BoundaryInfo->Constraint + max_baseline_bid;
            
            max_bid = max(max_bid, max_baseline_bid + elementsIn[k].BoundaryInfo->Constraint);
            
            if (elementsIn[k].BoundaryInfo->Left != NULL) {
                l = elementsIn[k].BoundaryInfo->Left->ElementIndex-1;
                elementsOut[cnt].BoundaryInfo->Left = &elementsOut[mesh.numberOfBulkElements*i+l];
            }
            if (elementsIn[k].BoundaryInfo->Right != NULL) {
                l = elementsIn[k].BoundaryInfo->Right->ElementIndex-1;
                elementsOut[cnt].BoundaryInfo->Right = &elementsOut[mesh.numberOfBulkElements*i+l];
            }
            
            if (elementsIn[k].Type.ElementCode >= 200) {
                elementsOut[cnt].NDOFs = 4;
                elementsOut[cnt].NodeIndexes = intvec(0, 3);
                elementsOut[cnt].sizeNodeIndexes = 4;
                
                ind[0] = elementsIn[k].NodeIndexes[0]+i*n;
                ind[1] = elementsIn[k].NodeIndexes[1]+i*n;
                ind[2] = elementsIn[k].NodeIndexes[1]+(i+1)*n;
                ind[3] = elementsIn[k].NodeIndexes[0]+(i+1)*n;
                memcpy(elementsOut[cnt].NodeIndexes, ind, 4*sizeof(int));
                elmType = NULL;
                elmType = [elementDescription getElementType:404 inMesh:mesh stabilization:NULL];
                elementsOut[cnt].Type = *elmType;
            } else {
                elementsOut[cnt].NDOFs = 1;
                elementsOut[cnt].NodeIndexes = intvec(0, elementsIn[k].sizeNodeIndexes-1);
                elementsOut[cnt].sizeNodeIndexes = elementsIn[k].sizeNodeIndexes;
                for (int l=0; l<elementsIn[k].sizeNodeIndexes; l++) {
                    elementsOut[cnt].NodeIndexes[l] = elementsIn[k].NodeIndexes[l]+i*n;
                }
                elementsOut[cnt].Type = elementsIn[k].Type;
            }
            elementsOut[cnt].ElementIndex = cnt + 1;
            elementsOut[cnt].DGDOFs = 0;
            elementsOut[cnt].DGIndexes = NULL;
            elementsOut[cnt].Pdefs = NULL;
            elementsOut[cnt].EdgeIndexes = NULL;
            elementsOut[cnt].FaceIndexes = NULL;
            elementsOut[cnt].BubbleIndexes = NULL;
            cnt++;
        }
    }
    
    // Take care of extra 101 elements
    if (cnt101 > 0) {
        for (int j=0; j<mesh.numberOfBoundaryElements; j++) {
            k = j + mesh.numberOfBulkElements;
            
            if (elementsIn[k].Type.ElementCode != 101) continue;
            
            memcpy(&elementsOut[cnt], &elementsIn[k], sizeof(Element_t));
            elementsOut[cnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
            memcpy(elementsOut[cnt].BoundaryInfo, elementsIn[k].BoundaryInfo, sizeof(BoundaryInfo_t));
            
            elementsOut[cnt].BoundaryInfo->Constraint = elementsOut[cnt].BoundaryInfo->Constraint + max_baseline_bid;
            
            max_bid = max(max_bid, max_baseline_bid + elementsIn[k].BoundaryInfo->Constraint);
            
            elementsOut[cnt].NDOFs = 1;
            elementsOut[cnt].NodeIndexes = intvec(0, 0);
            elementsOut[cnt].sizeNodeIndexes = 1;
            for (int l=0; l<elementsIn[k].sizeNodeIndexes; l++) {
                elementsOut[cnt].NodeIndexes[l] = elementsIn[k].NodeIndexes[l]+(inLevels+1)*n;
            }
            elementsOut[cnt].Type = elementsIn[k].Type;
            
            elementsOut[cnt].ElementIndex = cnt + 1;
            elementsOut[cnt].DGDOFs = 0;
            elementsOut[cnt].DGIndexes = NULL;
            elementsOut[cnt].Pdefs = NULL;
            elementsOut[cnt].EdgeIndexes = NULL;
            elementsOut[cnt].FaceIndexes = NULL;
            elementsOut[cnt].BubbleIndexes = NULL;
            cnt++;
        }
    }
    
    // TODO: Add support for parallel run
    
    fprintf(stdout, "FEMMeshUtils:extrudeMesh: first extruded BC set to: %d.\n", max_bid+1);
    
    max_body = 0;
    for (int i=0; i<mesh.numberOfBulkElements; i++) {
        max_body = max(max_body, elementsIn[i].BodyID);
    }
    // TODO: Add support for parallel run
    fprintf(stdout, "FEMMeshUtils:extrudeMesh: number of new BCs for layers: %d.\n", max_body);
    
    // Add bottom boundary
    for (int i=0; i<mesh.numberOfBulkElements; i++) {

        memcpy(&elementsOut[cnt], &elementsIn[i], sizeof(Element_t));
        
        ln = elementsIn[i].Type.NumberOfNodes;
        elementsOut[cnt].NDOFs = ln;
        
        elementsOut[cnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
        initBoundaryInfo(elementsOut[cnt].BoundaryInfo);
        elementsOut[cnt].BoundaryInfo->Left = &elementsOut[i];
        elementsOut[cnt].BoundaryInfo->Right = NULL;
        
        bcid = max_bid + elementsOut[cnt].BodyID;
        elementsOut[cnt].BoundaryInfo->Constraint = bcid;
        
        elementsOut[cnt].BodyID = 0;
        if (bcid <= model.numberOfBoundaryConditions) {
            boundaryConditionAtId = (model.boundaryConditions)[bcid-1];
            j = [listUtilities listGetInteger:model inArray:boundaryConditionAtId.valuesList forVariable:@"body id" info:&found minValue:NULL maxValue:NULL];
            if (found == YES) elementsOut[cnt].BodyID = j;
        }
        
        elementsOut[cnt].NodeIndexes = intvec(0, ln-1);
        elementsOut[cnt].sizeNodeIndexes = ln;
        memcpy(elementsOut[cnt].NodeIndexes, elementsIn[i].NodeIndexes, ln*sizeof(int));
        elementsOut[cnt].ElementIndex = cnt + 1;
        elementsOut[cnt].Type = elementsIn[i].Type;
        elementsOut[cnt].DGDOFs = 0;
        elementsOut[cnt].DGIndexes = NULL;
        elementsOut[cnt].Pdefs = NULL;
        elementsOut[cnt].EdgeIndexes = NULL;
        elementsOut[cnt].FaceIndexes = NULL;
        elementsOut[cnt].BubbleIndexes = NULL;
        cnt++;
    }
    
    // Add top boundary
    for (int i=0; i<mesh.numberOfBulkElements; i++) {

        memcpy(&elementsOut[cnt], &elementsIn[i], sizeof(Element_t));
        
        ln = elementsIn[i].Type.NumberOfNodes;
        elementsOut[cnt].NDOFs = ln;
        
        elementsOut[cnt].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
        initBoundaryInfo(elementsOut[cnt].BoundaryInfo);
        elementsOut[cnt].BoundaryInfo->Left = &elementsOut[inLevels*mesh.numberOfBulkElements+i];
        elementsOut[cnt].BoundaryInfo->Right = NULL;
        
        bcid = max_bid + elementsOut[cnt].BodyID + max_body;
        elementsOut[cnt].BoundaryInfo->Constraint = bcid;
        
        elementsOut[cnt].BodyID = 0;
        if (bcid <= model.numberOfBoundaryConditions) {
            boundaryConditionAtId = (model.boundaryConditions)[bcid-1];
            j = [listUtilities listGetInteger:model inArray:boundaryConditionAtId.valuesList forVariable:@"body id" info:&found minValue:NULL maxValue:NULL];
            if (found == YES) elementsOut[cnt].BodyID = j;
        }
        
        elementsOut[cnt].NodeIndexes = intvec(0, ln-1);
        elementsOut[cnt].sizeNodeIndexes = ln;
        for (j=0; j<ln; j++) {
            elementsOut[cnt].NodeIndexes[j] = elementsIn[i].NodeIndexes[j] + (inLevels+1)*n;
        }
        elementsOut[cnt].ElementIndex = cnt + 1;
        elementsOut[cnt].Type = elementsIn[i].Type;
        elementsOut[cnt].DGDOFs = 0;
        elementsOut[cnt].DGIndexes = NULL;
        elementsOut[cnt].Pdefs = NULL;
        elementsOut[cnt].EdgeIndexes = NULL;
        elementsOut[cnt].FaceIndexes = NULL;
        elementsOut[cnt].BubbleIndexes = NULL;
        cnt++;
    }
    
    [meshOut assignNodes:nodesOut];
    [meshOut assignElements:elementsOut];
    meshOut.numberOfBoundaryElements = cnt - meshOut.numberOfBulkElements;
    meshOut.name = mesh.name;
    meshOut.discontinuousMesh = mesh.discontinuousMesh;
    meshOut.maxElementDofs = meshOut.maxElementNodes;
    meshOut.dimension = 3;
    if (model.dimension != 3) model.dimension = 3;
    
    if (needEdges == YES) [self setEdgeFaceDofsMesh:meshOut edgeDofs:NULL faceDofs:NULL];
    [self setMaximumDofsMesh:meshOut];
    
    free_dvector(wTable, 0, inLevels+1);
    return meshOut;
}

/*********************************************************************************************************************************
 
    This method finds the structure of an extruded mesh even though it is given in an unstructured format. It may be used by some 
    special solutions that employ the special character of the mesh. The extrusion is found for a given direction and for each 
    node the corresponding up and down, and thereafter top and bottom node is computed.
 
    Returns the integer number of elements in the pointer structures
 
    Method corresponds to Elmer from git on October 27 2015

*********************************************************************************************************************************/
-(FEMVariable * __nullable)detectExtrudedStructureMesh:(FEMMesh * __nonnull)mesh solution:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model numberofNodes:(int)numberofNodes ifMask:(BOOL)ifMask mask:(variableArraysContainer * __nullable)mask isTopActive:(BOOL)isTopActive topNodePointer:(int * __nullable)topNodePointer isBottomActive:(BOOL)isBottomActive bottomNodePointer:(int * __nullable)bottomNodePointer isUpActive:(BOOL)isUpActive upNodePointer:(int * __nullable)upNodePointer isDownActive:(BOOL)isDownActive downNodePointer:(int * __nullable)downNodePointer numberOfLayers:(int * __nullable)numberOfLayers isMidNode:(BOOL)isMidNode midNodePointer:(int * __nullable)midNodePointer midLayerExists:(BOOL * __nullable)midLayerExists isNodeLayer:(BOOL)isNodeLayer nodeLayer:(int * __nullable)nodeLayer {
    
    int activeDirection, dim, ii, jj;
    double at0, at1, dotProduct, elemVector[3], eps, length, unitVector[3], *values = NULL, vector[3], vector2[3];
    NSString *coordTransform;
    FEMVariable *variable;
    variableArraysContainer *variableContainers = NULL;
    BOOL found;
    
    fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: determining extruded structure.\n");
    at0 = cputime();
    
    dim = mesh.dimension;
    if ((solution.solutionInfo)[@"active coordinate"] != nil) {
        activeDirection = [(solution.solutionInfo)[@"active coordinate"] intValue];
    }
    if (activeDirection < 1 || activeDirection > 3) {
        fatal("FEMMeshUtils:detectExtrudedStructureMesh", "Invalid value for active coordinate.");
    }
    memset( unitVector, 0.0, sizeof(unitVector) );
    unitVector[activeDirection-1] = 1.0;
    
    if ((solution.solutionInfo)[@"project to bottom"] != nil) {
        if ([(solution.solutionInfo)[@"project to bottom"] boolValue] == YES) {
            for (int i=0; i<3; i++) {
                unitVector[i] = -1.0 * unitVector[i];
            }
        }
    }
    for (int i=0; i<3; i++) {
        fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: unit vector of direction: %f\n", unitVector[i]);
    }
    
    // Set the dot product tolerance
    if ((solution.solutionInfo)[@"dot product tolerance"] != nil) {
        eps = [(solution.solutionInfo)[@"dot product tolerance"] doubleValue];
    } else eps = 1.0e-4;
    
    BOOL doCoordTransform = NO;
    if ((solution.solutionInfo)[@"mapping coordinate transformation"] != nil || ifMask == YES) {
        doCoordTransform = YES;
        coordTransform = (solution.solutionInfo)[@"mapping coordinate transformation"];
        values = doublevec(0, numberofNodes-1);
        memset( values, 0.0, numberofNodes*sizeof(double) );
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        variableArraysContainer *addContainers = allocateVariableContainer();
        addContainers->Values = values;
        addContainers->sizeValues = numberofNodes;
        if (ifMask == YES) {
            addContainers->Perm = mask->Perm;
            addContainers->sizePerm = numberofNodes;
            [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"extruded coordinate" dofs:1 container:addContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        } else {
            [utilities addVariableTo:mesh.variables mesh:mesh solution:solution name:@"extruded coordinate" dofs:1 container:addContainers component:NO ifOutput:NULL ifSecondary:NULL type:NULL];
        }
        variable = [utilities getVariableFrom:mesh.variables model:model name:@"extruded coordinate" onlySearch:NULL maskName:nil info:&found];
    } else if (activeDirection == 1) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        variable = [utilities getVariableFrom:mesh.variables model:model name:@"coordinate 1" onlySearch:NULL maskName:nil info:&found];
    } else if (activeDirection == 2) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        variable = [utilities getVariableFrom:mesh.variables model:model name:@"coordinate 2" onlySearch:NULL maskName:nil info:&found];
    } else if (activeDirection == 3) {
        FEMUtilities *utilities = [[FEMUtilities alloc] init];
        variable = [utilities getVariableFrom:mesh.variables model:model name:@"coordinate 3" onlySearch:NULL maskName:nil info:&found];
    }
    variableContainers = variable.getContainers;
    
    if (ifMask == YES || doCoordTransform == YES) {
        int j;
        Nodes_t *nodes = mesh.getNodes;
        for (int i=0; i<mesh.numberOfNodes; i++) {
            j = i;
            if (ifMask == YES) j = mask->Perm[i];
            vector[0] = nodes->x[i];
            vector[1] = nodes->y[i];
            vector[2] = nodes->z[i];
            if (doCoordTransform == YES) {
                [self FEMMeshUtils_coordinateTransformationNodalType:coordTransform vector:vector solution:solution];
            }
            values[j] = vector[activeDirection-1];
        }
    }
    
    // Check which direction is active
    BOOL upActive = NO, downActive = NO;
    upActive = (isUpActive == YES || isTopActive == YES) ? YES : NO;
    downActive = (isDownActive == YES || isBottomActive == YES) ? YES : NO;
    
    if (numberOfLayers != NULL || isNodeLayer == YES) {
        upActive = YES;
        downActive = YES;
    }
    
    if (!(upActive == YES || downActive == YES)) {
        fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: either up or down direction should be active.\n");
        return 0;
    }
    
    // Allocate pointers to top and bottom, and temporary pointers up and down
    if (upActive == YES) {
        for (int i=0; i<numberofNodes; i++) {
            topNodePointer[i] = i;
            upNodePointer[i] = i;
        }
    }
    if (downActive == YES) {
        for (int i=0; i<numberofNodes; i++) {
            bottomNodePointer[i] = i;
            downNodePointer[i] = i;
        }
    }
    
    // Determine the up and down pointers using dot product as criterian
    fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: determine the up and down pointers.\n");
    int n = mesh.maxElementNodes;
    Nodes_t *nodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    initNodes(nodes);
    nodes->x = doublevec(0, n-1);
    nodes->y = doublevec(0, n-1);
    nodes->z = doublevec(0, n-1);
    Element_t *meshElements = mesh.getElements;
    Nodes_t *meshNodes = mesh.getNodes;
    for (int t=0; t<mesh.numberOfBulkElements; t++) {
        n = meshElements[t].Type.NumberOfNodes;
        
        for (int i=0; i<n; i++) {
            nodes->x[i] = meshNodes->x[meshElements[t].NodeIndexes[i]];
            nodes->y[i] = meshNodes->y[meshElements[t].NodeIndexes[i]];
            nodes->z[i] = meshNodes->z[meshElements[t].NodeIndexes[i]];
        }
        
        // TODO: add support for parrallel run
        
        if (ifMask == YES) {
            BOOL any = NO;
            for (int i=0; i<n; i++) {
                if (mask->Perm[meshElements[t].NodeIndexes[i]] < 0) {
                    any = YES;
                    break;
                }
            }
            if (any == YES) continue;
        }
        
        for (int i=0; i<n; i++) {
            ii = meshElements[t].NodeIndexes[i];
            vector[0] = nodes->x[i];
            vector[1] = nodes->y[i];
            vector[2] = nodes->z[i];
            if (doCoordTransform == YES) {
                [self FEMMeshUtils_coordinateTransformationNodalType:coordTransform vector:vector solution:solution];
            }
            for (int j=i+1; j<n; j++) {
                jj = meshElements[t].NodeIndexes[j];
                vector2[0] = nodes->x[j];
                vector2[1] = nodes->y[j];
                vector2[2] = nodes->z[j];
                if (doCoordTransform == YES) {
                    [self FEMMeshUtils_coordinateTransformationNodalType:coordTransform vector:vector2 solution:solution];
                }

                vDSP_vsubD(vector, 1, vector2, 1, elemVector, 1, 3);
                length = sqrt(cblas_ddot(3, elemVector, 1, elemVector, 1));
                dotProduct = cblas_ddot(3, elemVector, 1, unitVector, 1) / length;
                
                if (dotProduct > 1.0 - eps) {
                    if (upActive == YES) upNodePointer[ii] = jj;
                    if (downActive == YES) downNodePointer[jj] = ii;
                } else if (dotProduct < eps - 1.0) {
                    if (downActive == YES) downNodePointer[ii] = jj;
                    if (upActive == YES) upNodePointer[jj] = ii;
                }
            }
        }
    }
    free_dvector(nodes->x, 0, mesh.maxElementNodes-1);
    free_dvector(nodes->y, 0, mesh.maxElementNodes-1);
    free_dvector(nodes->z, 0, mesh.maxElementNodes-1);
    free(nodes);
    
    // Pointer to top and bottom are found recursively using up and down
    fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: determine top and bottom pointers.\n");
    int j, downHit, rounds, upHit;
    for (rounds=1; rounds<=numberofNodes; rounds++) {
        downHit = 0;
        upHit = 0;
        for (int i=0; i<numberofNodes; i++) {
            if (ifMask == YES) {
                if (mask->Perm[i] < 0) continue;
            }
            if (upActive == YES) {
                j = upNodePointer[i];
                if (topNodePointer[i] != topNodePointer[j]) {
                    upHit++;
                    topNodePointer[i] = topNodePointer[j];
                }
            }
            if (downActive == YES) {
                j = downNodePointer[i];
                if (bottomNodePointer[i] != bottomNodePointer[j]) {
                    downHit++;
                    bottomNodePointer[i] = bottomNodePointer[j];
                }
            }
        }
        if (upHit == 0 && downHit == 0) break;
    }
    // The last round is always a check
    rounds--;
    
    fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: layered structure detected in %d cycles.\n", rounds);
    if (rounds == 0) {
        fprintf(stderr, "FEMMeshUtils:detectExtrudedStructureMesh: try to increase value for > dot product tolerance <.\n");
        fprintf(stderr, "FEMMeshUtils:detectExtrudedStructureMesh: zero rounds implies unsuccessfull operations.\n");
        fatal("FEMMeshUtils:detectExtrudedStructureMesh");
    }
    
    // Compute the number of layers. The rounds above may in some cases
    // be too small. Here just one layer is used to determine the number
    // of layers to save some time.
    if (numberOfLayers != NULL) {
        fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: compute the number of layers.\n");
        int i, k;
        for (i=0; i<numberofNodes; i++) {
            if (ifMask == YES) {
                if (mask->Perm < 0) continue;
            }
            break;
        }
        
        j = bottomNodePointer[i];
        *numberOfLayers = 0;
        while (1) {
            k = upNodePointer[j];
            if (k == j) {
                break;
            } else {
                *numberOfLayers = *numberOfLayers + 1;
                j = k;
            }
        }
        
        if (*numberOfLayers < rounds) {
            fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: there seems to be varying number of layers: %d vs. %d.\n", *numberOfLayers, rounds);
            *numberOfLayers = rounds;
        }
        fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: extruded structure layers: %d.\n", *numberOfLayers);
    }
    
    // Create layer index if requested
    if (isNodeLayer == YES) {
        fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: creating layer index.\n");
        int k;
        memset( nodeLayer, 1, numberofNodes*sizeof(int) );
        if (ifMask == YES) {
            for (int i=0; i<numberofNodes; i++) {
                if (mask->Perm[i] < 0) nodeLayer[i] = 0;
            }
        }
        for (int i=0; i<numberofNodes; i++) {
            if (ifMask == YES) {
                if (mask->Perm[i] < 0) continue;
            }
            rounds = 1;
            j = bottomNodePointer[i];
            nodeLayer[j] = rounds;
            while (1) {
                k = upNodePointer[j];
                if (k == j) break;
                rounds++;
                j = k;
                nodeLayer[j] = rounds;
            }
        }
        fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: layer range: ['%d', '%d'].\n", min_array(nodeLayer, numberofNodes), max_array(nodeLayer, numberofNodes));
    }
    
    if (isMidNode == YES) {
        FEMListUtilities *listUtilities = [FEMListUtilities sharedListUtilities];
        memset( midNodePointer, -1, numberofNodes*sizeof(int) );
        *midLayerExists = NO;
        for (int t=mesh.numberOfBulkElements; t<mesh.numberOfBulkElements+mesh.numberOfBoundaryElements; t++) {
            for (FEMBoundaryCondition *boundary in model.boundaryConditions) {
                if (meshElements[t].BoundaryInfo->Constraint == boundary.tag) {
                    if ([listUtilities listCheckPresentVariable:@"mid surface" inArray:boundary.valuesList] == YES) {
                        for (int i=0; i<meshElements[t].sizeNodeIndexes; i++) {
                            midNodePointer[meshElements[t].NodeIndexes[i]] = meshElements[t].NodeIndexes[i];
                            *midLayerExists = YES;
                        }
                    }
                    break;
                }
            }
        }
        
        if (*midLayerExists == YES) {
            fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: determine mid pointers.\n");
            for (rounds=1; rounds<=numberofNodes; rounds++) {
                downHit = 0;
                upHit = 0;
                for (int i=0; i<numberofNodes; i++) {
                    if (ifMask == YES) {
                        if (mask->Perm[i] < 0) continue;
                    }
                    if (midNodePointer[i] < 0) continue;
                    if (upActive == YES) {
                        j = upNodePointer[i];
                        if (midNodePointer[j] < 0) {
                            upHit++;
                            midNodePointer[j] = midNodePointer[i];
                        }
                    }
                    if (downActive == YES) {
                        j = downNodePointer[i];
                        if (midNodePointer[j] < 0) {
                            downHit++;
                            midNodePointer[j] = midNodePointer[i];
                        }
                    }
                }
                if (upHit == 0 && downHit == 0) break;
            }
            fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: mid layer structure detected in %d cycles. \n", rounds-1);
        }
    }
    
    // Count the numer of top and bottom nodes, for information only
    fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: counting top and bottom nodes.\n");
    int topNodes = 0, bottomNodes = 0;
    if (upActive == YES) {
        double minTop = HUGE_VAL;
        double maxTop = -HUGE_VAL;
        for (int i=0; i<numberofNodes; i++) {
            if (topNodePointer[i] == i) {
                minTop = min(minTop, variableContainers->Values[i]);
                maxTop = max(maxTop, variableContainers->Values[i]);
                topNodes++;
            }
        }
        fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: top range: %f %f.\n", minTop, maxTop);
    }
    if (downActive == YES) {
        double minBot = HUGE_VAL;
        double maxBot = -HUGE_VAL;
        for (int i=0; i<numberofNodes; i++) {
            if (bottomNodePointer[i] == i) {
                minBot = min(minBot, variableContainers->Values[i]);
                maxBot = max(maxBot, variableContainers->Values[i]);
                bottomNodes++;
            }
        }
        fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: bottom range: %f %f.\n", minBot, maxBot);
    }
    
    at1 = cputime();
    fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: top and bottom pointer init time: %f.\n", at1-at0);
    fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: top and bottom pointer init rounds: %d.\n", rounds);
    if (upActive == YES) fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: number of nodes at the top: %d.\n", topNodes);
    if (downActive == YES) fprintf(stdout, "FEMMeshUtils:detectExtrudedStructureMesh: number of nodes at the bottom: %d.\n.", bottomNodes);
        
    return variable;
}

/***************************************************************************************
 
    Given two interface meshes for nonconforming rotating boundaries, make
    a coordinate transformation to each node of the slave boundary (BMesh1) so that
    they hit the master boundary (BMesh2). In case of anti-periodic projector
    mark the nodes that need an odd number of periods.
 
    Method corresponds to Elmer from git on October 27 2015
 
***************************************************************************************/
-(void)preRotationalProjectorMesh1:(FEMMesh * __nonnull)bMesh1 mesh2:(FEMMesh * __nonnull)bMesh2 mirrorNode:(BOOL * __nullable)mirrorNode sizeMirrorNode:(int * __nullable)sizeMirrorNode {
    
    int nFii;
    double fii;
    
    BOOL antiPeriodic = (mirrorNode != NULL) ? YES : NO;
    
    Nodes_t *bMesh1Nodes = bMesh1.getNodes;
    Nodes_t *bMesh2Nodes = bMesh2.getNodes;
    
    double f2Min, f2Max;
    vDSP_minvD(bMesh2Nodes->x, 1, &f2Min, bMesh2.numberOfNodes);
    vDSP_maxvD(bMesh2Nodes->x, 1, &f2Max, bMesh2.numberOfNodes);
    
    double dFii2 = f2Max - f2Min;
    int sectorMax = ceil(360.0 / dFii2);
    
    fprintf(stdout, "FEMMeshUtils:preRotationalProjectorMesh1: maximum number of sectors: %d.\n", sectorMax);
    
    int *sectorCount = intvec(-sectorMax, sectorMax);
    for (int i=-sectorMax; i<=sectorMax; i++) {
        sectorCount[i] = 0;
    }
    
    for (int i=0; i<bMesh1.numberOfNodes; i++) {
        fii = bMesh1Nodes->x[i];
        nFii = floor((fii - f2Min) / dFii2);
        bMesh1Nodes->x[i] = bMesh1Nodes->x[i] - nFii * dFii2;
        sectorCount[nFii] = sectorCount[nFii] + 1;
        if (antiPeriodic == YES) {
            if (nFii % 2 != 0) {
                mirrorNode[i] = YES;
            }
        }
    }
    
    if (sectorCount[0] < bMesh1.numberOfNodes) {
        fprintf(stdout, "FEMMeshUtils:preRotationalProjectorMesh1: number of nodes by sectors: \n");
        for (int i=-sectorMax; i<=sectorMax; i++) {
            if (sectorCount[i] > 0) {
                fprintf(stdout, "FEMMeshUtils:preRotationalProjectorMesh1: sector: %d, nodes: %d\n", i, sectorCount[i]);
            }
        }
        if (antiPeriodic == YES) {
            int count = 0;
            for (int i=0; i<*sizeMirrorNode; i++) {
                if (mirrorNode[i] == YES) count++;
            }
            fprintf(stdout, "FEMMeshUtils:preRotationalProjectorMesh1: number of mirror nodes: %d.\n", count);
        }
    } else {
        fprintf(stdout, "FEMMeshUtils:preRotationalProjectorMesh1: no nodes needed mapping.\n");
    }
    free_ivector(sectorCount, -sectorMax, sectorMax);
}

/****************************************************************************
 
    Postprocess projector so that it changes the sign of the anti-periodic
    entries as assigned by the MirrorNode flag.
 
    Method corresponds to Elmer from git on October 27 2015
 
****************************************************************************/
-(void)postRotationalProjector:(FEMMatrix * __nonnull)projector mirrorNode:(BOOL * __nonnull)mirrorNode sizeMirrorNode:(int * __nonnull)sizeMirrorNode {
    
    if (mirrorNode == NULL) return;
    
    int count = 0;
    for (int i=0; i<*sizeMirrorNode; i++) {
        if (mirrorNode[i] == YES) count++;
    }
    if (count == 0) return;
    
    int n = projector.numberOfRows;
    matrixArraysContainer *matrixContainers = projector.getContainers;
    int *rows = matrixContainers->Rows;
    double *values = matrixContainers->Values;
    
    for (int i=0; i<n; i++) {
        if (mirrorNode[i] == YES) {
            for (int j=rows[i]; j<=rows[i+1]-1; j++) {
                values[j] = -values[j];
            }
        }
    }
}

/**************************************************************************

    Save projector, mainly a utility for debugging purposes
 
    Method corresponds to Elmer from git on October 27 2015
 
**************************************************************************/
-(void)saveProjector:(FEMMatrix * __nonnull)projector saveRowSum:(BOOL)saveRowSum prefix:(NSString * __nonnull)prefix invPerm:(int * __nullable)invPerm {
    
    double dia, rowsum;
    
    if (projector == nil) return;
    
    matrixArraysContainer *projectorContainers = projector.getContainers;
    
    int *intInvPerm = NULL;
    if (invPerm != NULL) {
        intInvPerm = invPerm;
    } else {
        intInvPerm = projectorContainers->InvPerm;
    }
    
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *fileName = [[fileManager currentDirectoryPath] stringByAppendingPathComponent:[prefix stringByAppendingPathExtension:@"dat"]];
    // TODO: Add support for parallel run
    
    FILE *f1 = fopen([fileName UTF8String], "w");
    for (int i=0; i<projector.numberOfRows; i++) {
        if (intInvPerm[i] < 0) continue;
        rowsum = 0.0;
        for (int j=projectorContainers->Rows[i]; j<=projectorContainers->Rows[i+1]-1; j++) {
            fprintf(f1, "%d %d %lf\n", intInvPerm[i], projectorContainers->Cols[j], projectorContainers->Values[j]);
            if (intInvPerm[i] == projectorContainers->Cols[j]) {
                dia = projectorContainers->Values[j];
            } else {
                rowsum = rowsum + projectorContainers->Values[j];
            }
        }
    }
    fclose(f1);
    
    if (saveRowSum == YES) {
        NSString *fileName = [[fileManager currentDirectoryPath] stringByAppendingPathComponent:[prefix stringByAppendingString:@"_rsum.dat"]];
        // TODO: Add support for parallel run
        f1 = fopen([fileName UTF8String], "w");
        for (int i=0; i<projector.numberOfRows; i++) {
            if (intInvPerm[i] < 0) continue;
            rowsum = 0.0;
            dia = 0.0;
            for (int j=projectorContainers->Rows[i]; j<=projectorContainers->Rows[i+1]-1; j++) {
                if (intInvPerm[i] == projectorContainers->Cols[j]) {
                    dia = projectorContainers->Values[j];
                }
                rowsum = rowsum + projectorContainers->Values[j];
            }
            
            fprintf(f1, "%d %d %lf %lf\n", intInvPerm[i], projectorContainers->Rows[i+1]-projectorContainers->Rows[i], dia, rowsum);
        }
        fclose(f1);
    }
}

/************************************************************************************************
 
        Color a mesh using the First Fit (FF) algorithm (also called greedy algorithm)
 
************************************************************************************************/
-(void)colorMesh:(FEMMesh * __nonnull)mesh {
    
    int i, j, numberOfColoredElements = 0, numberOfSameColors, numberOfElementsWithCurrentColor;
    BOOL isShared;
    RGBColors currentColor;
    
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    Element_t *elements = mesh.getElements;
    
    double at = cputime();
    
    // Color the elements
    while (numberOfColoredElements != mesh.numberOfBulkElements) {
        [utilities generateColor:&currentColor];
        numberOfElementsWithCurrentColor = 0;
        
        for (i=0; i<mesh.numberOfBulkElements; i++) {
            if (!elements[i].colored) {
                // Check if neighbors do not have current color
                numberOfSameColors = 0;
                for (j=0; j<mesh.numberOfBulkElements; j++) {
                    if (j == i) continue;
                    isShared = NO;
                    for (int k=0; k<elements[j].Type.NumberOfNodes; k++) {
                        for (int l=0; l<elements[i].Type.NumberOfNodes; l++) {
                            if (elements[j].NodeIndexes[k] == elements[i].NodeIndexes[l]) {
                                isShared = YES;
                                break;
                            }
                        }
                        if (isShared == YES) break;
                    }
                    
                    if (isShared == YES) {
                        if (elements[j].color.red == currentColor.red && elements[j].color.green == currentColor.green && elements[j].color.blue == currentColor.blue) numberOfSameColors++;
                    }
                }
                if (numberOfSameColors == 0) { // No neighbors with current color so color the element with current color
                    elements[i].color.red = currentColor.red;
                    elements[i].color.green = currentColor.green;
                    elements[i].color.blue = currentColor.blue;
                    elements[i].color.colorIndex = mesh.numberOfColors;
                    elements[i].colored = true;
                    numberOfColoredElements++;
                    numberOfElementsWithCurrentColor++;
                }
            }
        }
        NSMutableArray *color = [[NSMutableArray alloc] initWithObjects:@(numberOfElementsWithCurrentColor), @(mesh.numberOfColors), @(currentColor.red), @(currentColor.green), @(currentColor.blue), nil];
        [mesh.colors addObject:color];
        mesh.numberOfColors++;
    }
    
    fprintf(stdout, "FEMMeshUtils:colorMesh:Timing (s): %f\n", cputime() - at);
    
    fprintf(stdout, "FEMMeshUtils:FEMMesh_colorMesh: number of colors: %d\n", mesh.numberOfColors);
    fprintf(stdout, "FEMMeshUtils:FEMMesh_colorMesh: number of elements for each color set before optimization: \n");
    i = 1;
    for (NSMutableArray *color in mesh.colors) {
        fprintf(stdout, "color %d : number of elements: %d\n", i, [color[0] intValue]);
        i++;
    }
    
    // The FF algorithm produces an unbalanced distribution of elements per color set since the initial sets
    // are the first to get the elements. The following optimization redistributes elements from "rish" sets to "poor" sets
    int meanValue = mesh.numberOfBulkElements/mesh.numberOfColors;
    int kk = mesh.numberOfColors-1; // Start with the poorest set and then up to the richest set
    for (NSMutableArray *color in mesh.colors) {
        if ([color[0] intValue] > meanValue) {
            int k = [color[1] intValue];
            while ([color[0] intValue] > meanValue) {
                for (i=0; i<mesh.numberOfBulkElements; i++) {
                    if (elements[i].color.colorIndex == k) {
                        numberOfSameColors = 0;
                        for (j=0; j<mesh.numberOfBulkElements; j++) {
                            if (j == i) continue;
                            if (elements[j].color.colorIndex == kk) {
                                isShared = NO;
                                for (int l=0; l<elements[j].Type.NumberOfNodes; l++) {
                                    for (int ll=0; ll<elements[i].Type.NumberOfNodes; ll++) {
                                        if (elements[j].NodeIndexes[l] == elements[i].NodeIndexes[ll]) {
                                            isShared = YES;
                                            break;
                                        }
                                    }
                                    if (isShared == YES) break;
                                }
                                if (isShared == YES) numberOfSameColors++;
                            }
                        }
                        NSMutableArray *poorColor = mesh.colors[kk];
                        if (numberOfSameColors == 0) {
                            elements[i].color.red = [poorColor[2] doubleValue];
                            elements[i].color.green = [poorColor[3] doubleValue];
                            elements[i].color.blue = [poorColor[4] doubleValue];
                            elements[i].color.colorIndex = [poorColor[1] intValue];
                            color[0] = @([color[0] intValue]-1);
                            poorColor[0] = @([poorColor[0] intValue]+1);
                        }
                        if ([poorColor[0] intValue] >= (meanValue-DBL_EPSILON) && [poorColor[0] intValue] <= (meanValue+DBL_EPSILON)) {
                            kk--;
                            break;
                        }
                        if ([color[0] intValue] <= meanValue) break;
                    }
                }
                if (i == mesh.numberOfBulkElements) break; // Already done all elements for given color set so go to next set
            }
        }
    }
    
    fprintf(stdout, "FEMMeshUtils:colorMesh: number of elements for each color set after optimization: \n");
    i = 1;
    for (NSMutableArray *color in mesh.colors) {
        fprintf(stdout, "color %d : number of elements: %d\n", i, [color[0] intValue]);
        i++;
    }
    
    // Build the element color mapping
    int *colorMapping = intvec(0, mesh.numberOfBulkElements-1);
    int indx = 0;
    for (NSMutableArray *color in mesh.colors) {
        for (i=0; i<mesh.numberOfBulkElements; i++) {
            if (elements[i].color.colorIndex == [color[1] intValue]) {
                colorMapping[indx] = elements[i].ElementIndex-1;
                indx++;
            }
        }
    }
    [mesh assignColorMapping:colorMapping];
    
    // Make a last check to be sure...
    //    for (NSMutableArray *color in self.colors) {
    //         for (i=0; i<self.numberOfBulkElements; i++) {
    //             if (_elements[i].color.colorIndex == [color[1] intValue]) {
    //                 for (j=0; j<self.numberOfBulkElements; j++) {
    //                     if (j == i) continue;
    //                     isShared = NO;
    //                     for (int k=0; k<_elements[j].Type.NumberOfNodes; k++) {
    //                         for (int l=0; l<_elements[i].Type.NumberOfNodes; l++) {
    //                             if (_elements[j].NodeIndexes[k] == _elements[i].NodeIndexes[l]) {
    //                                 isShared = YES;
    //                                 break;
    //                             }
    //                         }
    //                         if (isShared == YES) break;
    //                     }
    //                     if (isShared == YES) {
    //                         if (_elements[i].color.colorIndex == _elements[j].color.colorIndex) {
    //                             NSLog(@"FEMMesh:FEMMesh_colorMesh: elements %d and %d have similar color but share nodes.\n", i, j);
    //                             errorfunct("FEMMesh:FEMMesh_colorMesh", "Programm terminating now...");
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //    }
}

-(void)saveColoredMesh:(FEMMesh * __nonnull)mesh meshdir:(NSString * __nonnull)dir meshName:(NSString * __nonnull)name elementsFileName:(NSString * __nonnull)elementsFileName saveAllElementData:(BOOL)saveAllElementData colorFileName:(NSString * __nonnull)colorFileName {
    
    NSFileHandle *fileHandle;
    
    FEMPost *post = [[FEMPost alloc] init];
    Element_t *elements = mesh.getElements;
    
    NSFileManager *fileManager = [NSFileManager defaultManager];
    NSString *elementFile = [[dir stringByAppendingPathComponent:name] stringByAppendingPathComponent:elementsFileName];
    if ([fileManager createFileAtPath:elementFile contents:nil attributes:nil] == YES) {
        fileHandle = [NSFileHandle fileHandleForWritingAtPath:elementFile];
    } else {
        fprintf(stderr, "FEMMeshUtils:saveColoredMesh: can't create the file %s.", [elementFile UTF8String]);
        fatal("FEMMeshUtils:saveColoredMesh");
    }
    
    char space = ' ';
    char newLine = '\n';
    NSData *spaceBuff = [NSMutableData dataWithBytes:&space length:sizeof(space)];
    NSData *newLineBuff = [NSMutableData dataWithBytes:&newLine length:sizeof(newLine)];
    
    fprintf(stdout, "FEMMeshUtils:saveColoredMesh: write file with colored elements: %s...\n", [elementFile UTF8String]);
    if (saveAllElementData == YES) {
        for (int i=0; i<mesh.numberOfBulkElements; i++) {
            [post writeInteger:elements[i].ElementIndex toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
            [post writeInteger:elements[i].BodyID toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
            [post writeInteger:elements[i].Type.ElementCode toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
            for (int j=0; j<elements[i].Type.NumberOfNodes; j++) {
                [post writeInteger:elements[i].NodeIndexes[j]+1 toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
            }
            [post writeInteger:elements[i].color.colorIndex+1 toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
            [fileHandle writeData:newLineBuff];
        }
    } else {
        for (int i=0; i<mesh.numberOfBulkElements; i++) {
            [post writeInteger:elements[i].ElementIndex toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
            [post writeInteger:elements[i].color.colorIndex+1 toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
            [fileHandle writeData:newLineBuff];
        }
    }
    [fileHandle closeFile];
    fprintf(stdout, "FEMMeshUtils:saveColoredMesh: done.\n");
    
    NSString *colorFile = [[dir stringByAppendingPathComponent:name] stringByAppendingPathComponent:colorFileName];
    if ([fileManager createFileAtPath:colorFile contents:nil attributes:nil] == YES) {
        fileHandle = [NSFileHandle fileHandleForWritingAtPath:colorFile];
    } else {
        fprintf(stderr, "FEMMeshUtils:saveColoredMesh: can't create the file %s.", [colorFile UTF8String]);
        fatal("FEMMeshUtils:saveColoredMesh");
    }
    
    fprintf(stdout, "FEMMeshUtils:saveColoredMesh: write color file: %s...\n", [colorFile UTF8String]);
    [post writeInteger:mesh.numberOfColors toFileHandle:fileHandle]; [fileHandle writeData:newLineBuff];
    for (NSMutableArray *color in mesh.colors) {
        [post writeInteger:[color[0] intValue] toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
        [post writeInteger:[color[1] intValue]+1 toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
        [post writeDouble:[color[2] doubleValue] toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
        [post writeDouble:[color[3] doubleValue] toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
        [post writeDouble:[color[4] doubleValue] toFileHandle:fileHandle]; [fileHandle writeData:spaceBuff];
        [fileHandle writeData:newLineBuff];
    }
    [fileHandle closeFile];
    fprintf(stdout, "FEMMeshUtils:saveColoredMesh: done.\n");
}

-(void)readColoredMesh:(FEMMesh * __nonnull)mesh name:(NSString * __nonnull)name directory:(NSString * __nonnull)dir readElementsFromFile:(BOOL)readElementsFromFile {
    
    int isLineBreak, lineCount = 0, nb;
    NSString *line;
    
    NSString *colorFile = [[dir stringByAppendingPathComponent:name] stringByAppendingPathComponent:@"mesh.colors"];
    
    FileReader * reader = [[FileReader alloc] initWithFilePath:colorFile];
    if (!reader) {
        fprintf(stderr, "FEMMeshUtils:readColoredMesh: file %s not found in mesh directory.\n", [colorFile UTF8String]);
        fatal("FEMMeshUtils:readColoredMesh");
    }
    
    // Used to separate strings and filter them from white spaces
    NSCharacterSet *whitespaces = [NSCharacterSet whitespaceCharacterSet];
    NSPredicate *noEmptyStrings = [NSPredicate predicateWithFormat:@"SELF != ''"];
    
    line = nil;
    while ((line = [reader readLine])) {
        lineCount++;
        fprintf(stdout, "FEMMeshUtils:readColoredMesh: %3.d: %s.\n", lineCount, [line UTF8String]);
        // Parse the line
        NSArray *stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
        NSArray *filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
        isLineBreak = 0;
        for (NSString *string in filteredArray) {
            if ([string isEqualToString:@"\n"] == YES) isLineBreak++;
        }
        if (lineCount == 1) {
            if (([filteredArray count]-isLineBreak) < 1 || ([filteredArray count]-isLineBreak) > 1) {
                fprintf(stderr, "FEMMeshUtils:readColoredMesh: not properlly formatted data in file %s at line %d.\n", [colorFile UTF8String], lineCount);
                fatal("FEMMeshUtils:readColoredMesh");
            }
        } else {
            if (([filteredArray count]-isLineBreak) < 5 || ([filteredArray count]-isLineBreak) > 5) {
                fprintf(stderr, "FEMMeshUtils:readColoredMesh: not properlly formatted data in file %s at line %d.\n", [colorFile UTF8String], lineCount);
                fatal("FEMMeshUtils:readColoredMesh");
            }
        }
        if (([filteredArray count]-isLineBreak) == 1) { // First line of color file
            mesh.numberOfColors = [filteredArray[0] intValue];
        } else if (([filteredArray count]-isLineBreak) == 5) {
            NSMutableArray *color = [[NSMutableArray alloc] initWithObjects:filteredArray[0], @([filteredArray[1] intValue]-1), filteredArray[2], filteredArray[3], filteredArray[4], nil];
            [mesh.colors addObject:color];
        }
    }
    
    [reader closeHandle];
    
    // Build the element color mapping
    Element_t *elements = NULL;
    if (readElementsFromFile == YES) {
        elements = (Element_t*) malloc( sizeof(Element_t) * mesh.numberOfBulkElements );
        initElements(elements, mesh.numberOfBulkElements );
        
        NSString *elementsFile = [[dir stringByAppendingPathComponent:name] stringByAppendingPathComponent:@"mesh.colored_elements"];
        
        FileReader * reader = [[FileReader alloc] initWithFilePath:colorFile];
        if (!reader) {
            fprintf(stderr, "FEMMeshUtils:readColoredMesh: file %s not found in mesh directory.\n", [elementsFile UTF8String]);
            fatal("FEMMeshUtils:readColoredMesh");
        }
        line = nil;
        lineCount = 0;
        nb = 0;
        while ((line = [reader readLine])) {
            lineCount++;
            fprintf(stdout, "FEMMeshUtils:readColoredMesh: %3.d: %s.\n", lineCount, [line UTF8String]);
            // Parse the line
            NSArray *stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
            NSArray *filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
            isLineBreak = 0;
            for (NSString *string in filteredArray) {
                if ([string isEqualToString:@"\n"] == YES) isLineBreak++;
            }
            if (([filteredArray count]-isLineBreak) < 2 || ([filteredArray count]-isLineBreak) > 2) {
                fprintf(stderr, "FEMMeshUtils:readColoredMesh: not properlly formatted data in file %s at line %d.\n", [elementsFile UTF8String], lineCount);
                fatal("FEMMeshUtils:readColoredMesh");
            }
            elements[nb].ElementIndex = [filteredArray[0] intValue];
            elements[nb].color.colorIndex = [filteredArray[1] intValue];
            nb++;
        }
        [reader closeHandle];
        
    } else elements = mesh.getElements;
    
    int indx;
    int *meshColorMapping = mesh.getColorMapping;
    if (meshColorMapping == NULL) {
        int *colorMapping = intvec(0, mesh.numberOfBulkElements-1);
        indx = 0;
        for (NSMutableArray *color in mesh.colors) {
            for (int i=0; i<mesh.numberOfBulkElements; i++) {
                if (elements[i].color.colorIndex-1 == [color[1] intValue]) {
                    colorMapping[indx] = elements[i].ElementIndex-1;
                    indx++;
                }
            }
        }
        [mesh assignColorMapping:colorMapping];
    }
    
    // Build the element permutation store. This is for the simplest form of getElementDofsSolution:model:forElement:atIndexes:disableDiscontinuousGalerkin:
    int *elementNodeIndexesStore = intvec(0, (mesh.numberOfBulkElements*mesh.maxElementDofs)-1);
    indx = 0;
    for (int t=0; t<mesh.numberOfBulkElements; t++) {
        for (int i=0; i<elements[t].NDOFs; i++) {
            elementNodeIndexesStore[indx] = elements[t].NodeIndexes[i];
            indx++;
        }
    }
    [mesh assignElementNodeIndexesStore:elementNodeIndexesStore];
    
    if (readElementsFromFile == YES) {
        free(elements);
    }
}

@end

