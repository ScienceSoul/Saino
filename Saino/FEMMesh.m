//
//  FEMMesh.m
//  Saino
//
//  Created by Seddik hakime on 15/04/11.
//  Copyright 2011 ScienceSoul. All rights reserved.
//

#import "FEMMesh.h"
#import "FEMSolution.h"
#import "FEMListUtilities.h"
#import "FEMElementDescription.h"
#import "FEMMeshUtils.h"
#import "FEMPElementMaps.h"
#import "FEMBoundaryCondition.h"
#import "FEMProjector.h"
#import "Utils.h"

@interface FEMMesh ()
-(void)FEMMesh_getMaxdefs:(FEMModel *)model element:(Element_t *)element elementDef:(NSString *)elementDef solverID:(int)solverID bodyID:(int)bodyID defDofs:(int *)defDofs;
-(void)FEMMesh_convertToACTetra:(Element_t *)tetra;
-(void)FEMMesh_deallocateQuadrantTree:(Quadrant_t *)root;
@end

@implementation FEMMesh

#pragma mark Private methods

-(void)FEMMesh_getMaxdefs:(FEMModel *)model element:(Element_t *)element elementDef:(NSString *)elementDef solverID:(int)solverID bodyID:(int)bodyID defDofs:(int *)defDofs {
    
    int i, j, l, n;
    double x, y, z, sum;
    FEMSolution *solution;
    solutionArraysContainer *solContainers = NULL;
    NSRange substr;
    
    solution = (model.solutions)[solverID];
    solContainers = solution.getContainers;
    
    if (solContainers->defDofs == NULL) {
        solContainers->defDofs = intmatrix(0, model.numberOfBodies-1, 0, 5);
        solContainers->size1DefDofs = model.numberOfBodies;
        solContainers->size2DefDofs = 6;
        for (i=0; i<model.numberOfBodies; i++) {
            for (j=0; j<6; j++) {
                solContainers->defDofs[i][j] = -1;
            }
        }
        for (i=0; i<model.numberOfBodies; i++) {
            solContainers->defDofs[i][0] = 1;
        }
    }
    
    substr = [elementDef rangeOfString:@"n:"];
    if (substr.location != NSNotFound) {
        l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr.location+2)]] intValue];
        solContainers->defDofs[bodyID][0] = l;
        defDofs[0] = max(defDofs[0], l);
    }
    
    substr = [elementDef rangeOfString:@"e:"];
    if (substr.location != NSNotFound) {
        l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr.location+2)]] intValue];
        solContainers->defDofs[bodyID][1] = l;
        defDofs[1] = max(defDofs[1], l);
    }
    
    substr = [elementDef rangeOfString:@"f:"];
    if (substr.location != NSNotFound) {
        l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr.location+2)]] intValue];
        solContainers->defDofs[bodyID][2] = l;
        defDofs[2] = max(defDofs[2], l);
    }
    
    substr = [elementDef rangeOfString:@"d:"];
    if (substr.location != NSNotFound) {
        l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr.location+2)]] intValue];
        solContainers->defDofs[bodyID][3] = l;
        defDofs[3] = max(defDofs[3], l);
    } else {
        if ([[solution.solutionInfo objectForKey:@"discountinuous galerkin"] boolValue] == YES) {
            solContainers->defDofs[bodyID][3] = l;
            defDofs[3] = max(defDofs[3], 0);
        }
    }
    
    substr = [elementDef rangeOfString:@"b:"];
    if (substr.location != NSNotFound) {
        l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr.location+2)]] intValue];
        solContainers->defDofs[bodyID][4] = l;
        defDofs[4] = max(defDofs[4], l);
    }
    
    substr = [elementDef rangeOfString:@"p:"];
    if (substr.location != NSNotFound) {
        if ([[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr.location+2)]] isEqualToString:@"%"]) {
            n = element->Type.NumberOfNodes;
            sum = 0.0;
            for (i=0; i<element->Type.NumberOfNodes; i++) {
                sum = sum + _globalNodes->x[element->NodeIndexes[i]];
            }
            x = sum / n;
            
            sum = 0.0;
            for (i=0; i<element->Type.NumberOfNodes; i++) {
                sum = sum + _globalNodes->y[element->NodeIndexes[i]];
            }
            y = sum / n;
            
            sum = 0.0;
            for (i=0; i<element->Type.NumberOfNodes; i++) {
                sum = sum + _globalNodes->z[element->NodeIndexes[i]];
            }
            z = sum / n;
            // TODO: The rest depends on a call to MATC. So we need to decide what to do here.
            
        } else {
            l = [[NSString stringWithFormat:@"%c", [elementDef characterAtIndex:(substr.location+2)]] intValue];
            solContainers->defDofs[bodyID][5] = l;
            defDofs[5] = max(defDofs[5], l);
        }
    }
}

/***********************************************************************
 
    Convert tetrahedral element to Ainsworth & Coyle type tetrahedron.

***********************************************************************/
-(void)FEMMesh_convertToACTetra:(Element_t *)tetra {
    
    int i, globalMin, globalMax, globalMinI;
    int *face, *globalFace;
    FEMPElementMaps *pMaps;
    
    face = intvec(0, 2);
    globalFace = intvec(0, 2);
    
    pMaps = [[FEMPElementMaps alloc] init];
    
    // Check if we can proceed any further
    if (tetra->Type.ElementCode != 504 || tetra->Pdefs == NULL) {
        warnfunct("FEMMesh_convertToACTetra", "Element to convert not a tetrahedron!");
        return;
    }
    
    // Find global min and max vertices
    globalMin = tetra->NodeIndexes[0];
    globalMinI = 0;
    globalMax = tetra->NodeIndexes[0];
    for (i=1; i<4; i++) {
        // Find min
        if (globalMin > tetra->NodeIndexes[i]) {
            globalMax = tetra->NodeIndexes[i];
            globalMinI = i;
        } else if (globalMax < tetra->NodeIndexes[i]) {
            globalMax = tetra->NodeIndexes[i];
        }
    }
    
    // Get face containing global min (either face 0 or 1)
    if (globalMinI == 3) {
        [pMaps getTetraFaceMap:face index:1 type:NULL];
    } else {
        [pMaps getTetraFaceMap:face index:0 type:NULL];
    }
    for (i=0; i<3; i++) {
        globalFace[i] = tetra->NodeIndexes[face[i]];
    }
    
    // Rotate face until first local index is min global
    while (1) {
        if (globalMin == globalFace[0]) break;
        cshift(globalFace, 3, 1);
    }
    // Assign new local numbering
    for (i=0; i<3; i++) {
        tetra->NodeIndexes[face[i]] = globalFace[i];
    }
    
    // Face 2 now contains global max
    [pMaps getTetraFaceMap:face index:2 type:NULL];
    for (i=0; i<3; i++) {
        globalFace[i] = tetra->NodeIndexes[face[i]];
    }
    // Rotate face until last local index is max global
    while (1) {
        if (globalMax == globalFace[2]) break;
        cshift(globalFace, 3, 1);
    }
    // Assign new local numbering
    for (i=0; i<3; i++) {
        tetra->NodeIndexes[face[i]] = globalFace[i];
    }
    
    // Set AC tetra type
    if (tetra->NodeIndexes[1] < tetra->NodeIndexes[2]) {
        tetra->Pdefs->TetraType = 1;
    } else if (tetra->NodeIndexes[2] < tetra->NodeIndexes[1]) {
        tetra->Pdefs->TetraType = 2;
    } else {
        errorfunct("FEMMesh_convertToACTetra", "Corrupted element type.");
    }
    
    free_ivector(face, 0, 2);
    free_ivector(globalFace, 0, 2);
    [pMaps deallocation];
}

-(void)FEMMesh_deallocateQuadrantTree:(Quadrant_t *)root {
    
    int i;
    
    if (root == NULL) return;
    
    if (root->elements != NULL) free_ivector(root->elements, 0, root->nElementsInQuadrant-1);
    root->elements = NULL;
    
    if (root->childQuadrants != NULL) {
        for (i=0; i<root->numberOfchildQuadrants; i++) {
            [self FEMMesh_deallocateQuadrantTree:root->childQuadrants[i].quadrant];
        }
        free(root->childQuadrants);
        root->childQuadrants = NULL;
    }
    free(root);
}

#pragma mark Public methods

@synthesize dimension = _dimension;
@synthesize numberOfNodes = _numberOfNodes;
@synthesize numberOfElements = _numberOfElements;
@synthesize numberOfBulkElements = _numberOfBulkElements;
@synthesize numberOfEdges = _numberOfEdges;
@synthesize numberOfFaces = _numberOfFaces;
@synthesize numberOfBoundaryElements = _numberOfBoundaryElements;
@synthesize numberOfViewFactors = _numberOfViewFactors;
@synthesize maxElementNodes = _maxElementNodes;
@synthesize maxElementDofs = _maxElementDofs;
@synthesize maxEdgeDofs = _maxEdgeDofs;
@synthesize maxFaceDofs = _maxFaceDofs;
@synthesize maxBdofs = _maxBdofs;
@synthesize numberOfPassiveBCs = _numberOfPassiveBCs;
@synthesize savesDone = _savesDone;
@synthesize outputActive = _outputActive;
@synthesize adaptiveMesh = _adaptiveMesh;
@synthesize changed = _changed;
@synthesize stabilize = _stabilize;
@synthesize name = _name;
@synthesize variables = _variables;
@synthesize projectors = _projectors;
@synthesize next = _next;
@synthesize parent = _parent;
@synthesize child = _child;

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        _savesDone = NO;
        _outputActive = NO;
        _adaptiveMesh = NO;
        _changed = NO;
        _stabilize = NO;
        
        _dimension = 0;
        _numberOfEdges = 0;
        _numberOfFaces = 0;
        _numberOfNodes = 0;
        _numberOfBulkElements = 0;
        _numberOfBoundaryElements = 0;
        
        _maxFaceDofs = 0;
        _maxEdgeDofs = 0;
        _maxBdofs = 0;
        _maxElementDofs = 0;
        _maxElementNodes = 0;
        _numberOfPassiveBCs = 0;
                
        _elements = NULL;
        _edges = NULL;
        _faces = NULL;
        _globalNodes = NULL;
        _rootQuadrant = NULL;
        _viewFactors = NULL;
        
        _variables = [[NSMutableArray alloc] init];
        _projectors = [[NSMutableArray alloc] init];
        
        _parent = nil;
        _child = nil;
    }
    
    return self;
}

-(void)allocatePDefinitionsForElement:(Element_t *)element {
    
    element->Pdefs = NULL;
    element->Pdefs = (PElementDefs_t*) malloc( sizeof(PElementDefs_t));
    if (element->Pdefs == NULL) {
        errorfunct("allocatePDefinitionsForElement", "Unable to allocate memory");
    }
    
    element->Pdefs->p = 0;
    element->Pdefs->TetraType = 0;
    element->Pdefs->isEdge = false;
    element->Pdefs->PyramidQuadEdge = false;
    element->Pdefs->LocalNumber = 0;
    element->Pdefs->GaussPoints = 0;
}

-(void)loadMeshForModel:(FEMModel *)model meshDirectory:(NSString *)dir meshName:(NSString *)name boundariesOnly:(BOOL)bd numberOfPartitions:(int *)numParts partitionID:(int *)partID definitions:(int *)defDofs {
    
    int i, j, k, n, body, type, bndry, left, right, tag;
    int addr1, addr2;
    int meshDim, saveDim;
    int minIndex, maxIndex, minEIndex, maxEIndex, dgIndex, bid, defaultTargetBC;
    int *countByType, *types, typeCount, *inDofs, *nodes;
    int *nodeTags, *localPerm, *localEPerm, *elementTags, *edgeDofs, *faceDofs;
    listBuffer coordMap = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer wrk = { NULL, NULL, NULL, NULL, 0, 0, 0};
    listBuffer bdList = { NULL, NULL, NULL, NULL, 0, 0, 0};
    double *cCoord, *coord, *coordScale;
    FEMMatrix *projector;
    FEMListUtilities *listUtil;
    SIOMeshIO *meshIO;
    FEMElementDescription *elmDescription;
    FEMPElementMaps *pMaps;
    FEMMeshUtils *meshUtils;
    ElementType_t *elmType;
    NSArray *bList;
    NSString *elementDef0, *elementDef;
    NSMutableString *str;
    BOOL parallel, found, gotIt, needEdges, any;
    
    parallel = NO;
    
    listUtil = [[FEMListUtilities alloc] init];
    
    countByType = intvec(0, 63);
    types = intvec(0, 63);
    
    if (numParts != NULL && partID != NULL) {
        if (*numParts > 1) {
            meshIO = [[SIOMeshIO alloc] initWithParallelNumberOfProcessors:*numParts processorID:*partID+1];
            parallel = YES;
        } else {
            meshIO = [[SIOMeshIO alloc] init];
        }
    } else {
        meshIO = [[SIOMeshIO alloc] init];
    }
    
    // Mesh
    [meshIO openMeshAtPath:name];
    if (meshIO.info != 0) {
        NSLog(@"Unable to load mesh: %@\n", name);
        errorfunct("loadMeshForModel", "Program terminating now...");
    }
    
    [meshIO getMeshDescriptionNodeCount:&_numberOfNodes elementCount:&_numberOfBulkElements boundaryElementCount:&_numberOfBoundaryElements usedElementTypes:&typeCount elementTypeTags:types elementCountByType:countByType];
    if (meshIO.info != 0) {
        NSLog(@"Unable to read mesh header for mesh: %@\n", name);
        errorfunct("loadMeshForModel", "Program terminating now...");
    }
    
    _globalNodes->numberOfNodes = self.numberOfNodes;
    
    if (bd == YES) self.numberOfBulkElements = 0;
    
    self.maxElementNodes = 0;
    self.maxElementDofs = 0;
    self.maxEdgeDofs = 0;
    self.maxFaceDofs = 0;
    self.maxBdofs = 0;
    for (i=0; i<typeCount; i++) {
        self.maxElementNodes = max( self.maxElementNodes, types[i]-100*(types[i]/100) );
    }
    
    _globalNodes = (Nodes_t*)malloc(sizeof(Nodes_t));
    if (_globalNodes == NULL) errorfunct("loadMeshForModel", "Failure to allocate nodes structure.");
    initNodes(_globalNodes);
    _globalNodes->x = doublevec(0, self.numberOfNodes-1);
    _globalNodes->y = doublevec(0, self.numberOfNodes-1);
    _globalNodes->z = doublevec(0, self.numberOfNodes-1);
    
    _elements = (Element_t*) malloc( sizeof(Element_t) * self.numberOfBulkElements+self.numberOfBoundaryElements );
    if (_elements == NULL) errorfunct("loadMeshForModel", "Failure to allocate elements structure.");
    initElements(_elements, self.numberOfBulkElements+self.numberOfBoundaryElements);
    self.numberOfElements = self.numberOfBulkElements+self.numberOfBoundaryElements;
    
    // Mesh nodes
    cCoord = doublevec(0, (3*self.numberOfNodes)-1);
    nodeTags = intvec(0, self.numberOfNodes-1);
    memset( cCoord, 0.0, (3*self.numberOfNodes)*sizeof(double) );
    memset( nodeTags, 0, self.numberOfNodes*sizeof(int) );
    
    [meshIO getMeshNodes:nodeTags coord:cCoord];
    
    found = [listUtil listGetIntegerArray:model inArray:model.simulation.valuesList forVariable:@"coordinate mapping" buffer:&coordMap];
    if (found == YES) {
        if (coordMap.m != 3) {
            warnfunct("loadMeshForModel", "Inconsistent coordinate mapping:");
            errorfunct("loadMeshForModel", "Coordinate mapping should be a permutation of 1, 2, 3.");
        }
        if (min_array(coordMap.ivector, 3) < 1 || max_array(coordMap.ivector, 3) > 3) {
            errorfunct("loadMeshForModel", "Coordinate mapping should be a permutation of 1, 2, 3.");
        }
        for (i=coordMap.ivector[0]; i<3*self.numberOfNodes; i+=3) {
            _globalNodes->x[i] = cCoord[i];
        }
        for (i=coordMap.ivector[1]; i<3*self.numberOfNodes; i+=3) {
            _globalNodes->y[i] = cCoord[i];
        }
        for (i=coordMap.ivector[2]; i<3*self.numberOfNodes; i+=3) {
            _globalNodes->z[i] = cCoord[i];
        }
    } else {
        for (i=0; i<3*self.numberOfNodes; i+=3) {
            _globalNodes->x[i] = cCoord[i];
        }
        for (i=1; i<3*self.numberOfNodes; i+=3) {
            _globalNodes->y[i] = cCoord[i];
        }
        for (i=2; i<3*self.numberOfNodes; i+=3) {
            _globalNodes->z[i] = cCoord[i];
        }
    }
    
    meshDim = 0;
    for (i=0; i<self.numberOfNodes; i++) {
        if (_globalNodes->x[i] != _globalNodes->x[0]) {
            meshDim++;
            break;
        }
    }
    for (i=0; i<self.numberOfNodes; i++) {
         if (_globalNodes->y[i] != _globalNodes->y[0]) {
             meshDim++;
             break;
         }
     }
     for (i=0; i<self.numberOfNodes; i++) {
         if (_globalNodes->z[i] != _globalNodes->z[0]) {
             meshDim++;
             break;
         }
     }
    
    saveDim = model.dimension;
    if (model.dimension <= 0) model.dimension = meshDim;
    
    self.dimension = model.dimension;
    
    free_dvector(cCoord, 0, (3*self.numberOfNodes)-1);
    if (coordMap.ivector != NULL) {
        free_ivector(coordMap.ivector, 0, coordMap.m-1);
        coordMap.ivector = NULL;
    }
    
    minIndex = min_array(nodeTags, self.numberOfNodes);
    maxIndex = max_array(nodeTags, self.numberOfNodes);
    
    localPerm = intvec(0, (maxIndex-minIndex+1)-1);
    memset( localPerm, 0, (maxIndex-minIndex+1)*sizeof(int) );
    for (i=0; i<self.numberOfNodes; i++) {
        localPerm[nodeTags[i]-minIndex] = i;
    }
    
    // Scaling of coordinates
    coordScale = doublevec(0, 2);

    found = [listUtil listGetConstRealArray:model inArray:model.simulation.valuesList forVariable:@"coordinate scaling" buffer:&wrk];
    if (found == YES) {
        for (i=0; i<3; i++) {
            coordScale[i] = 1.0;
        }
        for (i=0; i<meshDim; i++) {
            j = min(i, wrk.n);
            coordScale[i] = wrk.matrix[0][j];
        }
        NSLog(@"Scaling coordinates: %f %f %f\n", coordScale[0], coordScale[1], coordScale[2]);
        for (i=0; i<self.numberOfNodes; i++) {
            _globalNodes->x[i] = coordScale[0] * _globalNodes->x[i];
            if (meshDim > 1) _globalNodes->y[i] = coordScale[1] * _globalNodes->y[i];
            if (meshDim > 2) _globalNodes->z[i] = coordScale[2] * _globalNodes->z[i];
        }
    }
    
    free_dvector(coordScale, 0, 2);
    if (wrk.matrix != NULL) {
        free_dmatrix(wrk.matrix, 0, wrk.m-1, 0, wrk.n-1);
        wrk.matrix = NULL;
    }
    
    // Mesh elements
    elmDescription = [[FEMElementDescription alloc] init];
    elementTags = intvec(0, (self.numberOfBulkElements+1)-1);
    
    edgeDofs = NULL;
    edgeDofs = intvec(0, self.numberOfBulkElements-1);
    if (edgeDofs == NULL) errorfunct("loadMeshForModel", "Failure to allocate edge dofs.");
    
    faceDofs = NULL;
    faceDofs = intvec(0, self.numberOfBulkElements-1);
    if (faceDofs == NULL) errorfunct("loadMeshForModel", "Failure to allocate face dofs.");
    
    memset( elementTags, 0, (self.numberOfBulkElements+1)*sizeof(int) );
    dgIndex = 0;
    needEdges = NO;
    
    pMaps = [[FEMPElementMaps alloc] init];
    
    inDofs = intvec(0, 6);
    nodes = intvec(0, MAX_ELEMENT_NODES-1);
    for (i=0; i<self.numberOfBulkElements; i++) {
        memset( inDofs, 0, 7*sizeof(int) );
        [meshIO getMeshElementConnection:&elementTags[i] body:&body type:&type pdofs:inDofs nodes:nodes];
        if (meshIO.info != 0) break;
        
        if (defDofs != NULL) {
            for (j=0; j<6; j++) {
                if (inDofs[j] <= 0) inDofs[j] = defDofs[j];
            }
        }
        
        _elements[i].ElementIndex = i+1;
        //TODO: If parallel mesh, set the element partition index here
        _elements[i].PartIndex = *partID;
        
        
        elmType = NULL;
        _elements[i].BoundaryInfo = NULL;
        elmType = [elmDescription getElementType:type inMesh:self stabilization:NO];
        _elements[i].Type = *elmType;
        
        if (elmType != NULL) {
            n = _elements[i].Type.NumberOfNodes;
            _elements[i].NodeIndexes = intvec(0, n-1);
            _elements[i].sizeNodeIndexes = n;
            for (j=0; j<n; j++) {
                _elements[i].NodeIndexes[j] = localPerm[nodes[j]-minIndex];
            }
            
            _elements[i].BodyID = body;
            for (j=0; j<model.numberOfBodies; j++) {
                bList = [(model.bodies)[j] objectForKey:@"target bodies"];
                if (bList != nil) {
                    for (k=0; k<bList.count; k++) {
                        if (body == [bList[k] intValue]) _elements[i].BodyID = j+1;
                    }
                }
            }
            
            bid = _elements[i].BodyID;
            if ([(model.bodies)[bid-1] objectForKey:@"equation"] != nil) {
                j = [[(model.bodies)[bid-1] objectForKey:@"equation"] intValue];
                elementDef0 = [listUtil listGetString:model inArray:[(model.equations)[j-1] valuesList] forVariable:@"element" info:&found];
                k = 1;
                for (FEMSolution *solution in model.solutions) {
                    if (found == NO) elementDef0 = solution.solutionInfo[@"element"];
                    str = [NSMutableString stringWithString:@"element{"];
                    [str appendString:[[NSNumber numberWithInt:k] stringValue]];
                    [str appendString:@"}"];
                    elementDef = [listUtil listGetString:model inArray:[(model.equations)[j-1] valuesList] forVariable:str info:&gotIt];
                    if (gotIt == YES) {
                        [self FEMMesh_getMaxdefs:model element:&_elements[i] elementDef:elementDef solverID:k-1 bodyID:bid-1 defDofs:inDofs];
                    } else {
                        [self FEMMesh_getMaxdefs:model element:&_elements[i] elementDef:elementDef0 solverID:k-1 bodyID:bid-1 defDofs:inDofs];
                    }
                    k++;
                }
            }
            
            if (inDofs[0] != 0) {
                _elements[i].NDOFs = n;
            } else {
                _elements[i].NDOFs = 0;
            }
            
            edgeDofs[i] = max(0, inDofs[1]);
            faceDofs[i] = max(0, inDofs[2]);
            
            if (defDofs != NULL) {
                if (defDofs[3] == 0) inDofs[3] = n;
            }
            
            _elements[i].DGIndexes = NULL;
            if (inDofs[3] > 0) {
                _elements[i].DGIndexes = intvec(0, inDofs[3]-1);
                _elements[i].sizeDGIndexes = inDofs[3];
                for (j=0; j<inDofs[3]; j++) {
                    _elements[i].DGIndexes[j] = dgIndex;
                    dgIndex++;
                }
            } else {
                _elements[i].DGIndexes = NULL;
            }
            _elements[i].DGDOFs = max(0, inDofs[3]);
            any = NO;
            for (j=1; j<4; j++) {
                if (inDofs[j] > 0) {
                    any = YES;
                    break;
                }
            }
            needEdges = (needEdges == YES || any == YES) ? YES : NO;
            
            _elements[i].EdgeIndexes = NULL;
            _elements[i].FaceIndexes = NULL;
            _elements[i].BubbleIndexes = NULL;
            
            // Check if given element is a p element
            if (inDofs[5] > 0) {
                [self allocatePDefinitionsForElement:&_elements[i]];
                
                needEdges = YES;
                
                // Calculate element bubble dofs and set element p
                _elements[i].Pdefs->p = inDofs[5];
                if (inDofs[4] > 0) {
                    _elements[i].BDOFs = inDofs[4];
                } else {
                    _elements[i].BDOFs = [pMaps bubbleDofsForElement:&_elements[i] degree:_elements[i].Pdefs->p];
                }
                
                // All elements in actual mesh are not edges
                _elements[i].Pdefs->PyramidQuadEdge = false;
                _elements[i].Pdefs->isEdge = false;
                
                // If element is type tetrahedron and is a p element,
                // do the Ainsworth & Coyle trick
                if (_elements[i].Type.ElementCode == 504) [self FEMMesh_convertToACTetra:&_elements[i]];
                [pMaps getRefPElementNodesForElement:&_elements[i]
                                               nodeU:_elements[i].Type.NodeU
                                               nodeV:_elements[i].Type.NodeV
                                               nodeW:_elements[i].Type.NodeW];
            } else {
                // Clear P element definitions and set manual bubbles
                _elements[i].Pdefs = NULL;
                _elements[i].BDOFs = max(0, inDofs[4]);
            }
            
            self.maxElementNodes = max(self.maxElementNodes, _elements[i].Type.NumberOfNodes);
            
        } else {
            NSLog(@"loadMeshForModel: Unknown element type %d, ignoring element.\n", type);
            
        }
    }
    free_ivector(inDofs, 0, 6);
    free_ivector(nodes, 0, MAX_ELEMENT_NODES-1);
    
    minEIndex = min_array(elementTags, self.numberOfBulkElements);
    maxEIndex = max_array(elementTags, self.numberOfBulkElements);
    
    localEPerm = intvec(0, (maxEIndex-minEIndex+1)-1);
    memset( localEPerm, 0, (maxEIndex-minEIndex+1)*sizeof(int) );
    for (i=0; i<self.numberOfBulkElements; i++) {
        localEPerm[elementTags[i]-minEIndex] = i;
    }
    
    // TODO: Read element property file not supported yet. Is that really needed?
    
    //Mesh boundary elements
    coord = doublevec(0, self.maxElementNodes-1);
    memset( coord, 0.0, self.maxElementNodes*sizeof(double) );
    for (i=self.numberOfBulkElements; i<self.numberOfBulkElements+self.numberOfBoundaryElements; i++) {
        
        [meshIO getMeshBoundaryElement:&tag boundary:&bndry leftElement:&left rightElement:&right type:&type nodes:nodes coord:coord];
        if (meshIO.info != 0) {
            self.numberOfBoundaryElements = i - self.numberOfBulkElements;
            break;
        }
        
        if (left >= minEIndex && left <= maxEIndex) {
            left = localEPerm[left - minEIndex];
        } else if (left > 0) {
            NSLog(@"loadMeshForModel: %d Boundary parent out of range: %d %d\n", *partID, tag, left);
            left = 0;
        }
        
        if (right >= minEIndex && right <= maxEIndex) {
            right = localEPerm[right - minEIndex];
        } else if (right > 0) {
            NSLog(@"loadMeshForModel: %d Boundary parent out of range: %d %d\n", *partID, tag, right);
            right = 0;
        }
        
        _elements[i].ElementIndex = i+1;
        elmType = NULL;
        elmType = [elmDescription getElementType:type inMesh:self stabilization:NO];
        _elements[i].Type = *elmType;
        
        self.maxElementNodes = max(self.maxElementNodes, type-100*(type/100));
        
        if (elmType != NULL) {
            n = _elements[i].Type.NumberOfNodes;
            
            _elements[i].BoundaryInfo = (BoundaryInfo_t*) malloc( sizeof(BoundaryInfo_t));
            if (_elements[i].BoundaryInfo == NULL) errorfunct("loadMeshForModel", "Unable to allocate the boundary info.");
            initBoundaryInfo(_elements[i].BoundaryInfo);
            
            _elements[i].BoundaryInfo->Constraint = 0;
            
            for (j=0; j<model.numberOfBoundaries; j++) {
                if ([(model.boundaryID)[j] intValue] == bndry) {
                    addr1 = 1;
                    addr2 = model.numberOfBoundaryConditions;
                    _elements[i].BoundaryInfo->Constraint = [listUtil listGetInteger:model inArray:[(model.boundaries)[j] valuesList] forVariable:@"boundary condition" info:&found minValue:&addr1 maxValue:&addr2];
                    break;
                }
            }
            
            //TODO: If parallel mesh, set the element partition index here
            _elements[i].PartIndex = *partID;
            
            _elements[i].BodyID = 0;
        
            defaultTargetBC = 0;
            for (j=0; j<model.numberOfBoundaryConditions; j++) {
                if ([listUtil listGetLogical:model inArray:[(model.boundaryConditions)[j] valuesList] forVariable:@"default target" info:&found] == YES)
                    defaultTargetBC = j+1;
                
                gotIt = [listUtil listGetIntegerArray:model inArray:[(model.boundaryConditions)[j] valuesList] forVariable:@"target boundaries" buffer:&bdList];
                if (gotIt == YES) {
                    for (k=0; k<bdList.m; k++) {
                        if (bdList.ivector[k] == bndry) {
                            _elements[i].BoundaryInfo->Constraint = j+1;
                            addr1 = 1;
                            addr2 = model.numberOfBodies;
                            _elements[i].BodyID = [listUtil listGetInteger:model inArray:[(model.boundaryConditions)[j] valuesList] forVariable:@"body id" info:&found minValue:&addr1 maxValue:&addr2];
                        }
                    }
                }
            }
            
            j = _elements[i].BoundaryInfo->Constraint;
            if ((j<=0 || j>model.numberOfBoundaryConditions) && defaultTargetBC > 0); {
                _elements[i].BoundaryInfo->Constraint = defaultTargetBC;
                addr1 = 1;
                addr2 = model.numberOfBodies;
                _elements[i].BodyID = [listUtil listGetInteger:model inArray:[(model.boundaryConditions)[defaultTargetBC-1] valuesList] forVariable:@"body id" info:&found minValue:&addr1 maxValue:&addr2];
            }
            
            _elements[i].BoundaryInfo->Outbody = -1;
            j = _elements[i].BoundaryInfo->Constraint;
            if (j>0 && j<=model.numberOfBoundaryConditions) {
                addr1 = model.numberOfBodies;
                _elements[i].BoundaryInfo->Outbody = [listUtil listGetInteger:model inArray:[(model.boundaryConditions)[j-1] valuesList] forVariable:@"normal target body" info:&found minValue:NULL maxValue:&addr1];
            }
            
            _elements[i].NodeIndexes = intvec(0, n-1);
            _elements[i].sizeNodeIndexes = n;
            
            // Set local to global mapping for boundary element
            for (j=0; i<n; j++) {
                _elements[i].NodeIndexes[j] = localPerm[nodes[j] - minIndex];
            }
            
            _elements[i].EdgeIndexes = NULL;
            _elements[i].FaceIndexes = NULL;
            
            _elements[i].BoundaryInfo->Left = NULL;
            if (left >= 1 ) {
                _elements[i].BoundaryInfo->Left = &_elements[left];
                
            }
            
            _elements[i].BoundaryInfo->Right = NULL;
            if (right >= 1) {
                _elements[i].BoundaryInfo->Right = &_elements[right];
            }
            
            _elements[i].BoundaryInfo->GebhardtFactors = NULL;
            
            _elements[i].NDOFs = n;
            if (_elements[i].BoundaryInfo->Left != NULL) {
                if (_elements[i].BoundaryInfo->Left->NDOFs == 0) _elements[i].NDOFs = 0;
                
                if (_elements[i].Type.dimension == 1) {
                    _elements[i].BDOFs = edgeDofs[_elements[i].BoundaryInfo->Left->ElementIndex-1];
                } else {
                    _elements[i].BDOFs = faceDofs[_elements[i].BoundaryInfo->Left->ElementIndex-1];
                }
            }
            
            if (_elements[i].BoundaryInfo->Right != NULL) {
                if (_elements[i].BoundaryInfo->Right->NDOFs == 0) _elements[i].NDOFs = 0;
                
                if (_elements[i].Type.dimension == 1) {
                    _elements[i].BDOFs = edgeDofs[_elements[i].BoundaryInfo->Right->ElementIndex-1];
                } else {
                    _elements[i].BDOFs = faceDofs[_elements[i].BoundaryInfo->Right->ElementIndex-1];
                }
            }
            
            _elements[i].DGDOFs = 0;
            _elements[i].Pdefs = NULL;
            _elements[i].DGIndexes = NULL;
            _elements[i].EdgeIndexes = NULL;
            _elements[i].FaceIndexes = NULL;
            _elements[i].BubbleIndexes = NULL;
        } else {
            NSLog(@"loadMeshForModel: Unknown element type %d, ignoring element.\n", type);
        }
    }
    if (self.maxElementDofs <= 0) self.maxElementDofs = self.maxElementNodes;
    
    free_dvector(coord, 0, self.maxElementNodes-1);
    if (bdList.ivector != NULL) {
        free_ivector(bdList.ivector, 0, bdList.m-1);
        bdList.ivector = NULL;
    }
    free_ivector(localEPerm, 0, (maxEIndex-minEIndex+1)-1);
    
    meshUtils = [[FEMMeshUtils alloc] init];
    if (needEdges == YES) {
        [meshUtils findEdgesForMesh:self findEdges:NULL];
        
        // Set edge and face polynomial degree and degrees of freedom for all elements
        for (i=0; i<self.numberOfBulkElements; i++) {
            
            // Iterate in each edge of element
            for (j=0; j<_elements[i].Type.NumberOfEdges; j++) {
                
                // Set attributes of p element edges
                if (_elements[i].Pdefs != NULL) {
                    // Set edge polynomial degree and dofs
                    _edges[_elements[i].EdgeIndexes[j]].Pdefs->p = max(_elements[i].Pdefs->p, _edges[_elements[i].EdgeIndexes[j]].Pdefs->p);
                    _edges[_elements[i].EdgeIndexes[j]].BDOFs =
                            max(_edges[_elements[i].EdgeIndexes[j]].BDOFs, _edges[_elements[i].EdgeIndexes[j]].Pdefs->p - 1 );
                    _edges[_elements[i].EdgeIndexes[j]].Pdefs->isEdge = true;
                    // Get Gauss points for edge. If no dofs, two Gauss points are
                    // still needed for integration of linear equations.
                    _edges[_elements[i].EdgeIndexes[j]].Pdefs->GaussPoints =
                            pow( (_edges[_elements[i].EdgeIndexes[j]].BDOFs+2), _edges[_elements[i].EdgeIndexes[j]].Type.dimension );
                    if (_edges[_elements[i].EdgeIndexes[j]].BoundaryInfo->Left != NULL) {
                        [meshUtils assignLocalNumberToEdgeElement:&_edges[_elements[i].EdgeIndexes[j]] fromElement:_edges[_elements[i].EdgeIndexes[j]].BoundaryInfo->Left inMesh:self];
                    } else {
                        [meshUtils assignLocalNumberToEdgeElement:&_edges[_elements[i].EdgeIndexes[j]] fromElement:_edges[_elements[i].EdgeIndexes[j]].BoundaryInfo->Right inMesh:self];
                    }
                } else { // Other element types which need edge dofs
                    _edges[_elements[i].EdgeIndexes[j]].BDOFs = max(edgeDofs[i], _edges[_elements[i].EdgeIndexes[j]].BDOFs);
                }
                // Get maximum dof for edges
                self.maxEdgeDofs = max(_edges[_elements[i].EdgeIndexes[j]].BDOFs, self.maxEdgeDofs);
            }
            
            // Iterate in each face of element
            for (j=0; j<_elements[i].Type.NumberOfFaces; j++) {
            
                // Set attributes of element faces
                if (_elements[i].Pdefs != NULL) {
                    // Set face polynomial degree and dofs
                    _faces[_elements[i].FaceIndexes[j]].Pdefs->p = max(_elements[i].Pdefs->p, _faces[_elements[i].FaceIndexes[j]].Pdefs->p);
                    // Get number of face dofs
                    _faces[_elements[i].FaceIndexes[j]].BDOFs = max(_faces[_elements[i].FaceIndexes[j]].BDOFs, [pMaps getFaceDofsForElement:&_elements[i] polyDegree:_faces[_elements[i].FaceIndexes[j]].Pdefs->p faceNumber:j]);
                    _faces[_elements[i].FaceIndexes[j]].Pdefs->isEdge = true;
                    _faces[_elements[i].FaceIndexes[j]].Pdefs->GaussPoints = [pMaps getNumberOfGaussPointsForFace:&_faces[_elements[i].FaceIndexes[j]] inMesh:self];
                    if (_faces[_elements[i].FaceIndexes[j]].BoundaryInfo->Left != NULL) {
                        [meshUtils assignLocalNumberToEdgeElement:&_faces[_elements[i].FaceIndexes[j]] fromElement:_faces[_elements[i].FaceIndexes[j]].BoundaryInfo->Left inMesh:self];
                    } else {
                        [meshUtils assignLocalNumberToEdgeElement:&_faces[_elements[i].FaceIndexes[j]] fromElement:_faces[_elements[i].FaceIndexes[j]].BoundaryInfo->Right inMesh:self];
                    }
                } else {
                    _faces[_elements[i].FaceIndexes[j]].BDOFs  = max(faceDofs[i], _faces[_elements[i].FaceIndexes[j]].BDOFs);
                }
                
                // Get maximum dof for faces
                self.maxFaceDofs = max(_faces[_elements[i].FaceIndexes[j]].BDOFs, self.maxFaceDofs);
            }
        }
        
        // Set local edges for boundary elements
        for (i=self.numberOfBulkElements; i<self.numberOfBulkElements+self.numberOfBoundaryElements; i++) {
            
            // Here set local number and copy attributes to this boundary element for left parent.
            if (_elements[i].BoundaryInfo->Left != NULL) {
                // Local edges are only assigned for p elements
                if (_elements[i].BoundaryInfo->Left->Pdefs != NULL) {
                    [self allocatePDefinitionsForElement:&_elements[i]];
                    _elements[i].Pdefs->isEdge = true;
                    [meshUtils assignLocalNumberToEdgeElement:&_elements[i] fromElement:_elements[i].BoundaryInfo->Left inMesh:self];
                }
            }
            
            // Here set local number and copy attributes to this boundary element for right parent.
            if (_elements[i].BoundaryInfo->Right != NULL) {
                // Local edges are only assigned for p elements
                if (_elements[i].BoundaryInfo->Right->Pdefs != NULL) {
                    [self allocatePDefinitionsForElement:&_elements[i]];
                    _elements[i].Pdefs->isEdge = true;
                    [meshUtils assignLocalNumberToEdgeElement:&_elements[i] fromElement:_elements[i].BoundaryInfo->Right inMesh:self];
                }
            }
        }
    }
    
    // Set Gauss points for each element
    for (i=0; i<self.numberOfBulkElements; i++) {
        if (_elements[i].Pdefs != NULL) {
            _elements[i].Pdefs->GaussPoints = [pMaps getNumberOfGaussPointsForElement:&_elements[i] inMesh:self];
        }
        
        // Set max element dofs here (because element size may have changed when edges and faces have been set).
        // This is the absolute worst case. Element which has MaxElementDOFS may not even be present as a real
        // element.
        self.maxElementDofs = max(self.maxElementDofs, (_elements[i].Type.NumberOfNodes+_elements[i].Type.NumberOfEdges*self.maxEdgeDofs+_elements[i].Type.NumberOfFaces*self.maxFaceDofs+_elements[i].BDOFs), _elements[i].DGDOFs);
        
        self.maxBdofs = max(_elements[i].BDOFs, self.maxBdofs);
    }
    
    for (i=0; i<self.numberOfBulkElements; i++) {
        if (_elements[i].BDOFs > 0) {
            _elements[i].BubbleIndexes = intvec(0, _elements[i].BDOFs-1);
            _elements[i].sizeBubbleIndexes = _elements[i].BDOFs;
            for (j=0; j<_elements[i].BDOFs; j++) {
                _elements[i].BubbleIndexes[j] = self.maxBdofs*(i)+j;
            }
        }
    }
    if (edgeDofs != NULL) {
        free_ivector(edgeDofs, 0, self.numberOfBulkElements-1);
        free_ivector(faceDofs, 0, self.numberOfBulkElements-1);
        edgeDofs = NULL;
        faceDofs = NULL;
    }
    
    // Reallocate coordinate arrays for iso-parametric p-elements
    n = self.numberOfNodes + self.maxEdgeDofs * self.numberOfEdges + self.maxFaceDofs * self.numberOfFaces + self.maxBdofs * self.numberOfBulkElements;
    
    if (n>self.numberOfNodes) {
        cCoord = _globalNodes->x;
        _globalNodes->x = doublevec(0, n-1);
        memset( _globalNodes->x, 0.0, n*sizeof(double) );
        for (i=0; i<self.numberOfNodes; i++) {
            _globalNodes->x[i] = cCoord[i];
        }
        free_dvector(cCoord, 0, self.numberOfNodes-1);
        
        cCoord = _globalNodes->y;
        _globalNodes->y = doublevec(0, n-1);
        memset( _globalNodes->y, 0.0, n*sizeof(double) );
        for (i=0; i<self.numberOfNodes; i++) {
            _globalNodes->y[i] = cCoord[i];
        }
        free_dvector(cCoord, 0, self.numberOfNodes-1);
        
        cCoord = _globalNodes->z;
        _globalNodes->z = doublevec(0, n-1);
        memset( _globalNodes->z, 0.0, n*sizeof(double) );
        for (i=0; i<self.numberOfNodes; i++) {
            _globalNodes->z[i]= cCoord[i];
        }
        free_dvector(cCoord, 0, self.numberOfNodes-1);
    }
    
    if (parallel) {
        //TODO: Implementation here if we use a parallel mesh
    }
    
    free_ivector(localPerm, 0, (maxIndex-minIndex+1)-1);
    free_ivector(elementTags, 0, (self.numberOfBulkElements+1)-1);
    
    [meshIO closeMesh];
    [meshIO close];
    
    // If periodic BC given, compute boundary mesh projector    
    for (FEMBoundaryCondition *boundary in model.boundaryConditions) {
        boundary.pMatrix = nil;
    }
    
    i = 0;
    for (FEMBoundaryCondition *boundary in model.boundaryConditions) {
        addr1 = 1;
        addr2 = model.numberOfBoundaryConditions;
        k = [listUtil listGetInteger:model inArray:boundary.valuesList forVariable:@"periodic bc" info:&gotIt minValue:&addr1 maxValue:&addr2];
        projector = [meshUtils periodicProjectorInModel:model forMesh:self boundary:i target:k-1];
        if (projector != nil) boundary.pMatrix = projector;
        projector = nil;
        i++;
    }
    
    model.dimension = saveDim;
    
    [elmDescription deallocation];
    [pMaps deallocation];
    free_ivector(nodeTags, 0, self.numberOfNodes-1); // This should not deallocated if parallel mesh is supported
    free_ivector(countByType, 0, 63);
    free_ivector(types, 0, 63);
}

#pragma mark Nodes getter

-(Nodes_t *)getNodes {
    
    return _globalNodes;
}

#pragma mark Elements getter

-(Element_t *)getElements {
    
    return _elements;
}

#pragma mark Edges getter

-(Element_t *)getEdges {
    
    return _edges;
}

#pragma mark Faces getter

-(Element_t *)getFaces {
    
    return _faces;
}

#pragma mark Quadrant getter
-(Quadrant_t *)getQuadrant {
    
    return _rootQuadrant;
}

#pragma mark Test associativity

-(BOOL)isAssociatedEdges {
    
    if (_edges != NULL) {
        return YES;
    } else {
        return NO;
    }
}

-(BOOL)isAssociatedFaces {
    
    if (_faces != NULL) {
        return YES;
    } else {
        return NO;
    }
}

#pragma mark Deallocation

-(void)deallocationMeshVariables {
    
    variableArraysContainer *varContainers = NULL;
    
    for (FEMVariable *variable in self.variables) {
        varContainers = variable.getContainers;
        
        // Used to skip variables such as time, timestep, timestep size etc
        if (varContainers->sizeValues == variable.dofs) continue;
        
        if ([variable.name isEqualToString:@"coordinate 1"] || [variable.name isEqualToString:@"coordinate 2"] ||
                [variable.name isEqualToString:@"coordinate 3"]) continue;
        
        if (variable.secondary == YES) continue;
        
        if (varContainers->Perm != NULL) {
            free_ivector(varContainers->Perm, 0, varContainers->sizePerm-1);
            varContainers->Perm = NULL;
        }
        
        if (varContainers->Values != NULL) {
            free_dvector(varContainers->Values, 0, varContainers->sizeValues-1);
            varContainers->Values = NULL;
        }
        
        if (varContainers->PrevValues != NULL) {
            free_dmatrix(varContainers->PrevValues, 0, varContainers->size1PrevValues-1, 0, varContainers->size2PrevValues-1);
            varContainers->PrevValues = NULL;
        }
        
        if (varContainers->EigenValues != NULL) {
            free_cdvector(varContainers->EigenValues, 0, varContainers->sizeEigenValues-1);
            varContainers->EigenValues = NULL;
        }
        
        if (varContainers->EigenVectors != NULL) {
            free_cdmatrix(varContainers->EigenVectors, 0, varContainers->size1EigenVectors-1, 0, varContainers->size2EigenVectors-1);
            varContainers->EigenVectors = NULL;
        }
        
        if (varContainers->SteadyValues != NULL) {
            free_dvector(varContainers->SteadyValues, 0, varContainers->sizeSteadyValues-1);
            varContainers->SteadyValues = NULL;
        }
        
        if (varContainers->NonLinValues != NULL) {
            free_dvector(varContainers->NonLinValues, 0, varContainers->sizeNonLinValues-1);
            varContainers->NonLinValues = NULL;
        }
    }
    free(varContainers);
    self.variables = nil;
}

-(void)deallocateMeshEdgeTables {
    
    int i;
    
    if (_edges != NULL) {
        for (i=0; i<self.numberOfEdges; i++) {
            if (_edges[i].NodeIndexes != NULL) {
                free_ivector(_edges[i].NodeIndexes, 0, _edges[i].sizeNodeIndexes-1);
            }
            if (_edges[i].BoundaryInfo != NULL) {
                free(_edges[i].BoundaryInfo);
                _edges[i].BoundaryInfo = NULL;
            }
        }
        free(_edges);
    }
    _edges = NULL;
    self.numberOfEdges = 0;
    
    for (i=0; i<self.numberOfBulkElements; i++) {
        if (_elements[i].EdgeIndexes != NULL) {
            free_ivector(_elements[i].EdgeIndexes, 0, _elements[i].sizeEdgeIndexes-1);
            _elements[i].EdgeIndexes = NULL;
        }
    }
}

-(void)deallocateMeshFaceTables {
    
    int i;
    
    if (_faces != NULL) {
        for (i=0; i<self.numberOfFaces; i++) {
            if (_faces[i].NodeIndexes != NULL) {
                free_ivector(_faces[i].NodeIndexes, 0, _faces[i].sizeNodeIndexes-1);
            }
            if (_faces[i].BoundaryInfo != NULL) {
                free(_faces[i].BoundaryInfo);
                _faces[i].BoundaryInfo = NULL;
            }
        }
        free(_faces);
    }
    _faces = NULL;
    self.numberOfFaces = 0;
    
    for (i=0; i<self.numberOfBulkElements; i++) {
        if (_elements[i].FaceIndexes != NULL) {
            free_ivector(_elements[i].FaceIndexes, 0, _elements[i].sizeFaceIndexes-1);
            _elements[i].FaceIndexes = NULL;
        }
    }
}

-(void)deallocateMeshViewFactorTables {
    
    int i;
    
    if (_viewFactors != NULL) {
        for (i=0; i<self.numberOfViewFactors; i++) {
            if (_viewFactors[i].Factors != NULL) free_dvector(_viewFactors[i].Factors, 0, _viewFactors[i].sizeFactors-1);
            if (_viewFactors[i].Elements != NULL) free_ivector(_viewFactors[i].Elements, 0, _viewFactors[i].sizeElements-1);
        }
        free(_viewFactors);
    }
}

-(void)deallocateQuadrantTree {
    
    [self FEMMesh_deallocateQuadrantTree:_rootQuadrant];
}

-(void)deallocation {
    
    int i;
    
    //Deallocate mesh variables
    [self deallocationMeshVariables];
    
    // Deallocate mesh geometry (nodes, elements and edges)
    if (_globalNodes != NULL) {
        if (_globalNodes->x != NULL) free_dvector(_globalNodes->x, 0, self.numberOfNodes-1);
        if (_globalNodes->y != NULL) free_dvector(_globalNodes->y, 0, self.numberOfNodes-1);
        if (_globalNodes->z != NULL) free_dvector(_globalNodes->z, 0, self.numberOfNodes-1);
        
        //TODO: Add support for parallel runs
        
        free(_globalNodes);
    }
    _globalNodes = NULL;
    
    if (_edges != NULL) [self deallocateMeshEdgeTables];
    _edges = NULL;
    
    if (_faces != NULL) [self deallocateMeshFaceTables];
    _faces = NULL;
    
    if (_viewFactors != NULL) [self deallocateMeshViewFactorTables];
    _viewFactors = NULL;
    
    if (_elements != NULL) {
        for (i=0; i<self.numberOfBulkElements+self.numberOfBoundaryElements; i++) {
            // Boundary structure for boundary elements
            if (_elements[i].copy == true) continue;
            
            if (i >= self.numberOfBulkElements) {
                if (_elements[i].BoundaryInfo != NULL) {
                    if (_elements[i].BoundaryInfo->GebhardtFactors != NULL) {
                        if (_elements[i].BoundaryInfo->GebhardtFactors->Elements != NULL) {
                            free_ivector(_elements[i].BoundaryInfo->GebhardtFactors->Elements, 0, _elements[i].BoundaryInfo->GebhardtFactors->sizeElements-1);
                            free_dvector(_elements[i].BoundaryInfo->GebhardtFactors->Factors, 0, _elements[i].BoundaryInfo->GebhardtFactors->sizeFactors-1);
                        }
                        free(_elements[i].BoundaryInfo->GebhardtFactors);
                    }
                    free(_elements[i].BoundaryInfo);
                }
                _elements[i].BoundaryInfo = NULL;
            }
            
            if (_elements[i].NodeIndexes != NULL) free_ivector(_elements[i].NodeIndexes, 0, _elements[i].sizeNodeIndexes-1);
            _elements[i].NodeIndexes = NULL;
            
            if (_elements[i].EdgeIndexes != NULL) free_ivector(_elements[i].EdgeIndexes, 0, _elements[i].sizeEdgeIndexes-1);
            _elements[i].EdgeIndexes = NULL;
            
            if (_elements[i].FaceIndexes != NULL) free_ivector(_elements[i].FaceIndexes, 0, _elements[i].sizeFaceIndexes-1);
            _elements[i].FaceIndexes = NULL;
            
            if (_elements[i].DGIndexes != NULL) free_ivector(_elements[i].DGIndexes, 0, _elements[i].sizeDGIndexes-1);
            _elements[i].DGIndexes = NULL;
            
            if (_elements[i].BubbleIndexes != NULL) free_ivector(_elements[i].BubbleIndexes, 0, _elements[i].sizeBubbleIndexes-1);
            _elements[i].BubbleIndexes = NULL;
            
            if (_elements[i].Pdefs != NULL) free(_elements[i].Pdefs);
            _elements[i].Pdefs = NULL;
        }
        free(_elements);
    }
    _elements = NULL;
    
    // Deallocate mesh to mesh projector structures
    for (FEMProjector *projector in _projectors) {
        [projector.matrix deallocation];
        [projector.tMatrix deallocation];
    }
    _projectors = nil;

    // Deallocate quadrant tree (used in mesh to mesh interpolation)
    [self FEMMesh_deallocateQuadrantTree:_rootQuadrant];
    _rootQuadrant = NULL;
}

-(void)Simple2DMeshBorders:(double*)borders withSize:(int*) intervals elemetCode:(int) elementID {
    
    int i, j, k, l;
    int m, n, o, p;
    double minx, maxx, miny, maxy, incri, incrj;
    double meshsize1, meshsize2;
    int domainShape;
    int **side1, **side2, **side3, **side4;
    BOOL ct;
    
    minx = borders[0];
    maxx = borders[1];
    miny = borders[2];
    maxy = borders[3];
    
    // Rectangle or square?
    if ( (maxx-minx) > 0 && (maxy-miny) > 0 && (maxx-minx) == (maxy-miny) ) {
        domainShape = 1;                                                               // valid square
    } else  if ((maxx-minx) > 0 && (maxy-miny) > 0 && (maxx-minx) != (maxy-miny)) {
        domainShape = 2;                                                              // valid rectangle
    }
    else {
        domainShape = -1;                                                             // Something wrong with borders input
    }
    
    [self setDimension:2];
    
    // Build nodes
    switch (domainShape) {
        case 1:
            
            meshsize1 = (maxx-minx) / intervals[0];
            meshsize2 = meshsize1;
            intervals[1] = intervals[0];
            
            // Allocate memory for the GlobalNodes
            // TODO: [self AllocateNodes:0 :((intervals[0]+1)*(intervals[0]+1))-1];
            
            break;
            
        case 2:
            
            if (intervals[1] <= 0) {
                errorfunct("Simple2DMesh", "No intervals given for second direction discretization in rectangle domain.");
            }
            
            meshsize1 = (maxx-minx) / intervals[0];
            meshsize2 = (maxy-miny) / intervals[1];
            
            // Allocate memory for the GlobalNodes
            //TODO: [self AllocateNodes:0 :((intervals[0]+1)*(intervals[1]+1))-1];
            
            break;
            
        default:
            errorfunct("Simple2DMesh", "Failure in method Simple2DMesh. Cant initialize mesh.");
            break;
    }
    
    l = 1;
    k = 0;
    incri = minx;
    incrj = miny;
    
    for (i=0; i<=intervals[0]; i++) {
        for (j=0; j<=intervals[1]; j++) {
            
            _globalNodes->x[k] = incri;
            _globalNodes->y[k] = incrj;
            _globalNodes->z[k] = 0.0;
            incrj = incrj + meshsize2;
            l++;
            k++;
            
        }
        incrj = miny;
        incri = incri + meshsize1;
    }
    
    
    // Build elements and remember boundary elements
    switch (domainShape) {
        case 1:
            
            // Allocate memory for the Elements
            //Elements = ElementVec(0, (intervals[0]*intervals[0]));
            //Elements->NodeIndexes = intvec(0, 3);
            
            side1 = intmatrix(0, intervals[0], 0, 2);
            side2 = intmatrix(0, intervals[0], 0, 2);
            side3 = intmatrix(0, intervals[0], 0, 2);
            side4 = intmatrix(0, intervals[0], 0, 2);
            
            break;
            
        case 2:
            
            //Elements = ElementVec(0, (intervals[0]*intervals[1]));
            //Elements->NodeIndexes = intvec(0, 3);
            
            side1 = intmatrix(0, intervals[0], 0, 2);
            side2 = intmatrix(0, intervals[1], 0, 2);
            side3 = intmatrix(0, intervals[0], 0, 2);
            side4 = intmatrix(0, intervals[1], 0, 2);
            
            break;
            
    }
    
    k = 0;
    l = 1;
    
    m = 0;
    n = 0;
    o = 0;
    p = 0;
    
    for (i=0; i<((intervals[0]+1)*(intervals[1]+1))-(intervals[1]+2); i++) {
        
        // Top boundary elements
        for (j=1;j<(intervals[1]+1);j++) {
            if ( i == 0+(j*(intervals[1]+1)-1) ) {
                side3[o][2] = k-1;
                o++;
                ct = YES;
                break;
            }
        }
        
        if (ct == YES) {
            ct = NO;
            continue;
        }
        
        _elements[k].ElementIndex = l;
        _elements[k].BodyID = 1;
        _elements[k].Type.ElementCode = 404;
        
        // Left boundary elements
        if ( i < (intervals[1]+1) ) {
            side4[p][2] = k;
            p++;
        }
        
        // Bottom boundary elements
        for ( j=0;j<(intervals[1]+1);j++ ) {
            if ( i == 0+j*(intervals[1]+1) ) {
                side1[m][2] = k;
                m++;
                break;
            }
        }
        
        // Right boundary elements
        if ( i >= ((intervals[0]+1)*(intervals[1]+1))-(2*(intervals[1]+1)) ) {
            side2[n][2] = k;
            n++;
            break;
        }
        
        k++;
        l++;
    }
    
    // Build boundary elements
    
    switch (domainShape) {
        case 1:
            
            // Allocate memory for the boundary elements
            //BDElements = BDElementVec(0, (intervals[0]*4)-1);
            //BDElements->NodeIndexes = intvec(0, 1);
            
            break;
            
        case 2:
            
            //BDElements = BDElementVec(0, ((intervals[0]*2)+(intervals[1]*2))-1);
            //BDElements->NodeIndexes = intvec(0, 1);
            
            break;
            
    }
    
    k = 0;
    l = 1;
    
    // Boundaries are numbered anti-clockwise
    // First boundary -> bottom side
    for (i=0; i<m; i++) {
        
        //TODO:
    }
    
    // Second boundary -> right side
    for (i=0; i<n; i++) {
        
        //TODO:
    }
    
    // Third boundary -> top side
    for (i=0; i<o; i++) {
        
        //TODO:
    }
    
    // Fourth boundary -> left side
    for (i=0; i<p; i++) {
        
        //TODO:
    }
    
    switch (domainShape) {
        case 1:
            
            free_imatrix(side1, 0, intervals[0], 0, 2);
            free_imatrix(side2, 0, intervals[0], 0, 2);
            free_imatrix(side3, 0, intervals[0], 0, 2);
            free_imatrix(side4, 0, intervals[0], 0, 2);
            
            break;
            
        case 2:
            
            free_imatrix(side1, 0, intervals[0], 0, 2);
            free_imatrix(side2, 0, intervals[1], 0, 2);
            free_imatrix(side3, 0, intervals[0], 0, 2);
            free_imatrix(side4, 0, intervals[1], 0, 2);
            
            break;
    }
    
}


@end
