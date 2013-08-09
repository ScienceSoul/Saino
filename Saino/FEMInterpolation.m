//
//  FEMInterpolation.m
//  Saino
//
//  Created by Seddik hakime on 27/09/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMInterpolation.h"
#import "FEMPElementMaps.h"
#import "FEMElementDescription.h"
#import "Utils.h"

@interface FEMInterpolation ()
-(void)FEMInterpolation_findPointsQuadrantForPoint:(double *)point dimension:(int)dim motherQuadrant:(Quadrant_t *)mother;
@end

@implementation FEMInterpolation

#pragma mark Private methods

-(void)FEMInterpolation_findPointsQuadrantForPoint:(double *)point dimension:(int)dim motherQuadrant:(Quadrant_t *)mother {
    
    int i, k;
    double bbox[6], eps3, *valuesPtr, maxVal;
    double const eps2 = 0.0;
    Quadrant_t *childQuadrant;
    BOOL found;
    
    // Loop over child quadrants
    found = NO;
    for (i=0; i<(int)pow(2, dim); i++) {
        childQuadrant = mother->childQuadrants[i].quadrant;
        memcpy(bbox, childQuadrant->boundingBox, sizeof(childQuadrant->boundingBox));
        
        //********* Note: eps2 is set to 0 at the moment *********
        
        maxVal = -HUGE_VAL;
        k = 3;
        valuesPtr = &bbox[0];
        for (i=0; i<3; i++) {
            if ( (*(valuesPtr + k) - *(valuesPtr + i)) > maxVal ) {
                maxVal = *(valuesPtr + k) - *(valuesPtr + i);
            }
            k++;
        }
        eps3 = eps2 * maxVal;
        for (i=0; i<3; i++) {
            bbox[i] = bbox[i] - eps3;
        }
        for (i=3; i<6; i++) {
            bbox[i] = bbox[i] + eps3;
        }
        
        // Is the point in childQuadrant[i]
        if ( (point[0] >= bbox[0]) && (point[0] <= bbox[3]) && (point[1] >= bbox[1]) && (point[1] <= bbox[4]) &&
            (point[2] >= bbox[2]) && (point[2] <= bbox[5]) ) {
            found = YES;
            break;
        }
    }
    
    if ( found == NO) {
        mother = NULL;
        return;
    }
    
    mother = childQuadrant;
    
    // Are we already in the leaf quadrant?
    // If not, search for the next generation child quadrants for the point
    if (mother->childQuadrants != NULL) {
        [self FEMInterpolation_findPointsQuadrantForPoint:point dimension:dim motherQuadrant:mother];
    }
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

-(void)findLeafElementsForPoint:(double *)point dimension:(int)dim rootQuadrant:(Quadrant_t *)root leafQuadrant:(Quadrant_t *)leaf {
    
    leaf = root;
    
    // Find recursively the last generation quadrant the point belongs to
    [self FEMInterpolation_findPointsQuadrantForPoint:point dimension:dim motherQuadrant:leaf];
}

/*********************************************************************************************
 
    Check whether a given point belongs to a given bulk element.
    If if does, retuns the local coordinates in the bulk element.
 
    Element_t *element   ->  bulk element we are checking
    Nodes_t nodes        ->  the nodal points of the bulk element
    double *aPoint       ->  the point being searched
    double *localCoords  ->  local coordinates corresponding to the global ones
    double *globaleps    ->  required accuracy of global coordinates
    double *localeps     ->  required accuracy of local coordinates
    double *numericeps   ->  accuracy of numerical operations
    double *globaldist   ->  returns the distance from the element in global coordinates
    doubl *localdist     ->  returns the distance from the element in local coordinates

*********************************************************************************************/
-(BOOL)isPointInElement:(Element_t *)element elementNodes:(Nodes_t *)nodes point:(double *)aPoint localCoordinates:(double *)localCoords globalEpsilon:(double *)globaleps localEpsilon:(double *)localeps numericEpsilon:(double *)numericeps globalDistance:(double *)globaldist localDistance:(double *)localdist model:(FEMModel *)aModel {
    
    int i, n;
    double ug, vg, wg, xdist, ydist, zdist, sumdist, eps0, eps1, eps2, minx,
           maxx, miny, maxy, minz, maxz;
    FEMPElementMaps *elementMaps;
    FEMElementDescription *elementDescription;
    BOOL isInElement;
    
    elementMaps = [[FEMPElementMaps alloc] init];
    elementDescription = [FEMElementDescription sharedElementDescription];
    
    isInElement = NO;
    n = element->Type.NumberOfNodes;
    
    // The numeric precision
    if (numericeps != NULL) {
        eps0 = *numericeps;
    } else {
        eps0 = DBL_EPSILON;
    }
    
    // The rough check, used for global coordinates
    if (globaleps != NULL) {
        eps1 = *globaleps;
    } else {
        eps1 = 1.0e-4;
    }
    
    // The mode detailed condition used for local coordinates
    if (localeps != NULL) {
        eps2 = *localeps;
    } else {
        eps2 = 1.0e-10;
    }
    
    if (localdist != NULL) {
        *localdist = HUGE_VAL;
    }
    
    if (eps1 < 0.0) {
        // No ops
    } else if (globaldist != NULL) {
        // When distance to be computed all coordinate directions need to be checked
        
        minx = HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->x[i] < minx) {
                minx = nodes->x[i];
            }
        }
        maxx = -HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->x[i] > maxx) {
                maxx = nodes->x[i];
            }
        }
        
        miny = HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->y[i] < miny) {
                miny = nodes->y[i];
            }
        }
        maxy = -HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->y[i] > maxy) {
                maxy = nodes->y[i];
            }
        }
        
        minz = HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->z[i] < minz) {
                minz = nodes->z[i];
            }
        }
        maxz = -HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->z[i] > maxz) {
                maxz = nodes->z[i];
            }
        }
        
        xdist = max(max(aPoint[0] - maxx, 0.0), minx - aPoint[0]);
        ydist = max(max(aPoint[1] - maxy, 0.0), miny - aPoint[1]);
        zdist = max(max(aPoint[2] - maxz, 0.0), minz - aPoint[2]);
        
        *globaldist = sqrt(pow(xdist, 2.0) + pow(ydist, 2.0) + pow(zdist, 2.0));
        
        if (xdist > eps0 + eps1 * (maxx - minx) ) return isInElement;
        if (ydist > eps0 + eps1 * (maxy - miny) ) return isInElement;
        if (zdist > eps0 + eps1 * (maxz - minz) ) return isInElement;
    } else {
        // Otherwise make decision independently after each coordinate direction
        
        minx = HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->x[i] < minx) {
                minx = nodes->x[i];
            }
        }
        maxx = -HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->x[i] > maxx) {
                maxx = nodes->x[i];
            }
        }
        xdist = max(max(aPoint[0] - maxx, 0.0), minx - aPoint[0]);
        if (xdist > eps0 + eps1 * (maxx - minx)) return isInElement;
        
        
        miny = HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->y[i] < miny) {
                miny = nodes->y[i];
            }
        }
        maxy = -HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->y[i] > maxy) {
                maxy = nodes->y[i];
            }
        }
        ydist = max(max(aPoint[1] - maxy, 0.0), miny - aPoint[1]);
        if (ydist > eps0 + eps1 * (maxy - miny)) return isInElement;
        
        minz = HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->z[i] < minz) {
                minz = nodes->z[i];
            }
        }
        maxz = -HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->z[i] > maxz) {
                maxz = nodes->z[i];
            }
        }
        zdist = max(max(aPoint[2] - maxz, 0.0), minz - aPoint[2]);
        if (zdist > eps0 + eps1 * (maxz - minz)) return isInElement;
    }
    
    // Get element local coordinates from global coordinates of the point
    [elementDescription globalToLocalFromElement:element elementNodes:nodes localU:&ug localV:&vg localW:&wg x:aPoint[0] y:aPoint[1] z:aPoint[2] model:aModel];
    
    // Currently the eps of global coordinates is mixed with the eps of local
    // coordinates which is not totally satisfying. There could be a sloppier
    // global coordinate search and a more rigorous local coordinate search.
    sumdist = 0.0;
    switch (element->Type.ElementCode / 100) {
        case 1:
            sumdist = ug;
            break;
            
        case 2:
            sumdist = max(ug - 1.0, max(-ug - 1.0, 0.0));
            break;
            
        case 3:
            sumdist = max(-ug, 0.0) + max(-vg, 0.0);
            sumdist = sumdist + max(ug + vg - 1.0, 0.0);
            break;
            
        case 4:
            sumdist = max(ug - 1.0, max(-ug - 1.0, 0.0));
            sumdist = sumdist + max(vg - 1.0, max(-vg - 1.0, 0.0));
            break;
            
        case 5:
            sumdist = max(-ug, 0.0) + max(-vg, 0.0) + max(-wg, 0.0);
            sumdist = sumdist + max(ug + vg + wg - 1.0, 0.0);
            break;
            
        case 7:
            sumdist = max(-ug, 0.0) + max(-vg, 0.0);
            sumdist = sumdist + max(ug + vg - 1.0, 0.0);
            sumdist = sumdist + max(wg - 1.0, max(-wg - 1.0, 0.0));
            break;
            
        case 8:
            sumdist = max(ug - 1.0, max(-ug - 1.0, 0.0));
            sumdist = sumdist + max(vg - 1.0, max(-vg - 1.0, 0.0));
            sumdist = sumdist + max(wg - 1.0, max(-wg - 1.0, 0.0));
            break;
            
        default:
            NSLog(@"FEMInterpolation:isPointInElement: element code %d not implemented.", element->Type.ElementCode);
            break;
    }
    
    if (sumdist < eps0 + eps2) {
        isInElement = YES;
    }
    
    if (localdist != NULL) {
        *localdist = sumdist;
    }
    
    if ([elementMaps isPElement:element] == YES) {
        switch (element->Type.ElementCode / 100) {
            case 3:
                ug = 2.0 * ug + vg - 1.0;
                vg = sqrt(3.0) * vg;
                break;
                
            case 5:
                ug = 2.0 * ug + vg + wg - 1.0;
                vg = sqrt(3.0) * vg + 1/sqrt(3.0) * wg;
                wg = 2.0 * sqrt(2.0/3.0) * wg;
                break;
                
            case 6:
                wg = sqrt(2.0) * wg;
                break;
                
            case 7:
                ug = 2.0 * ug + vg - 1.0;
                vg = sqrt(3.0) * vg;
                break;
        }
    }
    
    localCoords[0] = ug;
    localCoords[1] = vg;
    localCoords[2] = wg;
    
    [elementMaps deallocation];
    return isInElement;
}

-(void)putElementsInChildQuadrants:(QuadrantPointer_t *)childQuadrants motherQuadrant:(Quadrant_t *)mother mesh:(FEMMesh *)aMesh dimension:(int)dim {
    
    int i, j, t, n;
    int elementList[(int)pow(2, dim)][mother->nElementsInQuadrant];
    double eps3, maxval, minval;
    double const eps2 = 2.5e-2;
    double bbox[6], xmin, xmax, ymin, ymax, zmin, zmax, elementSize;
    Element_t *elements;
    Nodes_t *nodes;
    BOOL elementInQuadrant;
    
    for (i=0; i<(int)pow(2, dim); i++) {
        childQuadrants[i].quadrant->nElementsInQuadrant = 0;
        childQuadrants[i].quadrant->minElementSize = 1.0e20;
    }
    
    elements = aMesh.getElements;
    nodes = aMesh.getNodes;
    
    for (t=0; t<mother->nElementsInQuadrant; t++) {
        n = elements[mother->elements[t]].Type.NumberOfNodes;
        
        // Get element coordinates
        minval = HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->x[elements[mother->elements[t]].NodeIndexes[i]] < minval ) {
                minval = nodes->x[elements[mother->elements[t]].NodeIndexes[i]];
            }
        }
        xmin = minval;
        
        maxval = -HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->x[elements[mother->elements[t]].NodeIndexes[i]] > maxval ) {
                maxval = nodes->x[elements[mother->elements[t]].NodeIndexes[i]];
            }
        }
        xmax = maxval;
        
        minval = HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->y[elements[mother->elements[t]].NodeIndexes[i]] < minval ) {
                minval = nodes->y[elements[mother->elements[t]].NodeIndexes[i]];
            }
        }
        ymin = minval;
        
        maxval = -HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->y[elements[mother->elements[t]].NodeIndexes[i]] > maxval ) {
                maxval = nodes->y[elements[mother->elements[t]].NodeIndexes[i]];
            }
        }
        ymax = maxval;
        
        minval = HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->z[elements[mother->elements[t]].NodeIndexes[i]] < minval ) {
                minval = nodes->z[elements[mother->elements[t]].NodeIndexes[i]];
            }
        }
        zmin = minval;

        maxval = -HUGE_VAL;
        for (i=0; i<n; i++) {
            if (nodes->z[elements[mother->elements[t]].NodeIndexes[i]] > maxval ) {
                maxval = nodes->z[elements[mother->elements[t]].NodeIndexes[i]];
            }
        }
        zmax = maxval;
        
        elementSize = max(max(xmax-xmin, ymax-ymin), zmax-zmin);
        
        // Is the element in one of the child quadrant?
        // Check whether element bounding box crosses any of the child
        // quadrant bounding box
        for (i=0; i<(int)pow(2, dim); i++) { // Loop over child quadrants
            
            memcpy(bbox, childQuadrants[i].quadrant->boundingBox, sizeof(childQuadrants[i].quadrant->boundingBox));
            
            eps3 = 0.0;
            eps3 = max(eps3, bbox[3]-bbox[0]);
            eps3 = max(eps3, bbox[4]-bbox[1]);
            eps3 = max(eps3, bbox[5]-bbox[2]);
            eps3 = eps2 * eps3;
            
            for (j=0; j<3; j++) {
                bbox[j] = bbox[j] - eps3;
            }
            for (j=3; j<6; j++) {
                bbox[j] = bbox[j] + eps3;
            }
            
            elementInQuadrant = YES;
            if (xmax < bbox[0] || xmin > bbox[3] || ymax < bbox[1] || ymin > bbox[4] ||
                zmax < bbox[2] || zmin > bbox[5]) elementInQuadrant = FALSE;
            
            if (elementInQuadrant == YES) {
                childQuadrants[i].quadrant->minElementSize = min(elementSize, childQuadrants[i].quadrant->minElementSize);
                
                // We allocate and store also the midlevels temporarily
                // (for the duration of the construction routine)
                elementList[i][childQuadrants[i].quadrant->nElementsInQuadrant] = mother->elements[t];
                childQuadrants[i].quadrant->nElementsInQuadrant++;
            }
        }
    }
    
    for (i=0; i<(int)pow(2, dim); i++) {
        if (childQuadrants[i].quadrant->nElementsInQuadrant != 0) {
            childQuadrants[i].quadrant->elements = intvec(0, childQuadrants[i].quadrant->nElementsInQuadrant-1);
            
            for (j=0; j<childQuadrants[i].quadrant->nElementsInQuadrant; j++) {
                childQuadrants[i].quadrant->elements[j] = elementList[i][j];
            }
        }
    }
}

-(void)createChildQuadrantFromMother:(Quadrant_t *)quadrant mesh:(FEMMesh *)aMesh dimension:(int)dim maxLeafElements:(int)maxLeafElems generation:(int *)gen {
    
    int i, n;
    QuadrantPointer_t childQuadrant[8];
    double xmin, xmax, ymin, ymax, zmin, zmax;
    
    // Create 2^dim child quadrants
    n = (int)pow(2, dim);
    quadrant->childQuadrants = (QuadrantPointer_t*) malloc( sizeof(QuadrantPointer_t) * n );
    quadrant->numberOfchildQuadrants = n;
    for (i=0; i<n; i++) {
        quadrant->childQuadrants[i].quadrant = (Quadrant_t*) malloc( sizeof(Quadrant_t) * 1);
        childQuadrant[i].quadrant = quadrant->childQuadrants[i].quadrant;
        childQuadrant[i].quadrant->nElementsInQuadrant = 0;
        childQuadrant[i].quadrant->elements = NULL;
        childQuadrant[i].quadrant->childQuadrants = NULL;
    }
    
    xmin = quadrant->boundingBox[0];
    ymin = quadrant->boundingBox[1];
    zmin = quadrant->boundingBox[2];
    
    xmax = quadrant->boundingBox[3];
    ymax = quadrant->boundingBox[4];
    zmax = quadrant->boundingBox[5];
    
    quadrant->size = max(max(xmax-xmin, ymax-ymin), zmax-zmin);
    
    childQuadrant[0].quadrant->boundingBox[0] = xmin;
    childQuadrant[0].quadrant->boundingBox[1] = ymin;
    childQuadrant[0].quadrant->boundingBox[2] = zmin;
    childQuadrant[0].quadrant->boundingBox[3] = (xmin + xmax) / 2.0;
    childQuadrant[0].quadrant->boundingBox[4] = (ymin + ymax) / 2.0;
    childQuadrant[0].quadrant->boundingBox[5] = (zmin + zmax) / 2.0;
    
    childQuadrant[1].quadrant->boundingBox[0] = (xmin + xmax) / 2.0;
    childQuadrant[1].quadrant->boundingBox[1] = ymin;
    childQuadrant[1].quadrant->boundingBox[2] = zmin;
    childQuadrant[1].quadrant->boundingBox[3] = xmax;
    childQuadrant[1].quadrant->boundingBox[4] = (ymin + ymax) / 2.0;
    childQuadrant[1].quadrant->boundingBox[5] = (zmin + zmax) / 2.0;
    
    if (dim >= 2) {
        childQuadrant[2].quadrant->boundingBox[0] = xmin;
        childQuadrant[2].quadrant->boundingBox[1] = (ymin + ymax) / 2.0;
        childQuadrant[2].quadrant->boundingBox[2] = zmin;
        childQuadrant[2].quadrant->boundingBox[3] = (xmin + xmax) / 2.0;
        childQuadrant[2].quadrant->boundingBox[4] = ymax;
        childQuadrant[2].quadrant->boundingBox[5] = (zmin + zmax) / 2.0;
        
        childQuadrant[3].quadrant->boundingBox[0] = (xmin + xmax) / 2.0;
        childQuadrant[3].quadrant->boundingBox[1] = (ymin + ymax) / 2.0;
        childQuadrant[3].quadrant->boundingBox[2] = zmin;
        childQuadrant[3].quadrant->boundingBox[3] = xmax;
        childQuadrant[3].quadrant->boundingBox[4] = ymax;
        childQuadrant[3].quadrant->boundingBox[5] = (zmin + zmax) / 2.0;
    }
    
    if (dim == 3) {
        childQuadrant[4].quadrant->boundingBox[0] = xmin;
        childQuadrant[4].quadrant->boundingBox[1] = ymin;
        childQuadrant[4].quadrant->boundingBox[2] = (zmin + zmax) / 2.0;
        childQuadrant[4].quadrant->boundingBox[3] = (xmin + xmax) / 2.0;
        childQuadrant[4].quadrant->boundingBox[4] = (ymin + ymax) / 2.0;
        childQuadrant[4].quadrant->boundingBox[5] = zmax;
        
        
        childQuadrant[5].quadrant->boundingBox[0] = (xmin + xmax) / 2.0;
        childQuadrant[5].quadrant->boundingBox[1] = ymin;
        childQuadrant[5].quadrant->boundingBox[2] = (zmin + zmax) / 2.0;
        childQuadrant[5].quadrant->boundingBox[3] = xmax;
        childQuadrant[5].quadrant->boundingBox[4] = (ymin + ymax) / 2.0;
        childQuadrant[5].quadrant->boundingBox[5] = zmax;
        
        childQuadrant[6].quadrant->boundingBox[0] = xmin;
        childQuadrant[6].quadrant->boundingBox[1] = (ymin + ymax) / 2.0;
        childQuadrant[6].quadrant->boundingBox[2] = (zmin + zmax) / 2.0;
        childQuadrant[6].quadrant->boundingBox[3] = (xmin + xmax) / 2.0;
        childQuadrant[6].quadrant->boundingBox[4] = ymax;
        childQuadrant[6].quadrant->boundingBox[5] = zmax;
        
        childQuadrant[7].quadrant->boundingBox[0] = (xmin + xmax) / 2.0;
        childQuadrant[7].quadrant->boundingBox[1] = (ymin + ymax) / 2.0;
        childQuadrant[7].quadrant->boundingBox[2] = (zmin + zmax) / 2.0;
        childQuadrant[7].quadrant->boundingBox[3] = xmax;
        childQuadrant[7].quadrant->boundingBox[4] = ymax;
        childQuadrant[7].quadrant->boundingBox[5] = zmax;
    }
    
    
    // Loop over all elements in the mother quadrant, placing them in
    // one of the 2^dim quadrants
    [self putElementsInChildQuadrants:childQuadrant motherQuadrant:quadrant mesh:aMesh dimension:dim];
    
    // Check whether we need to branch for the next level
    for (i=0; i<n; i++) {
        childQuadrant[i].quadrant->size = quadrant->size / 2;
        if (childQuadrant[i].quadrant->nElementsInQuadrant > maxLeafElems) {
            if (childQuadrant[i].quadrant->size > childQuadrant[i].quadrant->minElementSize) {
                if (*gen <= 8) {
                    *gen = *gen + 1;
                    [self createChildQuadrantFromMother:childQuadrant[i].quadrant mesh:aMesh dimension:dim maxLeafElements:maxLeafElems generation:gen];
                    *gen = *gen - 1;
                }
            }
        }
    }
    
    free_ivector(quadrant->elements, 0, quadrant->nElementsInQuadrant-1);
    quadrant->elements = NULL;
}

/*********************************************************************************
 
    Builds a tree hierarchy recursively bisectioning the geometry bounding box
    and partitioning the bulk element in the last level of the tree hierarchy
 
    FEMMesh *mesh  ->  finite element mesh
    double *box    ->  xmin, ymin, zmin, xmax, ymax, zmax
    Quandrant_t    ->  quadrant tree structure root
 
**********************************************************************************/
-(void)buildQuadrantTreeForMesh:(FEMMesh *)mesh model:(FEMModel *)aModel boundingBox:(double *)box rootQuadrant:(Quadrant_t *)quadrant {
    
    int dim, generation, i;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    Quadrant_t *motherQuadrant;
    int maxLeafElements;
    
    dim = aModel.dimension;
    
    if (dim == 3) {
        maxLeafElements = 16;
    } else {
        maxLeafElements = 8;
    }
    
    generation = 0;
    
    xmin = box[0];
    xmax = box[3];
    if (dim >= 2) {
        ymin = box[1];
        ymax = box[4];
    } else {
        ymin = 0.0;
        ymax = 0.0;
    }
    if (dim == 3) {
        zmin = box[2];
        zmax = box[5];
    } else {
        zmin = 0.0;
        zmax = 0.0;
    }
    
    // Create mother of all quadrants
     quadrant = (Quadrant_t*) malloc( sizeof(Quadrant_t) * 1 );
    
    quadrant->boundingBox[0] = xmin;
    quadrant->boundingBox[1] = ymin;
    quadrant->boundingBox[2] = zmin;
    quadrant->boundingBox[3] = xmax;
    quadrant->boundingBox[4] = ymax;
    quadrant->boundingBox[5] = zmax;
    
    quadrant->elements = intvec(0, mesh.numberOfBulkElements-1);
    quadrant->nElementsInQuadrant =  mesh.numberOfBulkElements;
    for (i=0; i<mesh.numberOfBulkElements; i++) {
        quadrant->elements[i] = i;
    }
    
    NSLog(@"FEMInterpolation:buildQuadrantTreeForMesh: Start...\n");
    motherQuadrant = quadrant;
    [self createChildQuadrantFromMother:motherQuadrant mesh:mesh dimension:dim maxLeafElements:maxLeafElements generation:&generation];
    NSLog(@"FEMInterpolation:buildQuadrantTreeForMesh: Ready.\n");
}

@end
