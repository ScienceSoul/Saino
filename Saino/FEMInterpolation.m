//===----------------------------------------------------------------------===//
//  FEMInterpolation.m
//  Saino
//
//  Created by Seddik hakime on 27/09/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

#import "FEMInterpolation.h"
#import "FEMPElementMaps.h"
#import "FEMElementDescription.h"
#import "Utils.h"

@interface FEMInterpolation ()
-(Quadrant_t * __nullable)FEMInterpolation_findPointsQuadrantForPoint:(double * __nonnull)point dimension:(int)dim motherQuadrant:(Quadrant_t * __nonnull)mother;
@end

@implementation FEMInterpolation

#pragma mark Private methods

/************************************************************
 
    Method corresponds to Elmer from git on October 27 2015
 
************************************************************/
-(Quadrant_t * __nullable)FEMInterpolation_findPointsQuadrantForPoint:(double * __nonnull)point dimension:(int)dim motherQuadrant:(Quadrant_t * __nonnull)mother {
    
    int i, k;
    double bbox[6], eps3, *valuesPtr = NULL, maxVal;
    double const eps2 = 0.0;
    Quadrant_t *childQuadrant = NULL, *returnedQuadrant = NULL;
    
    // Loop over child quadrants
    for (i=0; i<(int)pow(2.0, dim); i++) {
        childQuadrant = mother->childQuadrants[i].quadrant;
        memcpy(bbox, childQuadrant->boundingBox, sizeof(childQuadrant->boundingBox));
        
        //********* Note: eps2 is set to 0 at the moment *********
        
        maxVal = -HUGE_VAL;
        k = 3;
        valuesPtr = &bbox[0];
        for (int j=0; j<3; j++) {
            if ( (*(valuesPtr + k) - *(valuesPtr + j)) > maxVal ) {
                maxVal = *(valuesPtr + k) - *(valuesPtr + j);
            }
            k++;
        }
        eps3 = eps2 * maxVal;
        for (int j=0; j<3; j++) {
            bbox[j] = bbox[j] - eps3;
        }
        for (int j=3; j<6; j++) {
            bbox[j] = bbox[j] + eps3;
        }
        
        // Is the point in childQuadrant[i]
        if ( (point[0] >= bbox[0]) && (point[0] <= bbox[3]) && (point[1] >= bbox[1]) && (point[1] <= bbox[4]) &&
            (point[2] >= bbox[2]) && (point[2] <= bbox[5]) ) break;
    }
    
    if ( i >= (int)pow(2.0, dim)) {
        fprintf(stdout, "FEMInterpolation:FEMInterpolation_findPointsQuadrantForPoint: point not found in any of the quadrants?\n");
        return NULL;
    }
    
    mother = childQuadrant;
    returnedQuadrant = childQuadrant;
    
    // Are we already in the leaf quadrant?
    // If not, search for the next generation child quadrants for the point
    if (mother->childQuadrants != NULL) {
        returnedQuadrant = [self FEMInterpolation_findPointsQuadrantForPoint:point dimension:dim motherQuadrant:mother];
    }
    
    return returnedQuadrant;
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

-(Quadrant_t * __nullable)findLeafElementsForPoint:(double * __nonnull)point dimension:(int)dim rootQuadrant:(Quadrant_t * __nonnull)root {
    
    Quadrant_t *leaf = root;
    
    // Find recursively the last generation quadrant the point belongs to
    return [self FEMInterpolation_findPointsQuadrantForPoint:point dimension:dim motherQuadrant:leaf];
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
 
    Method corresponds to Elmer from git on October 27 2015

*********************************************************************************************/
-(BOOL)isPointInElement:(Element_t * __nonnull)element elementNodes:(Nodes_t * __nonnull)nodes point:(double * __nonnull)point localCoordinates:(double * __nonnull)localCoords globalEpsilon:(double * __nullable)globaleps localEpsilon:(double * __nullable)localeps numericEpsilon:(double * __nullable)numericeps globalDistance:(double * __nullable)globaldist localDistance:(double * __nullable)localdist model:(FEMModel * __nonnull)model elementDescription:(FEMElementDescription * __nonnull)elementDescription elementMaps:(FEMPElementMaps * __nonnull)elementMaps edgeBasis:(BOOL * __nullable)edgeBasis {
    
    int n;
    double ug, vg, wg, xdist, ydist, zdist, sumdist, eps0, eps1, eps2, minx,
           maxx, miny, maxy, minz, maxz;
    BOOL isInElement;
    
    // Initialize
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
        
        vDSP_minvD(nodes->x, 1, &minx, n);
        vDSP_maxvD(nodes->x, 1, &maxx, n);
        
        vDSP_minvD(nodes->y, 1, &miny, n);
        vDSP_maxvD(nodes->y, 1, &maxy, n);

        vDSP_minvD(nodes->z, 1, &minz, n);
        vDSP_maxvD(nodes->z, 1, &maxz, n);
        
        xdist = max(max(point[0] - maxx, 0.0), minx - point[0]);
        ydist = max(max(point[1] - maxy, 0.0), miny - point[1]);
        zdist = max(max(point[2] - maxz, 0.0), minz - point[2]);
        
        *globaldist = sqrt(pow(xdist, 2.0) + pow(ydist, 2.0) + pow(zdist, 2.0));
        
        if (xdist > eps0 + eps1 * (maxx - minx) ) return isInElement;
        if (ydist > eps0 + eps1 * (maxy - miny) ) return isInElement;
        if (zdist > eps0 + eps1 * (maxz - minz) ) return isInElement;
    } else {
        // Otherwise make decision independently after each coordinate direction
        vDSP_minvD(nodes->x, 1, &minx, n);
        vDSP_maxvD(nodes->x, 1, &maxx, n);
        xdist = max(max(point[0] - maxx, 0.0), minx - point[0]);
        if (xdist > eps0 + eps1 * (maxx - minx)) return isInElement;
        
        vDSP_minvD(nodes->y, 1, &miny, n);
        vDSP_maxvD(nodes->y, 1, &maxy, n);
        ydist = max(max(point[1] - maxy, 0.0), miny - point[1]);
        if (ydist > eps0 + eps1 * (maxy - miny)) return isInElement;
        
        vDSP_minvD(nodes->z, 1, &minz, n);
        vDSP_maxvD(nodes->z, 1, &maxz, n);
        zdist = max(max(point[2] - maxz, 0.0), minz - point[2]);
        if (zdist > eps0 + eps1 * (maxz - minz)) return isInElement;
    }
    
    // Get element local coordinates from global coordinates of the point
    [elementDescription globalToLocalFromElement:element elementNodes:nodes localU:&ug localV:&vg localW:&wg x:point[0] y:point[1] z:point[2] model:model];
    
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
            fprintf(stdout, "FEMInterpolation:isPointInElement: element code %d not implemented.", element->Type.ElementCode);
            break;
    }
    
    if (sumdist < eps0 + eps2) {
        isInElement = YES;
    }
    
    if (localdist != NULL) {
        *localdist = sumdist;
    }
    
    BOOL trans = (edgeBasis != NULL) ? YES : NO;
    if (trans == YES) trans = *edgeBasis;
    trans = (trans == YES || [elementMaps isPElement:element] == YES) ? YES : NO;
    
    if (trans == YES) {
        switch (element->Type.ElementCode / 100) {
            case 3:
                ug = 2.0 * ug + vg - 1.0;
                vg = sqrt(3.0) * vg;
                break;
                
            case 5:
                ug = 2.0 * ug + vg + wg - 1.0;
                vg = sqrt(3.0) * vg + 1.0/sqrt(3.0) * wg;
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
    
    return isInElement;
}

/***************************************************************
    Loop over all elements in the MotherQuadrant, placing them
    in one of the 2^dim child quadrants.
***************************************************************/
-(void)putElementsInChildQuadrants:(QuadrantPointer_t * __nonnull)childQuadrants motherQuadrant:(Quadrant_t * __nonnull)mother mesh:(FEMMesh * __nonnull)mesh dimension:(int)dim {
    
    int i, j, t, n;
    int elementList[(int)pow(2, dim)][mother->nElementsInQuadrant];
    double eps3, maxval, minval;
    double const eps2 = 2.5e-2;
    double bbox[6], xmin, xmax, ymin, ymax, zmin, zmax, elementSize;
    Element_t *elements = NULL;
    Nodes_t *nodes = NULL;
    BOOL elementInQuadrant;
    
    for (i=0; i<(int)pow(2.0, dim); i++) {
        childQuadrants[i].quadrant->nElementsInQuadrant = 0;
        childQuadrants[i].quadrant->minElementSize = 1.0e20;
    }
    
    elements = mesh.getElements;
    nodes = mesh.getNodes;
    
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
        for (i=0; i<(int)pow(2.0, dim); i++) { // Loop over child quadrants
            
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
                zmax < bbox[2] || zmin > bbox[5]) elementInQuadrant = NO;
            
            if (elementInQuadrant == YES) {
                childQuadrants[i].quadrant->minElementSize = min(elementSize, childQuadrants[i].quadrant->minElementSize);
                
                // We allocate and store also the midlevels temporarily
                // (for the duration of the construction routine)
                elementList[i][childQuadrants[i].quadrant->nElementsInQuadrant] = mother->elements[t];
                childQuadrants[i].quadrant->nElementsInQuadrant++;
            }
        }
    }
    for (i=0; i<(int)pow(2.0, dim); i++) {
        if (childQuadrants[i].quadrant->nElementsInQuadrant != 0) {
            childQuadrants[i].quadrant->elements = intvec(0, childQuadrants[i].quadrant->nElementsInQuadrant-1);
            
            for (j=0; j<childQuadrants[i].quadrant->nElementsInQuadrant; j++) {
                childQuadrants[i].quadrant->elements[j] = elementList[i][j];
            }
        }
    }
}

-(void)createChildQuadrantFromMother:(Quadrant_t * __nonnull)quadrant mesh:(FEMMesh * __nonnull)mesh dimension:(int)dim maxLeafElements:(int)maxLeafElems generation:(int * __nonnull)gen {
    
    int i, n;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    
    // Create 2^dim child quadrants
    n = (int)pow(2.0, dim);
    quadrant->childQuadrants = (QuadrantPointer_t*) malloc( sizeof(QuadrantPointer_t) * n );
    quadrant->numberOfchildQuadrants = n;
    for (i=0; i<n; i++) {
        quadrant->childQuadrants[i].quadrant = (Quadrant_t*) malloc( sizeof(Quadrant_t));
        quadrant->childQuadrants[i].quadrant->nElementsInQuadrant = 0;
        quadrant->childQuadrants[i].quadrant->elements = NULL;
        quadrant->childQuadrants[i].quadrant->childQuadrants = NULL;
    }
    
    xmin = quadrant->boundingBox[0];
    ymin = quadrant->boundingBox[1];
    zmin = quadrant->boundingBox[2];
    
    xmax = quadrant->boundingBox[3];
    ymax = quadrant->boundingBox[4];
    zmax = quadrant->boundingBox[5];
    
    quadrant->size = max(max(xmax-xmin, ymax-ymin), zmax-zmin);
    
    quadrant->childQuadrants[0].quadrant->boundingBox[0] = xmin;
    quadrant->childQuadrants[0].quadrant->boundingBox[1] = ymin;
    quadrant->childQuadrants[0].quadrant->boundingBox[2] = zmin;
    quadrant->childQuadrants[0].quadrant->boundingBox[3] = (xmin + xmax) / 2.0;
    quadrant->childQuadrants[0].quadrant->boundingBox[4] = (ymin + ymax) / 2.0;
    quadrant->childQuadrants[0].quadrant->boundingBox[5] = (zmin + zmax) / 2.0;
    
    quadrant->childQuadrants[1].quadrant->boundingBox[0] = (xmin + xmax) / 2.0;
    quadrant->childQuadrants[1].quadrant->boundingBox[1] = ymin;
    quadrant->childQuadrants[1].quadrant->boundingBox[2] = zmin;
    quadrant->childQuadrants[1].quadrant->boundingBox[3] = xmax;
    quadrant->childQuadrants[1].quadrant->boundingBox[4] = (ymin + ymax) / 2.0;
    quadrant->childQuadrants[1].quadrant->boundingBox[5] = (zmin + zmax) / 2.0;
    
    if (dim >= 2) {
        quadrant->childQuadrants[2].quadrant->boundingBox[0] = xmin;
        quadrant->childQuadrants[2].quadrant->boundingBox[1] = (ymin + ymax) / 2.0;
        quadrant->childQuadrants[2].quadrant->boundingBox[2] = zmin;
        quadrant->childQuadrants[2].quadrant->boundingBox[3] = (xmin + xmax) / 2.0;
        quadrant->childQuadrants[2].quadrant->boundingBox[4] = ymax;
        quadrant->childQuadrants[2].quadrant->boundingBox[5] = (zmin + zmax) / 2.0;
        
        quadrant->childQuadrants[3].quadrant->boundingBox[0] = (xmin + xmax) / 2.0;
        quadrant->childQuadrants[3].quadrant->boundingBox[1] = (ymin + ymax) / 2.0;
        quadrant->childQuadrants[3].quadrant->boundingBox[2] = zmin;
        quadrant->childQuadrants[3].quadrant->boundingBox[3] = xmax;
        quadrant->childQuadrants[3].quadrant->boundingBox[4] = ymax;
        quadrant->childQuadrants[3].quadrant->boundingBox[5] = (zmin + zmax) / 2.0;
    }
    
    if (dim == 3) {
        quadrant->childQuadrants[4].quadrant->boundingBox[0] = xmin;
        quadrant->childQuadrants[4].quadrant->boundingBox[1] = ymin;
        quadrant->childQuadrants[4].quadrant->boundingBox[2] = (zmin + zmax) / 2.0;
        quadrant->childQuadrants[4].quadrant->boundingBox[3] = (xmin + xmax) / 2.0;
        quadrant->childQuadrants[4].quadrant->boundingBox[4] = (ymin + ymax) / 2.0;
        quadrant->childQuadrants[4].quadrant->boundingBox[5] = zmax;
        
        quadrant->childQuadrants[5].quadrant->boundingBox[0] = (xmin + xmax) / 2.0;
        quadrant->childQuadrants[5].quadrant->boundingBox[1] = ymin;
        quadrant->childQuadrants[5].quadrant->boundingBox[2] = (zmin + zmax) / 2.0;
        quadrant->childQuadrants[5].quadrant->boundingBox[3] = xmax;
        quadrant->childQuadrants[5].quadrant->boundingBox[4] = (ymin + ymax) / 2.0;
        quadrant->childQuadrants[5].quadrant->boundingBox[5] = zmax;
        
        quadrant->childQuadrants[6].quadrant->boundingBox[0] = xmin;
        quadrant->childQuadrants[6].quadrant->boundingBox[1] = (ymin + ymax) / 2.0;
        quadrant->childQuadrants[6].quadrant->boundingBox[2] = (zmin + zmax) / 2.0;
        quadrant->childQuadrants[6].quadrant->boundingBox[3] = (xmin + xmax) / 2.0;
        quadrant->childQuadrants[6].quadrant->boundingBox[4] = ymax;
        quadrant->childQuadrants[6].quadrant->boundingBox[5] = zmax;
        
        quadrant->childQuadrants[7].quadrant->boundingBox[0] = (xmin + xmax) / 2.0;
        quadrant->childQuadrants[7].quadrant->boundingBox[1] = (ymin + ymax) / 2.0;
        quadrant->childQuadrants[7].quadrant->boundingBox[2] = (zmin + zmax) / 2.0;
        quadrant->childQuadrants[7].quadrant->boundingBox[3] = xmax;
        quadrant->childQuadrants[7].quadrant->boundingBox[4] = ymax;
        quadrant->childQuadrants[7].quadrant->boundingBox[5] = zmax;
    }
    
    // Loop over all elements in the mother quadrant, placing them in one of the 2^dim quadrants
    [self putElementsInChildQuadrants:quadrant->childQuadrants motherQuadrant:quadrant mesh:mesh dimension:dim];
    
    // Check whether we need to branch for the next level
    for (i=0; i<n; i++) {
        quadrant->childQuadrants[i].quadrant->size = quadrant->size / 2;
        if (quadrant->childQuadrants[i].quadrant->nElementsInQuadrant > maxLeafElems) {
            if (quadrant->childQuadrants[i].quadrant->size > quadrant->childQuadrants[i].quadrant->minElementSize) {
                if (*gen < 8) {
                    *gen = *gen + 1;
                    [self createChildQuadrantFromMother:quadrant->childQuadrants[i].quadrant mesh:mesh dimension:dim maxLeafElements:maxLeafElems generation:gen];
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
-(void)buildQuadrantTreeForMesh:(FEMMesh * __nonnull)mesh model:(FEMModel * __nonnull)model boundingBox:(double * __nonnull)box rootQuadrant:(Quadrant_t * __nullable)quadrant {
    
    int dim, generation, i;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    int maxLeafElements;
    
    dim = model.dimension;
    
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
    
    fprintf(stdout, "FEMInterpolation:buildQuadrantTreeForMesh: Start...\n");
    [self createChildQuadrantFromMother:quadrant mesh:mesh dimension:dim maxLeafElements:maxLeafElements generation:&generation];
    fprintf(stdout, "FEMInterpolation:buildQuadrantTreeForMesh: Ready.\n");
}

@end
