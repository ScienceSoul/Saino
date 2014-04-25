//
//  FEMRadiation.m
//  Saino
//
//  Created by Seddik hakime on 03/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMRadiation.h"

#import "FEMCore.h"
#import "FEMElementUtils.h"
#import "FEMListUtilities.h"
#import "FEMBoundaryCondition.h"
#import "Utils.h"

@implementation FEMRadiation

- (id)init
{
    self = [super init];
    if (self) {
        //TODO: Initialize here
        
    }
    
    return self;
}

-(double)computeRadiationLoadModel:(FEMModel *)model mesh:(FEMMesh *)mesh element:(Element_t *)elememt temperature:(double *)temperature reorder:(int *)reorder emissivity:(double)emissivity angleFraction:(double *)angleFraction {
    
    int i, j, *cols, n;
    double a1, a2, asum, emissivity1, sum, t, *vals;
    BOOL found;
    FEMCore *core;
    FEMElementUtils *elementUtils;
    FEMListUtilities *listUtilities;
    FEMBoundaryCondition *boundaryConditionAtID = nil;
    Element_t *currentElement = NULL, *parent = NULL;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    core = [FEMCore sharedCore];
    listUtilities = [[FEMListUtilities alloc] init];
    
    elementUtils = [[FEMElementUtils alloc] init];
    a1 = emissivity * [elementUtils elementArea:elememt numberOfNodes:elememt->Type.NumberOfNodes mesh:mesh nodel:model];
    
    cols = elememt->BoundaryInfo->GebhardtFactors->Elements;
    vals = elememt->BoundaryInfo->GebhardtFactors->Factors;
    
    t = 0.0;
    asum = 0.0;
    
    currentElement = mesh.getElements;
    for (i=0; i<elememt->BoundaryInfo->GebhardtFactors->NumberOfFactors; i++) {
        n = currentElement[cols[i]].Type.NumberOfNodes;
        
        boundaryConditionAtID = (model.boundaryConditions)[currentElement[cols[i]].BoundaryInfo->Constraint-1];
        found = [listUtilities listGetReal:model inArray:boundaryConditionAtID.valuesList forVariable:@"emissivity" numberOfNodes:n indexes:currentElement[cols[i]].NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
        sum = 0.0;
        if (found == YES) {
            for (j=0; j<buffer.m; j++) {
                sum = sum + buffer.vector[j];
            }
            emissivity1 = sum / n;
        } else {
            found = [core getParentMaterialProperty:@"emissivity" forElement:&currentElement[cols[i]] parentElement:parent model:model listUtilities:listUtilities buffer:&buffer];
            for (j=0; j<buffer.m; j++) {
                sum = sum + buffer.vector[j];
            }
            emissivity1 = sum / n;
        }
        
        a2 = emissivity1 * [elementUtils elementArea:&currentElement[cols[i]] numberOfNodes:currentElement[cols[i]].Type.NumberOfNodes mesh:mesh nodel:model];
        
        sum = 0.0;
        for (j=0; j<currentElement[cols[i]].Type.NumberOfNodes; j++) {
            sum = sum + temperature[reorder[currentElement[cols[i]].NodeIndexes[j]]] / n;
        }
        t = t + a2 * fabs(vals[i]) * pow(sum, 4.0);
        asum = asum + a2 * fabs(vals[i]);
    }
    
    t = pow((t/a1), (1.0/4.0));
    
    if (angleFraction != NULL) {
        *angleFraction = asum / a1;
    }
    
    if (buffer.vector != NULL) {
        free_dvector(buffer.vector, 0, buffer.m-1);
        buffer.vector = NULL;
    }

    return t;
}

-(double)computeRadiationCoeffModel:(FEMModel *)model mesh:(FEMMesh *)mesh element:(Element_t *)element index:(int)k {
    
    int i, n;
    double area, emissivity, sum, t;
    BOOL found;
    FEMCore *core;
    FEMElementUtils *elementUtils;
    FEMListUtilities *listUtilities;
    FEMBoundaryCondition *boundaryConditionAtID = nil;
    Element_t *currentElements = NULL, *parent = NULL;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    
    core = [FEMCore sharedCore];
    listUtilities = [[FEMListUtilities alloc] init];
    
    currentElements = model.getElements;
    n = currentElements[element->BoundaryInfo->GebhardtFactors->Elements[k]].Type.NumberOfNodes;
    
    boundaryConditionAtID = (model.boundaryConditions)[currentElements[element->BoundaryInfo->GebhardtFactors->Elements[k]].BoundaryInfo->Constraint-1];
    found = [listUtilities listGetReal:model inArray:boundaryConditionAtID.valuesList forVariable:@"emissivity" numberOfNodes:n indexes:currentElements[element->BoundaryInfo->GebhardtFactors->Elements[k]].NodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
    sum = 0.0;
    if (found == YES) {
        for (i=0; i<buffer.m; i++) {
            sum = sum + buffer.vector[i];
        }
        emissivity = sum / n;
    } else {
        found = [core getParentMaterialProperty:@"emissivity" forElement:&currentElements[element->BoundaryInfo->GebhardtFactors->Elements[k]] parentElement:parent model:model listUtilities:listUtilities buffer:&buffer];
        for (i=0; i<buffer.m; i++) {
            sum = sum + buffer.vector[i];
        }
        emissivity = sum / n;
    }
    
    elementUtils = [[FEMElementUtils alloc] init];
    area = emissivity * [elementUtils elementArea:&currentElements[element->BoundaryInfo->GebhardtFactors->Elements[k]] numberOfNodes:n mesh:mesh nodel:model];
    t = fabs(element->BoundaryInfo->GebhardtFactors->Factors[k]) * area;
    
    if (buffer.vector != NULL) {
        free_dvector(buffer.vector, 0, buffer.m-1);
        buffer.vector = NULL;
    }
    
    return t;
}

@end
