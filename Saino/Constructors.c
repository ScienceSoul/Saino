//===----------------------------------------------------------------------===//
//  Constructors.c
//  Saino
//
//  Created by Seddik hakime on 25/01/13.
//  Copyright (c) 2013 ScienceSoul. All rights reserved.
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

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include "Constructors.h"

void initNodes(Nodes_t * __nonnull nodes) {
    
    *nodes = (Nodes_t){.numberOfNodes=0, .x=NULL, .y=NULL, .z=NULL};
}

void initElements(Element_t * __nonnull elements, int n) {
    
    int i;
    
    for (i=0; i<n; i++) {
        elements[i] = (Element_t){.color.red=-1, .color.green=-1, .color.blue=-1, .color.colorIndex=-1, .ElementIndex=-1, .GElementIndex=-1,
                                  .colored=false, .copy=false, .BodyID=0, .Splitted=0, .NDOFs=-1, .BDOFs=-1, .DGDOFs=-1, .BoundaryInfo=NULL,
                                  .NodeIndexes=NULL, .EdgeIndexes=NULL, .FaceIndexes=NULL, .BubbleIndexes=NULL, .DGIndexes=NULL, .Pdefs=NULL};
    }
}

void initBoundaryInfo(BoundaryInfo_t * __nonnull boundaryInfo) {
    
    *boundaryInfo = (BoundaryInfo_t){.Constraint=0, .Outbody=-1, .GebhardtFactors=NULL, .Left=NULL, .Right=NULL};
}

variableArraysContainer * __nonnull allocateVariableContainer(void) {
    
    variableArraysContainer *varContainers = NULL;
    varContainers = (variableArraysContainer*)malloc(sizeof(variableArraysContainer));
    
    *varContainers = (variableArraysContainer){.Values = NULL, .SecondaryToValues = NULL, .ComponentValues = NULL, .ComponentSecondaryToValues = NULL, .Perm = NULL};
    
    return varContainers;
}

RungeKutta_t * __nonnull allocateRungeKutta(int n) {
    
    RungeKutta_t *rungeKutta = NULL;
    
    rungeKutta = (RungeKutta_t*)malloc(sizeof(RungeKutta_t) * n );
    *rungeKutta = (RungeKutta_t){.k1 = NULL, .k2 = NULL, .k3 = NULL, .k4 = NULL};
    
    return rungeKutta;
}
