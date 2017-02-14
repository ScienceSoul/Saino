//===----------------------------------------------------------------------===//
//  Lists.h
//  Saino
//
//  Created by Hakime Seddik on 07/03/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
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

#include "Constructors.h"

typedef struct Variable_t {
    
    struct Variable_t *Next;
    int NameLen;
    char *Name;
    
    int Valid, Output, ValuesChanged, Secondary;          // Those are logicals: 0 -> False; 1 -> True
    
    int Dofs;
    int *Perm;
    double Norm, PrevNorm, NonLinChange, SteadyChange;
    int NonLinConverged, SteadyConverged;
    int NonLinIter;
    double *Values, **PrevValues, *PValues, *NonLinValues, *SteadyValues;
    double complex *EigenValues, **EigenVectors;
    int sizePerm, sizeValues, sizePrevValues, sizePValues, sizeNonLinValues, sizeSteadyValues;
    
} Variable_t;

void addVariable(Variable_t * __nullable variables, char * __nonnull name, int dofs, double * __nonnull values, int * __nullable perm, int * __nullable output, int * __nullable secondary);
Variable_t * __nullable getVariable(Variable_t * __nullable variables, char * __nonnull name);
