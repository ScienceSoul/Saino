//===----------------------------------------------------------------------===//
//  Lists.c
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

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "Lists.h"

void addVariable(Variable_t * __nullable variables, char * __nonnull name, int dofs, double * __nonnull values, int * __nullable perm, int * __nullable output, int * __nullable secondary) {
    
    Variable_t *temp1, *temp2;
    
    temp1 = (Variable_t *)malloc(sizeof(Variable_t));
    
    temp1->Dofs = dofs;
    
    if (perm != NULL) {
        temp1->Perm = perm;
    } else {
        temp1->Perm = NULL;
    }
    
    temp1->NameLen = (int)strlen(name);
    strcpy(temp1->Name, name);
    
    temp1->Norm = 0.0;
    temp1->PrevNorm = 0.0;
    temp1->Values = values;
    temp1->PrevValues = NULL;
    temp1->EigenValues = NULL;
    temp1->EigenVectors = NULL;
    temp1->NonLinChange = 0.0;
    temp1->SteadyChange = 0.0;
    temp1->NonLinValues = NULL;
    temp1->SteadyValues = NULL;
    
    temp1->Valid = 1;
    temp1->Output = 1;
    temp1->Secondary = 0;
    temp1->ValuesChanged = 1;
    
    // Converged information undefined = -1, not = 0, yes = 1
    temp1->NonLinConverged = -1;
    temp1->SteadyConverged = -1;
    
    if (secondary != NULL) {
        temp1->Secondary = *secondary;
    }
    if (output != NULL) {
        temp1->Output = *output;
    }
    
    // Copying the head location into another node
    temp2 = variables;
    
    if (variables == NULL) {
        
        // If list is empty we create first node
        variables = temp1;
        variables->Next = NULL;
    } else {
        
        // Traverse down to end of the list
        while (temp2->Next != NULL) temp2 = temp2->Next;
        if (strcmp(temp2->Name, temp1->Name) == 0) {
            free(temp1);
            temp1 = NULL;
            return;
        }
        
        //Append at the end of the list
        temp1->Next = NULL;
        temp2->Next = temp1;
    }
    
}

Variable_t * __nullable getVariable(Variable_t * __nullable variables, char * __nonnull name) {
    
    int k;
    Variable_t *cur_ptr;
    Variable_t *var;
    
    cur_ptr = variables;
    var = NULL;
    
    k = (int)strlen(name);
    
    if (cur_ptr == NULL) {
        printf("List is empty.\n");
    } else {
        while (cur_ptr != NULL) {
            if (cur_ptr->NameLen == k) {
                if (strcmp(cur_ptr->Name, name) == 0) {
                    if (cur_ptr->Valid == 1) {
                        var = cur_ptr;
                        return var;
                    }
                }
            }
            cur_ptr= cur_ptr->Next;
        }
    }
    return var;
}
