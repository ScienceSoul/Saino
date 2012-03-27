//
//  Lists.c
//  Saino
//
//  Created by Hakime Seddik on 07/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "Lists.h"

void addVariable(Variable_t *variables, char *name, int dofs, double *values, int *perm, int *output, int *secondary) {
    
    Variable_t *temp1, *temp2;
    
    temp1 = (Variable_t *)malloc(sizeof(Variable_t));
    
    temp1->Dofs = dofs;
    
    if (perm != NULL) {
        temp1->Perm = perm;
    } else {
        temp1->Perm = NULL;
    }
    
    temp1->NameLen = strlen(name);
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
        
        //Apend at the end of the list
        temp1->Next = NULL;
        temp2->Next = temp1;
    }
    
}

Variable_t* getVariable(Variable_t *variables, char *name) {
    
    int k;
    Variable_t *cur_ptr;
    Variable_t *var;
    
    cur_ptr = variables;
    var = NULL;
    
    k = strlen(name);
    
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
