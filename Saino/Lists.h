//
//  Lists.h
//  Saino
//
//  Created by Hakime Seddik on 07/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

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
