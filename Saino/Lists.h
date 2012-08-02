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

void addVariable(Variable_t *variables, char *name, int dofs, double *values, int *perm, int *output, int *secondary);
Variable_t* getVariable(Variable_t *variables, char *name);
