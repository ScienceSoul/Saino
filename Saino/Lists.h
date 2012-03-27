//
//  Lists.h
//  Saino
//
//  Created by Hakime Seddik on 07/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#include "Constructors.h"

void addVariable(Variable_t *variables, char *name, int dofs, double *values, int *perm, int *output, int *secondary);
Variable_t* getVariable(Variable_t *variables, char *name);
