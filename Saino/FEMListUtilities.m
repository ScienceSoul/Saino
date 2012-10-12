//
//  FEMListUtilities.m
//  Saino
//
//  Created by Hakime Seddik on 18/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMListUtilities.h"

@implementation FEMListUtilities

-(void)listParseStrToValues:(FEMModel *)model: (NSString *)str: (int)ind: (NSString *)name: (double *)t: (int)count {

    int i, l, k1, n;
    FEMVariable *variable, *cVar;
    variableArraysContainer *varContainers, *cvarContainers;
    FEMUtilities *utilities;
    Element_t *element;
    BOOL found;
    
    utilities = [[FEMUtilities alloc] init];
    
    //TODO: This method only support one string passed at a time, not several strings separeted by a comma. 
    
    count = 0;
    
    if ([str isEqualToString:@"Coordinate"] == NO) {
        variable = [utilities getVariableFrom:model.variables model:model name:str onlySearch:NULL maskName:NULL info:&found];
        if (variable == nil) {
            warnfunct("listParseStrToValues", "Can't find indpendent variable:");
            printf("%s\n", [str UTF8String]);
            errorfunct("listParseStrToValues", "Abort...");
        }
    } else {
        variable = [utilities getVariableFrom:model.variables model:model name:@"Coordinate 1" onlySearch:NULL maskName:NULL info:&found];
    }
    
    varContainers = variable.getContainers;
    
    k1 = ind;
    if (variable.type == VARIABLE_ON_NODES_ON_ELEMENTS) {
        element = model.getCurrentElement;
        if (element != NULL) {
            if (element->DGIndexes != NULL) {
                n = element->Type.NumberOfNodes;
                if (element->sizeDGIndexes == n) {
                    for (i=0; i<n; i++) {
                        if (element->NodeIndexes[i] == ind) {
                            k1 = element->DGIndexes[i];
                        }
                    }
                }
            }
        }
    }
    if (varContainers->Perm != NULL) k1 = varContainers->Perm[k1];
    
    if (k1 > 0 && k1 <= varContainers->sizeValues) {
        if ([str isEqualToString:@"Coordinate"] == YES) {
            cVar = [utilities getVariableFrom:model.variables model:model name:@"Coordinate 1" onlySearch:NULL maskName:NULL info:&found];
            count++;
            cvarContainers = cVar.getContainers;
            t[0] = cvarContainers->Values[k1];
            
            cVar = [utilities getVariableFrom:model.variables model:model name:@"Coordinate 2" onlySearch:NULL maskName:NULL info:&found];
            count++;
            cvarContainers = cVar.getContainers;
            t[1] = cvarContainers->Values[k1];
            
            cVar = [utilities getVariableFrom:model.variables model:model name:@"Coordinate 3" onlySearch:NULL maskName:NULL info:&found];
            count++;
            cvarContainers = cVar.getContainers;
            t[2] = cvarContainers->Values[k1];
        }
        else {
            if (variable.dofs == 1) {
                t[count] = varContainers->Values[k1];
                count++;
            } else {
                for (l=0; l<variable.dofs; l++) {
                    t[count] = varContainers->Values[variable.dofs*(k1-1)+l];
                    count++;
                }
            }
        }
    } else {
        if (varContainers->Perm != NULL) {
            t[count] = HUGE_VAL;
        } else {
            t[count] = varContainers->Values[0];
        }
        count++;
    }
    
    variable = nil;
    cVar = nil;
    varContainers = NULL;
    cvarContainers = NULL;
}

-(NSString *)listGetString:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found {
    
    NSString *s;
    
    *found = NO;
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES) {
            s = [NSString stringWithString:list.cValue];
            break;
            *found = YES;
        }
    }
    
    return s;
}

-(BOOL)listGetReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes buffer:(listBuffer *)result minValue:(double *)minv maxValue:(double *)maxv {

    int i, j, k;
    double *buffer;
    double t[32];
    FEMUtilities *util;
    BOOL found;
    valueListArraysContainer *containers;
    
    found = NO;
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            if (list.type == LIST_TYPE_CONSTANT_SCALAR) {
                
                // TODO: implement the call of a user provided method if required
                 
                if (containers->fValues == NULL) {
                    warnfunct("listGetReal", "fValues not allocated in list:\n");
                    printf("%s\n", [varName UTF8String]);
                    errorfunct("listGetReal", "Program terminating now...");
                }
                // The buffer is allocated here, it's up to the caller to release this memory
                result->vector = doublevec(0, n-1);
                result->m = n;
                memset( result->vector, 0.0, (n*sizeof(result->vector)) );
                for (i=0; i<n; i++) {
                    result->vector[i] = containers->fValues[0][0][0];
                }
            } else if (list.type == LIST_TYPE_VARIABLE_SCALAR) {
                
                util = [[FEMUtilities alloc] init];
                buffer = doublevec(0, containers->sizeTValues-1);
                
                for (i=0; i<n; i++) {
                    buffer[i] = containers->fValues[0][0][i];
                }
                
                // The buffer is allocated here but it's up to the caller to release this memory
                result->vector = doublevec(0, n-1);
                result->m = n;
                memset( result->vector, 0.0, (n*sizeof(result->vector)) );
                for (i=0; i<n; i++) {
                    k = nodeIndexes[i];
                    [self listParseStrToValues:model :list.dependName :k :list.name :t :j];
                    
                    if (any(t, '=', HUGE_VAL, j) == 0) {
                        
                        // TODO: implement the call of a user provided method if required
                        
                        if (containers->fValues == NULL) {
                            warnfunct("listGetReal", "fValues not allocated in list:\n");
                            printf("%s\n", [varName UTF8String]);
                            errorfunct("listGetReal", "Program terminating now...");
                        }
                        result->vector[i] = [util interpolateCurve:containers->tValues :buffer :t[0] :containers->sizeTValues];
                    }
                }
                
                free_dvector(buffer, 0, containers->sizeTValues-1);
            } else if ([list type] == LIST_TYPE_CONSTANT_SCALAR_STR) {
                 // TODO: implement this case
            } else if ([list type] == LIST_TYPE_VARIABLE_SCALAR_STR) {
                 // TODO: implement this case
            }
            
            found = YES;
            containers = NULL;
            break;
        }
    }
    
    if (minv != NULL) {
        if (min_array(result->vector, n) < *minv) {
            warnfunct("listGetReal", "Value smaller than given value:");
            printf("Value: %f / Given value: %f / Property: %s\n", min_array(result->vector, n), *minv, [varName UTF8String]);
            errorfunct("listGetReal", "Program terminating now...");
        }
    }
    
    if (maxv != NULL) {
        if (max_array(result->vector, n) > *maxv) {
            warnfunct("listGetReal", "Value greater than given value:");
            printf("Value: %f / Given value: %f / Property: %s\n", max_array(result->vector, n), *minv, [varName UTF8String]);
            errorfunct("listGetReal", "Program terminating now...");
        }
    }
    
    return found;
}

-(BOOL)listGetRealArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes buffer:(listBuffer *)result {

    int i, j, k, n1, n2;
    BOOL found;
    listBuffer buffer;
    valueListArraysContainer *containers;
    
    found = NO;
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            n1 = containers->sizeFValues1;
            n2 = containers->sizeFValues2;
            
            if (list.type == LIST_TYPE_CONSTANT_TENSOR) {
                
                // The buffer is allocated here but it's up to the caller to release this memory
                result->tensor = d3tensor(0, n1-1, 0, n2-1, 0, n-1);
                result->m = n1;
                result->n = n2;
                result->p = n;
                for (i=0; i<n1; i++) {
                    for (j=0; j<n2; j++) {
                        for (k=0; k<n; k++) {
                            result->tensor[i][j][k] = 0.0;
                        }
                    }
                }
                
                for (i=0; i<n1; i++) {
                    for (j=0; j<n2; j++) {
                        for (k=0; k<n; k++) {
                            result->tensor[i][j][k] = containers->fValues[i][j][0];
                        }
                    }
                }
                
                // TODO: implement the call of a user provided method if required
            
            } else if (list.type == LIST_TYPE_VARIABLE_TENSOR || list.type == LIST_TYPE_VARIABLE_TENSOR_STR){
                
                // TODO: implement this case
                
            } else {
                     
                // The buffer is allocated here but it's up to the caller to release this memory
                result->tensor = d3tensor(0, n1-1, 0, n2-1, 0, n-1);
                result->m = n1;
                result->n = n2;
                result->p = n;
                for (i=0; i<n1; i++) {
                    for (j=0; j<n2; j++) {
                        for (k=0; k<n; k++) {
                            result->tensor[i][j][k] = 0.0;
                        }
                    }
                }
                
                for (i=0; i<n1; i++) {
                    for (j=0; j<n2; j++) {
                        
                        for (k=0; k<n; k++) {
                            result->tensor[i][j][k] = 0.0;
                        }
                    }
                }
                
                found = [self listGetReal:model inArray:array forVariable:varName numberOfNodes:n indexes:nodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
                for (i=0; i<n1; i++) {
                    for (k=0; k<n; k++) {
                        result->tensor[i][0][k] = buffer.vector[k];
                    }
                }
                
                free_dvector(buffer.vector, 0, buffer.m-1);
            }
            
            found = YES;
            containers = NULL;
            break; 
        }
    }
    
    return found;
}


-(double)listGetConstReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found minValue:(double *)minv maxValue:(double *)maxv {

    double f;
    valueListArraysContainer *containers;
    
    *found = NO;
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            
            if ([list type] >= 8) {
                
                // TODO: unimplemented case where type is greater than 8
                // In Elmer this calls matc. Figure out later what we do in Saino
                
            } else if (list.method == YES) {
                
                // TODO: implement the call of a user provided method
                
            } else {
                
                f = containers->fValues[0][0][0];
            }
            *found = YES;
            containers = NULL;
            break;
        }
    }
    
    if (minv != NULL) {
        if (f < *minv) {
            warnfunct("listGetConstReal", "Value smaller than given value:");
            printf("Value: %f / Given value: %f / Property: %s\n", f, *minv, [varName UTF8String]);
            errorfunct("listGetConstReal", "Program terminating now...");
        }
    }
    
    if (maxv != NULL) {
        if (f > *maxv) {
            warnfunct("listGetConstReal", "Value greater than given value:");
            printf("Value: %f / Given value: %f / Property: %s\n", f, *minv, [varName UTF8String]);
            errorfunct("listGetConstReal", "Program terminating now...");
        }
    }
    
    return f;
}

-(BOOL)listGetConstRealArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName buffer:(listBuffer *)result {

    int i, j, n1, n2;
    BOOL found;
    valueListArraysContainer *containers;
    
    found = NO;
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            n1 = containers->sizeFValues1;
            n2 = containers->sizeFValues2;
            
            // The buffer is allocated here but it's up to the caller to release this memory
            result->matrix = doublematrix(0, n1-1, 0, n2-1);
            result->m = n1;
            result->n = n2;
            for (i=0; i<n1; i++) {
                for (j=0; j<n2; j++) {
                    result->matrix[i][j] = 0.0;
                }
            }
            
            for (i=0; i<n1; i++) {
                for (j=0; j<n2; j++) {
                    result->matrix[i][j] = containers->fValues[i][j][0];
                }
            }
            
            // TODO: implement the call of a user provided method
            
            found = YES;
            containers = NULL;
            break;
        }
    }
    
    return found;
}

-(BOOL)listGetIntegerArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName buffer:(listBuffer *)result {
    
    int i, n;
    BOOL found;
    valueListArraysContainer *containers;
    
    found = NO;
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            if (containers->iValues == NULL) {
                warnfunct("listGetIntegerArray", "iValues not allocated in list:\n");
                printf("%s\n", [varName UTF8String]);
                errorfunct("listGetIntegerArray", "Program terminating now...");
            }
            
           // The buffer is allocated here but it's up to the caller to release this memory
            n = containers->sizeIValues;
            result->ivector = intvec(0, n-1);
            result->m = n;
            memset( result->ivector, 0, (n*sizeof(result->ivector)) );
            for (i=0; i<n; i++) {
                result->ivector[i] = containers->iValues[i];
            }
            
           // TODO: implement the call of a user provided method if required 
            
            found = YES;
            containers = NULL;
            break;
        }
    }
    
    return found;
}

-(int)listGetInteger:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found minValue:(int *)minv maxValue:(int *)maxv {
    
    int l;
    valueListArraysContainer *containers;
    
    *found = NO;
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            
            // TODO: implement the call of a user provided method if required
            
            if (containers->iValues == NULL) {
                warnfunct("listGetReal", "iValues not allocated in list:\n");
                printf("%s\n", [varName UTF8String]);
                errorfunct("listGetReal", "Program terminating now...");
            }
            
            l = containers->iValues[0];
            
            *found = YES;
            containers = NULL;
            break;
            
        }
    }
    
    if (minv != NULL) {
        if (l < *minv) {
            warnfunct("listGetInteger", "Value smaller than given value:");
            printf("Value: %d / Given value: %d / Property: %s\n", l, *minv, [varName UTF8String]);
            errorfunct("listGetInteger", "Program terminating now...");
        }
    }
    
    if (maxv != NULL) {
        if (l > *maxv) {
            warnfunct("listGetInteger", "Value greater than given value:");
            printf("Value: %d / Given value: %d / Property: %s\n", l, *minv, [varName UTF8String]);
            errorfunct("listGetInteger", "Program terminating now...");
        }
    }

    return l;
}

-(BOOL)listGetLogical:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found {
    
    BOOL l;
    
    *found = NO;
    l = NO;
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES) {
            l = list.isLvalue;
            break;
            *found = YES;
        }
    }
    
    return l;
}

-(FEMValueList *)listFindVariable:(NSString *)varName inArray:(NSArray *)array {

    for (FEMValueList *list in array) {
        if ([varName isEqualToString:[list name]] == YES) {
            return list;
        }
    }
    
    return nil;
}

-(BOOL)listCheckPresentVariable:(NSString *)varName inArray:(NSArray *)array {

    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES) {
            return YES;
        }
    }
    
    return NO;
}

-(void)addIntegerArrayInClassList:(id)className theVariable:(NSString *)varName withValues:(int *)values numberOfNodes:(int)n {
    
    int i;
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMValueList *newValueList;
    NSMutableArray *valuesArray;
    BOOL found;
    valueListArraysContainer *containers;
    
    found = NO;
    
    if ([className isKindOfClass:[FEMBodyForce class]]) {
        bodyForce = className;
        valuesArray = bodyForce.valuesList;
    } else if ([className isKindOfClass:[FEMBoundaryCondition class]]) {
        boundary = className;
        valuesArray = boundary.valuesList;
    } else if ([className isKindOfClass:[FEMSimulation class]]) {
        simulation = className;
        valuesArray = simulation.valuesList;
    }
    
    for (FEMValueList *list in valuesArray) {
        if ([varName isEqualToString:list.name] == YES) {
            found = YES;
            break;
        }
    }
    
    if (found == NO) {
        
        newValueList = [[FEMValueList alloc] init];
        newValueList.type = LIST_TYPE_CONSTANT_TENSOR;
        
        containers = newValueList.getContainers;
        containers->iValues = intvec(0, n-1);
        for (i=0; i<n; i++) {
            containers->iValues[i] = values[i];
        }
        
        newValueList.name = varName;
        [valuesArray addObject:newValueList];
        containers = NULL;
    }
    
    boundary = nil;
    bodyForce = nil;
    simulation = nil;
}

@end
