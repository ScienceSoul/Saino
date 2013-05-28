//
//  FEMListUtilities.m
//  Saino
//
//  Created by Hakime Seddik on 18/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMListUtilities.h"
#import "FEMBodyForce.h"
#import "FEMBoundaryCondition.h"
#import "FEMSimulation.h"
#import "FEMEquation.h"
#import "FEMUtilities.h"
#import "Utils.h"

@implementation FEMListUtilities {
    
    NSMutableString *_nameSpace;
}

- (id)init
{
    self = [super init];
    if (self) {
        _nameSpace = [NSMutableString stringWithString:@""];
    }
    
    return self;
}

-(void)listParseStrToValues:(FEMModel *)model string:(NSString *)str index:(int)ind name:(NSString *)name values:(double *)t count:(int *)count {

    int i, l, k1, n;
    FEMVariable *variable = nil, *cVar = nil;
    variableArraysContainer *varContainers = NULL, *cvarContainers = NULL;
    FEMUtilities *utilities;
    Element_t *element;
    BOOL found;
    
    utilities = [[FEMUtilities alloc] init];
    
    //TODO: This method only support one string passed at a time, not several strings separeted by a comma. 
    
    count = 0;
    
    if ([str isEqualToString:@"coordinate"] == NO) {
        variable = [utilities getVariableFrom:model.variables model:model name:str onlySearch:NULL maskName:NULL info:&found];
        if (variable == nil) {
            warnfunct("listParseStrToValues", "Can't find indpendent variable:");
            printf("%s\n", [str UTF8String]);
            errorfunct("listParseStrToValues", "Abort...");
        }
    } else {
        variable = [utilities getVariableFrom:model.variables model:model name:@"coordinate 1" onlySearch:NULL maskName:NULL info:&found];
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
        if ([str isEqualToString:@"coordinate"] == YES) {
            cVar = [utilities getVariableFrom:model.variables model:model name:@"coordinate 1" onlySearch:NULL maskName:NULL info:&found];
            count++;
            cvarContainers = cVar.getContainers;
            t[0] = cvarContainers->Values[k1];
            
            cVar = [utilities getVariableFrom:model.variables model:model name:@"coordinate 2" onlySearch:NULL maskName:NULL info:&found];
            count++;
            cvarContainers = cVar.getContainers;
            t[1] = cvarContainers->Values[k1];
            
            cVar = [utilities getVariableFrom:model.variables model:model name:@"coordinate 3" onlySearch:NULL maskName:NULL info:&found];
            count++;
            cvarContainers = cVar.getContainers;
            t[2] = cvarContainers->Values[k1];
        }
        else {
            if (variable.dofs == 1) {
                t[*count] = varContainers->Values[k1];
                count++;
            } else {
                for (l=0; l<variable.dofs; l++) {
                    t[*count] = varContainers->Values[variable.dofs*(k1-1)+l];
                    count++;
                }
            }
        }
    } else {
        if (varContainers->Perm != NULL) {
            t[*count] = HUGE_VAL;
        } else {
            t[*count] = varContainers->Values[0];
        }
        count++;
    }
}

-(void)listSetNameSpace:(NSString *)str {
    
    if (str != nil) {
        [_nameSpace setString:str];
    }
}

-(BOOL)listGetNameSpace:(NSMutableString *)str {
    
    BOOL l;

    if ([_nameSpace isEqual:@""] == YES) {
        str = nil;
        l = NO;
    } else {
        str = [NSMutableString stringWithString:_nameSpace];
        l = YES;
    }
    
    return l;
}

-(NSString *)listGetString:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found {
    
    NSString *s = nil;
    NSMutableString *strn;
    
    *found = NO;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES || [strn isEqualToString:list.name] == YES) {
            s = [NSString stringWithString:list.cValue];
            *found = YES;
            break;
        }
    }
    
    return s;
}

-(BOOL)listGetReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName numberOfNodes:(int)n indexes:(int *)nodeIndexes buffer:(listBuffer *)result minValue:(double *)minv maxValue:(double *)maxv {

    int i, j, k;
    double *buffer;
    double t[32];
    BOOL found;
    NSMutableString *strn;
    FEMUtilities *util;
    valueListArraysContainer *containers = NULL;
    
    found = NO;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES || [strn isEqualToString:list.name] == YES) {
            // The buffer is allocated here, it's up to the caller to release this memory
            // only once it doesn't need it anymore
            if (result->vector == NULL || result->m != n) {
                if (result->vector != NULL) free_dvector(result->vector, 0, result->m-1);
                result->vector = doublevec(0, n-1);
                result->m = n;
            }
            
            containers = list.getContainers;
            if (containers->fValues == NULL) {
                NSLog(@"listGetIntegerArray: fValues not allocated in list: %@\n", varName);
                errorfunct("listGetIntegerArray", "Program terminating now...");
            }
            
            if (list.type == LIST_TYPE_CONSTANT_SCALAR) {
                
                // TODO: implement the call of a user provided method if required
                 
                if (containers->fValues == NULL) {
                    warnfunct("listGetReal", "fValues not allocated in list:\n");
                    printf("%s\n", [varName UTF8String]);
                    errorfunct("listGetReal", "Program terminating now...");
                }
                memset( result->vector, 0.0, n*sizeof(double) );
                for (i=0; i<n; i++) {
                    result->vector[i] = containers->fValues[0][0][0];
                }
            } else if (list.type == LIST_TYPE_VARIABLE_SCALAR) {
                
                util = [[FEMUtilities alloc] init];
                buffer = doublevec(0, containers->sizeTValues-1);
                for (i=0; i<n; i++) {
                    buffer[i] = containers->fValues[0][0][i];
                }
                memset( result->vector, 0.0, n*sizeof(double) );
                for (i=0; i<n; i++) {
                    k = nodeIndexes[i];
                    [self listParseStrToValues:model string:list.dependName index:k name:list.name values:t count:&j];
                    
                    if (any(t, '=', HUGE_VAL, j) == false) {
                        
                        // TODO: implement the call of a user provided method if required
                        
                        if (containers->fValues == NULL) {
                            warnfunct("listGetReal", "fValues not allocated in list:\n");
                            printf("%s\n", [varName UTF8String]);
                            errorfunct("listGetReal", "Program terminating now...");
                        }
                        result->vector[i] = [util interpolateCurveTvalues:containers->tValues fValues:buffer value:t[0] sizeOfTValues:containers->sizeTValues cubicCoefficient:containers->cubicCoeff];
                    }
                }
                
                free_dvector(buffer, 0, containers->sizeTValues-1);
            } else if ([list type] == LIST_TYPE_CONSTANT_SCALAR_STR) {
                 // TODO: implement this case
            } else if ([list type] == LIST_TYPE_VARIABLE_SCALAR_STR) {
                 // TODO: implement this case
            }
            
            found = YES;
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
    NSMutableString *strn;
    listBuffer buffer;
    valueListArraysContainer *containers = NULL;
    
    found = NO;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES || [strn isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            if (containers->fValues == NULL) {
                NSLog(@"listGetIntegerArray: fValues not allocated in list: %@\n", varName);
                errorfunct("listGetIntegerArray", "Program terminating now...");
            }
            n1 = containers->sizeFValues1;
            n2 = containers->sizeFValues2;
            
            // The buffer is allocated here but it's up to the caller to release this memory
            // only once it doesn't need it anymore
            if (result->tensor == NULL || result->m != n1 || result->n != n2 || result->p != n) {
                if (result->tensor != NULL) free_d3tensor(result->tensor, 0, result->m-1, 0, result->n-1, 0, result->p-1);
                result->tensor = d3tensor(0, n1-1, 0, n2-1, 0, n-1);
                result->m = n1;
                result->n = n2;
                result->p = n;
            }
        
            if (list.type == LIST_TYPE_CONSTANT_TENSOR) {
                memset(**result->tensor, 0.0, (n1*n2*n)*sizeof(double) );
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
                memset(**result->tensor, 0.0, (n1*n2*n)*sizeof(double) );
                found = [self listGetReal:model inArray:array forVariable:varName numberOfNodes:n indexes:nodeIndexes buffer:&buffer minValue:NULL maxValue:NULL];
                for (i=0; i<n1; i++) {
                    for (k=0; k<n; k++) {
                        result->tensor[i][0][k] = buffer.vector[k];
                    }
                }
                
                free_dvector(buffer.vector, 0, buffer.m-1);
            }
            
            found = YES;
            break; 
        }
    }
    
    return found;
}


-(double)listGetConstReal:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found minValue:(double *)minv maxValue:(double *)maxv {

    double f = 0.0;
    NSMutableString *strn;
    valueListArraysContainer *containers = NULL;
    
    *found = NO;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];

    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES || [strn isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            if (containers->fValues == NULL) {
                NSLog(@"listGetIntegerArray: fValues not allocated in list: %@\n", varName);
                errorfunct("listGetIntegerArray", "Program terminating now...");
            }
            
            if ([list type] >= 8) {
                
                // TODO: unimplemented case where type is greater than 8
                // In Elmer this calls matc. Figure out later what we do in Saino
                
            } else if (list.method == YES) {
                
                // TODO: implement the call of a user provided method
                
            } else {
                
                f = containers->fValues[0][0][0];
            }
            *found = YES;
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
    NSMutableString *strn;
    valueListArraysContainer *containers = NULL;
    
    found = NO;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES || [strn isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            if (containers->fValues == NULL) {
                NSLog(@"listGetIntegerArray: fValues not allocated in list: %@\n", varName);
                errorfunct("listGetIntegerArray", "Program terminating now...");
            }
            n1 = containers->sizeFValues1;
            n2 = containers->sizeFValues2;
            
            // The buffer is allocated here, it's up to the caller to release this memory
            // only once it doesn't need it anymore
            if (result->matrix == NULL || result->m != n1 || result->n != n2) {
                if (result->matrix != NULL) free_dmatrix(result->matrix, 0, result->m-1, 0, result->n-1);
                result->matrix = doublematrix(0, n1-1, 0, n2-1);
                result->m = n1;
                result->n = n2;
            }
            memset( *result->matrix, 0.0, (n1*n2)*sizeof(double) );            
            for (i=0; i<n1; i++) {
                for (j=0; j<n2; j++) {
                    result->matrix[i][j] = containers->fValues[i][j][0];
                }
            }
            
            // TODO: implement the call of a user provided method
            
            found = YES;
            break;
        }
    }
    
    return found;
}

-(BOOL)listGetIntegerArray:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName buffer:(listBuffer *)result {
    
    int i, n;
    BOOL found;
    NSMutableString *strn;
    valueListArraysContainer *containers = NULL;
    
    found = NO;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES || [strn isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            if (containers->iValues == NULL) {
                NSLog(@"listGetIntegerArray: iValues not allocated in list: %@\n", varName);
                errorfunct("listGetIntegerArray", "Program terminating now...");
            }
            n = containers->sizeIValues;
            
            // The buffer is allocated here, it's up to the caller to release this memory
            // only once it doesn't need it anymore
            if (result->ivector == NULL || result->m != n) {
                if (result->ivector != NULL) free_ivector(result->ivector, 0, result->m-1);
                result->ivector = intvec(0, n-1);
                result->m = n;
            }
            memset( result->ivector, 0, n*sizeof(int) );
            for (i=0; i<n; i++) {
                result->ivector[i] = containers->iValues[i];
            }
            
           // TODO: implement the call of a user provided method if required 
            
            found = YES;
            break;
        }
    }
    
    return found;
}

-(int)listGetInteger:(FEMModel *)model inArray:(NSArray *)array forVariable:(NSString *)varName info:(BOOL *)found minValue:(int *)minv maxValue:(int *)maxv {
    
    int l = 0;
    NSMutableString *strn;
    valueListArraysContainer *containers = NULL;
    
    *found = NO;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES || [strn isEqualToString:list.name] == YES) {
            containers = list.getContainers;
            if (containers->iValues == NULL) {
                NSLog(@"listGetIntegerArray: iValues not allocated in list: %@\n", varName);
                errorfunct("listGetIntegerArray", "Program terminating now...");
            }

            // TODO: implement the call of a user provided method if required
            
            if (containers->iValues == NULL) {
                warnfunct("listGetReal", "iValues not allocated in list:\n");
                printf("%s\n", [varName UTF8String]);
                errorfunct("listGetReal", "Program terminating now...");
            }
            
            l = containers->iValues[0];
            
            *found = YES;
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
    
    BOOL l = NO;
    NSMutableString *strn;
    
    *found = NO;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES || [strn isEqualToString:list.name] == YES) {
            l = list.isLvalue;
            *found = YES;
            break;
        }
    }
    
    return l;
}

-(FEMValueList *)listFindVariable:(NSString *)varName inArray:(NSArray *)array {
    
     NSMutableString *strn;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];

    for (FEMValueList *list in array) {
        if ([varName isEqualToString:[list name]] == YES || [strn isEqualToString:list.name] == YES) {
            return list;
        }
    }
    
    return nil;
}

-(BOOL)listCheckPresentVariable:(NSString *)varName inArray:(NSArray *)array {
    
    NSMutableString *strn;
    
    if ([self listGetNameSpace:strn] == YES) {
        [strn appendString:@" "];
        [strn appendString:varName];
    } else strn = [NSMutableString stringWithString:@""];

    for (FEMValueList *list in array) {
        if ([varName isEqualToString:list.name] == YES || [strn isEqualToString:list.name] == YES) {
            return YES;
        }
    }
    
    return NO;
}

/**********************************************************************************
 Adds a string to a given class list of values
 **********************************************************************************/
-(void)addStringInClassList:(id)className theVariable:(NSString *)varName withValue:(NSString *)value {
    
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMValueList *newValueList;
    NSMutableArray *valuesArray;
    BOOL found;
    
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
        newValueList.type = LIST_TYPE_STRING;
        
        newValueList.cValue = [NSString stringWithString:value];
        
        newValueList.name = varName;
        [valuesArray addObject:newValueList];
    }
    
    boundary = nil;
    bodyForce = nil;
    simulation = nil;    
}

/**********************************************************************************
 Adds a logical entry to a given class list of values
 **********************************************************************************/
-(void)addLogicalInClassList:(id)className theVariable:(NSString *)varName withValue:(BOOL)value {
    
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMValueList *newValueList;
    NSMutableArray *valuesArray;
    BOOL found;
    
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
        newValueList.type = LIST_TYPE_LOGICAL;
        
        newValueList.lValue = value;
        
        newValueList.name = varName;
        [valuesArray addObject:newValueList];
    }

    boundary = nil;
    bodyForce = nil;
    simulation = nil;
}

/**********************************************************************************
 Adds an integer to a given class list of values
 **********************************************************************************/
-(void)addIntegerInClassList:(id)className theVariable:(NSString *)varName withValue:(int)value {
    
    FEMBoundaryCondition *boundary = nil;
    FEMBodyForce *bodyForce = nil;
    FEMSimulation *simulation = nil;
    FEMValueList *newValueList;
    NSMutableArray *valuesArray;
    BOOL found;
    valueListArraysContainer *containers = NULL;
    
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
        newValueList.type = LIST_TYPE_INTEGER;
        
        containers = newValueList.getContainers;
        containers->iValues = intvec(0, 0);
        containers->sizeIValues = 1;
        containers->iValues[0] = value;
        
        newValueList.name = varName;
        [valuesArray addObject:newValueList];
    }
}

/**********************************************************************************
    Adds an integer array to a given class list of values
**********************************************************************************/
-(void)addIntegerArrayInClassList:(id)className theVariable:(NSString *)varName withValues:(int *)values numberOfNodes:(int)n {
    
    int i;
    FEMBoundaryCondition *boundary = nil;
    FEMBodyForce *bodyForce = nil;
    FEMSimulation *simulation = nil;
    FEMValueList *newValueList;
    NSMutableArray *valuesArray;
    BOOL found;
    valueListArraysContainer *containers = NULL;
    
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
        containers->sizeIValues = n;
        for (i=0; i<n; i++) {
            containers->iValues[i] = values[i];
        }
        
        newValueList.name = varName;
        [valuesArray addObject:newValueList];
    }
}

/**********************************************************************************
 Adds an constant real value to a given class list of values
 **********************************************************************************/
-(void)addConstRealInClassList:(id)className theVariable:(NSString *)varName withValue:(double)value string:(NSString *)str {
    
    FEMBoundaryCondition *boundary = nil;
    FEMBodyForce *bodyForce = nil;
    FEMSimulation *simulation = nil;
    FEMValueList *newValueList;
    NSMutableArray *valuesArray;
    BOOL found;
    valueListArraysContainer *containers = NULL;
    
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
        newValueList.type = LIST_TYPE_CONSTANT_SCALAR;
        
        containers = newValueList.getContainers;
        containers->tValues = NULL;
        containers->fValues = d3tensor(0, 0, 0, 0, 0, 0);
        containers->sizeFValues1 = 1;
        containers->sizeFValues2 = 1;
        containers->sizeFValues3 = 1;
        containers->fValues[0][0][0] = value;
        
        if (str != nil) {
            newValueList.cValue = [NSString stringWithString:str];
            newValueList.type = LIST_TYPE_CONSTANT_SCALAR_STR;
        }
        
        newValueList.name = varName;
        [valuesArray addObject:newValueList];
    }
}

/*****************************************************************************************
    Check if given element belongs to a body for which a given equation should be solved
*****************************************************************************************/
-(BOOL)checkElementEquation:(FEMModel *)model forElement:(Element_t *)element andEquation:(NSString *)equation {
    
    int k, bodyId;
    BOOL flag, found;
    FEMEquation *equationAtId;
    
    flag = NO;
    bodyId = element->BodyID;
    if ( bodyId > 0 && bodyId <= model.numberOfBodies) {
        k = [(model.bodies)[bodyId-1][@"equation"] intValue];
        if (k > 0 && k <= model.numberOfEquations) {
            equationAtId = (model.equations)[k-1];
            flag = [self listGetLogical:model inArray:equationAtId.valuesList forVariable:equation info:&found];
        }
    }
    
    return flag;
}

/*****************************************************************************
    Check if the keyword is present in any boundary condition
*****************************************************************************/
-(BOOL)listCheckPresentAnyBoundaryCondition:(FEMModel *)model name:(NSString *)name {
    
    BOOL found;
    
    found = NO;
    for (FEMBoundaryCondition *boundaryCondition in model.boundaryConditions) {
        found = [self listCheckPresentVariable:name inArray:boundaryCondition.valuesList];
        if (found == YES) break;
    }
    
    return found;
}

/*****************************************************************************
 Check if the keyword is present in any body force
*****************************************************************************/
-(BOOL)listCheckPresentAnyBodyForce:(FEMModel *)model name:(NSString *)name {
    
    BOOL found;
    
    found = NO;
    for (FEMBodyForce *bodyForce in model.bodyForces) {
        found = [self listCheckPresentVariable:name inArray:bodyForce.valuesList];
        if (found == YES) break;
    }
    
    return found;
}

@end
