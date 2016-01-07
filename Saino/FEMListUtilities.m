//
//  FEMListUtilities.m
//  Saino
//
//  Created by Hakime Seddik on 18/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Accelerate/Accelerate.h>

#import "FEMListUtilities.h"
#import "FEMBodyForce.h"
#import "FEMBoundaryCondition.h"
#import "FEMSimulation.h"
#import "FEMEquation.h"
#import "FEMMaterial.h"
#import "FEMInitialConditions.h"
#import "FEMUtilities.h"
#import "Utils.h"
#import "TimeProfile.h"

@implementation FEMListUtilities {
    
    char *_nameSpace;
    char *_nameSpaceStr;
    BOOL _nameSpaceChanged;
    BOOL _timerPassive;
    BOOL _timerResults;
}

@synthesize timers = _timers;

- (id)init
{
    char emptyStr = '\0';
    
    self = [super init];
    if (self) {
        _nameSpace = &emptyStr;
        _nameSpaceStr = NULL;
        _nameSpaceChanged = YES;
        _timerPassive = NO;
        _timerResults = NO;
        _timers = [[NSMutableDictionary alloc] init];
    }
    
    return self;
}

-(void)listParseDependencies:(NSArray * __nonnull)dependencies index:(int)ind name:(NSString * __nonnull)name toValues:(double * __nonnull)t count:(int * __nonnull)count model:(FEMModel * __nonnull)model allGlobal:(BOOL * __nonnull)allGlobal {

    int k1, n;
    FEMVariable *variable = nil, *cVar = nil;
    variableArraysContainer *varContainers = NULL, *cvarContainers = NULL;
    Element_t *element;
    BOOL found;
    
    if (dependencies == nil) return;
    
    FEMUtilities *utilities = [[FEMUtilities alloc] init];
    
    *allGlobal = YES;
    *count = 0;
    
    for (NSString *variableName in dependencies) {
        if ([variableName isEqualToString:@"coordinate"] == NO) {
            variable = [utilities getVariableFrom:model.variables model:model name:variableName onlySearch:NULL maskName:NULL info:&found];
            if (variable == nil) {
                NSLog(@"FEMListUtilities:listParseStrToValues: can't find indpendent variable: %@ for dependent variable %@.\n", variableName, name);
                fatal("FEMListUtilities:listParseStrToValues", "Abort program now...");
            }
            varContainers = variable.getContainers;
            if (varContainers->sizeValues > 1) *allGlobal = NO;
        } else {
            *allGlobal = NO;
            variable = [utilities getVariableFrom:model.variables model:model name:@"coordinate 1" onlySearch:NULL maskName:NULL info:&found];
            varContainers = variable.getContainers;
        }
        
        k1 = ind;
        if (variable.type == VARIABLE_ON_NODES_ON_ELEMENTS) {
            element = model.getCurrentElement;
            if (element != NULL) {
                if (element->DGIndexes != NULL) {
                    n = element->Type.NumberOfNodes;
                    if (element->sizeDGIndexes == n) {
                        for (int i=0; i<n; i++) {
                            if (element->NodeIndexes[i] == ind) {
                                k1 = element->DGIndexes[i];
                                break;
                            }
                        }
                    }
                }
            }
        }
        if (varContainers->Perm != NULL) k1 = varContainers->Perm[k1];
        
        if (k1 >= 0 && k1 < varContainers->sizeValues) {
            if ([variableName isEqualToString:@"coordinate"] == YES) {
                cVar = [utilities getVariableFrom:model.variables model:model name:@"coordinate 1" onlySearch:NULL maskName:NULL info:&found];
                *count = *count + 1;
                cvarContainers = cVar.getContainers;
                t[0] = cvarContainers->Values[k1];
                
                cVar = [utilities getVariableFrom:model.variables model:model name:@"coordinate 2" onlySearch:NULL maskName:NULL info:&found];
                *count = *count + 1;
                cvarContainers = cVar.getContainers;
                t[1] = cvarContainers->Values[k1];
                
                cVar = [utilities getVariableFrom:model.variables model:model name:@"coordinate 3" onlySearch:NULL maskName:NULL info:&found];
                *count = *count + 1;
                cvarContainers = cVar.getContainers;
                t[2] = cvarContainers->Values[k1];
            } else {
                if (variable.dofs == 1) {
                    t[*count] = varContainers->Values[k1];
                    *count = *count + 1;
                } else {
                    for (int i=0; i<variable.dofs; i++) {
                        t[*count] = varContainers->Values[variable.dofs*k1+i];
                        *count = *count + 1;
                    }
                }
            }
        } else {
            if (varContainers->Perm != NULL) {
                t[*count] = HUGE_VAL;
            } else {
                t[*count] = varContainers->Values[0];
            }
            *count = *count + 1;
        }
    }
}

-(void)listSetNameSpace:(NSString * __nonnull)str {
    
    if (str != nil) {
        _nameSpace = (char *)[str UTF8String];
        _nameSpaceChanged = YES;
    }
}

-(char * __nonnull)listGetNameSpaceForVariable:(char * __nonnull)varName {
    
    char *str = NULL;
    char buffer[100];
    char emptyStr = '\0';

    if ( strcmp(_nameSpace, &emptyStr) == 0) {
        str = &emptyStr;
    } else {
        strncpy(buffer, _nameSpace, strlen(_nameSpace));
        strncat(buffer, " ", 1);
        strncat(buffer, varName, strlen(varName));
        str = buffer;
    }

    return str;
}

-(NSString * __nullable)listGetString:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName info:(BOOL * __nonnull)found {
    
    NSString *s = nil;
    char *nameStr = NULL;
    
    *found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            s = [NSString stringWithString:list.cValue];
            *found = YES;
            break;
        }
    }
    
    return s;
}

/*************************************************************************************
    Gets a real valued parameter in each node of an element.
************************************************************************************/
-(BOOL)listGetReal:(FEMModel * __nonnull)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName numberOfNodes:(int)n indexes:(int * __nonnull)nodeIndexes buffer:(listBuffer * __nonnull)result minValue:(double * __nullable)minv maxValue:(double * __nullable)maxv {

    int i, j, k;
    double *buffer;
    double t[32];
    BOOL allGlobal, found;
    char *nameStr = NULL;
    valueListArraysContainer *containers = NULL;
    
    found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            // The buffer is allocated here, it's up to the caller to release this memory
            // when it doesn't need it anymore
            if (result->vector == NULL || result->m != n) {
                if (result->vector != NULL) free_dvector(result->vector, 0, result->m-1);
                result->vector = doublevec(0, n-1);
                result->m = n;
            }
            
            containers = list.getContainers;
            if (list.type == LIST_TYPE_CONSTANT_SCALAR || list.type == LIST_TYPE_VARIABLE_SCALAR) {
                if (containers->fValues == NULL) {
                    NSLog(@"FEMListUtilities:listGetReal: fValues not allocated in list: %@.\n", varName);
                    fatal("FEMListUtilities:listGetReal", "Saino will abort the simulation now...");
                }
            }
            
            if (list.type == LIST_TYPE_CONSTANT_SCALAR) {
                memset( result->vector, 0.0, n*sizeof(double) );
                for (i=0; i<n; i++) {
                    result->vector[i] = containers->fValues[0][0][0];
                }
            } else if (list.type == LIST_TYPE_VARIABLE_SCALAR) {
                if (containers->tValues == NULL) {
                    NSLog(@"FEMListUtilities:listGetReal: tValues not allocated in list: %@\n", varName);
                    fatal("FEMListUtilities:listGetReal", "Saino will abort the simulation now...");
                }
                FEMUtilities *utilities = [[FEMUtilities alloc] init];
                buffer = doublevec(0, containers->sizeFValues3-1);
                for (i=0; i<n; i++) {
                    buffer[i] = containers->fValues[0][0][i];
                }
                memset( result->vector, 0.0, n*sizeof(double) );
                for (i=0; i<n; i++) {
                    k = nodeIndexes[i];
                    [self listParseDependencies:list.dependencies index:k name:list.name toValues:t count:&j model:model allGlobal:&allGlobal];
                    
                    if (any(t, '=', HUGE_VAL, j) == false) {
                        
                        // TODO: implement the call of a user function if given
                        result->vector[i] = [utilities interpolateCurveTvalues:containers->tValues fValues:buffer value:t[0] sizeOfTValues:containers->sizeTValues cubicCoefficient:containers->cubicCoeff];
                    }
                    if (allGlobal == YES) {
                        for (j=1; j<n; j++) {
                            result->vector[j] = result->vector[0];
                        }
                    }
                }
                free_dvector(buffer, 0, containers->sizeFValues3-1);
            } else if (list.type == List_TYPE_BLOCK) {
                for (i=0; i<n; i++) {
                    k = nodeIndexes[i];
                    [self listParseDependencies:list.dependencies index:k name:list.name toValues:t count:&j model:model allGlobal:&allGlobal];
                    if (any(t, '=', HUGE_VAL, j) == false) {
                        result->vector[i] = list.block(t);
                    }
                    if (allGlobal == YES) {
                         for (j=1; j<n; j++) {
                             result->vector[j] = result->vector[0];
                         }
                    }
              
                }
            }
            found = YES;
            if (minv != NULL) {
                double minVal;
                vDSP_minvD(result->vector, 1, &minVal, n);
                if (minVal < *minv) {
                    NSLog(@"FEMListUtilities:listGetReal: value smaller than given value: \n");
                    NSLog(@"Value: %f / Given value: %f / Property: %@.\n", minVal, *minv, varName);
                    fatal("FEMListUtilities:listGetReal", "Saino will abort the simulation now...");
                }
            }
            
            if (maxv != NULL) {
                double maxValue;
                vDSP_maxvD(result->vector, 1, &maxValue, n);
                if (maxValue > *maxv) {
                    NSLog(@"FEMListUtilities:listGetReal: value greater than given value: \n");
                    NSLog(@"Value: %f / Given value: %f / Property: %@.\n", maxValue, *maxv, varName);
                    fatal("FEMListUtilities:listGetReal", "Saino will abort the simulation now...");
                }
            }
            break;
        }
    }
    
    return found;
}

/*************************************************************************
    Gets a real valued parameter in one single point with value x.
    Note that this uses same logical in MDF as ListGetReal
    but the variable is just a dummy as the dependent function is
    assumed to be set inside the code. This should be used with caution
    is it sets some confusing limitations to the user.
**************************************************************************/
-(double)listGetValueParameter:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName value:(double)value info:(BOOL * __nonnull)found minValue:(double * __nullable)minv maxValue:(double * __nullable)maxv {
    
    double f, t[1];
    char *nameStr = NULL;
    valueListArraysContainer *containers = NULL;
    
    *found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    f = 0.0;
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            containers = list.getContainers;
            if (list.type == LIST_TYPE_VARIABLE_SCALAR) {
                
                // TODO: implement the call of a user function if given
                
                if (containers->tValues == NULL) {
                    NSLog(@"FEMListUtilities:listGetValueParameter: tValues not allocated in list: %@.\n", varName);
                    fatal("FEMListUtilities:listGetValueParameter", "Saino will abort the simulation now...");
                }
                FEMUtilities *utilities = [[FEMUtilities alloc] init];
                double *buffer = doublevec(0, containers->sizeFValues3-1);
                for (int i=0; i<containers->sizeFValues3; i++) {
                    buffer[i] = containers->fValues[0][0][i];
                }
                f = [utilities interpolateCurveTvalues:containers->tValues fValues:buffer value:value sizeOfTValues:containers->sizeTValues cubicCoefficient:containers->cubicCoeff];
                free_dvector(buffer, 0, containers->sizeFValues3-1);
            } else if (list.type == List_TYPE_BLOCK) {
                t[0] = value;
                f = list.block(t);
            } else {
                NSLog(@"FEMListUtilities:listGetValueParameter: list type not supported.\n");
                fatal("FEMListUtilities:listGetValueParameter", "Saino will abort the simulation now...");
            }
            *found = YES;
            if (minv != NULL) {
                if (f < *minv) {
                    NSLog(@"FEMListUtilities:listGetValueParameter: value smaller than given value: \n");
                    NSLog(@"Value: %f / Given value: %f / Property: %@.\n", f, *minv, varName);
                    fatal("FEMListUtilities:listGetValueParameter", "Saino will abort the simulation now...");
                }
            }
            
            if (maxv != NULL) {
                if (f > *maxv) {
                    NSLog(@"FEMListUtilities:listGetValueParameter: value greater than given value: \n");
                    NSLog(@"Value: %f / Given value: %f / Property: %@.\n", f, *maxv, varName);
                    fatal("FEMListUtilities:listGetValueParameter", "Saino will abort the simulation now...");
                }
            }
            break;
        }
    }
    return f;
}

/**********************************************************************************
    Gets a real array from the list by its name
**********************************************************************************/
-(BOOL)listGetRealArray:(FEMModel * __nonnull)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName numberOfNodes:(int)n indexes:(int * __nonnull)nodeIndexes buffer:(listBuffer * __nonnull)result {

    int i, j, k, n1, n2;
    double t[32];
    BOOL allGlobal=NO, found;
    char *nameStr = NULL;
    listBuffer buffer = { NULL, NULL, NULL, NULL, 0, 0, 0};
    valueListArraysContainer *containers = NULL;
    
    found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            containers = list.getContainers;
            
            if (list.type == LIST_TYPE_CONSTANT_TENSOR || list.type == LIST_TYPE_VARIABLE_TENSOR) {
                if (containers->fValues == NULL) {
                    NSLog(@"FEMListUtilities:listGetRealArray: fValues not allocated in list: %@.\n", varName);
                    fatal("FEMListUtilities:listGetRealArray", "Saino will abort the simulation now...");
                }
            }
            
            n1 = containers->sizeFValues1;
            n2 = containers->sizeFValues2;
            
            // The buffer is allocated here but it's up to the caller to release this memory
            // when it doesn't need it anymore
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
            } else if (list.type == List_TYPE_BLOCK || list.type == LIST_TYPE_VARIABLE_TENSOR){
                memset(**result->tensor, 0.0, (n1*n2*n)*sizeof(double) );
                for (i=0; i<n; i++) {
                    k = nodeIndexes[i];
                    [self listParseDependencies:list.dependencies index:k name:list.name toValues:t count:&j model:model allGlobal:&allGlobal];
                    if (any(t, '=', HUGE_VAL, j) == true) continue;

                    if (list.type == List_TYPE_BLOCK) {
                        for (j=0; j<n1; j++) {
                            for (k=0; k<n2; k++) {
                                result->tensor[j][k][i] = list.block(t);
                            }
                        }
                    } else {
                        FEMUtilities *utilities = [[FEMUtilities alloc] init];
                        double *buffer = doublevec(0, containers->sizeFValues3-1);
                        for (j=0; j<n1; j++) {
                            for (k=0; k<n2; k++) {
                                for (int l=0; l<n; l++) {
                                    buffer[l] = containers->fValues[j][k][l];
                                }
                                result->tensor[j][k][i] =  [utilities interpolateCurveTvalues:containers->tValues fValues:buffer value:t[0] sizeOfTValues:containers->sizeTValues cubicCoefficient:containers->cubicCoeff];
                            }
                        }
                        free_dvector(buffer, 0, containers->sizeFValues3-1);
                    }
                    if (allGlobal == YES) break;
                }
                if (allGlobal == YES) {
                    for (i=1; i<n; i++) {
                        for (j=0; j<n1; j++) {
                            for (k=0; k<n2; k++) {
                                result->tensor[j][k][i] = result->tensor[j][k][0];
                            }
                        }
                    }
                }
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

/**********************************************************************************
    Gets a constant real from the list by its name
**********************************************************************************/
-(double)listGetConstReal:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName info:(BOOL * __nonnull)found minValue:(double * __nullable)minv maxValue:(double * __nullable)maxv {

    double f = 0.0;
    char *nameStr = NULL;
    valueListArraysContainer *containers = NULL;
    
    *found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            containers = list.getContainers;
            if (containers->fValues == NULL) {
                NSLog(@"FEMListUtilities:listGetConstReal: fValues not allocated in list: %@.\n", varName);
                fatal("FEMListUtilities:listGetConstReal", "Saino will abort the simulation now...");
            }
            
            // In contrary to Elmer, we don't call any function or use MATC because even if we
            // use a block for evaluation, this is done when this variable is created since the
            // values are node-wise constant.
            f = containers->fValues[0][0][0];
            *found = YES;
            
            if (minv != NULL) {
                if (f < *minv) {
                    NSLog(@"FEMListUtilities:listGetConstReal: value smaller than given value: \n");
                    NSLog(@"Value: %f / Given value: %f / Property: %@.\n", f, *minv, varName);
                    fatal("FEMListUtilities:listGetConstReal", "Saino will abort the simulation now...");
                }
            }
            
            if (maxv != NULL) {
                if (f > *maxv) {
                    NSLog(@"FEMListUtilities:listGetConstReal: value greater than given value: \n");
                    NSLog(@"Value: %f / Given value: %f / Property: %@.\n", f, *maxv, varName);
                    fatal("FEMListUtilities:listGetConstReal", "Saino will abort the simulation now...");
                }
            }
            break;
        }
    }
    
    return f;
}

/**********************************************************************************
    Gets a constant real array from the list by its name
**********************************************************************************/
-(BOOL)listGetConstRealArray:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName buffer:(listBuffer * __nonnull)result {

    int i, j, n1, n2;
    BOOL found;
    char *nameStr = NULL;
    valueListArraysContainer *containers = NULL;
    
    found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            containers = list.getContainers;
            if (containers->fValues == NULL) {
                NSLog(@"FEMListUtilities:listGetConstRealArray: fValues not allocated in list: %@.\n", varName);
                fatal("FEMListUtilities:listGetConstRealArray", "Saino will abort the simulation now...");
            }
            n1 = containers->sizeFValues1;
            n2 = containers->sizeFValues2;
            
            // The buffer is allocated here, it's up to the caller to release this memory
            // when it doesn't need it anymore
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
            
            found = YES;
            break;
        }
    }
    
    return found;
}

/**********************************************************************************
    Gets an integer array from the list by its name
**********************************************************************************/
-(BOOL)listGetIntegerArray:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName buffer:(listBuffer * __nonnull)result {
    
    int i, n;
    BOOL found;
    char *nameStr = NULL;
    valueListArraysContainer *containers = NULL;
    
    found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            containers = list.getContainers;
            if (containers->iValues == NULL) {
                NSLog(@"FEMListUtilities:listGetIntegerArray: iValues not allocated in list: %@.\n", varName);
                fatal("FEMListUtilities:listGetIntegerArray", "Saino will abort the simulation now...");
            }
            n = containers->sizeIValues;
            
            // The buffer is allocated here, it's up to the caller to release this memory
            // when it doesn't need it anymore
            if (result->ivector == NULL || result->m != n) {
                if (result->ivector != NULL) free_ivector(result->ivector, 0, result->m-1);
                result->ivector = intvec(0, n-1);
                result->m = n;
            }
            memset( result->ivector, 0, n*sizeof(int) );
            for (i=0; i<n; i++) {
                result->ivector[i] = containers->iValues[i];
            }
            
            found = YES;
            break;
        }
    }
    
    return found;
}

-(int)listGetInteger:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName info:(BOOL * __nonnull)found minValue:(int * __nullable)minv maxValue:(int * __nullable)maxv {
    
    int l = 0;
    char *nameStr = NULL;
    valueListArraysContainer *containers = NULL;
    
    *found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            containers = list.getContainers;
            if (containers->iValues == NULL) {
                NSLog(@"FEMListUtilities:listGetInteger: iValues not allocated in list: %@.\n", varName);
                fatal("FEMListUtilities:listGetInteger", "Saino will abort the simulation now...");
            }
            l = containers->iValues[0];
            *found = YES;
            
            if (minv != NULL) {
                if (l < *minv) {
                    NSLog(@"FEMListUtilities:listGetInteger: value smaller than given value: \n");
                    NSLog(@"Value: %d / Given value: %d / Property: %@.\n", l, *minv, varName);
                    fatal("FEMListUtilities:listGetInteger", "Saino will abort the simulation now...");
                }
            }
            
            if (maxv != NULL) {
                if (l > *maxv) {
                    NSLog(@"FEMListUtilities:listGetInteger: value greater than given value: \n");
                    NSLog(@"Value: %d / Given value: %d / Property: %@.\n", l, *maxv, varName);
                    fatal("FEMListUtilities:listGetInteger", "Saino will abort the simulation now...");
                }
            }
            break;
        }
    }
    
    return l;
}

-(BOOL)listGetLogical:(FEMModel * __nullable)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName info:(BOOL * __nonnull)found {
    
    BOOL l = NO;
    char *nameStr = NULL;
    
    *found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            l = list.isLvalue;
            *found = YES;
            break;
        }
    }
    
    return l;
}

-(FEMValueList * __nullable)listFindVariable:(NSString * __nonnull)varName inArray:(NSArray * __nonnull)array {
    
    char *nameStr = NULL;
    
    if (varName == nil) return nil;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];

    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            return list;
        }
    }
    
    return nil;
}

/*******************************************************************************
    Finds an entry in the list by its name and returns a handle to it.
    This one just finds a keyword with the same start as specified by 'prefix'.
    
    Method corresponds to Elmer from git on October 27 2015
*******************************************************************************/
-(FEMValueList * __nullable)listFindPrefix:(NSString * __nonnull)prefix inArray:(NSArray * __nonnull)array info:(BOOL * __nonnull)found {
    
    int k, n;
    char *nameStr = NULL;
    FEMValueList *valueList;
    
    *found = NO;
    
    if (prefix == nil) return nil;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[prefix UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[prefix UTF8String];
    
    if (_nameSpaceStr != NULL) {
        k = (int)strlen(_nameSpaceStr);
        for (FEMValueList *list in array) {
            nameStr = (char *)[list.name UTF8String];
            n = (int)[list.name length];
            if (n >= k) {
                if (strncmp(_nameSpaceStr, nameStr, k) == 0) {
                    valueList = list;
                    break;
                }
            }
        }
    } else {
        k = (int)strlen(varStr);
        for (FEMValueList *list in array) {
            nameStr = (char *)[list.name UTF8String];
            n = (int)[list.name length];
             if (n >= k) {
                 if (strncmp(varStr, nameStr, k) == 0) {
                     valueList = list;
                     break;
                 }
             }
        }
    }
    
    if (valueList != nil) {
        *found = YES;
    } else {
        NSLog(@"FEMListUtilities:listFindPrefix: requested prefix [%@] not found.\n", prefix);
    }

    return valueList;
}

/************************************************************
    Just checks if a prefix is present in the list.
 
    Method corresponds to Elmer from git on October 27 2015
************************************************************/

-(BOOL)listCheckPrefix:(NSString * __nonnull)prefix inArray:(NSArray * __nonnull)array {
    
    FEMValueList *valueList;
    BOOL found;
    
    valueList = [self listFindPrefix:prefix inArray:array info:&found];
    return found;
}

-(BOOL)listGetDerivativeValue:(FEMModel * __nonnull)model inArray:(NSArray * __nonnull)array forVariable:(NSString * __nonnull)varName numberOfNodes:(int)n indexes:(int * __nonnull)nodeIndexes buffer:(listBuffer * __nonnull)result {
    
    int i, k;
    double *buffer, t;
    BOOL found, stat;
    char *nameStr = NULL;
    FEMVariable *variable;
    FEMUtilities *utilities;
    valueListArraysContainer *containers = NULL;
    variableArraysContainer *varContainers = NULL;
    
    found = NO;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];

    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            // The buffer is allocated here, it's up to the caller to release this memory
            // only once it doesn't need it anymore
            if (result->vector == NULL || result->m != n) {
                if (result->vector != NULL) free_dvector(result->vector, 0, result->m-1);
                result->vector = doublevec(0, n-1);
                result->m = n;
            }
            
            containers = list.getContainers;
            if (containers->fValues == NULL) {
                NSLog(@"FEMListUtilities:listGetDerivativeValue: fValues not allocated in list: %@.\n", varName);
                fatal("FEMListUtilities:listGetDerivativeValue", "Saino will abort the simulation now...");
            }
            
            if (list.type == LIST_TYPE_VARIABLE_SCALAR) {
                 if (containers->tValues == NULL) {
                     NSLog(@"FEMListUtilities:listGetDerivativeValue: tValues not allocated in list: %@.\n", varName);
                     fatal("FEMListUtilities:listGetDerivativeValue", "Saino will abort the simulation now...");
                 }
                utilities = [[FEMUtilities alloc] init];
                // Apparently this routine assumes that there is only one associated dependency
                NSString *dependent = list.dependencies[0];
                variable = [utilities getVariableFrom:model.variables model:model name:dependent onlySearch:NULL maskName:nil info:&stat];
                varContainers = variable.getContainers;
                memset( result->vector, 0.0, n*sizeof(double) );
                buffer = doublevec(0, containers->sizeTValues-1);
                for (i=0; i<n; i++) {
                    buffer[i] = containers->fValues[0][0][i];
                }
                for (i=0; i<n; i++) {
                    k = nodeIndexes[i];
                    if (varContainers->Perm != NULL) k = varContainers->Perm[k];
                    if (k >= 0) {
                        t = varContainers->Values[k];
                        result->vector[i] = [utilities derivateCurveTvalues:containers->tValues fValues:buffer value:t sizeOfTValues:containers->sizeTValues cubicCoefficient:containers->cubicCoeff];
                    }
                }
                free_dvector(buffer, 0, containers->sizeTValues-1);
            } else {
                NSLog(@"FEMListUtilities:listGetDerivativeValue: no automated derivation possible for %@.\n", varName);
            }
            
            found = YES;
            break;
        }
    }
    return found;
}

-(BOOL)listCheckPresentVariable:(NSString * __nonnull)varName inArray:(NSArray * __nonnull)array {
    
    char *nameStr = NULL;
    
    if (_nameSpaceChanged == YES) {
        _nameSpaceStr = [self listGetNameSpaceForVariable:(char *)[varName UTF8String]];
        _nameSpaceChanged = NO;
    }
    char *varStr = (char *)[varName UTF8String];
    
    for (FEMValueList *list in array) {
        nameStr = (char *)[list.name UTF8String];
        if (strcmp(varStr, nameStr) == 0 || strcmp(_nameSpaceStr, nameStr) == 0) {
            return YES;
        }
    }
    
    return NO;
}

/**********************************************************************************
    Adds a string to a given class list of values
**********************************************************************************/
-(void)addStringInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValue:(NSString * __nonnull)value {
    
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMConstants *constants;
    FEMEquation *equation;
    FEMMaterial *material;
    FEMInitialConditions *initialCondition;
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
    } else if ([className isKindOfClass:[FEMConstants class]]) {
        constants = className;
        valuesArray = constants.valuesList;
    } else if ([className isKindOfClass:[FEMEquation class]]) {
        equation = className;
        valuesArray = equation.valuesList;
    } else if ([className isKindOfClass:[FEMMaterial class]]) {
        material = className;
        valuesArray = material.valuesList;
    } else if ([className isKindOfClass:[FEMInitialConditions class]]) {
        initialCondition = className;
        valuesArray = initialCondition.valuesList;
    }
    
    for (FEMValueList *list in valuesArray) {
        if ([varName isEqualToString:list.name] == YES) {
            newValueList = list;
            found = YES;
            break;
        }
    }
    
    if (found == YES) { // if object exists, remove it and we will create one new
        [newValueList deallocation];
        [valuesArray removeObject:newValueList];
    }
        
    newValueList = [[FEMValueList alloc] init];
    newValueList.type = LIST_TYPE_STRING;
    
    newValueList.cValue = [NSString stringWithString:value];
    
    newValueList.name = varName;
    [valuesArray addObject:newValueList];

    
    boundary = nil;
    bodyForce = nil;
    simulation = nil;    
}

/**********************************************************************************
    Adds a logical entry to a given class list of values
**********************************************************************************/
-(void)addLogicalInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValue:(BOOL)value {
    
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMConstants *constants;
    FEMEquation *equation;
    FEMMaterial *material;
    FEMInitialConditions *initialCondition;
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
    } else if ([className isKindOfClass:[FEMConstants class]]) {
        constants = className;
        valuesArray = constants.valuesList;
    } else if ([className isKindOfClass:[FEMEquation class]]) {
        equation = className;
        valuesArray = equation.valuesList;
    } else if ([className isKindOfClass:[FEMMaterial class]]) {
        material = className;
        valuesArray = material.valuesList;
    } else if ([className isKindOfClass:[FEMInitialConditions class]]) {
        initialCondition = className;
        valuesArray = initialCondition.valuesList;
    }

    for (FEMValueList *list in valuesArray) {
        if ([varName isEqualToString:list.name] == YES) {
            newValueList = list;
            found = YES;
            break;
        }
    }
    
    if (found == YES) { // if object exists, remove it and we will create one new
        [newValueList deallocation];
        [valuesArray removeObject:newValueList];
    }
        
    newValueList = [[FEMValueList alloc] init];
    newValueList.type = LIST_TYPE_LOGICAL;
    
    newValueList.lValue = value;
    
    newValueList.name = varName;
    [valuesArray addObject:newValueList];
    
    boundary = nil;
    bodyForce = nil;
    simulation = nil;
}

/**********************************************************************************
    Adds an integer to a given class list of values
**********************************************************************************/
-(void)addIntegerInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValue:(int * __nullable)value orUsingBlock:(double (^ __nullable)())block {
    
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMConstants *constants;
    FEMEquation *equation;
    FEMMaterial *material;
    FEMInitialConditions *initialCondition;
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
    } else if ([className isKindOfClass:[FEMConstants class]]) {
        constants = className;
        valuesArray = constants.valuesList;
    } else if ([className isKindOfClass:[FEMEquation class]]) {
        equation = className;
        valuesArray = equation.valuesList;
    } else if ([className isKindOfClass:[FEMMaterial class]]) {
        material = className;
        valuesArray = material.valuesList;
    } else if ([className isKindOfClass:[FEMInitialConditions class]]) {
        initialCondition = className;
        valuesArray = initialCondition.valuesList;
    }
    
    for (FEMValueList *list in valuesArray) {
        if ([varName isEqualToString:list.name] == YES) {
            newValueList = list;
            found = YES;
            break;
        }
    }
    
    if (found == YES) { // if object exists, remove it and we will create one new
        [newValueList deallocation];
        [valuesArray removeObject:newValueList];
    }
        
    newValueList = [[FEMValueList alloc] init];
    newValueList.type = LIST_TYPE_INTEGER;
    
    containers = newValueList.getContainers;
    containers->iValues = intvec(0, 0);
    containers->sizeIValues = 1;
    
    if (value == NULL && block == nil) {
        NSLog(@"FEMListUtilities:addIntegerInClassList: no valid (value or block) input for variable %@.\n", varName);
        fatal("FEMListUtilities:addIntegerInClassList", "Saino will abort the simulation now...");
        
    }
    if (value != NULL && block != nil) {
        NSLog(@"FEMListUtilities:addIntegerInClassList: value and block are both non-null for variable %@.\n", varName);
        fatal("FEMListUtilities:addIntegerInClassList", "Saino will abort the simulation now...");
    }
    if (value != NULL) {
        containers->iValues[0] = *value;
    } else if (block != nil) {
        containers->iValues[0] = block();
    }
    
    newValueList.name = varName;
    [valuesArray addObject:newValueList];

}

/**********************************************************************************
    Adds an integer array to a given class list of values
**********************************************************************************/
-(void)addIntegerArrayInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValues:(int * __nullable)values size:(int)n orUsingBlock:(double (^ __nullable)())block {
    
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMConstants *constants;
    FEMEquation *equation;
    FEMMaterial *material;
    FEMInitialConditions *initialCondition;
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
    } else if ([className isKindOfClass:[FEMConstants class]]) {
        constants = className;
        valuesArray = constants.valuesList;
    } else if ([className isKindOfClass:[FEMEquation class]]) {
        equation = className;
        valuesArray = equation.valuesList;
    } else if ([className isKindOfClass:[FEMMaterial class]]) {
        material = className;
        valuesArray = material.valuesList;
    } else if ([className isKindOfClass:[FEMInitialConditions class]]) {
        initialCondition = className;
        valuesArray = initialCondition.valuesList;
    }

    for (FEMValueList *list in valuesArray) {
        if ([varName isEqualToString:list.name] == YES) {
            newValueList = list;
            found = YES;
            break;
        }
    }
    
    if (found == YES) { // if object exists, remove it and we will create one new
        [newValueList deallocation];
        [valuesArray removeObject:newValueList];
    }
        
    newValueList = [[FEMValueList alloc] init];
    newValueList.type = LIST_TYPE_CONSTANT_TENSOR;
    
    containers = newValueList.getContainers;
    containers->iValues = intvec(0, n-1);
    containers->sizeIValues = n;
    if (values == NULL && block == nil) {
        NSLog(@"FEMListUtilities:addIntegerArrayInClassList: no valid (value or block) input for variable %@.\n", varName);
        fatal("FEMListUtilities:addIntegerArrayInClassList", "Saino will abort the simulation now...");
        
    }
    if (values != NULL && block != nil) {
        NSLog(@"FEMListUtilities:addIntegerArrayInClassList: value and block are both non-null for variable %@.\n", varName);
        fatal("FEMListUtilities:addIntegerArrayInClassList", "Saino will abort the simulation now...");
    }
    if (values != NULL) {
        memcpy(containers->iValues, values, n*sizeof(int));
    } else if (block != nil) {
        for (int i=0; i<n; i++) {
            containers->iValues[i] = block();
        }
    }
    
    newValueList.name = varName;
    [valuesArray addObject:newValueList];
}

/**********************************************************************************
    Adds an constant real value to a given class list of values
**********************************************************************************/
-(void)addConstRealInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValue:(double * __nullable)value orUsingBlock:(double (^ __nullable)())block string:(NSString * __nullable)str {
    
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMConstants *constants;
    FEMEquation *equation;
    FEMMaterial *material;
    FEMInitialConditions *initialCondition;
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
    } else if ([className isKindOfClass:[FEMConstants class]]) {
        constants = className;
        valuesArray = constants.valuesList;
    } else if ([className isKindOfClass:[FEMEquation class]]) {
        equation = className;
        valuesArray = equation.valuesList;
    } else if ([className isKindOfClass:[FEMMaterial class]]) {
        material = className;
        valuesArray = material.valuesList;
    } else if ([className isKindOfClass:[FEMInitialConditions class]]) {
        initialCondition = className;
        valuesArray = initialCondition.valuesList;
    }
    
    for (FEMValueList *list in valuesArray) {
        if ([varName isEqualToString:list.name] == YES) {
            newValueList = list;
            found = YES;
            break;
        }
    }
    
    if (found == YES) { // if object exists, remove it and we will create one new
        [newValueList deallocation];
        [valuesArray removeObject:newValueList];
    }

    newValueList = [[FEMValueList alloc] init];
    newValueList.type = LIST_TYPE_CONSTANT_SCALAR;
    
    containers = newValueList.getContainers;
    containers->tValues = NULL;
    containers->fValues = d3tensor(0, 0, 0, 0, 0, 0);
    containers->sizeFValues1 = 1;
    containers->sizeFValues2 = 1;
    containers->sizeFValues3 = 1;
    
    if (value == NULL && block == nil) {
        NSLog(@"FEMListUtilities:addConstRealInClassList: no valid (value or block) input for variable %@.\n", varName);
        fatal("FEMListUtilities:addConstRealInClassList", "Saino will abort the simulation now...");

    }
    if (value != NULL && block != nil) {
        NSLog(@"FEMListUtilities:addConstRealInClassList: value and block are both non-null for variable %@.\n", varName);
        fatal("FEMListUtilities:addConstRealInClassList", "Saino will abort the simulation now...");
    }
    if (value != NULL) {
        containers->fValues[0][0][0] = *value;
    } else if (block != nil) {
        containers->fValues[0][0][0] = block();
    }
    
    if (str != nil) {
        newValueList.cValue = [NSString stringWithString:str];
        newValueList.type = LIST_TYPE_CONSTANT_SCALAR_STR;
    }
    
    newValueList.name = varName;
    [valuesArray addObject:newValueList];
}

-(void)addConstRealArrayInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName withValues:(double * __nullable * __nullable)fvalues size1:(int)m size2:(int)n orUsingBlock:(double (^ __nullable)())block string:(NSString * __nullable)str {
    
    int i, j;
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMConstants *constants;
    FEMEquation *equation;
    FEMMaterial *material;
    FEMInitialConditions *initialCondition;
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
    } else if ([className isKindOfClass:[FEMConstants class]]) {
        constants = className;
        valuesArray = constants.valuesList;
    } else if ([className isKindOfClass:[FEMEquation class]]) {
        equation = className;
        valuesArray = equation.valuesList;
    } else if ([className isKindOfClass:[FEMMaterial class]]) {
        material = className;
        valuesArray = material.valuesList;
    } else if ([className isKindOfClass:[FEMInitialConditions class]]) {
        initialCondition = className;
        valuesArray = initialCondition.valuesList;
    }

    for (FEMValueList *list in valuesArray) {
        if ([varName isEqualToString:list.name] == YES) {
            newValueList = list;
            found = YES;
            break;
        }
    }
    
    if (found == YES) { // if object exists, remove it and we will create one new
        [newValueList deallocation];
        [valuesArray removeObject:newValueList];
    }
        
    newValueList = [[FEMValueList alloc] init];
    newValueList.type = LIST_TYPE_CONSTANT_TENSOR;
    
    containers = newValueList.getContainers;
    containers->fValues = d3tensor(0, m-1, 0, n-1, 0, 0);
    containers->sizeFValues1 = m;
    containers->sizeFValues2 = n;
    containers->sizeFValues3 = 1;
    
    if (fvalues == NULL && block == nil) {
        NSLog(@"FEMListUtilities:addConstRealInClassList: no valid (value or block) input for variable %@.\n", varName);
        fatal("FEMListUtilities:addConstRealInClassList", "Saino will abort the simulation now...");
        
    }
    if (fvalues != NULL && block != nil) {
        NSLog(@"FEMListUtilities:addConstRealInClassList: value and block are both non-null for variable %@.\n", varName);
        fatal("FEMListUtilities:addConstRealInClassList", "Saino will abort the simulation now...");
    }
    if (fvalues != NULL) {
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                containers->fValues[i][j][0] = fvalues[i][j];
            }
        }
    } else if (block != nil) {
        for (i=0; i<m; i++) {
            for (j=0; j<n; j++) {
                containers->fValues[i][j][0] = block();
            }
        }
    }

    if (str != nil) {
        newValueList.cValue = [NSString stringWithString:str];
        newValueList.type = LIST_TYPE_CONSTANT_TENSIOR_STR;
    }
    
    newValueList.name = varName;
    [valuesArray addObject:newValueList];
}

-(void)addBlockInClassList:(id __nonnull)className theVariable:(NSString * __nonnull)varName usingBlock:(double (^ __nonnull)(double * __nullable variablesValues))block dependencies:(NSArray * __nullable)dependencies {
    
    FEMBoundaryCondition *boundary;
    FEMBodyForce *bodyForce;
    FEMSimulation *simulation;
    FEMConstants *constants;
    FEMEquation *equation;
    FEMMaterial *material;
    FEMInitialConditions *initialCondition;
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
    } else if ([className isKindOfClass:[FEMConstants class]]) {
        constants = className;
        valuesArray = constants.valuesList;
    } else if ([className isKindOfClass:[FEMEquation class]]) {
        equation = className;
        valuesArray = equation.valuesList;
    } else if ([className isKindOfClass:[FEMMaterial class]]) {
        material = className;
        valuesArray = material.valuesList;
    } else if ([className isKindOfClass:[FEMInitialConditions class]]) {
        initialCondition = className;
        valuesArray = initialCondition.valuesList;
    }
    
    for (FEMValueList *list in valuesArray) {
        if ([varName isEqualToString:list.name] == YES) {
            newValueList = list;
            found = YES;
            break;
        }
    }
    
    if (found == YES) { // if object exists, remove it and we will create one new
        [newValueList deallocation];
        [valuesArray removeObject:newValueList];
    }
    
    newValueList = [[FEMValueList alloc] init];
    newValueList.type = List_TYPE_BLOCK;
    
    newValueList.block = block;
    
    if (dependencies != nil) {
        newValueList.dependencies = [[NSArray alloc] initWithArray:dependencies];
    }
    
    newValueList.name = varName;
    [valuesArray addObject:newValueList];
}

/*****************************************************************************************
    Check if given element belongs to a body for which a given equation should be solved
*****************************************************************************************/
-(BOOL)checkElementEquation:(FEMModel * __nonnull)model forElement:(Element_t * __nonnull)element andEquation:(NSString * __nonnull)equation {
    
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
-(BOOL)listCheckPresentAnyBoundaryCondition:(FEMModel * __nonnull)model name:(NSString * __nonnull)name {
    
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
-(BOOL)listCheckPresentAnyBodyForce:(FEMModel * __nonnull)model name:(NSString * __nonnull)name {
    
    BOOL found;
    
    found = NO;
    for (FEMBodyForce *bodyForce in model.bodyForces) {
        found = [self listCheckPresentVariable:name inArray:bodyForce.valuesList];
        if (found == YES) break;
    }
    
    return found;
}

/************************************************************
    Check current time of the timer
    
    Method corresponds to Elmer from git on October 27 2015
************************************************************/
-(void)checkTimer:(NSString * __nonnull)timerName deleteTimer:(BOOL * __nullable)deleteTimer resetTimer:(BOOL * __nullable)resetTimer model:(FEMModel * __nonnull)model {
    
    if (_timerPassive == YES) return;
    
    double ct = 0.0;
    double rt = 0.0;
    
    NSString *name = [timerName stringByAppendingString:@" cpu time"];
    if ((self.timers)[name] != nil) {
        double ct0 = [(self.timers)[name] doubleValue];
        NSString *name = [timerName stringByAppendingString:@" real time"];
        double rt0 = [(self.timers)[name] doubleValue];
        ct = cputime() - ct0;
        rt = realtime() - rt0;
        
        NSLog(@"FEMListUtilities:checkTimer: elapsed time (CPU, REAL): %f %f.\n", ct, rt);
        if (_timerResults == YES) {
            NSString *string = [timerName stringByAppendingString:@" cpu time"];
            [self addConstRealInClassList:model.simulation theVariable:[@"res: " stringByAppendingString:string] withValue:&ct orUsingBlock:nil string:nil];
            string = [timerName stringByAppendingString:@" real time"];
            [self addConstRealInClassList:model.simulation theVariable:[@"res: " stringByAppendingString:string] withValue:&rt orUsingBlock:nil string:nil];
        }
    } else {
        NSLog(@"FEMListUtilities:checkTimer: requesting time from non-existing timer: %@.\n", timerName);
    }
    
    if (resetTimer != NULL) {
        if (*resetTimer == YES) {
            [self.timers setObject:@(ct) forKey:[timerName stringByAppendingString:@" cpu time"]];
            [self.timers setObject:@(rt) forKey:[timerName stringByAppendingString:@" real time"]];
        }
    }
    
    if (deleteTimer != NULL) {
        if (*deleteTimer == YES) [self deletTimer:timerName];
    }
}

/******************************************************************************
    A timer that uses a list structure to store the times making in
    generally applicable without any upper limit on the number of timers.
    This resets the timer.
 
    Method corresponds to Elmer from git on October 27 2015
******************************************************************************/
-(void)resetTimer:(NSString * __nonnull)timerName model:(FEMModel * __nonnull)model {
    
    double ct, rt;
    static BOOL firstTime = YES;
    BOOL found;
    
    if (firstTime == YES) {
        firstTime = NO;
        _timerPassive = [self listGetLogical:model inArray:model.simulation.valuesList forVariable:@"timer passive" info:&found];
        _timerResults = [self listGetLogical:model inArray:model.simulation.valuesList forVariable:@"timer results" info:&found];
    }
    
    if (_timerPassive == YES) return;
    
    ct = cputime();
    rt = realtime();
    
    [self.timers setObject:@(ct) forKey:[timerName stringByAppendingString:@" cpu time"]];
    [self.timers setObject:@(rt) forKey:[timerName stringByAppendingString:@" real time"]];
}

/*************************************************************
    Delete an existing timer
 
    Method corresponds to Elmer from git on October 27 2015
*************************************************************/
-(void)deletTimer:(NSString * __nonnull)timerName {
    
    if (_timerPassive == YES) return;
    
    [self.timers removeObjectForKey:[timerName stringByAppendingString:@" cpu time"]];
    [self.timers removeObjectForKey:[timerName stringByAppendingString:@" real time"]];
}



@end
