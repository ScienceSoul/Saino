//
//  FEMValueList.m
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMValueList.h"

static int LIST_TYPE_CONSTANT_SCALAR      =  1;
static int LIST_TYPE_CONSTANT_TENSOR      =  2;
static int LIST_TYPE_VARIABLE_SCALAR      =  3;
static int LIST_TYPE_VARIABLE_TENSOR      =  4;
static int LIST_TYPE_LOGICAL              =  5;
static int LIST_TYPE_STRING               =  6;
static int LIST_TYPE_INTEGER              =  7;
static int LIST_TYPE_CONSTANT_SCALAR_STR  =  8;
static int LIST_TYPE_CONSTANT_TENSIOR_STR =  9;
static int LIST_TYPE_VARIABLE_SCALAR_STR  =  10;
static int LIST_TYPE_VARIABLE_TENSOR_STR  =  11;


@implementation FEMValueList

-(int)sizeTValues {
    
    return sizeTValues;
}

-(int)sizeFValues {
    
    return sizeFValues;
}

-(int)type {
    
    return type;
}

-(NSString *)name {
    
    return name;
}

-(NSString *)dependName {
    
    return dependName;
}


-(void)setSizeTValues:(int)n {
    
    sizeTValues = n;
}

-(void)setSizeFValues:(int)n {
    
    sizeFValues = n;
}

-(void)setType:(int)n {
    
    type = n;
}

-(void)setName:(NSString *)string {
    
    name = [NSString stringWithString:string];
}

-(void)setDependName:(NSString *)string {
    
    dependName = [NSString stringWithString:string];
}

-(double *)returnPointerToTValues {
    
    return tValues;
}

-(double ***)returnPointerToFValues {
    
    return fValues;
}

-(int *)returnPointerToIValues {
    
    return iValues;
}


-(BOOL)listGetReal:(FEMModel *)model: (NSArray *)array: (NSString *)varName: (int)n: (int *)nodeIndexes: (double *)result {
    
    int i, j, k;
    double ***fVal;
    double *buffer;
    double t[32];
    FEMUtilities *util;
    BOOL found;
    
    found = NO;
    
    for (FEMValueList *list in array) {
        if ([varName isEqualToString:[list name]] == YES) {
            
            fVal = [list returnPointerToFValues];
            
            if ([list type] == LIST_TYPE_CONSTANT_SCALAR) {
                for (i=0; i<n; i++) {
                    result[i] = fVal[0][0][0];
                }
            }
            else if ([list type] == LIST_TYPE_VARIABLE_SCALAR) {
                
                util = [[FEMUtilities alloc] init];
                buffer = doublevec(0, [list sizeTValues]-1);
                
                for (i=0; i<n; i++) {
                    buffer[i] = fValues[0][0][i];
                }
                
                for (i=0; i<n; i++) {
                    k = nodeIndexes[i];
                    [self listParseStrToValues:model :[list dependName] :k :[list name] :t :j];
                    
                    if (any(t, '=', HUGE_VALF, j) == 0) {
                     
                        result[i] = [util interpolateCurve:[list returnPointerToTValues] :buffer :t[0] :[list sizeTValues]];
                    }
                }
                
                free_dvector(buffer, 0, [list sizeTValues]-1);
            }
            fVal = NULL;
            found = YES;
            break;
        }
    }
    
    return found;
    
}

-(void)listParseStrToValues:(FEMModel *)model: (NSString *)str: (int)ind: (NSString *)name: (double *)t: (int)count {
    
    int l, k1;
    char *theName;
    Variable_t *variable, *cVar;
    
    count = 0;
    
    theName = (char *)[str UTF8String];
    
    if ([str isEqualToString:@"Coordinate"] == NO) {
        variable = getVariable([model returnPointerToVariables], theName);
        if (variable == NULL) {
            warnfunct("listParseStrToValues", "Can't find indpendent variable:");
            printf("%s\n", theName);
            errorfunct("listParseStrToValues", "Abort...");
        }
    } else {
        theName = "Coordinate 1";
        variable = getVariable([model returnPointerToVariables], theName);
    }
    
    k1 = ind;
    if (variable->Perm != NULL) k1 = variable->Perm[k1];
    
    if (k1 > 0 && k1 <= variable->sizeValues) {
        if ([str isEqualToString:@"Coordinate"] == YES) {
            theName = "Coordinate 1";
            cVar = getVariable([model returnPointerToVariables], theName);
            count++;
            t[0] = cVar->Values[k1];
            
            theName = "Coordinate 2";
            cVar = getVariable([model returnPointerToVariables], theName);
            count++;
            t[1] = cVar->Values[k1];
            
            theName = "Coordinate 3";
            cVar = getVariable([model returnPointerToVariables], theName);
            count++;
            t[2] = cVar->Values[k1];
        }
        else {
            if (variable->Dofs == 1) {
                t[count] = variable->Values[k1];
                count++;
            } else {
                for (l=0; l<variable->Dofs; l++) {
                    t[count] = variable->Values[variable->Dofs*(k1-1)+l];
                    count++;
                }
            }
        }
    } else {
        if (variable->Perm != NULL) {
            t[count] = HUGE_VALF;
        } else {
            t[count] = variable->Values[0];
        }
    }
    
    variable = NULL;
    cVar = NULL;
    
}

@end
