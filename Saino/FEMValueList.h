//
//  FEMValueList.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMUtilities.h"
#import "FEMModel.h"
#import "Utils.h"
#import "Lists.h"

@interface FEMValueList : NSObject {
    
    double *tValues;
    double ***fValues;
    int sizeTValues, sizeFValues;
    
    int type;
    
    BOOL lValue;
    int *iValues;
    
    char *cValue;
    
    int nameLen, depNameLen;
    NSString *name, *dependName;
    
}

-(int)sizeTValues;
-(int)sizeFValues;
-(int)type;
-(NSString *)name;
-(NSString *)dependName;

-(void)setSizeTValues:(int)n;
-(void)setSizeFValues:(int)n;
-(void)setType:(int)n;
-(void)setName:(NSString *)string;
-(void)setDependName:(NSString *)string;

-(double *)returnPointerToTValues;
-(double ***)returnPointerToFValues;
-(int *)returnPointerToIValues;


-(BOOL)listGetReal:(FEMModel *)model: (NSArray *)array: (NSString *)varName: (int)n: (int *)nodeIndexes: (double *)result;
-(void)listParseStrToValues:(FEMModel *)model: (NSString *)str: (int)ind: (NSString *)name: (double *)t: (int)count;

@end
