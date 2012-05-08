//
//  FEMValueList.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMValueList : NSObject {
    
    double *tValues;
    double ***fValues;
    int sizeTValues, sizeFValues1, sizeFValues2, sizeFValues3;
    
    int type;
    BOOL method;
    
    BOOL lValue;
    int *iValues;
    int sizeIValues;
    
    char *cValue;
    
    int nameLen, depNameLen;
    NSString *name, *dependName;
    
}

-(double)tValues:(int)i;
-(double)fValues:(int)i: (int)j: (int)k;
-(int)sizeTValues;
-(int)sizeFValues1;
-(int)sizeFValues2;
-(int)sizeFValues3;
-(int)type;
-(BOOL)method;
-(BOOL)lValue;
-(int)iValues:(int)i;
-(int)sizeIValues;
-(NSString *)name;
-(NSString *)dependName;

-(void)setTValues:(int)i: (double)n;
-(void)setFValues:(int)i: (int)j: (int)k: (double)n;
-(void)setSizeTValues:(int)n;
-(void)setSizeFValues1:(int)n;
-(void)setSizeFValues2:(int)n;
-(void)setSizeFValues3:(int)n;
-(void)setType:(int)n;
-(void)setMethod:(BOOL)n;
-(void)setLValue:(BOOL)n;
-(void)setIValues:(int)i: (int)n;
-(void)setSizeIValues:(int)n;
-(void)setName:(NSString *)string;
-(void)setDependName:(NSString *)string;

-(double *)returnPointerToTValues;
-(double ***)returnPointerToFValues;
-(int *)returnPointerToIValues;

-(void)allocateIValues:(int)n;

@end
