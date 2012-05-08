//
//  FEMValueList.m
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMValueList.h"

#import "Constructors.h"
#import "memory.h"

@implementation FEMValueList

-(double)tValues:(int)i {
    
    return tValues[i];
}

-(double)fValues:(int)i :(int)j :(int)k {
    
    return fValues[i][j][k];
}

-(int)sizeTValues {
    
    return sizeTValues;
}

-(int)sizeFValues1 {
    
    return sizeFValues1;
}

-(int)sizeFValues2 {
    
    return sizeFValues2;
}

-(int)sizeFValues3 {
    
    return sizeFValues3;
}

-(int)type {
    
    return type;
}

-(BOOL)method {
    
    return method;
}

-(BOOL)lValue {
    
    return lValue;
}

-(int)iValues:(int)i {
    
    return iValues[i];
}

-(int)sizeIValues {
    
    return sizeIValues;
}

-(NSString *)name {
    
    return name;
}

-(NSString *)dependName {
    
    return dependName;
}

-(void)setTValues:(int)i :(double)n {
    
    tValues[i] = n;
}

-(void)setFValues:(int)i :(int)j :(int)k :(double)n {
    
    fValues[i][j][k] = n;
}

-(void)setSizeTValues:(int)n {
    
    sizeTValues = n;
}

-(void)setSizeFValues1:(int)n {
    
    sizeFValues1 = n;
}

-(void)setSizeFValues2:(int)n {
    
    sizeFValues2 = n;
}

-(void)setSizeFValues3:(int)n {
    
    sizeFValues3 = n;
}

-(void)setType:(int)n {
    
    type = n;
}

-(void)setMethod:(BOOL)n {
    
    method = n;
}

-(void)setLValue:(BOOL)n {
    
    lValue = n;
}

-(void)setIValues:(int)i :(int)n {
    
    iValues[i] = n;
}

-(void)setSizeIValues:(int)n {
    
    sizeIValues = n;
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

-(void)allocateIValues:(int)n {
    
    iValues = intvec(0, n-1);
}

@end
