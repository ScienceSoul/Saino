//
//  FEMListMatrix.h
//  Saino
//
//  Created by Seddik hakime on 07/01/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMatrix.h"
#import "Constructors.h"

@interface FEMListMatrix : NSObject {
    
    int _listMatrixGrowth;
}

@property(nonatomic, assign) int listMatrixGrowth;

-(ListMatrix_t * __nonnull)allocateMatrix:(int)n;
-(void)freeMatrix:(ListMatrix_t * __nonnull)list size:(int)n;
-(ListMatrix_t * __nonnull)enlargeMatrix:(ListMatrix_t * __nonnull)matrix toSize:(int)n;
-(ListMatrixEntry_t * __nonnull)getMatrixIndexInListMatrix:(ListMatrix_t * __nullable)list atIndex:(int)k1 andIndex:(int)k2;
-(void)addToMatrixElement:(ListMatrix_t * __nonnull)list atIndex:(int)k1 andIndex:(int)k2 value:(double)value setValue:(BOOL * __nullable)setValue;
-(void)convertToCRSMatrix:(FEMMatrix * __nonnull)matrix;

@end
