//
//  FEMListMatrix.h
//  Saino
//
//  Created by Seddik hakime on 07/01/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"

@interface FEMListMatrix : NSObject {
    
    int _listMatrixGrowth;
}

@property(nonatomic, assign) int listMatrixGrowth;

-(ListMatrix_t *)allocateMatrix:(int)n;
-(void)freeMatrix:(ListMatrix_t *)list size:(int)n;
-(ListMatrix_t *)enlargeMatrix:(ListMatrix_t *)matrix toSize:(int)n;
-(ListMatrixEntry_t *)getMatrixIndexInListMatrix:(ListMatrix_t *)list atIndex:(int)k1 andIndex:(int)k2;

@end
