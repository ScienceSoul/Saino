//
//  FEMLinearAlgebra.h
//  Saino
//
//  Created by Hakime Seddik on 01/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMLinearAlgebra : NSObject

-(void)LUDecompositionOfMatrix:(double * __nonnull * __nonnull)a ofSize:(int)n resultPivot:(int * __nonnull)pivot;
-(void)invertMatrix:(double * __nonnull * __nonnull)a ofSize:(int)n;

@end
