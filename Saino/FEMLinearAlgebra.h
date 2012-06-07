//
//  FEMLinearAlgebra.h
//  Saino
//
//  Created by Hakime Seddik on 01/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMLinearAlgebra : NSObject

-(void)LUDecompositionOfMatrix:(double **)a ofSize:(int)n resultPivot:(int *)pivot;
-(void)invertMatrix:(double **)a ofSize:(int)n;

@end
