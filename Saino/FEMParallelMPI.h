//
//  FEMParallelMPI.h
//  Saino
//
//  Created by Hakime Seddik on 09/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMParallelMPI : NSObject

-(double)parallelReductionOfValue:(double)r operArg:(int *)oper_arg;

@end
