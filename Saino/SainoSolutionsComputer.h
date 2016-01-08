//
//  SainoFieldSolutionsComputing.h
//  Saino
//
//  Created by Seddik hakime on 15/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"
#import "FEMSolution.h"

@protocol SainoSolutionsComputer <NSObject>
-(void)solutionComputer:(FEMSolution * __nonnull)solution model:(FEMModel * __nonnull)model timeStep:(int)timeStep transientSimulation:(BOOL)transient;
-(void)deallocation:(FEMSolution * __nonnull)solution;
@end