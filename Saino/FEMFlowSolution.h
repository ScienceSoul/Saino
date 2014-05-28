//
//  FEMFlowSolution.h
//  Saino
//
//  Created by Seddik hakime on 15/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "SainoSolutionsComputer.h"

// Class computing Navier-Stokes local matrices in general coordinate system
// (i.e. not cartesian, axisymmetric or cylindrically symmetric.

@interface FEMFlowSolution : NSObject <SainoSolutionsComputer>

@end
