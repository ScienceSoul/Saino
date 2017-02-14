//===----------------------------------------------------------------------===//
//  FEMFlowSolution.h
//  Saino
//
//  Created by Seddik hakime on 15/03/13.
//  Copyright (c) 2013 ScienceSoul. All rights reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//===----------------------------------------------------------------------===//

#import <Foundation/Foundation.h>
#import "SainoSolutionsComputer.h"

// Class computing Navier-Stokes local matrices in general coordinate system
// (i.e. not cartesian, axisymmetric or cylindrically symmetric.

// Class corresponds to Elmer from git on October 27 2015

@interface FEMFlowSolution : NSObject <SainoSolutionsComputer>

@end
