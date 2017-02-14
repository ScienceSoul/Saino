//===----------------------------------------------------------------------===//
//  FEMStructuredMeshMapper.h
//  Saino
//
//  Created by Seddik hakime on 13/10/2015.
//  Copyright © 2015 ScienceSoul. All rights reserved.
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
#import <SainoCore/SainoCore.h>

// This class implements a plug-in for mapping the mesh between given top and
// bottom surfaces. This solver assumes that the mesh is structural so that it
// could have been obtained by extrusion in the direction of interest. For the
// given direction the corresponding top and bottom node is computed for every
// node and this information is used to perform linear mapping in between.


@interface FEMStructuredMeshMapper : NSObject <SainoSolutionsComputer>

@end
