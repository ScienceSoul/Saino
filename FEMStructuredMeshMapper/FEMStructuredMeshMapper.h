//
//  FEMStructuredMeshMapper.h
//  Saino
//
//  Created by Seddik hakime on 13/10/2015.
//  Copyright © 2015 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import <SainoCore/SainoCore.h>

// This class implements a plug-in for mapping the mesh between given top and
// bottom surfaces. This solver assumes that the mesh is structural so that it
// could have been obtained by extrusion in the direction of interest. For the
// given direction the corresponding top and bottom node is computed for every
// node and this information is used to perform linear mapping in between.


@interface FEMStructuredMeshMapper : NSObject <SainoSolutionsComputer>

@end
