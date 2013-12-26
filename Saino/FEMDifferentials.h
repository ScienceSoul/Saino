//
//  FEMDifferentials.h
//  Saino
//
//  Created by Seddik hakime on 18/06/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"
#import "FEMMesh.h"
#import "FEMListUtilities.h"

/*******************************************************************************
    This class contains some built-in material laws, and also some
    vector utilities, curl, dot, cross, etc. Some of these may be 
    of no use currently.
*******************************************************************************/
@interface FEMDifferentials : NSObject

-(double)jouleHeatElement:(Element_t *)element nodes:(Nodes_t *)nodes numberOfNodes:(int)n integrationU:(double)u integrationV:(double)v integrationW:(double)w mesh:(FEMMesh *)mesh model:(FEMModel *)model listUtilities:(FEMListUtilities *)listUtil;

@end
