//
//  FEMFreeSurface.h
//  Saino
//
//  Created by Seddik hakime on 29/05/2014.
//  Copyright (c) 2014 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMModel.h"
#import "FEMNumericIntegration.h"

@interface FEMFreeSurface : NSObject

-(void)moveBoundaryModel:(FEMModel * __nonnull)model integration:(FEMNumericIntegration * __nonnull)integration relax:(double)relax;

@end
