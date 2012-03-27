//
//  FEMUtilities.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMUtilities : NSObject


-(double)interpolateCurve:(double *)tValues: (double *)fValues: (double)t: (int)n;

@end
