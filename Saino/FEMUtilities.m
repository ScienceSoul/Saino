//
//  FEMUtilities.m
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMUtilities.h"

@implementation FEMUtilities

-(double)interpolateCurve:(double *)tValues: (double *)fValues: (double)t: (int)n {
    
    int i;
    double f;
    
    for (i=0; i<n; i++) {
        if (tValues[i] >= t) break;
    }
    if (i > n-1) i = n-1;
    if (i < 1) i = 1;
    
    f = (t - tValues[i-1]) / (tValues[i] - tValues[i-1]);
    f = (1.0 - f)*fValues[i-1] + f*fValues[i];
    
    return f;
}

@end
