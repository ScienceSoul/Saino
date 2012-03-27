//
//  FEMBoundaryCondition.m
//  Saino
//
//  Created by Hakime Seddik on 27/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMBoundaryCondition.h"

@implementation FEMBoundaryCondition

-(int)tag {
    
    return tag;
}

-(void)setTag:(int)n {
    
    tag = n;
}

-(NSArray *)returnValuesList {
    
    return valuesList;
}


@end
