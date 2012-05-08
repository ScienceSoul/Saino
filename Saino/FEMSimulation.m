//
//  FEMSimulation.m
//  Saino
//
//  Created by Hakime Seddik on 26/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "FEMSimulation.h"

@implementation FEMSimulation

-(int)tag {
    
    return tag;
}

-(void)setTag:(int)n {
    
    tag = n;
}

-(NSMutableArray *)returnValuesList {
    
    return valuesList;
}


@end
