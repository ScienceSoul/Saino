//
//  FEMSimulation.h
//  Saino
//
//  Created by Hakime Seddik on 26/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMSimulation : NSObject {
    
    int tag;
    NSMutableArray *valuesList;  // Array of FEMValueList objects 

}

-(int)tag;
-(void)setTag:(int)n;
-(NSMutableArray *)returnValuesList;

@end
