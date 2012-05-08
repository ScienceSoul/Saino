//
//  FEMBodyForce.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMBodyForce : NSObject {
    
    NSMutableArray *valuesList;  // Array of FEMListValue objects 
}

-(NSMutableArray *)returnValuesList;

@end
