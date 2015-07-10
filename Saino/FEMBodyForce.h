//
//  FEMBodyForce.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMBodyForce : NSObject {
    
    NSMutableArray *_valuesList;  // Array of FEMValueList objects
}

@property(nonatomic, strong) NSMutableArray *valuesList;

@end
