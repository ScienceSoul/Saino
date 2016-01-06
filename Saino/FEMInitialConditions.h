//
//  FEMInitialConditions.h
//  Saino
//
//  Created by Seddik hakime on 13/03/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMInitialConditions : NSObject {
    
    NSMutableArray *_valuesList;  // Array of FEMValueList objects
}

@property(nonatomic, strong, nonnull) NSMutableArray *valuesList;

@end
