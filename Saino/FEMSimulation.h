//
//  FEMSimulation.h
//  Saino
//
//  Created by Hakime Seddik on 26/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMSimulation : NSObject {
    
    int _tag;
    NSMutableArray *_valuesList;  // Array of FEMValueList objects
}

@property(nonatomic, assign) int tag;
@property(nonatomic, strong, nonnull) NSMutableArray *valuesList;

@end
