//
//  FEMSimulation.h
//  Saino
//
//  Created by Hakime Seddik on 26/04/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMValueList.h"

@interface FEMSimulation : NSObject {
    
    int _tag;
    NSMutableArray *_valuesList;
}

@property(nonatomic, assign) int tag;
@property(nonatomic, strong, nonnull) NSMutableArray <FEMValueList *> *valuesList;

@end
