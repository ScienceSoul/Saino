//
//  FEMBodyForce.h
//  Saino
//
//  Created by Hakime Seddik on 16/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMValueList.h"

@interface FEMBodyForce : NSObject {
    
    NSMutableArray *_valuesList;
}

@property(nonatomic, strong, nonnull) NSMutableArray <FEMValueList *> *valuesList;

@end
