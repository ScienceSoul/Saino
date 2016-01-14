//
//  FEMEquation.h
//  Saino
//
//  Created by Seddik hakime on 13/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMValueList.h"

@interface FEMEquation : NSObject {
    
    NSMutableArray *_valuesList;
}

@property(nonatomic, strong, nonnull) NSMutableArray <FEMValueList *> *valuesList;

@end
