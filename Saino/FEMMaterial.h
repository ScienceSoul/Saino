//
//  FEMMaterial.h
//  Saino
//
//  Created by Seddik hakime on 17/05/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMValueList.h"

@interface FEMMaterial : NSObject {

    NSMutableArray *_valuesList;
}

@property(nonatomic, strong, nonnull) NSMutableArray <FEMValueList *> *valuesList;

@end

