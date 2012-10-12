//
//  FEMJob.h
//  Saino
//
//  Created by Seddik hakime on 09/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "FEMModel.h"

@interface FEMJob : NSObject {
    
    FEMModel *_model;
}

@property(nonatomic, strong) FEMModel *model;

@end
