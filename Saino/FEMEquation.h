//
//  FEMEquation.h
//  Saino
//
//  Created by Seddik hakime on 13/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMEquation : NSObject {
    
    NSMutableArray *_valuesList;  // Array of FEMValueList objects
}

@property(nonatomic, strong) NSMutableArray *valuesList;

@end
