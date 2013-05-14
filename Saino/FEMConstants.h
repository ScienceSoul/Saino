//
//  FEMConstants.h
//  Saino
//
//  Created by Seddik hakime on 08/05/13.
//  Copyright (c) 2013 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface FEMConstants : NSObject {

    NSMutableArray *_valuesList;  // Array of FEMValueList objects
}

@property(nonatomic, strong) NSMutableArray *valuesList;

@end
