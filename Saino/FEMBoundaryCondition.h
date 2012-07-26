//
//  FEMBoundaryCondition.h
//  Saino
//
//  Created by Hakime Seddik on 27/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMatrix.h"
#import "Constructors.h"

@interface FEMBoundaryCondition : NSObject {
    
    int tag;
    NSMutableArray *valuesList;  // Array of FEMValueList objects 
    FEMMatrix *pMatrix;

}

-(int)tag;
-(void)setTag:(int)n;
-(NSMutableArray *)returnValuesList;
-(FEMMatrix *)returnPMatrix;

@end
