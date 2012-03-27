//
//  FEMBoundaryCondition.h
//  Saino
//
//  Created by Hakime Seddik on 27/03/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "Constructors.h"

@interface FEMBoundaryCondition : NSObject {
    
    int tag;
    NSArray *valuesList;  // Array of FEMListValue objects 
    Matrix_t *pMatrix;

}

-(int)tag;
-(void)setTag:(int)n;
-(NSArray *)returnValuesList;

@end
