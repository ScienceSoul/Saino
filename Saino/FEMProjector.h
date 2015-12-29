//
//  FEMProjector.h
//  Saino
//
//  Created by Seddik hakime on 09/08/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

#import "FEMMatrix.h"
#import "FEMMesh.h"

@interface FEMProjector : NSObject {
    
    FEMMesh *_mesh;
    FEMMatrix *_matrix;
    FEMMatrix *_tMatrix;
    FEMProjector *_next;
}

@property(nonatomic, strong, nonnull) FEMMesh *mesh;
@property(nonatomic, strong, nonnull) FEMMatrix *matrix;
@property(nonatomic, strong, nonnull) FEMMatrix *tMatrix;
@property(nonatomic, strong, nullable) FEMProjector *next;

@end
