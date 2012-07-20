//
//  SIOInfoParallel.h
//  Saino
//
//  Created by Hakime Seddik on 09/07/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface SIOInfoParallel : NSObject {
@private
    
    BOOL _parallel;
    int _numProc;
    int _myProc;
}

@property (nonatomic, assign, getter=isParallel) BOOL parallel;
@property(nonatomic, assign) int numProc;
@property(nonatomic, assign) int myProc;

@end
