//
//  NSDataExtensions.h
//  Saino
//
//  Created by Hakime Seddik on 22/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//
//  Originally written by Dave DeLong, 
//  Source: http://stackoverflow.com/questions/3707427#3711079

// -----------------------------------------------------------------------------
// NSData additions.
// -----------------------------------------------------------------------------

#import <Foundation/Foundation.h>

@interface NSData (Additions)

- (NSRange)rangeOfData:(NSData*)dataToFind;
- (NSRange)rangeOfDataBackwardsSearch:(NSData*)dataToFind;
- (NSString*)stringValueWithEncoding:(NSStringEncoding)encoding;

@end


// -----------------------------------------------------------------------------
// NSMutableData additions.
// -----------------------------------------------------------------------------

@interface NSMutableData (Additions)

- (void)prepend:(NSData*)data;

@end
