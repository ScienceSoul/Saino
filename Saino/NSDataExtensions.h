//===----------------------------------------------------------------------===//
//  NSDataExtensions.h
//  Saino
//
//  Created by Hakime Seddik on 22/06/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
//  Original version written by Dave DeLong,
//  Source: http://stackoverflow.com/questions/3707427#3711079
//
//===----------------------------------------------------------------------===//
// NSData additions.
//===----------------------------------------------------------------------===//

#import <Foundation/Foundation.h>

@interface NSData (Additions)

- (NSRange)rangeOfData:(NSData * __nonnull)dataToFind;
- (NSRange)rangeOfDataBackwardsSearch:(NSData * __nonnull)dataToFind;
- (NSString * __nullable)stringValueWithEncoding:(NSStringEncoding)encoding;

@end


// -----------------------------------------------------------------------------
// NSMutableData additions.
// -----------------------------------------------------------------------------

@interface NSMutableData (Additions)

- (void)prepend:(NSData * __nonnull)data;

@end
