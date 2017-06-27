//===----------------------------------------------------------------------===//
//  NSDataExtensions.m
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
//===----------------------------------------------------------------------===//

#import "NSDataExtensions.h"

//===----------------------------------------------------------------------===//
// NSData additions.
//===----------------------------------------------------------------------===//

/**
 Extension of the NSData class. 
 Data can be found forwards or backwards. Further the extension supplies a function 
 to convert the contents to string for debugging purposes.
 @param Additions Category labeled Additions.
 @returns An initialized NSData object or nil if the object could not be created.
 */
@implementation NSData (Additions)

/**
 Returns a range of data.
 @param dataToFind Data object specifying the delimiter and encoding.
 @returns A range.
 */
- (NSRange)rangeOfData:(NSData * _Nonnull)dataToFind {
	
	const void* bytes = [self bytes];
	NSUInteger length = [self length];
	const void* searchBytes = [dataToFind bytes];
	NSUInteger searchLength = [dataToFind length];
	NSUInteger searchIndex = 0;
	
	NSRange foundRange = {NSNotFound, searchLength};
	for (NSUInteger index = 0; index < length; index++) {
		// The current character matches.
		if (((char*)bytes)[index] == ((char*)searchBytes)[searchIndex]) {
			// Store found location if not done earlier.
			if (foundRange.location == NSNotFound) {
				foundRange.location = index;
			}
			// Increment search character index to check for match.
			searchIndex++;
			// All search character match.
			// Break search routine and return found position.
			if (searchIndex >= searchLength) {
				return foundRange;
			}
		}
		// Match does not continue.
		// Return to the first search character.
		// Discard former found location.
		else {
			searchIndex = 0;
			foundRange.location = NSNotFound;
		}
	}
	return foundRange;
}


- (NSRange)rangeOfDataBackwardsSearch:(NSData * _Nonnull)dataToFind {
	
	const void* bytes = [self bytes];
	NSUInteger length = [self length];
	const void* searchBytes = [dataToFind bytes];
	NSUInteger searchLength = [dataToFind length];
	NSUInteger searchIndex = 0;
	
	NSRange foundRange = {NSNotFound, searchLength};
	if (length < searchLength) {
		return foundRange;
	}
	for (NSInteger index = length - searchLength; index >= 0;) {
        //		NSLog(@"%c == %c", ((char*)bytes)[index], ((char*)searchBytes)[searchIndex]); /* DEBUG LOG */
		if (((char*)bytes)[index] == ((char*)searchBytes)[searchIndex]) {
			// The current character matches.
			if (foundRange.location == NSNotFound) {
				foundRange.location = index;
			}
			index++;
			searchIndex++;
			if (searchIndex >= searchLength) {
				return foundRange;
			}
		}
		else {
			// Decrement to search backwards.
			if (foundRange.location == NSNotFound) {
				// Skip if first byte has been reached.
				if (index == 0) {
					foundRange.location = NSNotFound;
					return foundRange;
				}
				index--;
			}
			// Jump over the former found location
			// to avoid endless loop.
			else {
				index = index - 2;
			}
			searchIndex = 0;
			foundRange.location = NSNotFound;
		}
	}
	return foundRange;
}

- (NSString * _Nullable)stringValueWithEncoding:(NSStringEncoding)encoding {
	return [[NSString alloc] initWithData:self encoding:encoding];
}

@end


// -----------------------------------------------------------------------------
// NSMutableData additions.
// -----------------------------------------------------------------------------

/**
 Extension of the NSMutableData class. 
 Data can be prepended in addition to the append function of the framework.
 @param Additions Category labeled Additions.
 @returns An initialized NSMutableData object or nil if the object could not be created.
 */
@implementation NSMutableData (Additions)

/**
 Inserts the data before the data of the object.
 @param data Data to be prepended.
 */
- (void)prepend:(NSData * _Nonnull)data {
    
	
	NSMutableData* concat = [NSMutableData dataWithData:data];
	[concat appendData:self];
	[self setData:concat];
}

@end
