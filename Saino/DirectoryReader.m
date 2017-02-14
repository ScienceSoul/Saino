//===----------------------------------------------------------------------===//
//  DirectoryReader.m
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

#import "DirectoryReader.h"

/**
 A directory reader.
 */
@implementation DirectoryReader

#define kFileAttributeFileNameKey @"fileName"
/**
 Initializes a directory reader.
 @param path A directory path.
 @returns An initialized DirectoryReader object or nil if the object could not be created.
 */
- (id __nullable)initWithPath:(NSString * __nonnull)path {
    
	self = [super init];
	if (self != nil) {
		if (!path || [path length] <= 0) {
			return nil;
		}
		// Remove trailing slash if appended.
		NSMutableString* mutablePath = [NSMutableString stringWithString:path];
		[mutablePath replaceOccurrencesOfString:@"/" withString:@"" options:NSBackwardsSearch range:NSMakeRange([path length] - 1, 1)];
		m_path = mutablePath;
	}
	return self;
}

/**
 Reads the content of a directory.
 @param files Container for the listing.
 @returns YES if the directory was read otherwise NO.
 */
- (BOOL)readDirectory:(NSArray * __nonnull * __nullable)files {
    
	BOOL success = NO;
	
	NSArray* fileNames = [[NSFileManager defaultManager] contentsOfDirectoryAtPath:m_path error:nil];
	
	if (!fileNames || [fileNames count] <= 0) {
		return success;
	}
	
	NSMutableArray* fullNames = [NSMutableArray arrayWithCapacity:[fileNames count]];
	for (NSString* fileName in fileNames) {
		NSString* fullPath = [m_path stringByAppendingString:[NSString stringWithFormat:@"/%@", fileName]];
		[fullNames addObject:fullPath];
	}
	
	if (fullNames && [fullNames count] > 0 && files) {
		*files = fullNames;
		success = YES;
	}
	return success;
}

- (BOOL)readDirectoryWithFileAttributes:(NSArray * __nonnull * __nullable)fileAttributes {
    NSArray *fileNames;
    BOOL success = [self readDirectory:&fileNames];
    
    if (!success) // couldn't read directory
        return NO;
    
    NSError *fileIoError;
    // create array to gather attributes
    NSMutableArray *fileAttributesResults = [NSMutableArray arrayWithCapacity:[fileNames count]];
    for (NSString *fileName in fileNames){
        NSMutableDictionary *extendedDictionary = [NSMutableDictionary dictionaryWithDictionary:[[NSFileManager defaultManager] attributesOfItemAtPath:fileName error:&fileIoError]];
        extendedDictionary[kFileAttributeFileNameKey] = fileName;
        [fileAttributesResults addObject:extendedDictionary];
    }
    
    *fileAttributes = [NSArray arrayWithArray:fileAttributesResults];
    return YES;
}

- (BOOL)readDirectoryByFileModificationDateDesc:(NSArray * __nonnull * __nullable)files {
    NSArray *fileAttributes;
    BOOL success = [self readDirectoryWithFileAttributes:&fileAttributes];
    if (!success) // couldn't read directory
        return NO;
    NSSortDescriptor *sortDescriptor = [[NSSortDescriptor alloc] initWithKey:@"fileModificationDate" ascending:NO];
    NSArray *sortedAttributes = [fileAttributes sortedArrayUsingDescriptors:@[sortDescriptor]];
    NSMutableArray *fileNames = [NSMutableArray arrayWithCapacity:[fileAttributes count]];
    for( NSDictionary *fileAttributeDictionary in sortedAttributes){
        [fileNames addObject:fileAttributeDictionary[kFileAttributeFileNameKey]];
    }
    
    *files = [NSArray arrayWithArray:fileNames];
    return YES;
}

@end
