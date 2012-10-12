//
//  DirectoryReader.m
//  Saino
//
//  Created by Hakime Seddik on 22/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//
//  Originally written by Dave DeLong, 
//  Source: http://stackoverflow.com/questions/3707427#3711079

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
- (id)initWithPath:(NSString*)path {
    
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
- (BOOL)readDirectory:(NSArray**)files {
    
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

- (BOOL)readDirectoryWithFileAttributes:(NSArray**)fileAttributes{
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

- (BOOL)readDirectoryByFileModificationDateDesc:(NSArray**)files{
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
