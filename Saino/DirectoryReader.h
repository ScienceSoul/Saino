//
//  DirectoryReader.h
//  Saino
//
//  Created by Hakime Seddik on 22/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//
//  Originally written by Dave DeLong, 
//  Source: http://stackoverflow.com/questions/3707427#3711079

#import <Foundation/Foundation.h>

@interface DirectoryReader : NSObject {
    
	NSString * __nullable m_path;		/**< File path. */
}

- (id __nullable)initWithPath:(NSString * __nonnull)path;
- (BOOL)readDirectory:(NSArray * __nonnull * __nullable)files;
- (BOOL)readDirectoryWithFileAttributes:(NSArray * __nonnull * __nullable)fileAttributes;

/**
 Returns array of string paths for each file in the directory, sorted by file modification date, descending
 @param files output parameter containing the file paths
 @returns BOOL indicating whether the directory file list could be read
 */
- (BOOL)readDirectoryByFileModificationDateDesc:(NSArray * __nonnull * __nullable)files;
@end

