//===----------------------------------------------------------------------===//
//  DirectoryReader.h
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

#import <Foundation/Foundation.h>

@interface DirectoryReader : NSObject {
    
	NSString * _Nullable m_path;		/**< File path. */
}

- (id _Nullable)initWithPath:(NSString * _Nonnull)path;
- (BOOL)readDirectory:(NSArray * _Nonnull * _Nullable)files;
- (BOOL)readDirectoryWithFileAttributes:(NSArray * _Nonnull * _Nullable)fileAttributes;

/**
 Returns array of string paths for each file in the directory, sorted by file modification date, descending
 @param files output parameter containing the file paths
 @returns BOOL indicating whether the directory file list could be read
 */
- (BOOL)readDirectoryByFileModificationDateDesc:(NSArray * _Nonnull * _Nullable)files;
@end

