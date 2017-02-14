//===----------------------------------------------------------------------===//
//  FileReader.h
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


@interface FileReader : NSObject {
    
    NSString * __nonnull		m_filePath;				/**< File path. */
	NSFileHandle * __nonnull	m_fileHandle;			/**< File handle. */
	unsigned long long          m_currentOffset;		/**< Current offset is needed for forwards reading. */
	unsigned long long          m_currentInset;			/**< Current inset is needed for backwards reading. */
	NSRange                     m_prevDelimiterRange;	/**< Position and length of the last delimiter. */
	unsigned long long          m_totalFileLength;		/**< Total number of bytes in file. */
	NSString * __nonnull		m_lineDelimiter;		/**< Character for line break or page break. */
	NSUInteger              m_chunkSize;			/**< Standard block size. */
}

- (id __nullable)initWithFilePath:(NSString * __nonnull)filePathh;
- (void)closeHandle;
- (void)rewind;
- (NSString * __nullable)readLine;
- (NSString * __nullable)readLineBackwards;
- (NSString * __nullable)readTrimmedLine;

#if NS_BLOCKS_AVAILABLE
- (void)enumerateLinesUsingBlock:(void(^ __nonnull)(NSString * __nonnull, BOOL * __nonnull))block;
#endif

@end
