//
//  FileReader.h
//  Saino
//
//  Created by Hakime Seddik on 22/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//
//  Originally written by Dave DeLong, 
//  Source: http://stackoverflow.com/questions/3707427#3711079

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
