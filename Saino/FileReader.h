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
    
    NSString*			m_filePath;				/**< File path. */
	NSFileHandle*		m_fileHandle;			/**< File handle. */
	unsigned long long	m_currentOffset;		/**< Current offset is needed for forwards reading. */
	unsigned long long	m_currentInset;			/**< Current inset is needed for backwards reading. */
	NSRange				m_prevDelimiterRange;	/**< Position and length of the last delimiter. */	
	unsigned long long	m_totalFileLength;		/**< Total number of bytes in file. */
	NSString*			m_lineDelimiter;		/**< Character for line break or page break. */
	NSUInteger			m_chunkSize;			/**< Standard block size. */
}

- (id)initWithFilePath:(NSString*)filePath;
- (void)closeHandle;
- (void)rewind;
- (NSString*)readLine;
- (NSString*)readLineBackwards;
- (NSString*)readTrimmedLine;

#if NS_BLOCKS_AVAILABLE
- (void)enumerateLinesUsingBlock:(void(^)(NSString*, BOOL*))block;
#endif

@end
