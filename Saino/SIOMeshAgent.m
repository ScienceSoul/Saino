//
//  SIOMeshAgent.m
//  Saino
//
//  Created by Hakime Seddik on 21/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#import "SIOMeshAgent.h"

#import "memory.h"

static int step = 0;

static char **_extension = (char **)0;

enum {HEADER = 0, NODES, ELEMENTS, BOUNDARY, SHARED};

static char *_sequential_extensions[] = {
    "/mesh.header",
    "/mesh.nodes",
    "/mesh.elements",
    "/mesh.boundary"
};

static char *_parallel_extensions[] = {
    "%s/part.%d.header",
    "%s/part.%d.nodes",
    "%s/part.%d.elements",
    "%s/part.%d.boundary",
    "%s/part.%d.shared",
};

@interface SIOMeshAgent ()

-(void)SIOMeshAgent_makeFilename:(char * __nonnull)buf model:(const char * __nonnull)model suffix:(const char * __nonnull)suffix;
-(int)SIOMeshAgent_elementNodes:(const int)type;
-(void)SIOMeshAgent_cacheNodes;
-(int)SIOMeshAgent_copyCoords:(double * __nonnull)target address:(const int)address;
-(cacheNode *)SIOMeshAgent_searchNode:(const int)address;

@end

@implementation SIOMeshAgent {
    
    char _newdir[PATH_MAX];
    cacheNode *_clist;
}

@synthesize manager = _manager;
@synthesize meshFileStreams = _meshFileStreams;
@synthesize elementTypeTags = _elementTypeTags;
@synthesize elementTypeCount = _elementTypeCount;
@synthesize parts = _parts;
@synthesize me = _me;
@synthesize nodeCount = _nodeCount;
@synthesize elementCount = _elementCount;
@synthesize boundaryElementCount = _boundaryElementCount;
@synthesize elementTypes = _elementTypes;
@synthesize sharedNodeCount = _sharedNodeCount;
@synthesize borderElementCount = _borderElementCount;
@synthesize dim = _dim;
@synthesize parallel = _parallel;
@synthesize meshFiles = _meshFiles;

#pragma mark Private methods

-(void)SIOMeshAgent_makeFilename:(char * __nonnull)buf model:(const char * __nonnull)model suffix:(const char * __nonnull)suffix {
    
    buf[0] = '\0';
    strcat(buf, model);
    strcat(buf, suffix);
}

-(int)SIOMeshAgent_elementNodes:(const int)type {
    
    int cnt;
    
    cnt = type - 100*(type/100);
    return cnt;
}

-(void)SIOMeshAgent_cacheNodes {
    
    FileReader *reader;
    NSString *line;
    
    // Used to separate strings and filter them from white spaces
    NSCharacterSet *whitespaces = [NSCharacterSet whitespaceCharacterSet];
    NSPredicate *noEmptyStrings = [NSPredicate predicateWithFormat:@"SELF != ''"];
    
    if (!_clist) {
        _clist = (cacheNode*) malloc( sizeof(cacheNode) * self.nodeCount );
        
        line = nil;
        reader = (self.meshFileStreams)[NODES];
        for (int i=0; i<self.nodeCount; i++) {
            line = [reader readLine];
            NSArray *stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
            NSArray *filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
            if (self.parallel) { // assumes that everything is sorted by splitter 
                _clist[i].tag = [filteredArray[0] intValue];
                _clist[i].constraint = [filteredArray[1] intValue];
                _clist[i].x = [filteredArray[2] doubleValue];
                _clist[i].y = [filteredArray[3] doubleValue];
                _clist[i].z = [filteredArray[4] doubleValue];                
            } else {
                int tag = [filteredArray[0] intValue];
                _clist[tag-1].tag = tag;
                _clist[tag-1].constraint = [filteredArray[1] intValue];
                _clist[tag-1].x = [filteredArray[2] doubleValue];
                _clist[tag-1].y = [filteredArray[3] doubleValue];
                _clist[tag-1].z = [filteredArray[4] doubleValue];
            }
        }
        [reader rewind];
    }
}

-(cacheNode *)SIOMeshAgent_searchNode:(const int)address {
    
    cacheNode entry;
    cacheNode *retval;
    
    entry.tag = address;
    retval = (cacheNode *)bsearch((void *)&entry, (void *)_clist, self.nodeCount, sizeof(cacheNode), nodecomp);
    return retval;
}

-(int)SIOMeshAgent_copyCoords:(double * __nonnull)target address:(const int)address {
    
    int found = 1;
    if (self.parallel) {
        cacheNode *retval = [self SIOMeshAgent_searchNode:address];
        if (retval == NULL) {
            NSLog(@"SIOMeshAgent:SIOMeshAgent_copyCoords: partition error: something is going totally wrong. Address: %d.\n", address);
            found = 0;
        } else {
            target[0] = retval->x;
            target[1] = retval->y;
            target[2] = retval->z;
        }
    }
    else {
        int offset = address - 1;
        target[0] = _clist[offset].x;
        target[1] = _clist[offset].y;
        target[2] = _clist[offset].z;
    }
    return found;
}

#pragma mark Public methods

-(id __nullable)initWithManager:(SIOModelManager * __nonnull)mm split:(int)split part:(int)part
{
    self = [super init];
    if (self) {
        if (mm == nil) {
            return nil;
        }
        _manager = mm;
        _parts = split;
        _me = part;
        if (_me > 0) {
            _parallel = 1;
        } else {
            _parallel = 0;
        }
        
        _dim = 3;
        _clist = NULL;
        _elementTypeTags = NULL;
        _elementTypeCount = NULL;
        
        if (_parallel) {
            _meshFiles = 5;
            _extension = _parallel_extensions;
        } else {
            _meshFiles = 4;
            _extension = _sequential_extensions;
        }
        
        _meshFileStreams = [[NSMutableArray alloc] init];
        _elementTypeTags = [[NSMutableArray alloc] init];
        _elementTypeCount = [[NSMutableArray alloc] init];
    }
    
    return self;
}
    
-(int)createMesh:(NSString * __nonnull)dir {
    
    int i;
    char filename[PATH_MAX];
    NSFileHandle *file;
    
    for (i=0; i<self.meshFiles; i++) {
        [self SIOMeshAgent_makeFilename:filename model:[dir UTF8String] suffix:_extension[i]];
        [self.manager openStream:file name:@(filename) mode:@"write"];
        [self.meshFileStreams insertObject:file atIndex:i];
    }
    
    return 0;
}

-(int)openMesh:(NSString * __nonnull)dir {
    
    int i, j, isLineBreak, lineCount = 0;
    char filename[PATH_MAX];
    FileReader *reader;
    NSString *line;
    
    for (i=0; i<self.meshFiles; i++) {
        
        if (self.parallel) {
            snprintf(_newdir, sizeof(_newdir), "%s/partitioning.%d", [dir UTF8String], self.parts);
            snprintf(filename, sizeof(filename), _extension[i], _newdir, self.me);
            
        } else  [self SIOMeshAgent_makeFilename:filename model:[dir UTF8String] suffix:_extension[i]];
        
        reader = [[FileReader alloc] initWithFilePath:@(filename)];
        if (!reader) {
            return -1;
        } else {
            [self.meshFileStreams insertObject:reader atIndex:i];
        }
        
    }
    
    // Used to separate strings and filter them from white spaces
    NSCharacterSet *whitespaces = [NSCharacterSet whitespaceCharacterSet];
    NSPredicate *noEmptyStrings = [NSPredicate predicateWithFormat:@"SELF != ''"];
    
    // Read Header
    reader = (self.meshFileStreams)[HEADER];
    line = nil;
    j = 0;
    while ((line = [reader readLine])) {
        lineCount++;
        NSLog(@"SIOMeshAgent:openMesh: %3.d: %@", lineCount, line);
        // Parse the line
        NSArray *stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
        NSArray *filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
        isLineBreak = 0;
        for (NSString *string in filteredArray) {
            if ([string isEqualToString:@"\n"] == YES) isLineBreak++;
        }
        if (([filteredArray count]-isLineBreak) == 3) { // First line of header file, three elements
            self.nodeCount = [filteredArray[0] intValue];
            self.elementCount = [filteredArray[1] intValue];
            self.boundaryElementCount = [filteredArray[2] intValue];
        } else if (([filteredArray count]-isLineBreak) == 1) { // Second line, one element
            self.elementTypes = [filteredArray[0] intValue];
        } else if (([filteredArray count]-isLineBreak) == 2) { // The rest of the file, two elements
            [self.elementTypeTags addObject:filteredArray[0]];
            [self.elementTypeCount addObject:filteredArray[1]];
            j++;
        } else if (self.parallel && lineCount == (2+self.elementTypes)+1) { // In case of a parallel mesh
            self.sharedNodeCount = [filteredArray[0] intValue];
            self.borderElementCount = [filteredArray[1] intValue];
        }
    }
    
    step = 0;
    _clist = NULL;
    
    return 0;
}

-(int)closeMesh {
    
    int i;
    FileReader *reader;
    
    for (i=0; i<self.meshFiles; i++) {
        reader = (self.meshFileStreams)[i];
        [reader closeHandle];
    }

    free(_clist);
    _clist = NULL;
    
    return 0;
}

-(int)readDescriptorNode:(int * __nonnull)nodeC element:(int * __nonnull)elementC boundaryElement:(int * __nonnull)boundaryElementC usedElementTypes:(int * __nonnull)usedElementTypes usedElementTypeTags:(int * __nonnull)usedElementTypeTags usedElementTypeCount:(int * __nonnull)usedElementTypeCount {
    
    int i;
    
    *nodeC = self.nodeCount;
    *elementC = self.elementCount;
    *boundaryElementC = self.boundaryElementCount;
    *usedElementTypes = self.elementTypes;
    
    for (i=0; i<self.elementTypes; i++) {
        usedElementTypeTags[i] = [(self.elementTypeTags)[i] intValue];
        usedElementTypeCount[i] = [(self.elementTypeCount)[i] intValue];
    }
    
    return 0;
}

-(int)readPartDescriptor:(int * __nonnull)shared {
    
    *shared = self.sharedNodeCount;
    return 0;
}

-(int)readNextElementConnections:(int * __nonnull)tag part:(int * __nonnull)part body:(int * __nonnull)body type:(int * __nonnull)type pdofs:(int * __nonnull)pdofs nodes:(int * __nonnull)nodes colorIndex:(int * __nullable)colorIndex parallelAssembly:(BOOL * __nullable)parallelAssembly {
    
    int i, j, gotnodal;
    FileReader *reader;
    NSString *line, *tagstr, *typestr, *str1, *str2;
    NSRange subRange;
    char character[2];
    
    // Used to separate strings and filter them from white spaces
    NSCharacterSet *whitespaces = [NSCharacterSet whitespaceCharacterSet];
    NSPredicate *noEmptyStrings = [NSPredicate predicateWithFormat:@"SELF != ''"];

    line = nil;
    reader = (self.meshFileStreams)[ELEMENTS];
    if (step == self.elementCount) {
        [reader rewind];
        step = 0;
        return -1;
    }
    line = [reader readLine];
    NSArray *stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
    NSArray *filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
    tagstr = filteredArray[0];
    *body = [filteredArray[1] intValue];
    typestr = filteredArray[2];
    part = 0;
    subRange = [tagstr rangeOfString:@"/"];
    if (subRange.location != NSNotFound) {
        NSArray *components = [tagstr componentsSeparatedByString:@"/"];
        *tag = [components[0] intValue];
        *part = [components[1] intValue];        
    } else {
        *tag = [tagstr intValue];
    }
    
    // TODO: is this secure? What if pdofs is a smaller buffer?
    memset( pdofs, 0, 6*sizeof(int) );
    gotnodal = 0;
    for (i=0; i<[typestr length]; i++) {
        character[0] = [typestr characterAtIndex:i];
        str1 = @(character);
        if ([str1 isEqualToString:@"n"]) {
            character[0] = [typestr characterAtIndex:i+1];
            str2 = @(character);
            pdofs[0] = [str2 intValue];
            gotnodal = 1;
        } else if ([str1 isEqualToString:@"e"]) {
            character[0] = [typestr characterAtIndex:i+1];
            str2 = @(character);
            pdofs[1] = [str2 intValue];
        } else if ([str1 isEqualToString:@"f"]) {
            character[0] = [typestr characterAtIndex:i+1];
            str2 = @(character);
            pdofs[2] = [str2 intValue];
        } else if ([str1 isEqualToString:@"d"]) {
            character[0] = [typestr characterAtIndex:i+1];
            str2 = @(character);
            pdofs[3] = [str2 intValue];
        } else if ([str1 isEqualToString:@"b"]) {
            character[0] = [typestr characterAtIndex:i+1];
            str2 = @(character);
            pdofs[4] = [str2 intValue];
        } else if ([str1 isEqualToString:@"p"]) {
            character[0] = [typestr characterAtIndex:i+1];
            str2 = @(character);
            pdofs[5] = [str2 intValue];
        }
    }
    *type = [typestr intValue];
    
    int elNodes = [self SIOMeshAgent_elementNodes:*type];
    j = 3;
    for (i=0; i<elNodes; i++) {
        nodes[i] = [filteredArray[j] intValue];
        j++;
    }
    
    if (parallelAssembly != NULL) {
        if (*parallelAssembly == YES) {
            if (colorIndex != NULL) *colorIndex = [filteredArray[j] intValue];
        }
    }
    
    if (!gotnodal) pdofs[0] = 1;
    
    step++;

    return 0;
}

-(int)readNextElementCoordinates:(int * __nonnull)tag body:(int * __nonnull)body type:(int * __nonnull)type nodes:(int * __nonnull)nodes coord:(double * __nonnull)coord {
    
    int i, j;
    FileReader *reader;
    NSString *line;
    
    // Used to separate strings and filter them from white spaces
    NSCharacterSet *whitespaces = [NSCharacterSet whitespaceCharacterSet];
    NSPredicate *noEmptyStrings = [NSPredicate predicateWithFormat:@"SELF != ''"];

    line = nil;
    reader = (self.meshFileStreams)[ELEMENTS];
    if (step == self.elementCount) {
        
        [reader rewind];
        step = 0;
        return -1;

    } else if (step == 0) {
        [self SIOMeshAgent_cacheNodes];
    }
    
    line = [reader readLine];
    NSArray *stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
    NSArray *filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
    *tag = [filteredArray[0] intValue];
    *body = [filteredArray[1] intValue];
    *type = [filteredArray[2] intValue];

    int elNodes = [self SIOMeshAgent_elementNodes:*type];
    j = 3;
    for (i=0; i<elNodes; i++) {
        nodes[i] = [filteredArray[j] intValue];
        j++;
    }
    for (i=0; i<elNodes; i++) {
        if (![self SIOMeshAgent_copyCoords:coord+i*3 address:nodes[i]]) {
            NSLog(@"SIOMeshAgent:readNextElementCoordinates: an error occurred. Saino will abort the simulation now...");
            exit(14);
        }
    }
    
    step++;
    return 0;
}

-(int)readNextBoundaryElement:(int * __nonnull)tag part:(int * __nonnull)part boundary:(int * __nonnull)boundary leftElement:(int * __nonnull)leftElement rightElement:(int * __nonnull)rightElement type:(int * __nonnull)type nodes:(int * __nonnull)nodes coord:(double * __nonnull)coord {
    
    int i, j;
    FileReader *reader;
    NSString *line, *tagstr;
    NSRange subRange;
    
    // Used to separate strings and filter them from white spaces
    NSCharacterSet *whitespaces = [NSCharacterSet whitespaceCharacterSet];
    NSPredicate *noEmptyStrings = [NSPredicate predicateWithFormat:@"SELF != ''"];
    
    line = nil;
    reader = (self.meshFileStreams)[BOUNDARY];
    if (step == self.boundaryElementCount) {
        [reader rewind];
        step = 0;
        return -1;
    } else if (step == 0) {
        [self SIOMeshAgent_cacheNodes];
    }
    
    line = [reader readLine];
    NSArray *stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
    NSArray *filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
    tagstr = filteredArray[0];
    *boundary = [filteredArray[1] intValue];
    *leftElement = [filteredArray[2] intValue];
    *rightElement = [filteredArray[3] intValue];
    part = 0;
    subRange = [tagstr rangeOfString:@"/"];
    if (subRange.location != NSNotFound) {
        NSArray *components = [tagstr componentsSeparatedByString:@"/"];
        *tag = [components[0] intValue];
        *part = [components[1] intValue];        
    } else {
        *tag = [tagstr intValue];
    }
    
    *type = [filteredArray[4] intValue];
    int elNodes = [self SIOMeshAgent_elementNodes:*type];
    j = 5;
    for (i=0; i<elNodes; i++) {
        nodes[i] = [filteredArray[j] intValue];
        j++;
    }
    if (self.parallel) {
        int ok = 1;
        for (i=0; i<elNodes; i++) {
            if ([self SIOMeshAgent_searchNode:nodes[i]] == NULL) {
                ok = 0;
                break;
            }
        }
        if (!ok) {
            step++;
            return [self readNextBoundaryElement:tag part:part boundary:boundary leftElement:leftElement rightElement:rightElement type:type nodes:nodes coord:coord];
        }
    }
    
    for (i=0; i<elNodes; i++) {
        if (![self SIOMeshAgent_copyCoords:coord+i*3 address:nodes[i]]) {
            exit(14);
        }
    }
    
    step++;
    return 0;
}

-(int)readAllNodes:(int * __nonnull)tags coord:(double * __nonnull)coord {
    
    int i = 0;
    int pt = 0;
    
    [self SIOMeshAgent_cacheNodes];
    for (i=0; i<self.nodeCount; i++) {
        tags[i] = _clist[i].tag;
        coord[pt] = _clist[i].x;
        coord[pt+1] = _clist[i].y;
        coord[pt+2] = _clist[i].z;
        pt += 3;
    }
    
    return 0;
}

-(int)readSharedNode:(int * __nonnull)tag constraint:(int * __nonnull)constraint coord:(double * __nonnull)coord partCount:(int * __nonnull)partcount partitions:(int * __nonnull)partitions {
    
    int i, j;
    FileReader *reader;
    NSString *line;
    
    // Used to separate strings and filter them from white spaces
    NSCharacterSet *whitespaces = [NSCharacterSet whitespaceCharacterSet];
    NSPredicate *noEmptyStrings = [NSPredicate predicateWithFormat:@"SELF != ''"];
    
    line = nil;
    reader = (self.meshFileStreams)[SHARED];
    if (step == self.sharedNodeCount) {
        
        [reader rewind];
        step = 0;
        return -1;
    } else if (step == 0) {
        [self SIOMeshAgent_cacheNodes];
    }
    
    line = [reader readLine];
    NSArray *stringParts = [line componentsSeparatedByCharactersInSet:whitespaces];
    NSArray *filteredArray = [stringParts filteredArrayUsingPredicate:noEmptyStrings];
    *tag = [filteredArray[0] intValue];
    *partcount = [filteredArray[1] intValue];
    j = 2;
    for (i=0; i<*partcount; i++) {
        partitions[i] = [filteredArray[j] intValue];
        j++;
    }
    
    cacheNode *retval = [self SIOMeshAgent_searchNode:*tag];
    if (retval == NULL) {
        NSLog(@"SIOMeshAgent:readSharedNode: partition error: something is going totally wrong. Address: %d.\n", *tag);
        exit(23);
    } else {
        *constraint = retval->constraint;
        coord[0] = retval->x;
        coord[1] = retval->y;
        coord[2] = retval->z;
    }
    
    step++;
    return 0;
}

@end
