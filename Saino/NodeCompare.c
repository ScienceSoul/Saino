//
//  File.c
//  Saino
//
//  Created by Hakime Seddik on 26/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#include <stdio.h>
#include "NodeCompare.h"

int __attribute__ ((cdecl)) nodecomp(const void * __nonnull a, const void * __nonnull b) {
    
    cacheNode *aptr = (cacheNode *)a;
    cacheNode *bptr = (cacheNode *)b;
    
    if (aptr->tag < bptr->tag) return -1;
    else if (aptr->tag > bptr->tag) return 1;
    return 0;
}