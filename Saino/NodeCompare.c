//===----------------------------------------------------------------------===//
//  File.c
//  Saino
//
//  Created by Hakime Seddik on 26/06/12.
//  Copyright (c) 2012 ScienceSoul. All rights reserved.
//  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland.
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
//===----------------------------------------------------------------------===//

#include <stdio.h>
#include "NodeCompare.h"

int __attribute__ ((cdecl)) nodecomp(const void * _Nonnull a, const void * _Nonnull b) {
    
    cacheNode *aptr = (cacheNode *)a;
    cacheNode *bptr = (cacheNode *)b;
    
    if (aptr->tag < bptr->tag) return -1;
    else if (aptr->tag > bptr->tag) return 1;
    return 0;
}
