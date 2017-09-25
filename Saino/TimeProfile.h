//===----------------------------------------------------------------------===//
//  TimeProfile.h
//  Saino
//
//  Created by Hakime Seddik on 08/02/12.
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
//===----------------------------------------------------------------------===//

#include <time.h>
#include <mach/mach_time.h>

double cputime(void);
double realtime(void);
double cpumemory(void);
double machcore(uint64_t endTime, uint64_t startTime);
