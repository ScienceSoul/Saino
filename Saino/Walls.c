//===----------------------------------------------------------------------===//
//  Walls.c
//  Saino
//
//  Created by Seddik hakime on 23/05/2014.
//  Copyright (c) 2014 ScienceSoul. All rights reserved.
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
//
//===----------------------------------------------------------------------===//
//  Contains functions for the wall boundary conditions of
//  the k-epsilon turbulence model on walls.
//===----------------------------------------------------------------------===//

#include "Walls.h"

extern inline double wall_law(double ufric, double ut, double densit, double viscos, double dist, double rough);
extern inline double d_wall_law(double ufric, double ut, double densit, double viscos, double dist, double rough);
extern inline void solve_ufric(double densit, double viscos, double dist, double rough, double ut, double * _Nonnull ufric, double * _Nonnull dfx);
extern inline void kewall(double * _Nonnull tk, double * _Nonnull teps, double * _Nonnull tomg, double ut, double dist, double rough, double viscos, double densit);
