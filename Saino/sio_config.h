//
//  sio_config.h
//  Saino
//
//  Created by Hakime Seddik on 21/06/12.
//  Copyright (c) 2012 Institute of Low Temperature Science. All rights reserved.
//

#ifndef Saino_sio_config_h
#define Saino_sio_config_h

#define PATH_MAX 1024

typedef struct {
    int tag;
    int constraint;
    double x, y, z;
} cacheNode;

#endif
