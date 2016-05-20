//
//  GPUUtils.c
//  Saino
//
//  Created by Hakime Seddik on 13/05/2016.
//  Copyright © 2016 Institute of Low Temperature Science. All rights reserved.
//

#include "GPUUtils.h"

void setPrecision(bool single) {
    
    if (single == true) {
        precisionMode = 1;
    }
}

int precision(void) {
    return precisionMode;
}

