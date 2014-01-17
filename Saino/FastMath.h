//
//  FastMath.h
//  Saino
//
//  Created by Seddik hakime on 15/01/2014.
//  Copyright (c) 2014 Institute of Low Temperature Science. All rights reserved.
//

inline double fastPow(double a, double b) {
    
    union {
        double d;
        int x[2];
    } u = { a };
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    
    return u.d;
}

// should be much more precise with large b
inline double fastPrecisePow(double a, double b) {
    
    // calculate approximation with fraction of the exponent
    int e = (int) b;
    union {
        double d;
        int x[2];
    } u = { a };
    u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
    u.x[0] = 0;
    
    // exponentiation by squaring with the exponent's integer part
    // double r = u.d makes everything much slower, not sure why
    double r = 1.0;
    while (e) {
        if (e & 1) {
            r *= a;
        }
        a *= a;
        e >>= 1;
    }
    
    return r * u.d;
}