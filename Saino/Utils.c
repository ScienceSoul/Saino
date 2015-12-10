//
//  Utils.c
//  Saino
//
//  Created by Seddik hakime on 05/07/11.
//  Copyright 2011 Institute of Low Temperature Science. All rights reserved.
//

#include "Utils.h"
#include "TimeProfile.h"

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

char timeFormat[20];
static double advanceTime1;
static double advanceTime2;

int __attribute__((overloadable)) max(int x, int y) {
    
    return (x > y) ? x : y;
}

int __attribute__((overloadable)) max(int x, int y, int z) {
    
    int m = (x > y) ? x : y;
    return (z > m) ? z : m;
}

int __attribute__((overloadable)) max(int w, int x, int y, int z) {
    
    int m = (w > x) ? w : x;
    int n = (y > m) ? y : m;
    return (z > n) ? z : n;
}

double __attribute__((overloadable)) max(double x, double y) {
    
    return (x > y) ? x : y;
}

double __attribute__((overloadable)) max(double x, double y, double z) {
    
    double m = (x > y) ? x : y;
    return (z > m) ? z : m;
}

int __attribute__((overloadable)) min(int x, int y) {
    
    return (x < y) ? x : y;
}

int __attribute__((overloadable)) min(int x, int y, int z) {
    
    int m = (x < y) ? x : y;
    return (z < m) ? z : m;
}

int __attribute__((overloadable)) min(int w, int x, int y, int z) {
 
    int m = (w < x) ? w : x;
    int n = (y < m) ? y : m;
    return (z < n) ? z : n;
}

double __attribute__((overloadable)) min(double x, double y) {
    
    return (x < y) ? x : y;
}

double __attribute__((overloadable)) min(double x, double y, double z) {
    
    double m = (x < y) ? x : y;
    return (z < m) ? z : m;
}

int __attribute__((overloadable)) max_array(int *a, int num_elements)
{
    int i, max = -INT_MAX;
    for (i=0; i<num_elements; i++) {
        if (a[i]>max) {
            max = a[i];
        }
    }
    return max;
}

double __attribute__((overloadable)) max_array(double *a, int num_elements) {
    
    int i;
    double max = -HUGE_VAL;
    for (i=0; i<num_elements; i++) {
        if (a[i]>max) {
            max = a[i];
        }
    }
    return max;
    
}

int __attribute__((overloadable)) min_array(int *a, int num_elements) {
    
    int i, min = INT_MAX;
    for (i=0; i<num_elements; i++) {
        if (a[i]<min) {
            min = a[i];
        }
    }
    return min;
    
}

double __attribute__((overloadable)) min_array(double *a, int num_elements) {
    
    int i;
    double min = HUGE_VAL;
    for (i=0; i<num_elements; i++) {
        if (a[i]<min) {
            min = a[i];
        }
    }
    return min;
}

bool __attribute__((overloadable)) all(int *v, char mask, int val, int range) {
/*************************************************************************************************
    Test if all values in v satisfies the mask for value val for the range of data range
    Work on an array of int values
*************************************************************************************************/

    int i;
    bool rlt = true;
    
    switch (mask) {
        case '=':
            for (i=0; i<range; i++) {
                if (v[i] == val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
        case '>':
            for (i=0; i<range; i++) {
                if (v[i] > val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
        case '<':
            for (i=0; i<range; i++) {
                if (v[i] < val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
    }
    
    return rlt;
}

bool __attribute__((overloadable)) all(float *v, char mask, float val, int range) {
/*************************************************************************************************
    Test if all values in v satisfies the mask for value val for the range of data range
    Work on an array of float values
*************************************************************************************************/
    
    int i;
    bool rlt = true;
    
    switch (mask) {
        case '=':
            for (i=0; i<range; i++) {
                if (v[i] == val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
        case '>':
            for (i=0; i<range; i++) {
                if (v[i] > val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
        case '<':
            for (i=0; i<range; i++) {
                if (v[i] < val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
    }
    
    return rlt;
}

bool __attribute__((overloadable)) all(double *v, char mask, double val, int range) {
/*************************************************************************************************
    Test if all values in v satisfies the mask for value val for the range of data range
    Work on an array of float values
*************************************************************************************************/
    
    int i;
    bool rlt = true;
    
    switch (mask) {
        case '=':
            for (i=0; i<range; i++) {
                if (v[i] == val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
        case '>':
            for (i=0; i<range; i++) {
                if (v[i] > val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
        case '<':
            for (i=0; i<range; i++) {
                if (v[i] < val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
    }
    
    return rlt;
}

bool __attribute__((overloadable)) all_in_range(int *v, char mask, int val, int start, int range) {
/*************************************************************************************************
    Test if all values in v satisfies the mask for value val for the range of data range 
    starting at start
*************************************************************************************************/
    
    int i;
    bool rlt = true;
    
    switch (mask) {
        case '=':
            for (i=start; i<range; i++) {
                if (v[i] == val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
        case '>':
            for (i=start; i<range; i++) {
                if (v[i] > val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
        case '<':
            for (i=start; i<range; i++) {
                if (v[i] < val) { 
                    continue;
                } else {
                    rlt = false;
                    return rlt;
                }
            }
            break;
    }
    
    return rlt;
    
}


bool __attribute__((overloadable)) any(int *v, char mask, int val, int range) {
    
    int i;
    bool rlt = false;
    
    switch (mask) {
        case '=':
            for (i=0; i<range; i++) {
                if (v[i] == val) {
                    rlt = true;
                    return rlt;
                } else {
                    continue;
                }
            }
            break;
        case '>':
            for (i=0; i<range; i++) {
                if (v[i] > val) { 
                    rlt = true;
                    return rlt;
                } else {
                     continue;
                }
            }
            break;
        case '<':
            for (i=0; i<range; i++) {
                if (v[i] < val) { 
                    rlt = true;
                    return rlt;
                } else {
                    continue;
                }
            }
            break;
    }
    
    return rlt;
    
}

bool __attribute__((overloadable)) any(float *v, char mask, float val, int range) {
    
    int i;
    bool rlt = false;
    
    switch (mask) {
        case '=':
            for (i=0; i<range; i++) {
                if (v[i] == val) {
                    rlt = true;
                    return rlt;
                } else {
                    continue;
                }
            }
            break;
        case '>':
            for (i=0; i<range; i++) {
                if (v[i] > val) { 
                    rlt = true;
                    return rlt;
                } else {
                    continue;
                }
            }
            break;
        case '<':
            for (i=0; i<range; i++) {
                if (v[i] < val) { 
                    rlt = true;
                    return rlt;
                } else {
                    continue;
                }
            }
            break;
    }
    
    return rlt;
    
}

bool __attribute__((overloadable)) any(double *v, char mask, double val, int range) {
    
    int i;
    bool rlt = false;
    
    switch (mask) {
        case '=':
            for (i=0; i<range; i++) {
                if (v[i] == val) {
                    rlt = true;
                    return rlt;
                } else {
                    continue;
                }
            }
            break;
        case '>':
            for (i=0; i<range; i++) {
                if (v[i] > val) { 
                    rlt = true;
                    return rlt;
                } else {
                    continue;
                }
            }
            break;
        case '<':
            for (i=0; i<range; i++) {
                if (v[i] < val) { 
                    rlt = true;
                    return rlt;
                } else {
                    continue;
                }
            }
            break;
    }
    
    return rlt;
    
}

int __attribute__((overloadable)) count(int *v, char mask, int val, int range) {
    
    int i;
    int rlt = 0;
    
    switch (mask) {
        case '=':
            for (i=0; i<range; i++) {
                if (v[i] == val) {
                    rlt++;
                } else {
                    continue;
                }
            }
            break;
        case '>':
            for (i=0; i<range; i++) {
                if (v[i] > val) { 
                    rlt++;
                } else {
                    continue;
                }
            }
            break;
        case '<':
            for (i=0; i<range; i++) {
                if (v[i] < val) { 
                    rlt++;
                } else {
                    continue;
                }
            }
            break;
    }
    
    return rlt;
    
}

int __attribute__((overloadable)) count(float *v, char mask, float val, int range) {
    
    int i;
    int rlt = 0;
    
    switch (mask) {
        case '=':
            for (i=0; i<range; i++) {
                if (v[i] == val) {
                    rlt++;
                } else {
                    continue;
                }
            }
            break;
        case '>':
            for (i=0; i<range; i++) {
                if (v[i] > val) { 
                    rlt++;
                } else {
                    continue;
                }
            }
            break;
        case '<':
            for (i=0; i<range; i++) {
                if (v[i] < val) { 
                    rlt++;
                } else {
                    continue;
                }
            }
            break;
    }
    
    return rlt;
    
}

int __attribute__((overloadable)) count(double *v, char mask, double val, int range) {
    
    int i;
    int rlt = 0;
    
    switch (mask) {
        case '=':
            for (i=0; i<range; i++) {
                if (v[i] == val) {
                    rlt++;
                } else {
                    continue;
                }
            }
            break;
        case '>':
            for (i=0; i<range; i++) {
                if (v[i] > val) { 
                    rlt++;
                } else {
                    continue;
                }
            }
            break;
        case '<':
            for (i=0; i<range; i++) {
                if (v[i] < val) { 
                    rlt++;
                } else {
                    continue;
                }
            }
            break;
    }
    
    return rlt;
    
}

void __attribute__((overloadable)) sort(unsigned long n, int *arr) {
/*******************************************************************************************
 Sorts an array arr[1...n] into ascending order using Quicksort algorithm. n is imput;
 arr is replaced on output by its sorted rearrangment
*******************************************************************************************/
    
    unsigned long i, ir=n, j, k, l=1, *istack;
    int jstack=0;
    int a, temp;
    
    istack = ulongvec(1, NSTACK);
    
    for (;;) {
        if (ir-l < M) {
            for (j=l+1; j<=ir; j++) {
                a = arr[j];
                for (i=j-1; i>=l; i--) {
                    if (arr[i] <= a) break;
                    arr[i+1] = arr[i];
                }
                arr[i+1] = a;
            }
            if (jstack == 0) break;
            ir = istack[jstack--];
            l = istack[jstack--];
        } else {
            k = (l+ir) >> 1;
            SWAP(arr[k], arr[l+1]);
            if (arr[l] > arr[ir]) {
                SWAP(arr[l], arr[ir]);
            }
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1], arr[ir]);
            }
            if (arr[l] > arr[l+1]) {
                SWAP(arr[l], arr[l+1]);
            }
            i = l+1;
            j = ir;
            a = arr[l+1];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i], arr[j]);
            }
            arr[l+1] = arr[j];
            arr[j] = a;
            jstack += 2;
            if (jstack > NSTACK) fatal("Sort", "NSTACK too small in sort routine.");
            if (ir-i+1 >= j-l) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack] = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
    
    free_ulvector(istack, 1, NSTACK);
}

void __attribute__((overloadable)) sort(unsigned long n, int *arr, int *brr) {
/*******************************************************************************************
 Sorts an array arr[1...n] into ascending order using Quicksort, while making the 
 corresponding rearrangment of array brr[1...n]
*******************************************************************************************/
    
    unsigned long i, ir=n, j, k, l=1, *istack;
    int jstack=0;
    int a, b, temp;
    
    istack = ulongvec(1, NSTACK);
    
    for (;;) {
        if (ir-l < M) {
            for (j=l+1; j<=ir; j++) {
                a = arr[j];
                b = brr[j];
                for (i=j-1; i>=l; i--) {
                    if (arr[i] <= a) break;
                    arr[i+1] = arr[i];
                    brr[i+1] = brr[i];
                }
                arr[i+1] = a;
                brr[i+1] = b;
            }
            if (!jstack) {
                free_ulvector(istack, 1, NSTACK);
                return;
            }
            ir = istack[jstack];
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir) >> 1;
            SWAP(arr[k], arr[l+1]);
            SWAP(brr[k], brr[l+1]);
            if (arr[l] > arr[ir]) {
                SWAP(arr[l], arr[ir]);
                SWAP(brr[l], brr[ir]);
            }
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1], arr[ir]);
                SWAP(brr[l+1], brr[ir]);
            }
            if (arr[l] > arr[l+1]) {
                SWAP(arr[l], arr[l+1]);
                SWAP(brr[l], brr[l+1]);
            }
            i = l+1;
            j = ir;
            a = arr[l+1];
            b = brr[l+1];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i], arr[j]);
                SWAP(brr[i], brr[j]);
            }
            arr[l+1] = arr[j];
            arr[j] = a;
            brr[l+1] = brr[j];
            brr[j] = b;
            jstack += 2;
            if (jstack > NSTACK) fatal("Sort", "NSTACK too small in sort routine.");
            if (ir-i+1 >= j-l) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack] = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}

void __attribute__((overloadable)) sort(unsigned long n, int *arr, float *brr) {
    
    unsigned long i, ir=n, j, k, l=1, *istack;
    int jstack=0;
    int a, temp;
    float b;
    
    istack = ulongvec(1, NSTACK);
    
    for (;;) {
        if (ir-l < M) {
            for (j=l+1; j<=ir; j++) {
                a = arr[j];
                b = brr[j];
                for (i=j-1; i>=l; i--) {
                    if (arr[i] <= a) break;
                    arr[i+1] = arr[i];
                    brr[i+1] = brr[i];
                }
                arr[i+1] = a;
                brr[i+1] = b;
            }
            if (!jstack) {
                free_ulvector(istack, 1, NSTACK);
                return;
            }
            ir = istack[jstack];
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir) >> 1;
            SWAP(arr[k], arr[l+1]);
            SWAP(brr[k], brr[l+1]);
            if (arr[l] > arr[ir]) {
                SWAP(arr[l], arr[ir]);
                SWAP(brr[l], brr[ir]);
            }
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1], arr[ir]);
                SWAP(brr[l+1], brr[ir]);
            }
            if (arr[l] > arr[l+1]) {
                SWAP(arr[l], arr[l+1]);
                SWAP(brr[l], brr[l+1]);
            }
            i = l+1;
            j = ir;
            a = arr[l+1];
            b = brr[l+1];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i], arr[j]);
                SWAP(brr[i], brr[j]);
            }
            arr[l+1] = arr[j];
            arr[j] = a;
            brr[l+1] = brr[j];
            brr[j] = b;
            jstack += 2;
            if (jstack > NSTACK) fatal("Sort", "NSTACK too small in sort routine.");
            if (ir-i+1 >= j-l) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack] = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}

void __attribute__((overloadable)) sort(unsigned long n, int *arr, double *brr) {
    
    unsigned long i, ir=n, j, k, l=1, *istack;
    int jstack=0;
    int a, temp;
    double b;
    
    istack = ulongvec(1, NSTACK);
    
    for (;;) {
        if (ir-l < M) {
            for (j=l+1; j<=ir; j++) {
                a = arr[j];
                b = brr[j];
                for (i=j-1; i>=l; i--) {
                    if (arr[i] <= a) break;
                    arr[i+1] = arr[i];
                    brr[i+1] = brr[i];
                }
                arr[i+1] = a;
                brr[i+1] = b;
            }
            if (!jstack) {
                free_ulvector(istack, 1, NSTACK);
                return;
            }
            ir = istack[jstack];
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir) >> 1;
            SWAP(arr[k], arr[l+1]);
            SWAP(brr[k], brr[l+1]);
            if (arr[l] > arr[ir]) {
                SWAP(arr[l], arr[ir]);
                SWAP(brr[l], brr[ir]);
            }
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1], arr[ir]);
                SWAP(brr[l+1], brr[ir]);
            }
            if (arr[l] > arr[l+1]) {
                SWAP(arr[l], arr[l+1]);
                SWAP(brr[l], brr[l+1]);
            }
            i = l+1;
            j = ir;
            a = arr[l+1];
            b = brr[l+1];
            for (;;) {
                do i++; while (arr[i] < a);
                do j--; while (arr[j] > a);
                if (j < i) break;
                SWAP(arr[i], arr[j]);
                SWAP(brr[i], brr[j]);
            }
            arr[l+1] = arr[j];
            arr[j] = a;
            brr[l+1] = brr[j];
            brr[j] = b;
            jstack += 2;
            if (jstack > NSTACK) fatal("Sort", "NSTACK too small in sort routine.");
            if (ir-i+1 >= j-l) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack] = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }    
}

void __attribute((overloadable)) reverse(int *arr, size_t narr) {
    
    size_t i;
    
    for (i=0; i < narr / 2; ++i) {
        int tmp = arr[i];
        arr[i] = arr[narr-i-1];
        arr[narr-i-1] = tmp;
    }
}

/*******************************************************************************************
 Circular shift an array. Only works for shift > 0
*******************************************************************************************/
void __attribute((overloadable)) cshift(int *arr, size_t narr, unsigned long shift) {
    reverse(arr, shift);
    reverse(arr + shift, narr - shift);
    reverse(arr, narr);
}

/*******************************************************************************************
 Retunrs current data and time
*******************************************************************************************/
char *dateAndTime(void){
    time_t currtime;
    time(&currtime);
    strftime(timeFormat,sizeof(timeFormat),"%Y/%m/%d %H:%M:%S",localtime(&currtime));
    return (char *)timeFormat;
}

void startAdvanceOutput(char *solverName, char *outputType) {
    
    advanceTime1 = realtime();
    advanceTime2 = realtime();
    printf("%s: %s\n", solverName, outputType);
}

void advanceOutput(int t, int n, double *dot_t, double *percent_t) {
    
    int i;
    double d_t, p_t;
    static bool newLine = true;
    
    d_t = 1.0;
    p_t = 20.0;
    if (dot_t != NULL) d_t = *dot_t;
    if (percent_t != NULL) p_t = *percent_t;
    
    if (realtime() - advanceTime1 > d_t) {
        if (newLine == true) {
            printf(": ");
            newLine = false;
        }
        printf(".");
        
        if (realtime() -  advanceTime2 > p_t) {
            i = (int)round(t*100.0/n);
            printf("%d%%\n", i);
            newLine = true;
            advanceTime2 = realtime();
        }
        advanceTime1 = realtime();
    }
    if (t == n-1) {
        printf("\n");
        newLine = true;
    }
}