/* @(#)globals.h
 */

#ifndef _GLOBALS_H
#define _GLOBALS_H 1

// need this to work around bugi in gcc 4.7.2 when using stdint.h
#define __STDC_LIMIT_MACROS 1

// this is the number of sites that the HMM will run on
#define NUMSITES 1024
#define WN 16 // word num = NUMSITES / 64

#define WORDSHIFT 6
// input wordsize is 64 because we are working with uint64_t
#define WORDSIZE 64
#define WORDMOD 63

#define UINT32T_SIZE 32

// need to wait with defining norm until nvcc understands c++11
// #define NORM powf(FLT_MIN, 2.0f / 3.0f);

// definition of codebook specific globals
// this always needs to be a power of 2!!!!
#define BITSPERCODE 4

#endif /* _GLOBALS_H */

