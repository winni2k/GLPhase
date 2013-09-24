/* @(#)kMeans.h
 */

#ifndef _KMEANS_H
#define _KMEANS_H 1

#include        <gsl/gsl_rng.h>
#include        <vector>
#include        <cassert>
#include <iostream>
#include        "haplotype.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

class KMeans {

private:
    bool m_bInitialized = false;

    unsigned m_uNumWordsPerHap = 0;
    unsigned m_uNumSites = 0;
    unsigned m_uNumHaps = 0;

    double 
public:
    
}

#endif /* _KMEANS_H */

