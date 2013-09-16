/* @(#)kMedoids.h
 */

#ifndef _KMEDOIDS_H
#define _KMEDOIDS_H 1

#include        <gsl/gsl_rng.h>
#include        <vector>
#include        <cassert>
#include <iostream>
#include        "haplotype.h"
#include <float.h>


//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

class KMedoids {

private:

    bool m_bInitialized = false;
    
    // need these for clustering
    std::vector< unsigned > m_vuHapMedNum; // keeps track of which medioid each haplotype is closest to

    // keeps track of which haplotype each medoid is
    // m_vuMedoidHapNum.size() is number of clusters
    std::vector< unsigned > m_vuMedoidHapNum; 
    std::vector< unsigned > m_vuHapHammingDist; // keeps track of hamming distance between each hap and its medoid
//    std::vector< Haplotype > m_uClusterHaps;

    // vars for clustering
    unsigned m_uNumWordsPerHap = 0;
    unsigned m_uNumSites = 0;
    unsigned m_uNumHaps = 0;

    // 0 = simple (slow) kMedoids
    // 1 = clustering according to Park and Jun 2008, 10% subsample
    unsigned m_uClusterType;
    double m_dDelta = DBL_MIN * 10;
    
    double MedoidLoss(const std::vector< uint64_t > * pvuHaplotypes, double dPower ) {
        std::vector< unsigned > vuMedHaps;
        return MedoidLoss(pvuHaplotypes, vuMedHaps, dPower);
    };
    double MedoidLoss(const std::vector< uint64_t > * pvuHaplotypes ) {
        std::vector< unsigned > vuMedHaps;
        return MedoidLoss(pvuHaplotypes, vuMedHaps);
    };
    double MedoidLoss(const std::vector< uint64_t > * pvuHaplotypes, const std::vector< unsigned > & vuMedHaps )
        { return MedoidLoss(pvuHaplotypes,vuMedHaps, 2.0); };
    double MedoidLoss(const std::vector< uint64_t > * pvuHaplotypes, const std::vector< unsigned > & vuMedHaps, double dPower );
    double UpdateMedoidPAM(const std::vector< uint64_t > * pvuHaplotypes, double dBestLoss, unsigned uMedNum );
    double UpdateMedoidParkJun(const std::vector< uint64_t > * pvuHaplotypes, unsigned uMedNum );
    void InputTesting(const std::vector< uint64_t > * pvuHaplotypes);
    void AssignHapsToBestMedoids(const std::vector< uint64_t > * pvuHaplotypes );
    
public:

    KMedoids(unsigned uClusterType, unsigned uNumClust, const std::vector< uint64_t > * pvuHaplotypes,
             unsigned uNumWordsPerHap, unsigned uNumSites, gsl_rng *rng)
        : m_uClusterType( uClusterType ) {
        m_vuMedoidHapNum.resize(uNumClust);
        init(pvuHaplotypes, uNumWordsPerHap, uNumSites, rng);
    };

    KMedoids(unsigned uClusterType, unsigned uNumClust)
        : m_uClusterType( uClusterType ) {
        m_vuMedoidHapNum.resize(uNumClust);
    };

    void init(const std::vector< uint64_t > * pvuHaplotypes, unsigned uNumWordsPerHap, unsigned uNumSites, gsl_rng *rng);
    
    // returns a haplotype sampled using the relationship graph
    unsigned SampleHap(unsigned uInd, gsl_rng *rng);

    // update the medoids
    void UpdateMedoids(const std::vector< uint64_t > * pvuHaplotypes );

};

#endif /* _KMEDOIDS_H */

