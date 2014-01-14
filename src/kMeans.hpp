/* @(#)kMeans.hpp
 */

<<<<<<< HEAD
#ifndef _KNN_H
#define _KNN_H 1
=======
#ifndef _KMEANS_H
#define _KMEANS_H 1
>>>>>>> refs/heads/synced/master

#include        <gsl/gsl_rng.h>
#include        <vector>
#include        <cassert>
#include        <iostream>
#include        <math.h>
#include        <algorithm>
#include        "haplotype.hpp"


//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

<<<<<<< HEAD
class KNN {
=======
class KMeans {
>>>>>>> refs/heads/synced/master

  private:
    bool m_bInitialized = false;
    unsigned m_uNumWordsPerHap = 0;
    unsigned m_uNumSites = 0;
    unsigned m_uNumHaps = 0;
    unsigned m_uNumClusters = 0;

    // scaffold data
    double m_dFreqCutoff = -2;
    std::vector<double> m_vdVarAfs;
    std::vector<unsigned> m_vuCommonSiteNums;
    std::vector < std::vector <unsigned> > m_vvuNeighbors; // sample then hap

    void CalculateVarAfs(const std::vector < uint64_t > & vuScaffoldHaps);
    void CutDownHaps();
    void AssignHap(Haplotype &hHap, const uint64_t *oa);
    bool UsingScaffold() {
        return (m_dFreqCutoff >= 0 && m_dFreqCutoff <= 1);
    }
    
  public:

<<<<<<< HEAD
    KNN(unsigned uNumClust)
=======
    KMeans(unsigned uNumClust)
>>>>>>> refs/heads/synced/master
        : m_uNumClusters(uNumClust) {};

    void init(const std::vector< uint64_t > & vuHaplotypes,
              unsigned uNumWordsPerHap, unsigned uNumSites,
              double dFreqCutoff);

    // returns a haplotype sampled using the relationship graph
    unsigned SampleHap(unsigned uInd, gsl_rng *rng);

    // fills returnHapNums with haplotypes closest to the sample that uHapNum is a part of
    void Neighbors(unsigned uHapNum, std::vector< unsigned > & returnHapNums){
        assert(m_bInitialized == true);
        for(auto iter : m_vvuNeighbors[floor(uHapNum/2)])
            returnHapNums.push_back(iter);
    };

    // fills returnSiteNums with sites that are used for clustering
    void ClusterSites(std::vector<unsigned> & returnSiteNums){
        assert(m_bInitialized);
        returnSiteNums = m_vuCommonSiteNums;
    }

    // return the calculated variant allele frequencies
    void VarAfs(std::vector<double> & returnVarAfs){
        assert(m_bInitialized);
        returnVarAfs = m_vdVarAfs;
    }
    
    void ClusterHaps(const std::vector< uint64_t > & vuHaplotypes);


};

#endif /* _KMEANS_H */


