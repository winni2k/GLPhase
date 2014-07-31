/* @(#)kNN.hpp
 */

#ifndef _KNN_H
#define _KNN_H 1

#include <gsl/gsl_rng.h>
#include <vector>
#include <cassert>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include "haplotype.hpp"
#include "sampler.hpp"
#include <omp.h>
#include <string>

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

enum class kNNDistT { hamming, tracLen };

namespace KNNhelper {

// returns the neighbors of haplotype hapNum
struct dist {
  unsigned idx;
  unsigned distance;
};

struct by_dist {
  bool operator()(dist const &a, dist const &b) {
    return a.distance < b.distance;
  }
};
}

class KNN : public Sampler {

private:
  unsigned m_numWordsPerHap = 0;
  unsigned m_numSites = 0;
  unsigned m_numClusters = 0;
  kNNDistT m_distMetric;

  // scaffold data
  const double m_freqLB; // lower bound
  const double m_freqUB; // upper bound
  const bool m_usingMAF;
  std::vector<double> m_varAFs;
  std::vector<unsigned> m_vuCommonSiteNums;
  const bool m_usingScaffold = false;

  // stores the list of haplotypes that are closest to each sample
  // sample then hap
  std::vector<std::vector<unsigned> > m_neighbors;
  //  std::vector<std::vector<KNNhelper::dist> > m_distsHapA;
  //  std::vector<std::vector<KNNhelper::dist> > m_distsHapB;

  void CalculateVarAfs(const std::vector<uint64_t> &vuScaffoldHaps);
  void CutDownHaps();
  void AssignHap(Haplotype &hHap, const uint64_t *oa);
  bool UsingScaffold() { return m_usingScaffold; }

public:
  KNN(unsigned numClust, const std::vector<uint64_t> &vuHaplotypes,
      unsigned numWordsPerHap, unsigned numSites, double freqLB, double freqUB,
      bool usingMAF, gsl_rng *rng, kNNDistT distMetric = kNNDistT::hamming);

  // returns a haplotype sampled uniformly from k-nearest neighbors
  unsigned SampleHap(unsigned uInd) override;

  // update proposal distribution based on the result of an MCMC proposal
  // (input)
  void UpdatePropDistProp(const std::vector<unsigned> &, unsigned, bool,
                          float) override{};

  // update proposal distribution based on the input haplotype set
  void UpdatePropDistHaps(const std::vector<uint64_t> &) override{};

  // fills returnHapNums with haplotypes closest to the sample that uHapNum is a
  // part of
  std::vector<unsigned> Neighbors(unsigned uHapNum);

  // fills returnSiteNums with sites that are used for clustering
  void ClusterSites(std::vector<unsigned> &returnSiteNums) {
    returnSiteNums = m_vuCommonSiteNums;
  }

  // return the calculated variant allele frequencies
  void VarAfs(std::vector<double> &returnVarAfs) { returnVarAfs = m_varAFs; }

  void ClusterHaps(const std::vector<uint64_t> &vuHaplotypes);

  void Save(const std::string &fileName,
            const std::vector<std::string> &sampNames) override;
};

#endif /* _KNN_H */
