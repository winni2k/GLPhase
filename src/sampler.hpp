#ifndef _SAMPLER_H
#define _SAMPLER_H 1

#include <gsl/gsl_rng.h>
#include <vector>
#include <cstdint>
#include <string>
#include "utils.hpp"

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

class Sampler {

protected:
  gsl_rng *m_rng;
  unsigned m_numSamples;
  unsigned m_numHaps;

public:
  virtual ~Sampler() {};
  Sampler(gsl_rng *rng, unsigned numSamples, unsigned numHaps);
    
  // returns a haplotype sampled using the sampler, but only from the
  // reference haplotypes if onlyFromRef = true
  unsigned SampleHap(unsigned excludeInd, bool onlyFromRef);

  // returns a haplotype sampled using the sampler
  virtual unsigned SampleHap(unsigned excludeInd) = 0;

  // update proposal distribution based on the result of an MCMC proposal
  // (input)
  virtual void UpdatePropDistProp(const std::vector<unsigned> &propHaps,
                                  unsigned updateIndNum, bool accepted,
                                  float penRatio = 1) = 0;

  // update proposal distribution based on the input haplotype set
  virtual void UpdatePropDistHaps(const std::vector<uint64_t> &haplotypes) = 0;

  virtual void Save(std::string outputFile,
                    std::vector<std::string> &sampNames) {
    throw myException("Attempted to save sampler state of " +
                      sutils::uint2str(sampNames.size()) +
                      " samples to file: " + outputFile +
                      "\nHowever, no state worth saving was encountered");
  };
};

// Derived class for uniform sampling
class UnifSampler : public Sampler {

public:
  UnifSampler(gsl_rng *rng, unsigned numSamples, unsigned numHaps)
      : Sampler(rng, numSamples, numHaps) {};

  // perform uniform sampling to get hap
  unsigned SampleHap(unsigned excludeInd);

  // don't do any updating of the proposal distribution
  void UpdatePropDistProp(const std::vector<unsigned> &, unsigned, bool,
                          float) {};

  // don't do any updating of the proposal distribution
  void UpdatePropDistHaps(const std::vector<uint64_t> &) {};
};

/*
  This is a sampler based on a relationship graph.
  The graph is updated based on the success of proposal haplotype sets.
  There are two possible graphs corresponding to row and columns:
  1. Sample/Sample
  2. Sample/Hap

*/

enum class GraphSampT { sampSamp, sampHap };
class GraphSampler : public Sampler {

private:
  // 1 based number of rows and columns
  const unsigned m_rows;
  unsigned m_cols = 0;

  /*
    numerator and denominator of relationship matrix
    numerator = number of accepted proposals
    denominator = number of total proposals
    first index is individual for which proposal was made
    second index is individual from which was copied
  */
  std::vector<std::vector<float> > m_2dRelationshipMatNum;
  std::vector<std::vector<float> > m_2dRelationshipMatDen;

  GraphSampT m_graphType;
  bool m_usingHaps = false;
  unsigned m_colHapFactor = 0;

  // hap to column index converter
  unsigned Hap2Col(unsigned uHap) {
    return uHap / m_colHapFactor;
  };

  // column to hap index converter
  unsigned Col2Hap(unsigned uCol) {
    return uCol * m_colHapFactor;
  };

public:
  GraphSampler(gsl_rng *rng, unsigned numSamples, unsigned numHaps,
               GraphSampT graphType);

  // returns a haplotype sampled using the sampler
  unsigned SampleHap(unsigned excludeInd);

  /*
    update proposal distribution based on the result of an MCMC (input)

    Takes sample num and haplotype num as well as graph type.
    Samples can be updated, while haplotypes can only be copied from.
    Every sample has two haplotypes.
    Therefore: any haplotypes that don't match a sample are reference haps
  */

  void UpdatePropDistProp(const std::vector<unsigned> &propHaps,
                          unsigned updateIndNum, bool accepted,
                          float penRatio = 1);

  // update proposal distribution based on the input haplotype set
  void UpdatePropDistHaps(const std::vector<uint64_t> &) {
    return;
  };

  void Save(std::string fileName, const std::vector<std::string> &sampNames);
};

#endif
