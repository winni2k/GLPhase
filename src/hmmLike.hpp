#ifndef _HMMLIKE_H
#define _HMMLIKE_H 1

/*
  This is a CUDA based version of hmm_like in SNPTools.
  This library updates multiple samples in parallel,
  thereby breaking the Markov property in exchange for
  a speed boost through parallelization.
*/

#include <vector>
#include <cstdint>

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

class HMMLike {
private:
  unsigned m_numSites;
  unsigned m_totalNumHaps; // this may be more than total number of samples
  unsigned m_totalNumSamps;

public:
  /*
    HMMLike constructor.  takes input to a 2*NxnumWordsPerHap matrix.
    It is assumed that the word length is 64 bit.
  */
  HMMLike(const uint64_t &hapPanel, unsigned numSites, unsigned numHaps,
          unsigned numSamps, const std::vector<float> &GLs,
          unsigned sampleStride);

  /*
    Returns an NxS matrix of S haplotype numbers
    that have been sampled using C HMM iterations.
  */
  std::vector<std::vector<unsigned> > SampleHaps();
}

#endif
