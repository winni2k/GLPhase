#ifndef _HMMLIKE_HPP
#define _HMMLIKE_HPP 1

/*
  This is a CUDA based version of hmm_like in SNPTools.
  This library updates multiple samples in parallel,
  thereby breaking the Markov property in exchange for
  a speed boost through parallelization.
*/

#include <memory>
#include <vector>
#include <cstdint>
#include <cmath>
#include <cuda_runtime.h>
#include <sstream>
#include <gsl/gsl_rng.h>
#include <exception>
#include <utility>
#include <thrust/host_vector.h>
#include "sampler.hpp"
#include "utils.hpp"
#include "glPack.hpp"
#include "hmmLike.h" // for defines and extern functions

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");


static_assert(NUMSITES % WORDSIZE == 0,
              "Numsites is not evenly divisible by wordsize");

static_assert(
    sizeof(unsigned int) >= 4,
    "Size of unsigned int is < 4. We might run into indexing issues otherwise");

namespace HMMLikeHelper {
std::vector<unsigned> transposeHapIdxs(const std::vector<unsigned> &hapIdxs);
}

class HMMLike {
private:
  const std::vector<uint64_t> &m_inHapPanel;

  // this may be more than total number of samples in the case of a reference
  // panel
  const unsigned m_totalNumHaps;

  GLPack &m_glPack;
  const unsigned m_totalNumSamps;

  const unsigned m_numCycles;

  // transition matrix
  const std::vector<float> &m_tran;

  // mutation matrix.  pc[4][4] in impute.cpp
  const float (*m_mutationMat)[4][4];

  std::shared_ptr<Sampler> &m_sampler;
  unsigned m_nextSampIdx;
  gsl_rng &m_rng;

  /*
    Pulls out GLs for next run and repackages them ready to call cuda code
  */
  void CheckDevice() const { HMMLikeCUDA::CheckDevice(); }
  void CopyTranToDevice() const;
  void CopyMutationMatToDevice() const;

public:
  /*
    HMMLike constructor.  takes input to a 2*NxnumWordsPerHap matrix.
    It is assumed that the word length is 64 bit.
  */
  HMMLike(const std::vector<uint64_t> &hapPanel, unsigned numHaps,
          GLPack &glPack, unsigned numCycles, const std::vector<float> &tran,
          const float (*mutationMat)[4][4], std::shared_ptr<Sampler> &sampler,
          gsl_rng &rng);

  ~HMMLike() { HMMLikeCUDA::Cleanup(); }
  /*
    Replaces the sampler to be used to create prop haps
    Prop haps are produced when RunHMMOnSamples is called
   */
  void UpdateSampler(std::shared_ptr<Sampler> &sampler) {
    if (sampler)
      m_sampler = sampler;
    else
      throw std::runtime_error(
          "HMMLike::UpdateSampler was called with null pointer");
  };

  /*
    Fills sampldHaps with m_sampleStride sets of NUMHAPS haplotype numbers
    that have been sampled using numCycles HMM iterations.
    Also Fills hapIdxs with the m_sampleStride indices of the samples that the
    hapNumber sets were sampled for.
  */
  std::vector<unsigned> RunHMMOnSamples(unsigned &firstSampIdx,
                                        unsigned &lastSampIdx);
};

#endif
