#ifndef _HMMLIKE_HPP
#define _HMMLIKE_HPP 1

/*
  This is a CUDA based version of hmm_like in SNPTools.
  This library updates multiple samples in parallel,
  thereby breaking the Markov property in exchange for
  a speed boost through parallelization.
*/

#include <vector>
#include <cstdint>
#include <cmath>
#include <cuda_runtime.h>
#include <sstream>
#include <gsl/gsl_rng.h>
#include "sampler.hpp"
#include "utils.hpp"
#include "glPack.hpp"

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

// this is the number of sites that the HMM will run on
#define NUMSITES 512
// input wordsize is 64 because we are working with uint64_t
#define WORDSIZE 64
static_assert(NUMSITES % WORDSIZE == 0,
              "Numsites is not evenly divisible by wordsize");

// this is the number of haplotypes we use in the HMM
// That is, the HMM has NUMSITES * NUMHAPS states
#define NUMHAPS 4
static_assert(
    sizeof(unsigned int) >= 4,
    "Size of unsigned int is < 4. We might run into indexing issues otherwise");

// these functions are implemented in HMMLike.cu
namespace HMMLikeCUDA {
extern "C" void CheckDevice();
extern "C" cudaError_t CopyTranToDevice(const std::vector<float> &tran);
extern "C" cudaError_t CopyMutationMatToDevice(const float (*mutMat)[4][4]);
extern "C" void RunHMMOnDevice(const std::vector<char> &packedGLs,
                               const std::vector<uint64_t> &hapPanel,
                               const std::vector<unsigned> &extraPropHaps,
                               unsigned numSites, unsigned numSamples,
                               unsigned numCycles,
                               std::vector<unsigned> &hapIdxs,
                               unsigned long seed);
}

class HMMLike {
private:
  const std::vector<uint64_t> &m_inHapPanel;
  const unsigned m_numSites = NUMSITES;

  // this may be more than total number of samples in the case of a reference
  // panel
  const unsigned m_totalNumHaps;
  const unsigned m_totalNumSamps;

  const unsigned m_numCycles;

  // transition matrix
  const std::vector<float> &m_tran;

  // mutation matrix.  pc[4][4] in impute.cpp
  const float (*m_mutationMat)[4][4];

  Sampler &m_sampler;
  unsigned m_nextSampIdx;
  GLPack m_glPack;
  gsl_rng &m_rng;

  /*
    Pulls out GLs for next run and repackages them ready to call cuda code
  */
  void CheckDevice() const;
  void CopyTranToDevice() const;
  void CopyMutationMatToDevice() const;
  void CopyPackedGLsToDevice(const std::vector<char> &packedGLs) const;

public:
  /*
    HMMLike constructor.  takes input to a 2*NxnumWordsPerHap matrix.
    It is assumed that the word length is 64 bit.
  */
  HMMLike(const std::vector<uint64_t> &hapPanel, unsigned numHaps,
          const std::vector<float> &GLs, unsigned numSamps,
          unsigned sampleStride, unsigned numCycles,
          const std::vector<float> &tran, const float (*mutationMat)[4][4],
          Sampler &sampler, gsl_rng &rng);

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
