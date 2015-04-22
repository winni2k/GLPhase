#include "hmmLike.hpp"

using namespace std;

namespace HMMLikeHelper {
vector<unsigned> transposeHapIdxs(const std::vector<unsigned> &hapIdxs) {
  const unsigned hapsPerSamp = 4;
  assert(hapIdxs.size() % hapsPerSamp == 0);
  const unsigned numSamps = hapIdxs.size() / hapsPerSamp;

  // this is the transposed vector to return
  vector<unsigned> retHapIdxs(hapIdxs.size());

  // transpose
  for (unsigned sampNum = 0; sampNum < numSamps; ++sampNum)
    for (unsigned hapNum = 0; hapNum < hapsPerSamp; ++hapNum)
      retHapIdxs[hapNum + sampNum * hapsPerSamp] =
          hapIdxs[sampNum + hapNum * numSamps];

  // return by reference due to c++11 magic (we hope)
  return retHapIdxs;
}
}

/*
  Constructor
  Checks that a CUDA enabled device exists
  Copies the transition matrix to the device
  Copies the mutation matrix to the device
  Stores everything else for copying later
 */
HMMLike::HMMLike(const vector<uint64_t> &hapPanel, unsigned numHaps,
                 GLPack &glPack, unsigned numCycles, const vector<float> &tran,
                 const float (*mutationMat)[4][4], shared_ptr<Sampler> &sampler,
                 gsl_rng &rng)
    : m_inHapPanel(hapPanel), m_totalNumHaps(numHaps), m_glPack(glPack),
      m_totalNumSamps(m_glPack.GetNumSamps()), m_numCycles(numCycles),
      m_tran(tran), m_mutationMat(mutationMat), m_sampler(sampler), m_rng(rng) {

  // Checking expectations.
  assert(m_numCycles > 0);
  assert(m_inHapPanel.size() <= m_totalNumHaps * NUMSITES / WORDSIZE);
  assert(NUMSITES >= m_glPack.GetNumSites());
  assert(m_glPack.GetCodeBook().size() == 1 << BITSPERCODE);
  assert(m_glPack.GetNumBitsPerGL() == BITSPERCODE);

  // make sure we have a K20 or better installed
  // also define some constants
  CheckDevice();

  // copy the transition matrix to constant memory on device
  HMMLikeCUDA::CopyTranToDevice(m_tran);

  // copy the mutation matrixt to constant memory on the device
  HMMLikeCUDA::CopyMutationMatToDevice(m_mutationMat);

  // copy code book for VQ to device
  HMMLikeCUDA::CopyCodeBookToDevice(glPack.GetCodeBook());
  
  // copy all strided GLs across if the sample stride is equal to the number of
  // samples
  if (m_totalNumSamps == m_glPack.GetSampleStride())
      HMMLikeCUDA::CopyPackedGLsToDevice(m_glPack.GetPackedGLs());
}

// returns range of samples sampled, and vector of four unsigned hap indices per
// sample
vector<unsigned> HMMLike::RunHMMOnSamples(unsigned &firstSampIdx,
                                          unsigned &lastSampIdx) {

  size_t numRunSamps = m_glPack.GetSampleStride();

  // if we are running on complete sample set
  if (m_glPack.GetSampleStride() == m_totalNumSamps) {
    firstSampIdx = 0;
    lastSampIdx = numRunSamps - 1;
  } else {

    // get next set of GLs
    firstSampIdx = m_glPack.GetNextSampIdx();
    HMMLikeCUDA::CopyPackedGLsToDevice(m_glPack.GetPackedGLs());
    lastSampIdx = m_glPack.GetLastSampIdx();
  }

  // generate initial list of four hap nums for kernel to use
  // generate list of numCycles haplotype nums for the kernel to choose from
  vector<unsigned> hapIdxs(4 * numRunSamps);
  vector<unsigned> extraPropHaps(m_glPack.GetSampleStride() * m_numCycles);
  unsigned sampIdx = firstSampIdx;
  for (unsigned sampNum = 0; sampNum < numRunSamps; ++sampNum, ++sampIdx) {

    // fill the initial haps
    for (unsigned propHapIdx = 0; propHapIdx < 4; ++propHapIdx)
        hapIdxs.at(sampNum + propHapIdx * numRunSamps) =
          m_sampler->SampleHap(sampIdx);

    // fill the proposal haps
    for (unsigned cycleIdx = 0; cycleIdx < m_numCycles; ++cycleIdx)
      extraPropHaps[cycleIdx * m_glPack.GetSampleStride() + sampNum] =
          m_sampler->SampleHap(sampIdx);
  }

  // run kernel
  HMMLikeCUDA::RunHMMOnDevice(
      m_inHapPanel, extraPropHaps, m_glPack.GetNumSites(),
      m_glPack.GetSampleStride(), m_numCycles, hapIdxs, gsl_rng_get(&m_rng));

  // return
  return HMMLikeHelper::transposeHapIdxs(hapIdxs);
}


