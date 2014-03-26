#include "hmmLike.hpp"

using namespace std;

HMMLike::HMMLike(const vector<uint64_t> &hapPanel, unsigned numHaps,
                 const vector<float> &GLs, unsigned numSamps,
                 unsigned sampleStride, unsigned numCycles,
                 const vector<float> &tran, const float (*mutationMat)[4][4],
                 Sampler &sampler)
    : m_inHapPanel(hapPanel), m_totalNumHaps{ numHaps },
      m_totalNumSamps{ numSamps }, m_numCycles{ numCycles }, m_tran(tran),
      m_mutationMat(mutationMat), m_sampler(sampler),
      m_glPack(GLs, numSamps, sampleStride) {

  // Checking expectations.
  assert(m_inHapPanel.size() == m_numSites / WORDSIZE * m_totalNumHaps);
  assert(GLs.size() == m_numSites * 3 * m_totalNumSamps);

  // make sure we have a K20 or better installed
  CheckDevice();

  // copy the transition matrix to constant memory on device
  CopyTranToDevice();

  // copy the mutation matrixt to constant memory on the device
  CopyMutationMatToDevice();
}

// returns range of samples sampled, and vector of four unsigned hap indices per
// sample
vector<unsigned> HMMLike::RunHMMOnSamples(unsigned &firstSampIdx,
                                          unsigned &lastSampIdx) {

  // get next set of GLs
  firstSampIdx = m_glPack.GetNextSampIdx();
  vector<char> packedGLs = m_glPack.GetPackedGLs();
  lastSampIdx = m_glPack.GetLastSampIdx();

  // generate initial list of four hap nums for kernel to use
  // generate list of numCycles haplotype nums for the kernel to choose from
  vector<unsigned> hapIdxs(4 * m_totalNumSamps);
  assert(hapIdxs.size() == (lastSampIdx - firstSampIdx + 1) * 4);
  vector<unsigned> extraPropHaps(m_totalNumSamps * m_numCycles);
  for (unsigned sampIdx = firstSampIdx; sampIdx < lastSampIdx; ++sampIdx) {

    // fill the initial haps
    for (unsigned propHapIdx = 0; propHapIdx < 4; ++propHapIdx)
      hapIdxs[sampIdx * 4 + propHapIdx] = m_sampler.SampleHap(sampIdx);

    // fill the proposal haps
    for (unsigned cycleIdx = 0; cycleIdx < m_numCycles; ++cycleIdx)
      extraPropHaps[cycleIdx * m_glPack.GetSampleStride() + sampIdx] =
          m_sampler.SampleHap(sampIdx);
  }

  // run kernel
  HMMLikeCUDA::RunHMMOnDevice(packedGLs, m_inHapPanel, extraPropHaps,
                              m_glPack.GetNumSites(), m_totalNumSamps, m_numCycles,
                              firstSampIdx, hapIdxs);

  // unpack results

  // return
  return hapIdxs;
}

void HMMLike::CheckDevice() const { HMMLikeCUDA::CheckDevice(); }
void HMMLike::CopyTranToDevice() const {

  cudaError_t err = HMMLikeCUDA::CopyTranToDevice(m_tran);
  if (err != cudaSuccess) {
    stringstream outerr(
        "Could not copy transition matrix to device with error: ");
    outerr << err;
    throw myException(outerr.str());
  }
}
void HMMLike::CopyMutationMatToDevice() const {

  cudaError_t err = HMMLikeCUDA::CopyMutationMatToDevice(m_mutationMat);
  if (err != cudaSuccess) {
    stringstream outerr(
        "Could not copy mutation matrix to device with error: ");
    outerr << err;
    throw myException(outerr.str());
  }
}
