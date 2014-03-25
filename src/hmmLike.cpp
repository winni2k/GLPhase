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

  vector<unsigned> hapIdxs;

  firstSampIdx = m_glPack.GetNextSampIdx();
  vector<char> packedGLs = m_glPack.GetPackedGLs();
  lastSampIdx = m_glPack.GetLastSampIdx();

  // run kernel

  // unpack results

  // return
  assert(hapIdxs.size() == (lastSampIdx - firstSampIdx + 1) * 4);
  return hapIdxs;
}

void HMMLike::CheckDevice() { HMMLikeCUDA::CheckDevice(); }
void HMMLike::CopyTranToDevice() {

  cudaError_t err = HMMLikeCUDA::CopyTranToDevice(m_tran);
  if (err != cudaSuccess) {
    stringstream outerr(
        "Could not copy transition matrix to device with error: ");
    outerr << err;
    throw myException(outerr.str());
  }
}
void HMMLike::CopyMutationMatToDevice() {

  cudaError_t err = HMMLikeCUDA::CopyMutationMatToDevice(m_mutationMat);
  if (err != cudaSuccess) {
    stringstream outerr(
        "Could not copy mutation matrix to device with error: ");
    outerr << err;
    throw myException(outerr.str());
  }
}
