#include "hmmLike.hpp"

using namespace std;

HMMLike::HMMLike(const vector<uint64_t> &hapPanel, unsigned numHaps,
                 GLPack &glPack, unsigned numCycles, const vector<float> &tran,
                 const float (*mutationMat)[4][4], Sampler &sampler,
                 gsl_rng &rng)
    : m_inHapPanel(hapPanel), m_totalNumHaps(numHaps), m_glPack(glPack),
      m_totalNumSamps(m_glPack.GetNumSamps()), m_numCycles(numCycles),
      m_tran(tran), m_mutationMat(mutationMat), m_sampler(sampler), m_rng(rng) {

  // Checking expectations.
  assert(m_numCycles > 0);
  assert(m_inHapPanel.size() == m_totalNumHaps * NUMSITES / WORDSIZE);
  assert(NUMSITES == m_glPack.GetNumSites());

  // make sure we have a K20 or better installed
  // also define some constants
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

  unsigned numRunSamps = lastSampIdx - firstSampIdx + 1;
  // generate initial list of four hap nums for kernel to use
  // generate list of numCycles haplotype nums for the kernel to choose from
  vector<unsigned> hapIdxs(4 * numRunSamps);
  vector<unsigned> extraPropHaps(m_glPack.GetSampleStride() * m_numCycles);
  unsigned sampIdx = firstSampIdx;
  for (unsigned sampNum = 0; sampNum < numRunSamps; ++sampNum, ++sampIdx) {

    // fill the initial haps
    for (unsigned propHapIdx = 0; propHapIdx < 4; ++propHapIdx)
      hapIdxs[sampNum + propHapIdx * numRunSamps] =
          m_sampler.SampleHap(sampIdx);

    // fill the proposal haps
    for (unsigned cycleIdx = 0; cycleIdx < m_numCycles; ++cycleIdx)
      extraPropHaps[cycleIdx * m_glPack.GetSampleStride() + sampNum] =
          m_sampler.SampleHap(sampIdx);
  }

  // run kernel
  HMMLikeCUDA::RunHMMOnDevice(
      packedGLs, m_inHapPanel, extraPropHaps, m_glPack.GetNumSites(),
      m_glPack.GetSampleStride(), m_numCycles, hapIdxs, gsl_rng_get(&m_rng));

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
