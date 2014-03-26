
#include "gtest/gtest.h"
#include "hmmLike.hpp"
#include "sampler.hpp"
#include <cuda_runtime.h>
#include <gsl/gsl_rng.h>

using namespace std;

namespace HMMLikeCUDATest {
extern "C" bool UnpackGLs(char GLset, float *GLs);
extern "C" cudaError_t CopyTranToHost(vector<float> &tran);
extern "C" cudaError_t CopyMutMatToHost(vector<float> &mutMat);
}

TEST(FindDevice, FoundDevice) { HMMLikeCUDA::CheckDevice(); }

TEST(CopyToTransitionMat, CopySuccess) {

  vector<float> tran(NUMSITES * 3);
  for (int i = 0; i < NUMSITES * 3; ++i)
    tran[i] = i;
  ASSERT_NEAR(tran[513], 513, 0.001);

  HMMLikeCUDA::CopyTranToDevice(tran);

  vector<float> postDTran(NUMSITES * 3);
  ASSERT_EQ(HMMLikeCUDATest::CopyTranToHost(postDTran),
            0); // 0 equals cudaSuccess
  ASSERT_EQ(postDTran.size(), NUMSITES * 3);
  for (unsigned i = 0; i < postDTran.size(); ++i)
    EXPECT_FLOAT_EQ(i, postDTran[i]);
}

TEST(CopyToMutMat, CopySuccess) {

  float pc[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      pc[i][j] = i + 4 * j;

  ASSERT_FLOAT_EQ(pc[2][3], 2 + 4 * 3);
  HMMLikeCUDA::CopyMutationMatToDevice(&pc);

  vector<float> postDMutMat(4 * 4);
  ASSERT_EQ(HMMLikeCUDATest::CopyMutMatToHost(postDMutMat),
            0); // 0 equals cudaSuccess
  ASSERT_EQ(4 * 4, postDMutMat.size());
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      EXPECT_FLOAT_EQ(i + 4 * j, postDMutMat[i + 4 * j]);
}

// testing UnpackGLs
TEST(UnpackGLs, UnpackOk) {

  // testing silly vaules
  char num = 255;
  float GLs[3];
  for (int i = 0; i < 4; ++i)
    GLs[i] = 0;
  ASSERT_TRUE(HMMLikeCUDATest::UnpackGLs(num, GLs));
  EXPECT_FLOAT_EQ(15.5f / 16, GLs[0]);
  EXPECT_FLOAT_EQ(15.5f / 16, GLs[1]);
  EXPECT_FLOAT_EQ(0, GLs[2]);

  // testing realistic values
  num = 17;
  for (int i = 0; i < 4; ++i)
    GLs[i] = 0;
  ASSERT_TRUE(HMMLikeCUDATest::UnpackGLs(num, GLs));
  EXPECT_FLOAT_EQ(1.5f / 16, GLs[0]);
  EXPECT_FLOAT_EQ(1.5f / 16, GLs[1]);
  EXPECT_FLOAT_EQ(13.0f / 16, GLs[2]);
}

TEST(HMMLike, createsOK) {

  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, time(NULL));

  const unsigned numSamps = 2;
  const unsigned numHaps = 4;
  const unsigned numSites = 512;
  vector<uint64_t> hapPanel(numSites * numHaps / 64);
  vector<float> GLs(3 * numSites * numSamps);
  unsigned sampleStride = 2;
  unsigned numCycles = 10;
  vector<float> tran(numSites * 3);
  float mutMat[4][4];
  UnifSampler sampler(rng, numSamps, numHaps);

  // initialize mutation matrix
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      mutMat[i][j] = gsl_rng_uniform(rng);

  for (auto &GL : GLs)
    GL = gsl_rng_uniform(rng);

  HMMLike hmmLike(hapPanel, numHaps, GLs, numSamps, sampleStride, numCycles,
                  tran, &mutMat, sampler);

  unsigned firstSampIdx;
  unsigned lastSampIdx;
  vector<unsigned> hapIdxs =
      hmmLike.RunHMMOnSamples(firstSampIdx, lastSampIdx);
  ASSERT_EQ(0, firstSampIdx);
  ASSERT_EQ(1, lastSampIdx);
  ASSERT_EQ(numSamps * 4, hapIdxs.size());
  for (int i = 0; i < 8; ++i)
    EXPECT_EQ(i, hapIdxs[i]);
}
