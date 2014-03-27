
#include "gtest/gtest.h"
#include "hmmLike.hpp"
#include "sampler.hpp"
#include <cuda_runtime.h>
#include <gsl/gsl_rng.h>

using namespace std;

namespace HMMLikeCUDATest {
extern "C" bool UnpackGLs(char GLset, float *GLs);
extern "C" void FillEmit(const vector<float> &GLs, vector<float> &emit);
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

TEST(HMMLike, FillEmitFillsOK) {

  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, time(NULL));

  float mutMat[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      mutMat[i][j] = 1.0f;
  HMMLikeCUDA::CopyMutationMatToDevice(&mutMat);

  vector<float> GLs(3);
  for (auto &GL : GLs)
    GL = gsl_rng_uniform(rng);
  vector<float> emit(4);
  HMMLikeCUDATest::FillEmit(GLs, emit);
  for (int i = 0; i < 4; ++i)
    EXPECT_FLOAT_EQ((GLs[0] + 2 * GLs[1] + GLs[2]), emit[i]);
}

TEST(HMMLike, createsOK) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, time(NULL));

  const unsigned numSamps = 2;
  const unsigned numHaps = 4;
  const unsigned numSites = 512;
  const unsigned numWords = numSites / 64;
  vector<uint64_t> hapPanel(numWords * numHaps);
  vector<float> GLs(3 * numSites * numSamps);
  unsigned sampleStride = 2;
  unsigned numCycles = 1000;
  vector<float> tran(numSites * 3);
  float mutMat[4][4];
  UnifSampler sampler(rng, numSamps, numHaps);

  // initialize mutation matrix
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      mutMat[i][j] = gsl_rng_uniform(rng);

  for (auto &GL : GLs)
    GL = gsl_rng_uniform(rng);

  // test hmmLike function first

  HMMLike hmmLike(hapPanel, numHaps, GLs, numSamps, sampleStride, numCycles,
                  tran, &mutMat, sampler, *rng);

  unsigned firstSampIdx;
  unsigned lastSampIdx;
  vector<unsigned> hapIdxs = hmmLike.RunHMMOnSamples(firstSampIdx, lastSampIdx);
  ASSERT_EQ(0, firstSampIdx);
  ASSERT_EQ(1, lastSampIdx);
  ASSERT_EQ(numSamps * 4, hapIdxs.size());
  for (int i = 0; i < 8; i += 2) {
    EXPECT_GT(4, hapIdxs[i]);
    EXPECT_LT(1, hapIdxs[i]);
  }
  for (int i = 1; i < 8; i += 2)
    EXPECT_GT(2, hapIdxs[i]);

  // ok, let's try this again with a larger data set

  // this should be 40 haplotypes
  const unsigned bigNumHaps = numHaps * 10;
  hapPanel.resize(numWords * bigNumHaps);
  hapPanel[5 * numWords + 5] = ~hapPanel[5 * numWords + 5];
  hapPanel[6 * numWords + 5] = ~hapPanel[6 * numWords + 5];
  hapPanel[7 * numWords + 7] = ~hapPanel[7 * numWords + 7];
  hapPanel[8 * numWords + 7] = ~hapPanel[8 * numWords + 7];

  // initialize site mutation probability matrix
  // diagonal is chance of no mutation
  // the diagonal "rotated by 90 degrees" is the chance of both positions
  // mutating
  // all other entries are chance of just one mutation
  const float mu = 0.001;
  mutMat[0][0] = mutMat[1][1] = mutMat[2][2] = mutMat[3][3] =
      (1 - mu) * (1 - mu); //	  probability of mutating no positions for each
  // parents haplotype
  mutMat[0][1] = mutMat[0][2] = mutMat[1][0] = mutMat[1][3] = mutMat[2][0] =
      mutMat[2][3] = mutMat[3][1] = mutMat[3][2] =
          mu * (1 - mu); //	  probability of mutating one position
  // for each parental haplotype
  mutMat[0][3] = mutMat[1][2] = mutMat[2][1] = mutMat[3][0] =
      mu * mu; //	  probability of mutating both positions for each
  // parental haplotype

  unsigned offset = 3 * numSites * 5 + numWords * 5 * 64;
  for (int i = 0; i < 64; ++i) {
    GLs[i + offset + 2] = 10;                      // sample 5, 5th word
    GLs[i + offset + 2 + numSites] = 10;           // sample 6, 5th word
    GLs[i + offset + 2 + 2 * numSites + 128] = 10; // sample 7, 7th word
    GLs[i + offset + 2 + 3 * numSites + 128] = 10; // sample 8, 7th word
  }
  UnifSampler sampler2(rng, bigNumHaps / 2, bigNumHaps);
  HMMLike hmmLike2(hapPanel, bigNumHaps, GLs, numSamps, sampleStride, numCycles,
                   tran, &mutMat, sampler2, *rng);

  firstSampIdx = 0;
  lastSampIdx = 0;
  vector<unsigned> hapIdxs2 =
      hmmLike2.RunHMMOnSamples(firstSampIdx, lastSampIdx);
  ASSERT_EQ(0, firstSampIdx);
  ASSERT_EQ(1, lastSampIdx);
  ASSERT_EQ(numSamps * 4, hapIdxs2.size());
  for (int i = 0; i < 8; i += 2) {
    EXPECT_GT(7, hapIdxs[i]);
    EXPECT_LT(4, hapIdxs[i]);
  }
  for (int i = 1; i < 8; i += 2) {
    EXPECT_GT(9, hapIdxs[i]);
    EXPECT_LT(6, hapIdxs[i]);
  }
}
