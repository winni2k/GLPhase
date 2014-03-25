
#include "gtest/gtest.h"
#include "hmmLike.hpp"
#include <cuda_runtime.h>

using namespace std;

namespace HMMLikeCUDATest {
extern bool UnpackGLs(char GLset, float *GLs);
extern cudaError_t CopyTranToHost(vector<float> &tran);
extern cudaError_t CopyMutMatToHost(vector<float> &mutMat);
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

/*
TEST(HMMLike, createsOK) {

  hapPanel panel1;
  vector<vector<char> > haps;
  vector<snp> sites = { 0, 1, 2, 3 };
  vector<string> sampIDs = { "s1", "s2", "s3", "s4" };
  for 
  haps.push_back(
  panel1.Init();
}
*/
