
#include "gtest/gtest.h"
#include "hmmLike.hpp"
#include "sampler.hpp"
#include "glPack.hpp"
#include <cuda_runtime.h>
#include <gsl/gsl_rng.h>
#include <memory>
#include <limits>
#include <utility>
#include <thrust/host_vector.h>

using namespace std;

namespace HMMLikeCUDATest {
extern bool UnpackGLs(unsigned char GLset, float *GLs);
extern void FillEmit(const vector<float> &GLs, vector<float> &emit);
extern cudaError_t CopyTranToHost(vector<float> &tran);
extern cudaError_t CopyMutMatToHost(vector<float> &mutMat);
extern float CallHMMLike(unsigned idx, const unsigned (*hapIdxs)[4],
                         unsigned packedGLStride,
                         const vector<uint64_t> &h_hapPanel);
extern void UnpackGLsWithCodeBook(uint32_t GLcodes, vector<float> &GLs,
                                  unsigned char glIdx);
extern void FillRNs(thrust::host_vector<unsigned> &h_rns, size_t numRNs);
}

TEST(FindDevice, FoundDevice) { HMMLikeCUDA::CheckDevice(); }

TEST(CopyToTransitionMat, CopySuccess) {

  vector<float> tran(NUMSITES * 3);
  for (int i = 0; i < NUMSITES * 3; ++i)
    tran[i] = 0.25;
  ASSERT_NEAR(tran[513], 0.25, 0.001);

  HMMLikeCUDA::CopyTranToDevice(tran);

  vector<float> postDTran(NUMSITES * 3);
  ASSERT_EQ(HMMLikeCUDATest::CopyTranToHost(postDTran),
            0); // 0 equals cudaSuccess
  ASSERT_EQ(postDTran.size(), NUMSITES * 3);
  for (unsigned i = 0; i < postDTran.size(); ++i)
    EXPECT_FLOAT_EQ(0.25, postDTran[i]);
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
TEST(UnpackGLsWithCodeBook, UnpackOk) {

  // testing silly values
  uint32_t num = 1;
  vector<pair<float, float> > codeBook(exp2(BITSPERCODE));
  for (int i = 0; i < codeBook.size(); ++i) {
    codeBook[i].first = static_cast<float>(i) / 10;
    codeBook[i].second = static_cast<float>(i) * 2 / 100;
  }
  HMMLikeCUDA::CopyCodeBookToDevice(codeBook);

  vector<float> GLs(3, 0);
  HMMLikeCUDATest::UnpackGLsWithCodeBook(num, GLs,
                                         (UINT32T_SIZE / BITSPERCODE) - 1);
  EXPECT_FLOAT_EQ(0.1, GLs[0]);
  EXPECT_FLOAT_EQ(0.02, GLs[1]);
  EXPECT_FLOAT_EQ(0.88, GLs[2]);

  // testing realistic values
  num = 19088743;
  for (int i = 0; i < 4; ++i)
    GLs[i] = 0;
  for (size_t i = 0; i < UINT32T_SIZE / BITSPERCODE; ++i) {
    HMMLikeCUDATest::UnpackGLsWithCodeBook(num, GLs, i);
    EXPECT_FLOAT_EQ(static_cast<float>(i) / 10, GLs[0]);
    EXPECT_FLOAT_EQ(static_cast<float>(i) * 2 / 100, GLs[1]);
    EXPECT_FLOAT_EQ(1.0f - static_cast<float>(i) / 10 -
                        static_cast<float>(i) * 2 / 100,
                    GLs[2]);
  }
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

TEST(HMMLike, rngTest) {

  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, 112);

  const unsigned numSamps = 2;
  const unsigned numHaps = 4;
  const unsigned numSites = 1024;
  const unsigned wordSize = 64;
  const unsigned numWords = numSites / wordSize;
  vector<uint64_t> hapPanel(numWords * numHaps);
  const unsigned sampleStride = numSamps;
  const unsigned numCycles = 100;

  // initialize transition matrix
  vector<float> tran(numSites * 3);
  float rho = 2 * 10e-8;
  for (unsigned m = numSites - 1; m; m--) {
    float rhoTDist = rho * 100; // approximate distance between SNPs
    float r = rhoTDist / (rhoTDist + numHaps);
    tran[m * 3] = (1 - r) * (1 - r);
    tran[m * 3 + 1] = r * (1 - r);
    tran[m * 3 + 2] = r * r; // for each position, transition.  r= alternative,
                             // 1-r= refrence? 4 state HMM with three
                             // transitions at each position
  }

  // initialize mutation matrix
  float mutMat[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      mutMat[i][j] = gsl_rng_uniform(rng);

  // initialize GLs
  vector<float> GLs(3 * numSites * numSamps, 0);
  for (auto &GL : GLs)
    GL = gsl_rng_uniform(rng);
  assert(GLs.size() == 1024 * 3 * 2);

  GLPackHelper::Init init(GLs, *rng);
  init.numSamps = numSamps;
  init.sampleStride = sampleStride;
  init.useVQ = true;
  init.numBitsPerGL = BITSPERCODE;
  GLPack glPack1(init);

  vector<uint32_t> packedGLs = glPack1.GetPackedGLs();
  ASSERT_EQ(numSites * BITSPERCODE / UINT32T_SIZE * numSamps, packedGLs.size());
  ASSERT_EQ(1 << BITSPERCODE, glPack1.GetCodeBook().size());

  shared_ptr<Sampler> sampler =
      make_shared<UnifSampler>(rng, numSamps, numHaps);

  for (unsigned i = 0; i < numHaps * 10; ++i) {
    unsigned val = sampler->SampleHap(0);
    ASSERT_LT(1, val);
    ASSERT_GT(numHaps, val);
  }
  {
    // now create hmmLike so we can test the RNG
    HMMLike hmmLike(hapPanel, numHaps, glPack1, numCycles, tran, &mutMat,
                    sampler, *rng);

    // getting 32 random numbers
    thrust::host_vector<unsigned> rns;
    HMMLikeCUDATest::FillRNs(rns, 32);

    cout << "32 random numbers from MTRNG:";
    for (auto rn : rns)
      cout << "\t" << rn;
    cout << endl;

    cout << "hex of those numbers:";
    for (auto rn : rns)
      cout << "\t" << std::hex << rn << std::dec;
    cout << endl;
  }
  {
    // and do it again...
    HMMLike hmmLike(hapPanel, numHaps, glPack1, numCycles, tran, &mutMat,
                    sampler, *rng);

    // getting 32 random numbers
    thrust::host_vector<unsigned> rns;
    HMMLikeCUDATest::FillRNs(rns, 32);

    cout << "32 random numbers from MTRNG:";
    for (auto rn : rns)
      cout << "\t" << rn;
    cout << endl;

    cout << "hex for those numbers:";
    for (auto rn : rns)
      cout << "\t" << std::hex << rn << std::dec;
    cout << endl;
  }
}


TEST(HMMLike, createsOK) {
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, 112);

  const unsigned numSamps = 2;
  const unsigned numHaps = 4;
  const unsigned numSites = 1024;
  const unsigned wordSize = 64;
  const unsigned numWords = numSites / wordSize;
  vector<uint64_t> hapPanel(numWords * numHaps);
  const unsigned sampleStride = numSamps;
  const unsigned numCycles = 100;

  // initialize transition matrix
  vector<float> tran(numSites * 3);
  float rho = 2 * 10e-8;
  for (unsigned m = numSites - 1; m; m--) {
    float rhoTDist = rho * 100; // approximate distance between SNPs
    float r = rhoTDist / (rhoTDist + numHaps);
    tran[m * 3] = (1 - r) * (1 - r);
    tran[m * 3 + 1] = r * (1 - r);
    tran[m * 3 + 2] = r * r; // for each position, transition.  r= alternative,
                             // 1-r= refrence? 4 state HMM with three
                             // transitions at each position
  }
  /* debugging
  cout << "printing tran sums: ";
  for (unsigned m = 0; m < tran.size(); m += 3)
    cout << tran[m] + 2 * tran[m + 1] + tran[m + 2] << " ";
  cout << endl;
  */

  // initialize mutation matrix
  float mutMat[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      mutMat[i][j] = gsl_rng_uniform(rng);

  vector<float> GLs(3 * numSites * numSamps, 0);
  for (auto &GL : GLs)
    GL = gsl_rng_uniform(rng);
  assert(GLs.size() == 1024 * 3 * 2);

  GLPackHelper::Init init(GLs, *rng);
  init.numSamps = numSamps;
  init.sampleStride = sampleStride;
  init.useVQ = true;
  init.numBitsPerGL = BITSPERCODE;
  GLPack glPack1(init);

  vector<uint32_t> packedGLs = glPack1.GetPackedGLs();
  ASSERT_EQ(numSites * BITSPERCODE / UINT32T_SIZE * numSamps, packedGLs.size());
  ASSERT_EQ(1 << BITSPERCODE, glPack1.GetCodeBook().size());

  shared_ptr<Sampler> sampler =
      make_shared<UnifSampler>(rng, numSamps, numHaps);

  for (unsigned i = 0; i < numHaps * 10; ++i) {
    unsigned val = sampler->SampleHap(0);
    ASSERT_LT(1, val);
    ASSERT_GT(numHaps, val);
  }
  {
    // now test HMMLike functionality
    HMMLike hmmLike(hapPanel, numHaps, glPack1, numCycles, tran, &mutMat,
                    sampler, *rng);

    unsigned firstSampIdx;
    unsigned lastSampIdx;
    {
      vector<unsigned> hapIdxs =
          hmmLike.RunHMMOnSamples(firstSampIdx, lastSampIdx);
      ASSERT_EQ(0, firstSampIdx);
      ASSERT_EQ(1, lastSampIdx);
      ASSERT_EQ(numSamps * 4, hapIdxs.size());
      for (int i = 0; i < 4; ++i) {
        EXPECT_GT(4, hapIdxs[i]);
        EXPECT_LT(1, hapIdxs[i]);
      }
      for (int i = 4; i < 8; ++i)
        EXPECT_GT(2, hapIdxs[i]);
    }

    {
      // test hmmLike function
      GLPackHelper::Init init2(GLs, *rng);
      init2.numSamps = numSamps;
      init2.sampleStride = sampleStride;
      init2.useVQ = true;
      init2.numBitsPerGL = BITSPERCODE;

      GLPack glPack0(init2);
      auto packedGLs = glPack0.GetPackedGLs();
      unsigned sampIdx = 0;
      unsigned fixedHapIdxs[4];
      for (int i = 0; i < 4; ++i)
        fixedHapIdxs[i] = 2;
      HMMLikeCUDA::CopyPackedGLsToDevice(packedGLs);
      HMMLikeCUDA::CopyCodeBookToDevice(glPack0.GetCodeBook());
      float like = HMMLikeCUDATest::CallHMMLike(
          sampIdx, &fixedHapIdxs, glPack0.GetSampleStride(), hapPanel);
      ASSERT_GE(1, like);
    }
  }
  //  cout << "Likelihood of Model: " << like << endl << endl;

  // ok, let's try this again with a larger data set

  // this should be 40 haplotypes
  const unsigned bigNumHaps = 12;
  hapPanel.resize(numWords * bigNumHaps);
  hapPanel.at(5 *numWords) = ~0;
  hapPanel.at(6 *numWords) = ~0;
  hapPanel.at(5 *numWords + 5) = ~0;
  hapPanel.at(6 *numWords + 5) = ~0;

  hapPanel.at(7 *numWords + 1) = ~0;
  hapPanel.at(8 *numWords + 1) = ~0;
  hapPanel.at(7 *numWords + 7) = ~0;
  hapPanel.at(8 *numWords + 7) = ~0;

  // initialize site mutation probability matrix
  // diagonal is chance of no mutation
  // the diagonal "rotated by 90 degrees" is the chance of both positions
  // mutating
  // all other entries are chance of just one mutation
  const unsigned avgSiteDist = 1000;
  float mu = 0;
  for (unsigned i = 1; i < bigNumHaps; i++)
    mu += 1.0 / i;
  mu = 1 / mu;
  rho = 0.5 * mu * (numSites - 1) / (numSites * avgSiteDist) / 1;
  mu =
      mu / (bigNumHaps + mu); // rho is recombination rate?  mu is mutation rate

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

  // now try creating tran based on mu
  for (unsigned m = numSites - 1; m; m--) {
    float rhoTDist = rho * avgSiteDist; // approximate distance between SNPs
    float r = rhoTDist / (rhoTDist + numHaps);
    tran[m * 3] = (1 - r) * (1 - r);
    tran[m * 3 + 1] = r * (1 - r);
    tran[m * 3 + 2] = r * r; // for each position, transition.  r= alternative,
                             // 1-r= refrence? 4 state HMM with three
                             // transitions at each position
  }

  // reset all GLs as HOM REF
  /*
    for (int i = 0; i < GLs.size(); i += 3) {
      GLs[i] = 1;
      GLs[i + 1] = 0;
      GLs[i + 2] = 0;
    }
  */
  {
    const float highExp = 100;
    unsigned offset = wordSize * 5 * 3;
    for (unsigned i = 0; i < wordSize * 3; i += 3) {
      // sample 0, 1st word set to ALT/ALT
      GLs.at(i + 2) = highExp;
      // sample 0, 6th word set to ALT/ALT
      GLs.at(offset + i + 2) = highExp;
      // sample 1, 2nd word set to ALT/ALT
      GLs.at(numSites * 3 + wordSize * 3 + i + 2) = highExp;
      // sample 1, 8th word set to ALT/ALT
      GLs.at(offset + i + 2 + numSites * 3 + wordSize * 2 * 3) = highExp;
    }
    assert(GLs.size() == numSites * numSamps * 3);
  }

  shared_ptr<Sampler> sampler2 =
      make_shared<UnifSampler>(rng, bigNumHaps / 2, bigNumHaps);

  // make sure sampler is giving us what we expect
  vector<unsigned> sampledHaps;
  for (unsigned i = 0; i < bigNumHaps * 10; ++i)
    sampledHaps.push_back(sampler2->SampleHap(0));
  ASSERT_EQ(bigNumHaps - 1,
            *std::max_element(sampledHaps.begin(), sampledHaps.end()));

  {
    gsl_rng_set(rng, 1);
    GLPackHelper::Init init(GLs, *rng);
    init.numSamps = numSamps;
    init.sampleStride = sampleStride;
    init.useVQ = true;
    init.numBitsPerGL = BITSPERCODE;

    GLPack glPack2(init);
    HMMLike hmmLike2(hapPanel, bigNumHaps, glPack2, numCycles, tran, &mutMat,
                     sampler2, *rng);

    unsigned firstSampIdx = 0;
    unsigned lastSampIdx = 0;
    {
      vector<unsigned> hapIdxs2 =
          hmmLike2.RunHMMOnSamples(firstSampIdx, lastSampIdx);
      ASSERT_EQ(0, firstSampIdx);
      ASSERT_EQ(1, lastSampIdx);
      ASSERT_EQ(numSamps * 4, hapIdxs2.size());

      // only one of the father and mother pairs needs to be correct
      for (unsigned i = 0; i < 4; i += 2) {
        cout << "Hap idx " << i << " is " << hapIdxs2[i] << endl;
        cout << "Hap idx " << i + 1 << " is " << hapIdxs2[i + 1] << endl;
        EXPECT_TRUE((7 > hapIdxs2.at(i) && 4 < hapIdxs2.at(i)) ||
                    (7 > hapIdxs2.at(i + 1) && 4 < hapIdxs2.at(i + 1)));
      }

      for (unsigned i = 4; i < 8; i += 2) {
        cout << "Hap idx " << i << " is " << hapIdxs2[i] << endl;
        cout << "Hap idx " << i + 1 << " is " << hapIdxs2[i + 1] << endl;
        EXPECT_TRUE((9 > hapIdxs2.at(i) && 6 < hapIdxs2.at(i)) ||
                    (9 > hapIdxs2.at(i + 1) && 6 < hapIdxs2.at(i + 1)));
      }
    }
  }

  /*
    Now let's create a data set where both father and mother pairs are defined
  */
  // bigNumHaps should be 12
  hapPanel.at(10 *numWords + 2) = ~0;
  hapPanel.at(11 *numWords + 3) = ~0;

  // debugging
  {
    vector<unsigned> haps = { 5, 6, 7, 8, 10, 11 };
    for (auto hapNum : haps) {
      cout << endl << "HapNum: " << hapNum << endl;
      for (unsigned wordNum = 0; wordNum < numWords; ++wordNum) {
        cout << hapPanel.at(numWords * hapNum + wordNum) << ' ';
      }
    }
    cout << endl;
  }

  {
    const float highExp = 100;
    for (unsigned i = 0; i < wordSize * 3; i += 3) {
      // sample 0, 3rd word set to ALT/ALT
      GLs.at(i + 2 + wordSize * 3 * 2) = highExp;
      // sample 1, 4th word set to ALT/ALT
      GLs.at(i + 2 + wordSize * 3 * 3 + numSites * 3) = highExp;
    }
  }

  // print GL states for debugging
  {
    unsigned sampNum = 0;
    for (unsigned siteNum = 0; siteNum < numSites * numSamps;
         siteNum += wordSize) {

      if (siteNum % numSites == 0) {
        cout << endl << "SampNum: " << sampNum << endl;
        ++sampNum;
      }
      cout << GLs[siteNum * 3] << ',' << GLs[siteNum * 3 + 1] << ','
           << GLs[siteNum * 3 + 2] << "\t";
    }
    cout << endl;
  }

  {
    {
      GLPackHelper::Init init(GLs, *rng);
      init.numSamps = numSamps;
      init.sampleStride = sampleStride;
      init.useVQ = true;
      init.numBitsPerGL = BITSPERCODE;

      // let's test the hmmLike function first
      GLPack glPack4(init);
      auto packedGLs = glPack4.GetPackedGLs();
      unsigned sampIdx = 0;
      unsigned fixedHapIdxs[4];
      for (int i = 0; i < 4; ++i)
        fixedHapIdxs[i] = 2;
      HMMLikeCUDA::CopyPackedGLsToDevice(packedGLs);
      HMMLikeCUDA::CopyCodeBookToDevice(glPack4.GetCodeBook());
      float badLike = HMMLikeCUDATest::CallHMMLike(
          sampIdx, &fixedHapIdxs, glPack4.GetSampleStride(), hapPanel);
      ASSERT_GE(1, badLike);

      GLPack glPack5(init);
      auto packedGLs2 = glPack4.GetPackedGLs();
      unsigned fixedHapIdxs2[4] = { 5, 10, 6, 10 };
      HMMLikeCUDA::CopyPackedGLsToDevice(packedGLs2);
      HMMLikeCUDA::CopyCodeBookToDevice(glPack5.GetCodeBook());
      float goodLike = HMMLikeCUDATest::CallHMMLike(
          sampIdx, &fixedHapIdxs2, glPack5.GetSampleStride(), hapPanel);

      ASSERT_GT(goodLike, badLike);
    }

    gsl_rng_set(rng, 114);
    const unsigned numCycles3 = 100;
    GLPackHelper::Init init(GLs, *rng);
    init.numSamps = numSamps;
    init.sampleStride = sampleStride;
    init.useVQ = true;
    init.numBitsPerGL = BITSPERCODE;

    GLPack glPack3(init);
    HMMLike hmmLike3(hapPanel, bigNumHaps, glPack3, numCycles3, tran, &mutMat,
                     sampler2, *rng);

    // getting 32 random numbers
    thrust::host_vector<unsigned> rns;
    HMMLikeCUDATest::FillRNs(rns, 32);

    cout << "\n32 random numbers from MTRNG:";
    for (auto rn : rns)
      cout << "\t" << rn;
    cout << endl;

    cout << "hex for those numbers:";
    for (auto rn : rns)
      cout << "\t" << std::hex << rn << std::dec;
    cout << endl << endl;

    // printing out codebook
    auto codeBook = glPack3.GetCodeBook();
    cout << "CodeBook:" << endl;

    for (unsigned i = 0; i < codeBook.size(); ++i)
      cout << i << ": " << codeBook[i].first << ", " << codeBook[i].second
           << endl;

    // printing out VQ GLs
    cout << "VQ packed GLS:\n";
    auto pGLs = glPack3.GetPackedGLs();
    for (unsigned i = 0; i < pGLs.size(); ++i) {
      if (i % numWords == 0)
        cout << "Word: " << i / numWords << "\t\t";
      cout << std::hex << pGLs.at(i) << std::dec << "\t";
      if (i % numWords == 4) {
        cout << "..." << endl;
        i += numWords - 5;
      }
    }
    cout << endl;

    unsigned firstSampIdx = 0;
    unsigned lastSampIdx = 0;

    gsl_rng_set(rng, 111);

    // this has been moved to hmmLike constructor
    //    HMMLikeCUDA::SetUpRNGs(glPack3.GetSampleStride(), gsl_rng_get(rng));

    vector<unsigned> hapIdxs3 =
        hmmLike3.RunHMMOnSamples(firstSampIdx, lastSampIdx);
    ASSERT_EQ(0, firstSampIdx);
    ASSERT_EQ(1, lastSampIdx);
    ASSERT_EQ(numSamps * 4, hapIdxs3.size());

    // both of the father and mother pairs needs to be correct
    cout << "Hap indices:\n";
    for (auto h : hapIdxs3)
      cout << h << " ";
    cout << endl << endl;
    for (unsigned i = 0; i != 4; i += 2) {
      unsigned hap1 = hapIdxs3[i];
      unsigned hap2 = hapIdxs3[i + 1];
      if (hap1 > hap2)
        swap(hap1, hap2);
      EXPECT_EQ(10, hap2);
      EXPECT_LE(5, hap1);
      EXPECT_GE(6, hap1);

      // these are haps for second sample
      unsigned hap3 = hapIdxs3[i + 4];
      unsigned hap4 = hapIdxs3[i + 4 + 1];
      if (hap3 > hap4)
        swap(hap3, hap4);
      EXPECT_EQ(11, hap4);
      EXPECT_LE(7, hap3);
      EXPECT_GE(8, hap3);
    }
  }
  /*
    Now let's add two more samples to the GLs and try again
  */
  const unsigned numSamps2 = numSamps * 2;
  GLs.resize(3 * numSites * numSamps2);
  for (unsigned i = 3 * numSites * numSamps; i < 3 * numSites * numSamps2; ++i)
    GLs.at(i) = gsl_rng_uniform(rng);
  {
    const float highExp = 100;
    for (unsigned i = 0; i < wordSize * 3; i += 3) {
      // sample 2, 1st,4th and 6th word set to ALT/ALT
      GLs.at(i + 2 + 2 *numSites * 3) = highExp;
      GLs.at(i + 2 + 2 *numSites * 3 + wordSize * 3 * 3) = highExp;
      GLs.at(i + 2 + 2 *numSites * 3 + wordSize * 3 * 5) = highExp;

      // sample 3, 2nd,3rd and 8th word set to ALT/ALT
      GLs.at(i + 2 + 3 *numSites * 3 + wordSize * 3 * 1) = highExp;
      GLs.at(i + 2 + 3 *numSites * 3 + wordSize * 3 * 2) = highExp;
      GLs.at(i + 2 + 3 *numSites * 3 + wordSize * 3 * 7) = highExp;
    }
  }

  /* print GL state
     for debugging...
  // print GL states for debugging
  {
    unsigned sampNum = 0;
    for (unsigned siteNum = 0; siteNum < numSites * numSamps2;
         siteNum += wordSize) {

      if (siteNum % wordSize == 0)
        cout << endl;
      if (siteNum % numSites == 0) {
        cout << "SampNum: " << sampNum << endl;
        ++sampNum;
      }
      cout << GLs[siteNum * 3] << ',' << GLs[siteNum * 3 + 1] << ','
           << GLs[siteNum * 3 + 2] << endl;
    }
    cout << endl;
  }
  */
  {
    gsl_rng_set(rng, 112);
    const unsigned numCycles3 = 100;
    ASSERT_EQ(3 * numSites * numSamps2, GLs.size());

    GLPackHelper::Init init(GLs, *rng);
    init.numSamps = numSamps2;
    init.sampleStride = sampleStride;
    init.useVQ = true;
    init.numBitsPerGL = BITSPERCODE;

    GLPack glPack3(init);
    HMMLike hmmLike3(hapPanel, bigNumHaps, glPack3, numCycles3, tran, &mutMat,
                     sampler2, *rng);

    unsigned firstSampIdx = 0;
    unsigned lastSampIdx = 0;
    {

      vector<unsigned> hapIdxs3 =
          hmmLike3.RunHMMOnSamples(firstSampIdx, lastSampIdx);
      ASSERT_EQ(0, firstSampIdx);
      ASSERT_EQ(1, lastSampIdx);
      ASSERT_EQ(numSamps * 4, hapIdxs3.size());

      // both of father and mother pairs need to be correct
      for (unsigned i = 0; i < 4; i += 2) {
        unsigned hap1 = hapIdxs3[i];
        unsigned hap2 = hapIdxs3[i + 1];
        if (hap1 > hap2)
          swap(hap1, hap2);
        EXPECT_EQ(10, hap2);
        EXPECT_LE(5, hap1);
        EXPECT_GE(6, hap1);

        // these are haps for second sample
        unsigned hap3 = hapIdxs3[i + 4];
        unsigned hap4 = hapIdxs3[i + 1 + 4];
        if (hap3 > hap4)
          swap(hap3, hap4);
        EXPECT_EQ(11, hap4);
        EXPECT_LE(7, hap3);
        EXPECT_GE(8, hap3);
      }
    }
    {

      // debugging

      vector<unsigned> haps = { 5, 6, 7, 8, 10, 11 };
      for (auto hapNum : haps) {
        cout << endl << "HapNum: " << hapNum << endl;
        for (unsigned wordNum = 0; wordNum < numWords; ++wordNum) {
          cout << hapPanel.at(numWords * hapNum + wordNum) << ' ';
        }
      }
      cout << endl;

      // printing out codebook
      auto codeBook = glPack3.GetCodeBook();
      cout << "CodeBook:" << endl;

      for (unsigned i = 0; i < codeBook.size(); ++i)
        cout << i << ": " << codeBook[i].first << ", " << codeBook[i].second
             << endl;

      // printing out VQ GLs
      cout << "VQ packed GLS:\n";
      auto pGLs = glPack3.GetPackedGLs();
      for (unsigned i = 0; i < pGLs.size(); ++i) {
        if (i % numWords == 0)
          cout << "Word: " << i / numWords << "\t\t";
        cout << std::hex << (pGLs.at(i) >> 28) << std::dec << "\t";
        if (i % numWords == 4) {
          cout << "..." << endl;
          i += numWords - 5;
        }
      }
      cout << endl;

      // pulling out GLs again to cycle back round to the right set
      glPack3.GetPackedGLs();
    }
    {
      gsl_rng_set(rng, 114);

      // this has been moved to hmmLike constructor
      //      HMMLikeCUDA::SetUpRNGs(glPack3.GetSampleStride(),
      // gsl_rng_get(rng));
      vector<unsigned> hapIdxs3 =
          hmmLike3.RunHMMOnSamples(firstSampIdx, lastSampIdx);
      ASSERT_EQ(2, firstSampIdx);
      ASSERT_EQ(3, lastSampIdx);
      ASSERT_EQ(numSamps * 4, hapIdxs3.size());

      cout << "Hap indices:\n";
      for (auto h : hapIdxs3)
        cout << h << " ";
      cout << endl << endl;

      // both of father and mother pairs need to be correct
      for (unsigned i = 0; i < 4; i += 2) {
        unsigned hap1 = hapIdxs3[i];
        unsigned hap2 = hapIdxs3[i + 1];
        if (hap1 > hap2)
          swap(hap1, hap2);
        EXPECT_EQ(11, hap2);
        EXPECT_LE(5, hap1);
        EXPECT_GE(6, hap1);

        // these are haps for second sample
        unsigned hap3 = hapIdxs3[i + 4];
        unsigned hap4 = hapIdxs3[i + 1 + 4];
        if (hap3 > hap4)
          swap(hap3, hap4);
        EXPECT_EQ(10, hap4);
        EXPECT_LE(7, hap3);
        EXPECT_GE(8, hap3);
      }
    }
  }
}
