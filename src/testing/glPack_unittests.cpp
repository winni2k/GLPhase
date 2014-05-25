
#include "gtest/gtest.h"
#include "glPack.hpp"
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <utility>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

// function class
class MatchPair {
public:
  MatchPair(float x, float y) : RR(x), Het(y) {}
  bool operator()(const std::pair<float, float> &entry) {
    return abs(RR - entry.first) < 0.01 && abs(Het - entry.second) < 0.01;
  };

private:
  float RR, Het;
};

TEST(GLPack, packGLsOK) {

  unsigned numSamps = 5;
  unsigned stride = 3;
  vector<float> GLs{ 0,        1.0 / 16, 15.0 / 16, // samp1
                     2.0 / 16, 3.0 / 16, 11.0 / 16, // samp2
                     4.0 / 16, 5.0 / 16, 7.0 / 16,  // samp3
                     0,        0,        1,         // samp4
                     0,        0,        1          // samp5
  };

  GLPackHelper::Init init(GLs, *rng);
  init.numSamps = numSamps;
  init.sampleStride = stride;
  GLPack glPack(init);

  for (int i = 0; i < 2; ++i) {
    ASSERT_EQ(0, glPack.GetNextSampIdx());
    vector<uint32_t> packedGLs = glPack.GetPackedGLs();
    ASSERT_EQ(stride - 1, glPack.GetLastSampIdx());
    ASSERT_EQ(stride, packedGLs.size());

    // make sure gls got packed correctly
    const unsigned char mask = (1 << CHAR_BIT) - 1;
    ASSERT_EQ(3, glPack.GetNextSampIdx());
    EXPECT_EQ(0, static_cast<unsigned>(
                     mask & (packedGLs[0] >>
                             (UINT32T_SIZE - glPack.GetNumBitsPerGL()))));
    EXPECT_EQ(18, static_cast<unsigned>(
                      mask & (packedGLs[1] >>
                              (UINT32T_SIZE - glPack.GetNumBitsPerGL()))));
    EXPECT_EQ(52, static_cast<unsigned>(
                      mask & (packedGLs[2] >>
                              (UINT32T_SIZE - glPack.GetNumBitsPerGL()))));

    vector<uint32_t> packedGLs2 = glPack.GetPackedGLs();
    ASSERT_EQ(numSamps - 1, glPack.GetLastSampIdx());
    ASSERT_EQ(2, packedGLs2.size());
  }
}
TEST(GLPack, packGLsMultiSiteOK) {
  gsl_rng_set(rng, 111);
  unsigned numSamps = 3;
  unsigned stride = 2;
  unsigned numSites = 2;
  vector<float> GLs{
    0,        1.0 / 16, 15.0 / 16, 2.0 / 16, 3.0 / 16, 11.0 / 16, // sample1
    4.0 / 16, 5.0 / 16, 7.0 / 16,  0.9 / 16, 3.5 / 16, 11.6 / 16, // sample2
    0,        0,        1,         0.5,      0.5,      0          // sample3
  };

  GLPackHelper::Init init(GLs, *rng);
  init.numSamps = numSamps;
  init.sampleStride = stride;
  GLPack glPack(init);

  ASSERT_EQ(numSites, glPack.GetNumSites());
  const unsigned char mask = (1 << CHAR_BIT) - 1;
  ASSERT_EQ(8, glPack.GetNumBitsPerGL());
  for (int i = 0; i < 2; ++i) {
    ASSERT_EQ(0, glPack.GetNextSampIdx());
    vector<uint32_t> packedGLs = glPack.GetPackedGLs();
    ASSERT_EQ(stride - 1, glPack.GetLastSampIdx());
    ASSERT_EQ(stride, packedGLs.size());

    // make sure gls got packed correctly
    ASSERT_EQ(stride, glPack.GetNextSampIdx());
    EXPECT_EQ(0, static_cast<unsigned>(
                     mask & (packedGLs.at(0) >>
                             (UINT32T_SIZE - glPack.GetNumBitsPerGL()))));
    EXPECT_EQ(52, static_cast<unsigned>(
                      mask & (packedGLs.at(1) >>
                              (UINT32T_SIZE - glPack.GetNumBitsPerGL()))));
    EXPECT_EQ(18, static_cast<unsigned>(
                      mask & (packedGLs.at(0) >>
                              (UINT32T_SIZE - 2 * glPack.GetNumBitsPerGL()))));
    EXPECT_EQ(3, static_cast<unsigned>(
                     mask & (packedGLs.at(1) >>
                             (UINT32T_SIZE - 2 * glPack.GetNumBitsPerGL()))));

    vector<uint32_t> packedGLs2 = glPack.GetPackedGLs();
    ASSERT_EQ(numSamps - 1, glPack.GetLastSampIdx());
    ASSERT_EQ(1, packedGLs2.size());
    EXPECT_EQ(0, static_cast<unsigned>(
                     mask & (packedGLs2.at(0) >>
                             (UINT32T_SIZE - glPack.GetNumBitsPerGL()))));
    EXPECT_EQ(119, static_cast<unsigned>(
                       mask & (packedGLs2.at(0) >>
                               (UINT32T_SIZE - 2 * glPack.GetNumBitsPerGL()))));
  }
}

TEST(GLVQHelper, EuclidianDist) {

  GLVQHelper::EuclidianDist eDist(0.5, 1);
  const std::pair<float, float> point1(0, 0);
  const std::pair<float, float> point2(2, 2);
  ASSERT_TRUE(eDist(point1, point2));
  ASSERT_FALSE(eDist(point2, point1));
}

TEST(GLVQ, FindGLCode) {
  unsigned numSamps = 5;
  vector<float> GLs{ 0,        1.0 / 16, 15.0 / 16, // samp1
                     2.0 / 16, 3.0 / 16, 11.0 / 16, // samp2
                     4.0 / 16, 5.0 / 16, 7.0 / 16,  // samp3
                     0,        0,        1,         // samp4
                     0,        0,        1          // samp5
  };
  const size_t codeBookSize = 16;

  GLVQ vq(GLs, *rng, codeBookSize);

  auto codeBook = vq.GetCodeBook();
  cout << "CodeBook:" << endl;
  unsigned codeNum = 0;
  for (auto c : codeBook)
    cout << codeNum++ << "\t" << c.first << "\t" << c.second << endl;
  // make sure each set of GLs has a cluster
  // and that the cluster found is the correct one
  cout << "samp num : GLs RR, Het : code : code RR, Het" << endl;
  for (unsigned sampNum = 0; sampNum < numSamps; ++sampNum) {
    const size_t GLidx = sampNum * 3;
    cout << sampNum << " : " << GLs[GLidx] << ", " << GLs[GLidx + 1] << " : ";
    MatchPair samp(GLs[GLidx], GLs[GLidx + 1]);
    EXPECT_EQ(true, std::any_of(codeBook.begin(), codeBook.end(), samp));

    auto sampIter = std::find_if(codeBook.begin(), codeBook.end(), samp);
    ASSERT_NE(codeBook.end(), sampIter);
    const unsigned char handGLCode = std::distance(codeBook.begin(), sampIter);

    const unsigned char glCode = vq.FindGLCode(GLs[GLidx], GLs[GLidx + 1]);
    ASSERT_EQ(handGLCode, glCode);
    cout << static_cast<unsigned>(glCode) << " : " << codeBook[glCode].first
         << ", " << codeBook[glCode].second << endl;
  }
}

TEST(GLPack, packGLsVQ) {
  unsigned numSamps = 5;
  unsigned stride = 3;
  vector<float> GLs{ 0,        1.0 / 16, 15.0 / 16, // samp1
                     2.0 / 16, 3.0 / 16, 11.0 / 16, // samp2
                     4.0 / 16, 5.0 / 16, 7.0 / 16,  // samp3
                     0,        0,        1,         // samp4
                     0,        0,        1          // samp5
  };

  GLPackHelper::Init init(GLs, *rng);
  init.numSamps = numSamps;
  init.sampleStride = stride;
  init.useVQ = true;
  init.numBitsPerGL = 4;
  //  init.VQcoolingRate =
  GLPack glPack(init);

  auto codeBook = glPack.GetCodeBook();
  cout << "CodeBook:" << endl;
  unsigned codeNum = 0;
  for (auto c : codeBook)
    cout << codeNum++ << "\t" << c.first << "\t" << c.second << endl;

  for (int i = 0; i < 2; ++i) {
    ASSERT_EQ(0, glPack.GetNextSampIdx());
    vector<uint32_t> packedGLs = glPack.GetPackedGLs();
    ASSERT_EQ(stride - 1, glPack.GetLastSampIdx());
    ASSERT_EQ(stride, packedGLs.size());

    auto codeBook = glPack.GetCodeBook();
    ASSERT_EQ(exp2(init.numBitsPerGL), codeBook.size());

    // make sure gls got packed correctly
    ASSERT_EQ(3, glPack.GetNextSampIdx());

    cout << "samp num : GLs RR, Het : code : code RR, Het" << endl;
    const unsigned siteNum = 0;
    const unsigned numSites = 1;
    for (unsigned sampNum = 0; sampNum < numSamps; ++sampNum) {

      const size_t GLidx = sampNum * 3 * numSites + siteNum * 3;
      cout << sampNum << " : " << GLs[GLidx] << ", " << GLs[GLidx + 1] << endl;

      MatchPair samp(GLs[GLidx], GLs[GLidx + 1]);
      EXPECT_EQ(true, std::any_of(codeBook.begin(), codeBook.end(), samp));
    }

    vector<uint32_t> packedGLs2 = glPack.GetPackedGLs();
    ASSERT_EQ(numSamps - 1, glPack.GetLastSampIdx());
    ASSERT_EQ(2, packedGLs2.size());
  }
}

//  NEED TO MAKE TEST SET WITH LOTS OF SITES AND MORE THAN ONE SAMPLE
TEST(GLPack, packGLsVQmultisite) {
  unsigned numSamps = 5;
  unsigned stride = 5;
  size_t numSites = 1024;
  vector<float> GLs;
  GLs.resize(numSamps * numSites * 3, 0);
  for (size_t sampIdx = 0; sampIdx < numSamps; ++sampIdx)
    for (size_t i = 0; i < numSites; ++i)
        GLs[ 3 * i + numSites * 3 * sampIdx] = i + sampIdx * numSites;

  GLPackHelper::Init init(GLs, *rng);
  init.numSamps = numSamps;
  init.sampleStride = stride;
  init.useVQ = true;
  init.numBitsPerGL = BITSPERCODE;
  //  init.VQcoolingRate =
  GLPack glPack(init);

  auto codeBook = glPack.GetCodeBook();
  cout << "CodeBook:" << endl;
  unsigned codeNum = 0;
  for (auto c : codeBook)
    cout << codeNum++ << "\t" << c.first << "\t" << c.second << endl;

  for (int i = 0; i < 2; ++i) {
    ASSERT_EQ(0, glPack.GetNextSampIdx());
    vector<uint32_t> packedGLs = glPack.GetPackedGLs();
    ASSERT_EQ(numSamps * numSites * BITSPERCODE / UINT32T_SIZE,
              packedGLs.size());
  }
}

TEST(GLVQ, FindGLCodeMultisite) {
  gsl_rng_set(rng, 111);
  unsigned numSamps = 3;
  unsigned numSites = 2;
  vector<float> GLs{
    0,        1.0 / 16, 15.0 / 16, 2.0 / 16, 3.0 / 16, 11.0 / 16, // sample1
    4.0 / 16, 5.0 / 16, 7.0 / 16,  0.9 / 16, 3.5 / 16, 11.6 / 16, // sample2
    0,        0,        1,         0.5,      0.5,      0          // sample3
  };
  const size_t codeBookSize = 16;
  GLVQ vq(GLs, *rng, codeBookSize);
  auto codeBook = vq.GetCodeBook();
  cout << "CodeBook:" << endl;
  unsigned codeNum = 0;
  for (auto c : codeBook)
    cout << codeNum++ << "\t" << c.first << "\t" << c.second << endl;
  // make sure each set of GLs has a cluster
  // and that the cluster found is the correct one
  cout << "samp num : GLs RR, Het : code : code RR, Het" << endl;
  for (unsigned sampNum = 0; sampNum < numSamps; ++sampNum) {
    for (size_t siteNum = 0; siteNum < numSites; ++siteNum) {
      const size_t GLidx = sampNum * 3 * numSites + siteNum * 3;
      cout << sampNum << " : " << GLs[GLidx] << ", " << GLs[GLidx + 1] << " : ";

      MatchPair samp(GLs[GLidx], GLs[GLidx + 1]);
      EXPECT_EQ(true, std::any_of(codeBook.begin(), codeBook.end(), samp));

      auto sampIter = std::find_if(codeBook.begin(), codeBook.end(), samp);
      ASSERT_NE(codeBook.end(), sampIter);
      const unsigned char handGLCode =
          std::distance(codeBook.begin(), sampIter);

      const unsigned char glCode = vq.FindGLCode(GLs[GLidx], GLs[GLidx + 1]);
      ASSERT_EQ(handGLCode, glCode);
      cout << static_cast<unsigned>(glCode) << " : " << codeBook[glCode].first
           << ", " << codeBook[glCode].second << endl;
    }
  }
}

TEST(GLVQ, FindGLCodeMultisiteOneSamp) {
  gsl_rng_set(rng, 111);
  const unsigned numSamps = 1;
  const unsigned numSites = 12;
  vector<float> GLs{ 0,         1.0 / 16,  15.0 / 16, 2.0 / 16,  3.0 / 16,
                     11.0 / 16, 4.0 / 16,  5.0 / 16,  7.0 / 16,  0.9 / 16,
                     3.5 / 16,  11.6 / 16, 0,         0,         1,
                     0.5,       0.5,       0,         0,         1.0 / 16,
                     15.0 / 16, 2.0 / 16,  3.0 / 16,  11.0 / 16, 4.0 / 16,
                     5.0 / 16,  7.0 / 16,  0.9 / 16,  3.5 / 16,  11.6 / 16,
                     0,         0,         1,         0.5,       0.5,
                     0 };
  // add some random points
  for (int i = 0; i < 20; ++i) {
    const size_t siteNum = gsl_rng_uniform_int(rng, numSites);
    float RR = abs(GLs[siteNum * 3] + (gsl_rng_uniform(rng) - 0.5) * 0.05);
    float Het = abs(GLs[siteNum * 3 + 1] + (gsl_rng_uniform(rng) - 0.5) * 0.05);
    GLs.push_back(RR);
    GLs.push_back(Het);
    GLs.push_back(max(0.0f, 1 - RR - Het));
  }
  const size_t codeBookSize = 16;
  GLVQ vq(GLs, *rng, codeBookSize);
  auto codeBook = vq.GetCodeBook();
  cout << "CodeBook:" << endl;
  unsigned codeNum = 0;
  for (auto c : codeBook)
    cout << codeNum++ << "\t" << c.first << "\t" << c.second << endl;
  // make sure each set of GLs has a cluster
  // and that the cluster found is the correct one
  cout << "samp num : GLs RR, Het : code : code RR, Het" << endl;
  const unsigned sampNum = numSamps - 1;
  for (size_t siteNum = 0; siteNum < numSites; ++siteNum) {
    const size_t GLidx = sampNum * 3 * numSites + siteNum * 3;
    cout << sampNum << " : " << GLs[GLidx] << ", " << GLs[GLidx + 1] << " : ";

    MatchPair samp(GLs[GLidx], GLs[GLidx + 1]);
    EXPECT_EQ(true, std::any_of(codeBook.begin(), codeBook.end(), samp));

    auto sampIter = std::find_if(codeBook.begin(), codeBook.end(), samp);
    ASSERT_NE(codeBook.end(), sampIter);
    const unsigned char handGLCode = std::distance(codeBook.begin(), sampIter);

    const unsigned char glCode = vq.FindGLCode(GLs[GLidx], GLs[GLidx + 1]);
    ASSERT_EQ(handGLCode, glCode);
    cout << static_cast<unsigned>(glCode) << " : " << codeBook[glCode].first
         << ", " << codeBook[glCode].second << endl;
  }
}
