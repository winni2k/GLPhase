
#include "gtest/gtest.h"
#include "haplotype.hpp"
#include "kNN.hpp"
#include <algorithm>
#include <gsl/gsl_rng.h>

using namespace std;

string sampleDir = "../../samples";
string brokenDir = sampleDir + "/brokenFiles";
string sampleLegend =
    sampleDir + "/20_011976121_012173018.bin.onlyThree.legend";
string sampleHap = sampleDir + "/20_011976121_012173018.bin.onlyThree.hap";
string sampleBin = sampleDir + "/20_011976121_012173018.bin.onlyThree.bin";
string refHap =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.hap";
string refLegend =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.legend";
string refHaps =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.haps";

string scaffHapLegSampSample = sampleDir + "/onlyThree.hapLegSamp.sample";
string scaffHapsSampSample = sampleDir + "/onlyThree.hapsSample.sample";
string scaffHapsSampUnorderedSample =
    sampleDir + "/onlyThree.hapsSample.unordered.sample";

string scaffoldUnorderedHaps =
    sampleDir +
    "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.unordered.haps";
string scaffoldHap =
    sampleDir + "/20_011976121_012173018.bin.onlyThree.scaffold50.hap";
string scaffoldHaps =
    sampleDir + "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.haps";
string scaffoldHapsWrongReg =
    sampleDir +
    "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.wrongRegion.haps";
string scaffoldLegend =
    sampleDir + "/20_011976121_012173018.bin.onlyThree.scaffold50.legend";

string brokenHapLegSampSample =
    brokenDir + "/onlyThree.hapLegSample.extraLine.sample";
string brokenHapsSampSample =
    brokenDir + "/onlyThree.hapsSample.extraLine.sample";
string unsortedRefHaps =
    sampleDir +
    "/20_0_62000000.011976121_012173018.paste.onlyThree.unsorted.haps";

gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

TEST(Haplotype, StoresOK) {

  // testing to see if init and testing works ok
  Haplotype simpleA(4);

  for (unsigned i = 0; i < 4; i++)
    EXPECT_FALSE(simpleA.TestSite(i));

  EXPECT_DEATH(simpleA.TestSite(4), "uSite < m_uNumAlleles");

  Haplotype simpleB(4);
  simpleB.Set(0, 1);
  simpleB.Set(3, 1);
  EXPECT_TRUE(simpleB.TestSite(0));
  EXPECT_TRUE(simpleB.TestSite(3));
  EXPECT_FALSE(simpleB.TestSite(2));

  // test hamming distance
  EXPECT_EQ(2, simpleA.HammingDist(simpleB));
  EXPECT_EQ(0, simpleA.HammingDist(simpleA));

  Haplotype longA(128);
  Haplotype longB(128);

  EXPECT_DEATH(simpleA.Set(128, 1), "uSite < m_uNumAlleles");
  EXPECT_DEATH(simpleA.TestSite(128), "uSite < m_uNumAlleles");

  longA.Set(127, 1);
  longA.Set(1, 1);
  EXPECT_EQ(2, longA.HammingDist(longB));
  EXPECT_EQ(2, longB.HammingDist(longA));
  longB.Set(120, 1);
  EXPECT_EQ(3, longB.HammingDist(longA));

  vector<uint64_t> hapWords;
  hapWords.push_back(longA.GetWord(0));
  ASSERT_TRUE(longB.TestSite(1, hapWords.data()));
}

TEST(KNN, clustersOK) {

  gsl_rng_set(rng, time(NULL));
  std::srand(1);

  // testing to see if nearest neighbor clustering works ok
  unsigned numClusters = 4;
  unsigned numHaps = 16;
  unsigned numSites = numHaps + numClusters * 2;

  // create a set of test haplotypes to cluster and sample from
  vector<Haplotype> haplotypes;
  vector<unsigned> shuffledIndexes(numHaps);

  for (unsigned i = 0; i < numHaps; i++) {
    shuffledIndexes[i] = i;
    Haplotype temp(numSites);

    for (unsigned j = 0; j < numClusters + 2; j++)
      temp.Set(i + j, 1);

    haplotypes.push_back(temp);
  }

  // make the second haplotype match hap number 16 and 15 (0 based) closest
  // but not as well as hap 0 would match hap 4
  for (unsigned i = numSites - 6; i < numSites; i++) {
    haplotypes[1].Set(i - numSites + 7, false);
    haplotypes[1].Set(i, true);
  }

  // shuffle the haplotypes except for the first two
  std::random_shuffle(shuffledIndexes.begin() + 2, shuffledIndexes.end());

  vector<uint64_t> passHaps;

  for (unsigned j = 0; j < shuffledIndexes.size(); j++)
    passHaps.push_back(haplotypes[shuffledIndexes[j]].GetWord(0));

/*  for (unsigned idx = 0; idx != haplotypes.size(); idx++) {
    cerr << idx << ": ";

    for (unsigned i = 0; i < numSites; i++)
      cerr << haplotypes[shuffledIndexes[idx]].TestSite(i);

    cerr << endl;
    }*/

  KNN kNN(numClusters);
  kNN.init(passHaps, 1, numSites, 0);

  // check to make sure kNN has the haps stored correctly
  vector<unsigned> neighborHapNums = kNN.Neighbors(0);

  EXPECT_EQ(numClusters, neighborHapNums.size());
  EXPECT_EQ(2, shuffledIndexes[neighborHapNums[0]]);
  EXPECT_EQ(15, shuffledIndexes[neighborHapNums[1]]);
  EXPECT_EQ(3, shuffledIndexes[neighborHapNums[2]]);
  EXPECT_EQ(14, shuffledIndexes[neighborHapNums[3]]);

  // testing sampling
  for (unsigned i = 0; i < 10; i++) {
    unsigned sampHap = kNN.SampleHap(0, rng);
    EXPECT_LT(shuffledIndexes[sampHap], numHaps);
    EXPECT_GT(shuffledIndexes[sampHap], 1);

    for (unsigned i = 4; i < numHaps - 2; i++)
      EXPECT_NE(shuffledIndexes[sampHap], i);
  }

  // undo change of hap 1
  for (unsigned i = numSites - 6; i < numSites; i++) {
    haplotypes[1].Set(i - numSites + 7, true);
    haplotypes[1].Set(i, false);
  }

  // now test the thresholding option
  // first site above threshold
  haplotypes[1].Set(0, true);
  haplotypes[5].Set(0, true);
  haplotypes[6].Set(0, true);
  haplotypes[8].Set(0, true);
  haplotypes[9].Set(0, true);
  haplotypes[10].Set(0, true);

  // second site above threshold
  haplotypes[5].Set(2, true);
  haplotypes[6].Set(2, true);
  haplotypes[8].Set(2, true);
  haplotypes[11].Set(2, true);
  haplotypes[12].Set(2, true);

  haplotypes[0].Set(numSites - 1, true);
  haplotypes[1].Set(numSites - 1, true);
  haplotypes[5].Set(numSites - 1, true);
  haplotypes[7].Set(numSites - 1, true);
  haplotypes[15].Set(numSites - 1, true);
  haplotypes[14].Set(numSites - 1, true);
  haplotypes[13].Set(numSites - 1, true);

  haplotypes[0].Set(numSites - 2, true);
  haplotypes[1].Set(numSites - 2, true);
  haplotypes[2].Set(numSites - 2, true);
  haplotypes[3].Set(numSites - 2, true);
  haplotypes[4].Set(numSites - 2, true);
  haplotypes[5].Set(numSites - 2, true);
  haplotypes[8].Set(numSites - 2, true);

  /*    for(auto hap : haplotypes){
          for(unsigned i = 0; i < numSites; i++)
              cerr << hap.TestSite(i);
          cerr << endl;
          }*/

  passHaps.clear();

  for (unsigned j = 0; j < shuffledIndexes.size(); j++)
    passHaps.push_back(haplotypes[shuffledIndexes[j]].GetWord(0));

  KNN kNN2(numClusters);
  kNN2.init(passHaps, 1, numSites, 0.4375);

  // check to make sure variant allele freqs are calculated correctly
  vector<double> varAfs;
  kNN2.VarAfs(varAfs);
  EXPECT_EQ(numSites, varAfs.size());
  EXPECT_EQ(0.4375, varAfs[0]);
  EXPECT_EQ(0.125, varAfs[1]);
  EXPECT_EQ(0.5, varAfs[2]);

  // check to make sure kNN is thresholding the correct sites
  vector<unsigned> commonSites;
  kNN2.ClusterSites(commonSites);
  EXPECT_EQ(4, commonSites.size());

  // check to make sure kNN has the haps stored correctly
  neighborHapNums = kNN2.Neighbors(0);

  EXPECT_EQ(numClusters, neighborHapNums.size());
  EXPECT_EQ(5, shuffledIndexes[neighborHapNums[0]]);
  EXPECT_EQ(6, shuffledIndexes[neighborHapNums[1]]);
  EXPECT_EQ(8, shuffledIndexes[neighborHapNums[2]]);
  EXPECT_EQ(2, shuffledIndexes[neighborHapNums[3]]);

  // testing sampling
  for (unsigned i = 0; i < 10; i++) {
    unsigned sampHap = kNN2.SampleHap(0, rng);
    EXPECT_LT(shuffledIndexes[sampHap], 9);
    EXPECT_GT(shuffledIndexes[sampHap], 1);
    EXPECT_NE(shuffledIndexes[sampHap], 3);
    EXPECT_NE(shuffledIndexes[sampHap], 4);
    EXPECT_NE(shuffledIndexes[sampHap], 7);
  }
}
