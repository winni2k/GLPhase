
#include "gtest/gtest.h"
#include "insti.hpp"
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
string refVCFGZ =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.vcf.gz";

string scaffHapLegSampSample = sampleDir + "/onlyThree.hapLegSamp.sample";
string scaffHapsSampSample = sampleDir + "/onlyThree.hapsSample.sample";
string scaffHapsSampUnorderedSample =
    sampleDir + "/onlyThree.hapsSample.unordered.sample";

string scaffoldUnorderedHaps =
    sampleDir +
    "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.unordered.haps";
string scaffoldUnorderedVCFGZ =
    sampleDir +
    "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.unordered.vcf.gz";
string scaffoldHap =
    sampleDir + "/20_011976121_012173018.bin.onlyThree.scaffold50.hap";
string scaffoldHaps =
    sampleDir + "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.haps";
string scaffoldTabHaps =
    sampleDir +
    "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.tabhaps.gz";
string scaffoldVCFGZ =
    sampleDir +
    "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.vcf.gz";
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

string geneticMap =
    sampleDir + "/geneticMap/genetic_map_chr20_combined_b37.txt";

gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

TEST(Insti, loadBin) {

  InstiHelper::Init init;
  init.geneticMap = geneticMap;
  Insti lp(init);
  lp.load_bin(sampleBin);

  ASSERT_EQ(3, lp.in);

  // testing sites
  EXPECT_EQ(1024, lp.m_sites.size());

  // chr
  EXPECT_EQ("20", lp.m_sites.chrom());
  EXPECT_EQ("20", lp.m_sites.at(0)->chr);
  EXPECT_EQ("20", lp.m_sites.at(5)->chr);
  EXPECT_EQ("20", lp.m_sites.at(1023)->chr);
  EXPECT_NE("16", lp.m_sites.at(1023)->chr);

  // posi
  EXPECT_EQ(11976121, lp.m_sites.at(0)->pos);
  EXPECT_EQ(11977230, lp.m_sites.at(5)->pos);
  EXPECT_EQ(12173018, lp.m_sites.at(1023)->pos);

  // all
  EXPECT_EQ("T", lp.m_sites.at(0)->ref);
  EXPECT_EQ("C", lp.m_sites.at(0)->alt);
  EXPECT_EQ("A", lp.m_sites.at(5)->ref);
  EXPECT_EQ("G", lp.m_sites.at(5)->alt);
  EXPECT_EQ("G", lp.m_sites.at(1023)->ref);
  EXPECT_EQ("T", lp.m_sites.at(1023)->alt);

  // prob
  // making sure prob size is correct
  EXPECT_EQ(1024 * 2 * lp.in, lp.prob.size());

  EXPECT_EQ(0, lp.prob[0]);
  EXPECT_EQ(1, lp.prob[1]);
  EXPECT_EQ(0.2f, lp.prob[2]);
  EXPECT_EQ(0, lp.prob[3]);
  EXPECT_EQ(0, lp.prob[4]);
  EXPECT_EQ(0, lp.prob[5]);
  EXPECT_EQ(0, lp.prob[6]);
  EXPECT_EQ(1, lp.prob[7]);
  EXPECT_EQ(1, lp.prob[lp.in * 6 - 4]);

  //    cerr << "BLOINC1\n";
  // now initialize lp and see if probs still make sense
  lp.initialize();

  EXPECT_EQ(1024 * 3 * lp.in, lp.prob.size());
  EXPECT_EQ(0, lp.prob[0]);
  EXPECT_EQ(0, lp.prob[1]);
  EXPECT_EQ(1, lp.prob[2]);

  // 3 = lp.pn
  EXPECT_EQ(0.8f, lp.prob[0 + lp.mn * 3]);
  EXPECT_EQ(0.2f, lp.prob[1 + lp.mn * 3]);
  EXPECT_EQ(0, lp.prob[2 + lp.mn * 3]);

  EXPECT_EQ(1, lp.prob[0 + lp.mn * 2 * 3]);
  EXPECT_EQ(0, lp.prob[1 + lp.mn * 2 * 3]);
  EXPECT_EQ(0, lp.prob[2 + lp.mn * 2 * 3]);

  EXPECT_EQ(0, lp.prob[4]);
  EXPECT_EQ(1, lp.prob[5]);

  //  cerr << "BLOINC2\n";
  // now test refpanel loading
  lp.LoadHapLegSamp(refLegend, refHap, "", InstiPanelType::REFERENCE);

  for (unsigned i = 0; i != 601; i++) {
    EXPECT_EQ(0, lp.TestRefHap(0, i));
    EXPECT_EQ(0, lp.TestRefHap(2, i));
    EXPECT_EQ(1, lp.TestRefHap(3, i));
  }

  for (unsigned i = 601; i != 1024; i++) {
    EXPECT_EQ(1, lp.TestRefHap(0, i));
    EXPECT_EQ(0, lp.TestRefHap(1, i));
    EXPECT_EQ(0, lp.TestRefHap(2, i));
    EXPECT_EQ(1, lp.TestRefHap(3, i));
  }

  // test the scaffold loading not implemented yet
  //    lp.LoadHapLegSamp(refLegend, refHap, scaffHapLegSampSample,
  // InstiPanelType::SCAFFOLD);
}

TEST(Insti, loadHapsSamp) {

  InstiHelper::Init init;
  init.geneticMap = geneticMap;
  Insti lp(init);
  lp.load_bin(sampleBin);
  lp.initialize();
  const string nothing = "";
  lp.LoadHapsSamp(refHaps, nothing, InstiPanelType::REFERENCE);

  for (unsigned i = 0; i != 601; i++) {
    EXPECT_EQ(0, lp.TestRefHap(0, i));
    EXPECT_EQ(0, lp.TestRefHap(2, i));
    EXPECT_EQ(1, lp.TestRefHap(3, i));
  }

  for (unsigned i = 601; i != 1024; i++) {
    EXPECT_EQ(1, lp.TestRefHap(0, i));
    EXPECT_EQ(0, lp.TestRefHap(1, i));
    EXPECT_EQ(0, lp.TestRefHap(2, i));
    EXPECT_EQ(1, lp.TestRefHap(3, i));
  }

  // test the scaffold loading
  lp.LoadHapsSamp(scaffoldHaps, scaffHapsSampSample, InstiPanelType::SCAFFOLD);
  ASSERT_EQ("samp1", lp.GetScaffoldID(0));
  ASSERT_EQ("samp2", lp.GetScaffoldID(1));
  ASSERT_EQ("samp3", lp.GetScaffoldID(2));

  ASSERT_EQ(6, lp.GetScaffoldNumHaps());
  ASSERT_EQ(1, lp.GetScaffoldNumWordsPerHap());
  ASSERT_EQ(50, lp.GetScaffoldNumSites());

  for (unsigned i = 0; i != 26; i++) {
    EXPECT_EQ(0, lp.TestScaffoldSite(0, i));
    EXPECT_EQ(1, lp.TestScaffoldSite(1, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(2, i));
    EXPECT_EQ(1, lp.TestScaffoldSite(3, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(4, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(5, i));
  }

  for (unsigned i = 26; i != 50; i++) {
    EXPECT_EQ(1, lp.TestScaffoldSite(0, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(1, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(2, i));
    EXPECT_EQ(1, lp.TestScaffoldSite(3, i));
    EXPECT_EQ(1, lp.TestScaffoldSite(4, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(5, i));
  }

  // unordered haps test
  Insti lp2(init);
  lp2.load_bin(sampleBin);
  lp2.initialize();

  // test the scaffold loading
  lp2.LoadHapsSamp(scaffoldHaps, scaffHapsSampUnorderedSample,
                   InstiPanelType::SCAFFOLD);
  ASSERT_EQ("samp1", lp2.GetScaffoldID(0));
  ASSERT_EQ("samp2", lp2.GetScaffoldID(1));
  ASSERT_EQ("samp3", lp2.GetScaffoldID(2));

  for (unsigned i = 0; i != 26; i++) {
    EXPECT_EQ(0, lp2.TestScaffoldSite(4, i));
    EXPECT_EQ(1, lp2.TestScaffoldSite(5, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(0, i));
    EXPECT_EQ(1, lp2.TestScaffoldSite(1, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(2, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(3, i));
  }

  for (unsigned i = 26; i != 50; i++) {
    EXPECT_EQ(1, lp2.TestScaffoldSite(4, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(5, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(0, i));
    EXPECT_EQ(1, lp2.TestScaffoldSite(1, i));
    EXPECT_EQ(1, lp2.TestScaffoldSite(2, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(3, i));
  }

  // try some quirky input
  Insti lp3(init);
  lp3.load_bin(sampleBin);
  lp3.initialize();
  lp3.LoadHapsSamp(refHaps, nothing, InstiPanelType::REFERENCE);

  // test the scaffold loading
  ASSERT_EXIT(lp3.LoadHapsSamp(scaffoldHapsWrongReg, scaffHapsSampSample,
                               InstiPanelType::SCAFFOLD),
              ::testing::ExitedWithCode(1),
              "Number of haplotypes in haps file is 0.  Haps file empty?");

  // test the scaffold loading using tabhaps
  // unordered haps test
  Insti lp4(init);
  lp4.load_bin(sampleBin);
  lp4.initialize();

  // test the scaffold loading
  lp4.LoadHapsSamp(scaffoldTabHaps, scaffHapsSampUnorderedSample,
                   InstiPanelType::SCAFFOLD);
  ASSERT_EQ("samp1", lp4.GetScaffoldID(0));
  ASSERT_EQ("samp2", lp4.GetScaffoldID(1));
  ASSERT_EQ("samp3", lp4.GetScaffoldID(2));

  for (unsigned i = 0; i != 26; i++) {
    EXPECT_EQ(0, lp4.TestScaffoldSite(4, i));
    EXPECT_EQ(1, lp4.TestScaffoldSite(5, i));
    EXPECT_EQ(0, lp4.TestScaffoldSite(0, i));
    EXPECT_EQ(1, lp4.TestScaffoldSite(1, i));
    EXPECT_EQ(0, lp4.TestScaffoldSite(2, i));
    EXPECT_EQ(0, lp4.TestScaffoldSite(3, i));
  }

  for (unsigned i = 26; i != 50; i++) {
    EXPECT_EQ(1, lp4.TestScaffoldSite(4, i));
    EXPECT_EQ(0, lp4.TestScaffoldSite(5, i));
    EXPECT_EQ(0, lp4.TestScaffoldSite(0, i));
    EXPECT_EQ(1, lp4.TestScaffoldSite(1, i));
    EXPECT_EQ(1, lp4.TestScaffoldSite(2, i));
    EXPECT_EQ(0, lp4.TestScaffoldSite(3, i));
  }
}

TEST(Insti, loadHapLegSampErrors) {

  InstiHelper::Init init;
  init.geneticMap = geneticMap;
  Insti lp(init);
  lp.load_bin(sampleBin);
  lp.initialize();

  //    cerr << "BLOINC1\n";
  ASSERT_EXIT(lp.LoadHapLegSamp("", sampleHap, "", InstiPanelType::REFERENCE),
              ::testing::ExitedWithCode(1),
              "Need to define a legend file if defining a hap file");
  ASSERT_EXIT(
      lp.LoadHapLegSamp(sampleLegend, "", "", InstiPanelType::REFERENCE),
      ::testing::ExitedWithCode(1),
      "Need to define a hap file if defining a legend file");
  EXPECT_ANY_THROW(
      lp.LoadHapsSamp(refHaps, brokenHapsSampSample, InstiPanelType::SCAFFOLD));

  /*
                  "Input haplotypes file " + unsortedRefHaps +
                  " needs to be sorted by position");
  */
  EXPECT_ANY_THROW(lp.LoadHapsSamp(
      unsortedRefHaps, scaffHapsSampUnorderedSample, InstiPanelType::SCAFFOLD));
}

TEST(Insti, initializingHapsFromScaffold) {
  InstiHelper::Init init;
  init.geneticMap = geneticMap;
  init.scaffoldHapsFile = scaffoldHaps;
  init.initPhaseFromScaffold = true;
  init.scaffoldSampleFile = scaffHapsSampSample;
  Insti lp(init);

  lp.load_bin(sampleBin);
  // test the scaffold loading
  lp.initialize();

  // test to see if main haps were initialized correctly
  vector<unsigned> siteIdxs = {6, 316, 576};
  for (auto i : siteIdxs) {
    EXPECT_EQ(0, lp.TestMainHap_(0, i));
    EXPECT_EQ(1, lp.TestMainHap_(1, i));
    EXPECT_EQ(0, lp.TestMainHap_(2, i));
    EXPECT_EQ(1, lp.TestMainHap_(3, i));
    EXPECT_EQ(0, lp.TestMainHap_(4, i));
    EXPECT_EQ(0, lp.TestMainHap_(5, i));
  }

  siteIdxs = {635, 739, 1021};
  for (auto i : siteIdxs) {
    EXPECT_EQ(1, lp.TestMainHap_(0, i));
    EXPECT_EQ(0, lp.TestMainHap_(1, i));
    EXPECT_EQ(0, lp.TestMainHap_(2, i));
    EXPECT_EQ(1, lp.TestMainHap_(3, i));
    EXPECT_EQ(1, lp.TestMainHap_(4, i));
    EXPECT_EQ(0, lp.TestMainHap_(5, i));
  }
}

TEST(Insti, LoadVCFGZ) {

  string region = "20:43281-61958963";

  InstiHelper::Init init;
  init.geneticMap = geneticMap;
  Insti lp(init);
  lp.load_bin(sampleBin);
  lp.initialize();
  lp.LoadVCFGZ(refVCFGZ, InstiPanelType::REFERENCE, region);

  for (unsigned i = 0; i != 601; i++) {
    EXPECT_EQ(0, lp.TestRefHap(0, i));
    EXPECT_EQ(0, lp.TestRefHap(2, i));
    EXPECT_EQ(1, lp.TestRefHap(3, i));
  }

  for (unsigned i = 601; i != 1024; i++) {
    EXPECT_EQ(1, lp.TestRefHap(0, i));
    EXPECT_EQ(0, lp.TestRefHap(1, i));
    EXPECT_EQ(0, lp.TestRefHap(2, i));
    EXPECT_EQ(1, lp.TestRefHap(3, i));
  }

  // test the scaffold loading
  lp.LoadVCFGZ(scaffoldVCFGZ, InstiPanelType::SCAFFOLD, region);
  ASSERT_EQ("samp1", lp.GetScaffoldID(0));
  ASSERT_EQ("samp2", lp.GetScaffoldID(1));
  ASSERT_EQ("samp3", lp.GetScaffoldID(2));

  ASSERT_EQ(6, lp.GetScaffoldNumHaps());
  ASSERT_EQ(1, lp.GetScaffoldNumWordsPerHap());
  ASSERT_EQ(50, lp.GetScaffoldNumSites());

  for (unsigned i = 0; i != 26; i++) {
    EXPECT_EQ(0, lp.TestScaffoldSite(0, i));
    EXPECT_EQ(1, lp.TestScaffoldSite(1, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(2, i));
    EXPECT_EQ(1, lp.TestScaffoldSite(3, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(4, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(5, i));
  }

  for (unsigned i = 26; i != 50; i++) {
    EXPECT_EQ(1, lp.TestScaffoldSite(0, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(1, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(2, i));
    EXPECT_EQ(1, lp.TestScaffoldSite(3, i));
    EXPECT_EQ(1, lp.TestScaffoldSite(4, i));
    EXPECT_EQ(0, lp.TestScaffoldSite(5, i));
  }

  // unordered haps test
  Insti lp2(init);
  lp2.load_bin(sampleBin);
  lp2.initialize();

  // test the scaffold loading
  lp2.LoadVCFGZ(scaffoldUnorderedVCFGZ, InstiPanelType::SCAFFOLD, region);
  ASSERT_EQ("samp1", lp2.GetScaffoldID(0));
  ASSERT_EQ("samp2", lp2.GetScaffoldID(1));
  ASSERT_EQ("samp3", lp2.GetScaffoldID(2));

  for (unsigned i = 0; i != 26; i++) {
    EXPECT_EQ(0, lp2.TestScaffoldSite(4, i));
    EXPECT_EQ(1, lp2.TestScaffoldSite(5, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(0, i));
    EXPECT_EQ(1, lp2.TestScaffoldSite(1, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(2, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(3, i));
  }

  for (unsigned i = 26; i != 50; i++) {
    EXPECT_EQ(1, lp2.TestScaffoldSite(4, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(5, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(0, i));
    EXPECT_EQ(1, lp2.TestScaffoldSite(1, i));
    EXPECT_EQ(1, lp2.TestScaffoldSite(2, i));
    EXPECT_EQ(0, lp2.TestScaffoldSite(3, i));
  }
}
