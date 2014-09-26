
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
string sampleVCF = sampleDir + "/20_011976121_012173018.bin.onlyThree.vcf.gz";
string sampleVCFExtraSites =
    sampleDir + "/20_011976121_012173018.bin.onlyThree.extraSites.vcf.gz";
string refHap =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.hap";
string refLegend =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.legend";
string refSamp =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.hap.indv";
string refSample =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.sample";
string refHaps =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.haps";
string refVCFGZ =
    sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.vcf.gz";
string chrom("20");
Bio::Region refRegion(chrom, 43281, 61958963);

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
  init.inputGLFile = sampleBin;
  Insti lp(init);

  ASSERT_EQ(3, lp.in);

  // testing sites
  EXPECT_EQ(1024, lp.m_glSites.size());

  // chr
  EXPECT_EQ("20", lp.m_glSites[0].chr);
  EXPECT_EQ("20", lp.m_glSites[5].chr);
  EXPECT_EQ("20", lp.m_glSites[1023].chr);
  EXPECT_NE("16", lp.m_glSites[1023].chr);

  // posi
  EXPECT_EQ(11976121, lp.m_glSites[0].pos);
  EXPECT_EQ(11977230, lp.m_glSites[5].pos);
  EXPECT_EQ(12173018, lp.m_glSites[1023].pos);

  // all
  EXPECT_EQ("T", lp.m_glSites[0].ref);
  EXPECT_EQ("A", lp.m_glSites[5].ref);
  EXPECT_EQ("G", lp.m_glSites[1023].ref);
  EXPECT_EQ("C", lp.m_glSites[0].alt);
  EXPECT_EQ("G", lp.m_glSites[5].alt);
  EXPECT_EQ("T", lp.m_glSites[1023].alt);

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
  lp.LoadHapLegSamp(refLegend, refHap, refSamp, InstiPanelType::REFERENCE,
                    Bio::Region());

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

TEST(Insti, loadGLVCF) {

  float absErr = 0.00001;
  InstiHelper::Init init;
  init.geneticMap = geneticMap;
  init.inputGLFile = sampleVCF;
  init.inputGLFileType = "bcf";
  Insti lp(init);

  ASSERT_EQ(3, lp.in);

  // testing sites
  EXPECT_EQ(1024, lp.m_glSites.size());

  // chr
  EXPECT_EQ("20", lp.m_glSites[0].chr);
  EXPECT_EQ("20", lp.m_glSites[5].chr);
  EXPECT_EQ("20", lp.m_glSites[1023].chr);
  EXPECT_NE("16", lp.m_glSites[1023].chr);

  // posi
  EXPECT_EQ(11976121, lp.m_glSites[0].pos);
  EXPECT_EQ(11977230, lp.m_glSites[5].pos);
  EXPECT_EQ(12173018, lp.m_glSites[1023].pos);

  // all
  EXPECT_EQ("T", lp.m_glSites[0].ref);
  EXPECT_EQ("A", lp.m_glSites[5].ref);
  EXPECT_EQ("G", lp.m_glSites[1023].ref);
  EXPECT_EQ("C", lp.m_glSites[0].alt);
  EXPECT_EQ("G", lp.m_glSites[5].alt);
  EXPECT_EQ("T", lp.m_glSites[1023].alt);

  // prob
  // making sure prob size is correct
  EXPECT_EQ(1024 * 2 * lp.in, lp.prob.size());

  EXPECT_NEAR(0, lp.prob[0], absErr);
  EXPECT_NEAR(1, lp.prob[1], absErr);
  EXPECT_NEAR(0.2f, lp.prob[2], absErr);
  EXPECT_NEAR(0, lp.prob[3], absErr);
  EXPECT_NEAR(0, lp.prob[4], absErr);
  EXPECT_NEAR(0, lp.prob[5], absErr);
  EXPECT_NEAR(0, lp.prob[6], absErr);
  EXPECT_NEAR(1, lp.prob[7], absErr);
  EXPECT_NEAR(1, lp.prob[lp.in * 6 - 4], absErr);

  //    cerr << "BLOINC1\n";
  // now initialize lp and see if probs still make sense
  lp.initialize();

  EXPECT_EQ(1024 * 3 * lp.in, lp.prob.size());
  EXPECT_NEAR(0, lp.prob[0], absErr);
  EXPECT_NEAR(0, lp.prob[1], absErr);
  EXPECT_NEAR(1, lp.prob[2], absErr);

  // 3 = lp.pn
  EXPECT_NEAR(0.8f, lp.prob[0 + lp.mn * 3], absErr);
  EXPECT_NEAR(0.2f, lp.prob[1 + lp.mn * 3], absErr);
  EXPECT_NEAR(0, lp.prob[2 + lp.mn * 3], absErr);

  EXPECT_NEAR(1, lp.prob[0 + lp.mn * 2 * 3], absErr);
  EXPECT_NEAR(0, lp.prob[1 + lp.mn * 2 * 3], absErr);
  EXPECT_NEAR(0, lp.prob[2 + lp.mn * 2 * 3], absErr);

  EXPECT_NEAR(0, lp.prob[4], absErr);
  EXPECT_NEAR(1, lp.prob[5], absErr);

  //  cerr << "BLOINC2\n";
  // now test refpanel loading
  lp.LoadHapLegSamp(refLegend, refHap, refSamp, InstiPanelType::REFERENCE,
                    Bio::Region());

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

TEST(Insti, loadGLVCFExtraSites) {

  float absErr = 0.00001;
  InstiHelper::Init init;
  init.geneticMap = geneticMap;
  init.inputGLFile = sampleVCF;
  init.inputGLFileType = "bcf";
  init.inputGLRegion = "20:11976121-12173018";
  Insti lp(init);

  ASSERT_EQ(3, lp.in);

  // testing sites
  EXPECT_EQ(1024, lp.m_glSites.size());

  // chr
  EXPECT_EQ("20", lp.m_glSites[0].chr);
  EXPECT_EQ("20", lp.m_glSites[5].chr);
  EXPECT_EQ("20", lp.m_glSites[1023].chr);
  EXPECT_NE("16", lp.m_glSites[1023].chr);

  // posi
  EXPECT_EQ(11976121, lp.m_glSites[0].pos);
  EXPECT_EQ(11977230, lp.m_glSites[5].pos);
  EXPECT_EQ(12173018, lp.m_glSites[1023].pos);

  // all
  EXPECT_EQ("T", lp.m_glSites[0].ref);
  EXPECT_EQ("A", lp.m_glSites[5].ref);
  EXPECT_EQ("G", lp.m_glSites[1023].ref);
  EXPECT_EQ("C", lp.m_glSites[0].alt);
  EXPECT_EQ("G", lp.m_glSites[5].alt);
  EXPECT_EQ("T", lp.m_glSites[1023].alt);

  // prob
  // making sure prob size is correct
  EXPECT_EQ(1024 * 2 * lp.in, lp.prob.size());

  EXPECT_NEAR(0, lp.prob[0], absErr);
  EXPECT_NEAR(1, lp.prob[1], absErr);
  EXPECT_NEAR(0.2f, lp.prob[2], absErr);
  EXPECT_NEAR(0, lp.prob[3], absErr);
  EXPECT_NEAR(0, lp.prob[4], absErr);
  EXPECT_NEAR(0, lp.prob[5], absErr);
  EXPECT_NEAR(0, lp.prob[6], absErr);
  EXPECT_NEAR(1, lp.prob[7], absErr);
  EXPECT_NEAR(1, lp.prob[lp.in * 6 - 4], absErr);

  //    cerr << "BLOINC1\n";
  // now initialize lp and see if probs still make sense
  lp.initialize();

  EXPECT_EQ(1024 * 3 * lp.in, lp.prob.size());
  EXPECT_NEAR(0, lp.prob[0], absErr);
  EXPECT_NEAR(0, lp.prob[1], absErr);
  EXPECT_NEAR(1, lp.prob[2], absErr);

  // 3 = lp.pn
  EXPECT_NEAR(0.8f, lp.prob[0 + lp.mn * 3], absErr);
  EXPECT_NEAR(0.2f, lp.prob[1 + lp.mn * 3], absErr);
  EXPECT_NEAR(0, lp.prob[2 + lp.mn * 3], absErr);

  EXPECT_NEAR(1, lp.prob[0 + lp.mn * 2 * 3], absErr);
  EXPECT_NEAR(0, lp.prob[1 + lp.mn * 2 * 3], absErr);
  EXPECT_NEAR(0, lp.prob[2 + lp.mn * 2 * 3], absErr);

  EXPECT_NEAR(0, lp.prob[4], absErr);
  EXPECT_NEAR(1, lp.prob[5], absErr);

  //  cerr << "BLOINC2\n";
  // now test refpanel loading
  lp.LoadHapLegSamp(refLegend, refHap, refSamp, InstiPanelType::REFERENCE,
                    Bio::Region());

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
  init.inputGLFile = sampleBin;
  Insti lp(init);

  lp.initialize();
  lp.LoadHapsSamp(refHaps, refSample, InstiPanelType::REFERENCE, Bio::Region());

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
  lp.LoadHapsSamp(scaffoldHaps, scaffHapsSampSample, InstiPanelType::SCAFFOLD,
                  refRegion);
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

  lp2.initialize();

  // test the scaffold loading
  lp2.LoadHapsSamp(scaffoldHaps, scaffHapsSampUnorderedSample,
                   InstiPanelType::SCAFFOLD, refRegion);
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

  lp3.initialize();
  lp3.LoadHapsSamp(refHaps, refSample, InstiPanelType::REFERENCE,
                   Bio::Region());

  // test the scaffold loading
  ASSERT_EXIT(lp3.LoadHapsSamp(scaffoldHapsWrongReg, scaffHapsSampSample,
                               InstiPanelType::SCAFFOLD, refRegion),
              ::testing::ExitedWithCode(1),
              "Number of haplotypes in haps file is 0.  Haps file empty?");

  // test the scaffold loading using tabhaps
  // unordered haps test
  Insti lp4(init);

  lp4.initialize();

  // test the scaffold loading
  lp4.LoadHapsSamp(scaffoldTabHaps, scaffHapsSampUnorderedSample,
                   InstiPanelType::SCAFFOLD, refRegion);
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
  init.inputGLFile = sampleBin;
  Insti lp(init);

  lp.initialize();

  //    cerr << "BLOINC1\n";
  ASSERT_EXIT(lp.LoadHapLegSamp("", sampleHap, "", InstiPanelType::REFERENCE,
                                Bio::Region()),
              ::testing::ExitedWithCode(1),
              "Need to define a legend file if defining a hap file");
  ASSERT_EXIT(lp.LoadHapLegSamp(sampleLegend, "", "", InstiPanelType::REFERENCE,
                                Bio::Region()),
              ::testing::ExitedWithCode(1),
              "Need to define a hap file if defining a legend file");
  ASSERT_EXIT(lp.LoadHapsSamp(refHaps, brokenHapsSampSample,
                              InstiPanelType::SCAFFOLD, refRegion),
              ::testing::ExitedWithCode(2), "Error in sample file "
                                            "../../samples/brokenFiles/"
                                            "onlyThree.hapsSample.extraLine."
                                            "sample: empty lines detected.");

  /*
                  "Input haplotypes file " + unsortedRefHaps +
                  " needs to be sorted by position");
  */
  ASSERT_EXIT(lp.LoadHapsSamp(unsortedRefHaps, scaffHapsSampUnorderedSample,
                              InstiPanelType::SCAFFOLD, refRegion),
              ::testing::ExitedWithCode(2), "Input haplotypes file "
                                            "../../samples/"
                                            "20_0_62000000.011976121_012173018."
                                            "paste.onlyThree.unsorted.haps "
                                            "needs to be sorted by position");
}

TEST(Insti, initializingHapsFromScaffold) {
  InstiHelper::Init init;
  init.geneticMap = geneticMap;
  init.scaffoldFiles["h"] = scaffoldHaps;
  init.scaffoldFiles["s"] = scaffHapsSampSample;
  init.initPhaseFromScaffold = true;
  init.inputGLFile = sampleBin;
  Insti lp(init);

  // test the scaffold loading
  lp.initialize();

  // test to see if main haps were initialized correctly
  vector<unsigned> siteIdxs = { 6, 316, 576 };
  for (auto i : siteIdxs) {
    EXPECT_EQ(0, lp.TestMainHap_(0, i));
    EXPECT_EQ(1, lp.TestMainHap_(1, i));
    EXPECT_EQ(0, lp.TestMainHap_(2, i));
    EXPECT_EQ(1, lp.TestMainHap_(3, i));
    EXPECT_EQ(0, lp.TestMainHap_(4, i));
    EXPECT_EQ(0, lp.TestMainHap_(5, i));
  }

  siteIdxs = { 635, 739, 1021 };
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

  InstiHelper::Init init;
  init.geneticMap = geneticMap;
  init.inputGLFile = sampleBin;
  Insti lp(init);

  lp.initialize();
  lp.LoadVCFGZ(refVCFGZ, InstiPanelType::REFERENCE, Bio::Region());

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
  lp.LoadVCFGZ(scaffoldVCFGZ, InstiPanelType::SCAFFOLD, refRegion);
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

  lp2.initialize();

  // test the scaffold loading
  lp2.LoadVCFGZ(scaffoldUnorderedVCFGZ, InstiPanelType::SCAFFOLD, refRegion);
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
