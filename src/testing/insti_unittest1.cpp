
#include "gtest/gtest.h"
#include "insti.hpp"
#include "haplotype.hpp"
#include "kNN.hpp"
#include "relationship.hpp"
#include <algorithm>
#include <gsl/gsl_rng.h>

string sampleDir = "../../samples";
string brokenDir = sampleDir + "/brokenFiles";
string sampleLegend = sampleDir +
                      "/20_011976121_012173018.bin.onlyThree.legend";
string sampleHap = sampleDir + "/20_011976121_012173018.bin.onlyThree.hap";
string sampleBin = sampleDir + "/20_011976121_012173018.bin.onlyThree.bin";
string refHap =  sampleDir +
                 "/20_0_62000000.011976121_012173018.paste.onlyThree.hap";
string refLegend =  sampleDir +
                    "/20_0_62000000.011976121_012173018.paste.onlyThree.legend";
string refHaps = sampleDir +
                 "/20_0_62000000.011976121_012173018.paste.onlyThree.haps";

string scaffHapLegSampSample = sampleDir + "/onlyThree.hapLegSamp.sample";
string scaffHapsSampSample = sampleDir + "/onlyThree.hapsSample.sample";
string scaffHapsSampUnorderedSample = sampleDir +
                                      "/onlyThree.hapsSample.unordered.sample";

string scaffoldUnorderedHaps = sampleDir +
    "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.unordered.haps";
string scaffoldHap = sampleDir +
                     "/20_011976121_012173018.bin.onlyThree.scaffold50.hap";
string scaffoldHaps = sampleDir +
                      "/20_011976121_012173018.bin.onlyThree.scaffold50.sorted.haps";
string scaffoldLegend = sampleDir +
                        "/20_011976121_012173018.bin.onlyThree.scaffold50.legend";

string brokenHapLegSampSample = brokenDir +
                                "/onlyThree.hapLegSample.extraLine.sample";
string brokenHapsSampSample = brokenDir +
                              "/onlyThree.hapsSample.extraLine.sample";
string unsortedRefHaps = sampleDir +
                 "/20_0_62000000.011976121_012173018.paste.onlyThree.unsorted.haps";


gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

TEST(Insti, loadBin)
{

    Insti lp;
    lp.load_bin(sampleBin.c_str());

    ASSERT_EQ(3, lp.in);

    // testing sites
    EXPECT_EQ(1024, lp.site.size());

    // chr
    EXPECT_EQ("20", lp.site[0].chr);
    EXPECT_EQ("20", lp.site[5].chr);
    EXPECT_EQ("20", lp.site[1023].chr);
    EXPECT_NE("16", lp.site[1023].chr);

    // posi
    EXPECT_EQ(11976121, lp.site[0].pos);
    EXPECT_EQ(11977230, lp.site[5].pos);
    EXPECT_EQ(12173018, lp.site[1023].pos);

    // all
    EXPECT_EQ("TC", lp.site[0].all);
    EXPECT_EQ("AG", lp.site[5].all);
    EXPECT_EQ("GT", lp.site[1023].all);

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
    lp.LoadHapLegSamp(refLegend, refHap, "", PanelType::REFERENCE);

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
    //    lp.LoadHapLegSamp(refLegend, refHap, scaffHapLegSampSample, PanelType::SCAFFOLD);

}

TEST(Insti, loadHapsSamp)
{

    Insti lp;
    lp.load_bin(sampleBin.c_str());
    lp.initialize();
    lp.LoadHapsSamp(refHaps, "",  PanelType::REFERENCE);

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
    lp.LoadHapsSamp(scaffoldHaps, scaffHapsSampSample, PanelType::SCAFFOLD);
    ASSERT_EQ("samp1", lp.GetScaffoldID(0));
    ASSERT_EQ("samp2", lp.GetScaffoldID(1));
    ASSERT_EQ("samp3", lp.GetScaffoldID(2));

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
    Insti lp2;
    lp2.load_bin(sampleBin.c_str());
    lp2.initialize();

    // test the scaffold loading
    lp2.LoadHapsSamp(scaffoldHaps, scaffHapsSampUnorderedSample, PanelType::SCAFFOLD);
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


TEST(Insti, loadHapLegSampErrors)
{

    Insti lp;
    lp.load_bin(sampleBin.c_str());
    lp.initialize();

    //    cerr << "BLOINC1\n";
    ASSERT_EXIT(lp.LoadHapLegSamp("", sampleHap, "", PanelType::REFERENCE),
                ::testing::ExitedWithCode(1),
                "Need to define a legend file if defining a hap file");
    ASSERT_EXIT(lp.LoadHapLegSamp(sampleLegend, "", "", PanelType::REFERENCE),
                ::testing::ExitedWithCode(1),
                "Need to define a hap file if defining a legend file");
    ASSERT_EXIT(lp.LoadHapsSamp(refHaps, brokenHapsSampSample,
                                PanelType::SCAFFOLD), ::testing::ExitedWithCode(1),
                "Error in sample file " + brokenHapsSampSample +
                ": empty lines detected.");

    ASSERT_EXIT(lp.LoadHapsSamp(unsortedRefHaps,  scaffHapsSampUnorderedSample,
                                PanelType::SCAFFOLD), ::testing::ExitedWithCode(1),
                "Input haplotypes file " + unsortedRefHaps + " needs to be sorted by position");

}

TEST(Haplotype, StoresOK)
{

    // testing to see if init and testing works ok
    Haplotype simpleA(4);

    for (unsigned i = 0 ; i < 4; i++)
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

TEST(KNN, clustersOK)
{

    gsl_rng_set(rng, time(NULL));
    std::srand(1);

    // testing to see if nearest neighbor clustering works ok
    unsigned numClusters = 4;
    unsigned numHaps = 16;
    unsigned numSites = numHaps + numClusters * 2;

    // create a set of test haplotypes to cluster and sample from
    vector< Haplotype > haplotypes;
    vector< unsigned > shuffledIndexes(numHaps);

    for (unsigned i = 0; i < numHaps; i++) {
        shuffledIndexes[i] = i;
        Haplotype temp(numSites);

        for (unsigned j = 0; j < numClusters + 2; j++)
            temp.Set(i + j, 1);

        haplotypes.push_back(temp);
    }

    // shuffle the haplotypes except for the first two
    std::random_shuffle(shuffledIndexes.begin() + 2, shuffledIndexes.end());

    vector < uint64_t > passHaps;

    for (unsigned j = 0; j < shuffledIndexes.size(); j++)
        passHaps.push_back(haplotypes[shuffledIndexes[j]].GetWord(0));
/*
    for(unsigned idx = 0; idx != haplotypes.size(); idx++){
        cerr << idx << ": ";
        for(unsigned i = 0; i < numSites; i++)
            cerr << haplotypes[shuffledIndexes[idx]].TestSite(i);
        cerr << endl;
    }
*/

    KNN kNN(numClusters);
    kNN.init(passHaps, 1, numSites, 0);
    
    // check to make sure kNN has the haps stored correctly
    vector< unsigned > neighborHapNums = kNN.Neighbors(0);

    EXPECT_EQ(numClusters, neighborHapNums.size());
    EXPECT_EQ(2, shuffledIndexes[neighborHapNums[0]]);
    EXPECT_EQ(4, shuffledIndexes[neighborHapNums[1]]);
    EXPECT_EQ(3, shuffledIndexes[neighborHapNums[2]]);
    EXPECT_EQ(5, shuffledIndexes[neighborHapNums[3]]);

    // testing sampling
    for (unsigned i = 0; i < 10; i++) {
        unsigned sampHap = kNN.SampleHap(0, rng);
        EXPECT_LT(shuffledIndexes[sampHap], numClusters + 2);
        EXPECT_GT(shuffledIndexes[sampHap], 1);
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








