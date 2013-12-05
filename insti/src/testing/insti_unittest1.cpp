
#include "gtest/gtest.h"
#include "insti.hpp"
#include "haplotype.hpp"
#include "kMeans.hpp"
#include "relationship.hpp"
#include <algorithm>
#include <gsl/gsl_rng.h>

string sampleDir = "../../samples";
string sampleLegend = sampleDir +
                      "/20_011976121_012173018.bin.onlyThree.legend";
string sampleHap = sampleDir + "/20_011976121_012173018.bin.onlyThree.hap";
string sampleBin = sampleDir + "/20_011976121_012173018.bin.onlyThree.bin";
string refHap =  sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.hap";
string refLegend =  sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.legend";
string refHaps = sampleDir + "/20_0_62000000.011976121_012173018.paste.onlyThree.haps";

gsl_rng *rng = gsl_rng_alloc(gsl_rng_default); 

TEST(Insti, loadBin) {

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

    // test the scaffold loading

}

TEST(Insti, loadHapsSamp) {

    Insti lp;
    lp.load_bin(sampleBin.c_str());
    lp.initialize();
    lp.LoadHapsSamp("", refHaps, PanelType::REFERENCE);

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

}


TEST(Insti, loadHapLegSampErrors) {

    Insti lp;
    lp.load_bin(sampleBin.c_str());

    //    cerr << "BLOINC1\n";
    ASSERT_EXIT(lp.LoadHapLegSamp("", sampleHap, "", PanelType::REFERENCE), ::testing::ExitedWithCode(1),
                "Need to define a legend file if defining a hap file");
    ASSERT_EXIT(lp.LoadHapLegSamp(sampleLegend, "", "", PanelType::REFERENCE), ::testing::ExitedWithCode(1),
                "Need to define a hap file if defining a legend file");

}

TEST(Haplotype, StoresOK) {

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
    ASSERT_TRUE(longB.TestSite(1,hapWords.data()));

}

TEST(KMeans, clustersOK) {

    gsl_rng_set(rng, time(NULL));
    std::srand(1);

    // testing to see if nearest neighbor clustering works ok
    unsigned numClusters = 3;
    unsigned numHaps = 16;
    unsigned numSites = numHaps + numClusters*2;

    // create a set of test haplotypes to cluster and sample from
    vector< Haplotype > haplotypes;
    vector< unsigned > shuffledIndexes(numHaps);
    for (unsigned i = 0; i < numHaps; i++) {
        shuffledIndexes[i] = i;
        Haplotype temp(numSites);

        for (unsigned j = 0; j < numClusters+2; j++)
            temp.Set(i + j, 1);

        haplotypes.push_back(temp);
    }

    // shuffle the haplotypes except for the first two
    std::random_shuffle(shuffledIndexes.begin() + 2, shuffledIndexes.end());

    vector < uint64_t > passHaps;
    for (unsigned j = 0; j < shuffledIndexes.size(); j++)
        passHaps.push_back(haplotypes[shuffledIndexes[j]].GetWord(0));

    KMeans kMeans(numClusters);
    kMeans.init(passHaps, 1, numSites, 0);

    // check to make sure kMeans has the haps stored correctly        
    vector< unsigned > neighborHapNums;
    kMeans.Neighbors(0, neighborHapNums);
    
    EXPECT_EQ(numClusters, neighborHapNums.size());
    EXPECT_EQ(2,shuffledIndexes[neighborHapNums[0]]);
    EXPECT_EQ(3,shuffledIndexes[neighborHapNums[1]]);
    EXPECT_EQ(4,shuffledIndexes[neighborHapNums[2]]);

    // testing sampling
    for(unsigned i = 0; i <10; i++){
        unsigned sampHap = kMeans.SampleHap(0, rng);
        EXPECT_LT( shuffledIndexes[sampHap], 5);
        EXPECT_GT( shuffledIndexes[sampHap], 1);
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

    haplotypes[0].Set(numSites -1, true);
    haplotypes[1].Set(numSites -1, true);
    haplotypes[5].Set(numSites -1, true);
    haplotypes[15].Set(numSites -1, true);
    haplotypes[14].Set(numSites -1, true);
    haplotypes[13].Set(numSites -1, true);


    haplotypes[0].Set(numSites -2, true);
    haplotypes[1].Set(numSites -2, true);
    haplotypes[3].Set(numSites -2, true);
    haplotypes[4].Set(numSites -2, true);
    haplotypes[5].Set(numSites -2, true);
    haplotypes[8].Set(numSites -2, true);

/*    
    for(auto hap : haplotypes){
        for(unsigned i = 0; i < numSites; i++)
            cerr << hap.TestSite(i);
        cerr << endl;
    }
*/
    
    passHaps.clear();
    for (unsigned j = 0; j < shuffledIndexes.size(); j++)
        passHaps.push_back(haplotypes[shuffledIndexes[j]].GetWord(0));

    KMeans kMeans2(numClusters);
    kMeans2.init(passHaps, 1, numSites, 0.375);

    // check to make sure variant allele freqs are calculated correctly
    vector<double> varAfs;
    kMeans2.VarAfs(varAfs);
    EXPECT_EQ(numSites, varAfs.size());
    EXPECT_EQ(0.4375, varAfs[0]);
    EXPECT_EQ(0.125, varAfs[1]);
    EXPECT_EQ(0.5, varAfs[2]);
    
    // check to make sure kMeans is thresholding the correct sites
    vector<unsigned> commonSites;
    kMeans2.ClusterSites(commonSites);
    EXPECT_EQ(4, commonSites.size());
    
    // check to make sure kMeans has the haps stored correctly        
    neighborHapNums.clear();
    kMeans2.Neighbors(0, neighborHapNums);
    
    EXPECT_EQ(numClusters, neighborHapNums.size());
    EXPECT_EQ(5,shuffledIndexes[neighborHapNums[0]]);
    EXPECT_EQ(8,shuffledIndexes[neighborHapNums[1]]);
    EXPECT_EQ(6,shuffledIndexes[neighborHapNums[2]]);

    // testing sampling
    for(unsigned i = 0; i <10; i++){
        unsigned sampHap = kMeans2.SampleHap(0, rng);
        EXPECT_LT( shuffledIndexes[sampHap], 9);
        EXPECT_GT( shuffledIndexes[sampHap], 4);
        EXPECT_NE( shuffledIndexes[sampHap], 7);
    }


}



