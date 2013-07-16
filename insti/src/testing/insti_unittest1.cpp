
#include "gtest/gtest.h"
#include "insti.h"

string sampleDir = "../../samples";

TEST(Insti, loadBin){
    Insti lp;
    string sampleBin = sampleDir + "/20_011976121_012173018.bin.onlyThree.bin";
    lp.load_bin( sampleBin.c_str());
    
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
    EXPECT_EQ(0, lp.prob[0]);
    EXPECT_EQ(1, lp.prob[1]);
    EXPECT_EQ(1, lp.prob[2]);
    EXPECT_EQ(0, lp.prob[3]);
    EXPECT_EQ(0, lp.prob[4]);
    EXPECT_EQ(0, lp.prob[5]);

    EXPECT_EQ(0, lp.prob[6]);
    EXPECT_EQ(1, lp.prob[7]);
    EXPECT_EQ(1, lp.prob[lp.in*6-4]);

}

TEST(Insti, loadHaps){

    Insti lp;
    string sampleLegend = sampleDir + "/20_011976121_012173018.bin.onlyThree.legend";
    string sampleHaps = sampleDir + "/20_011976121_012173018.bin.onlyThree.haps";
    string sampleBin = sampleDir + "/20_011976121_012173018.bin.onlyThree.bin";
    lp.load_bin( sampleBin.c_str());

    ASSERT_EXIT(lp.load_refPanel( "", sampleHaps), ::testing::ExitedWithCode(1),"Need to define a legend file if defining a reference haplotypes file");
    ASSERT_EXIT(lp.load_refPanel( sampleLegend, ""), ::testing::ExitedWithCode(1),"Need to define a reference haplotypes file if defining a legend file");
    
    lp.load_refPanel( sampleLegend, sampleHaps);

    for(unsigned i = 0; i != 601; i++){
        EXPECT_EQ(0,lp.TestRefHap(0,i));
        EXPECT_EQ(1,lp.TestRefHap(1,i));
        EXPECT_EQ(0,lp.TestRefHap(2,i));
        EXPECT_EQ(1,lp.TestRefHap(3,i));
    }

    for(unsigned i = 601; i != 1024; i++){
        EXPECT_EQ(1,lp.TestRefHap(0,i));
        EXPECT_EQ(0,lp.TestRefHap(1,i));
        EXPECT_EQ(0,lp.TestRefHap(2,i));
        EXPECT_EQ(1,lp.TestRefHap(3,i));
    }

}


