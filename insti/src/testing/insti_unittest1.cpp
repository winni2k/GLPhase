
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
    // making sure prob size is correct
    EXPECT_EQ(1024*2*lp.in, lp.prob.size());
        
    EXPECT_EQ(0, lp.prob[0]);
    EXPECT_EQ(1, lp.prob[1]);
    EXPECT_EQ(0.2f, lp.prob[2]);
    EXPECT_EQ(0, lp.prob[3]);
    EXPECT_EQ(0, lp.prob[4]);
    EXPECT_EQ(0, lp.prob[5]);
    EXPECT_EQ(0, lp.prob[6]);
    EXPECT_EQ(1, lp.prob[7]);
    EXPECT_EQ(1, lp.prob[lp.in*6-4]);

    // now initialize lp and see if probs still make sense
    lp.initialize();

    EXPECT_EQ(1024*3*lp.in, lp.prob.size());
    EXPECT_EQ(0, lp.prob[0]);
    EXPECT_EQ(0, lp.prob[1]);
    EXPECT_EQ(1, lp.prob[2]);

    // 3 = lp.pn
    EXPECT_EQ(0.8f, lp.prob[0 + lp.mn * 3]);
    EXPECT_EQ(0.2f, lp.prob[1 + lp.mn * 3]);
    EXPECT_EQ(0, lp.prob[2 + lp.mn * 3]);

    EXPECT_EQ(0, lp.prob[0 + lp.mn *2 * 3]);
    EXPECT_EQ(0, lp.prob[1 + lp.mn *2 * 3]);
    EXPECT_EQ(0, lp.prob[2 + lp.mn *2 * 3]);

    EXPECT_EQ(0, lp.prob[4]);
    EXPECT_EQ(1, lp.prob[5]);
}

TEST(Insti, loadHaps){

    Insti lp;
    string sampleLegend = sampleDir + "/20_011976121_012173018.bin.onlyThree.legend";
    string sampleHaps = sampleDir + "/20_011976121_012173018.bin.onlyThree.haps";
    string sampleBin = sampleDir + "/20_011976121_012173018.bin.onlyThree.bin";
    lp.load_bin( sampleBin.c_str());

//    cerr << "BLOINC1\n";
    ASSERT_EXIT(lp.load_refPanel( "", sampleHaps), ::testing::ExitedWithCode(1),"Need to define a legend file if defining a reference haplotypes file");
    ASSERT_EXIT(lp.load_refPanel( sampleLegend, ""), ::testing::ExitedWithCode(1),"Need to define a reference haplotypes file if defining a legend file");
    
    lp.load_refPanel( sampleLegend, sampleHaps);

//    cerr << "BLOINC2\n";
    for(unsigned i = 0; i != 601; i++){
//        cerr << "\tbloincing0\t" << i << endl;
//        cerr << lp.TestRefHap(0,i);
        EXPECT_EQ(0,lp.TestRefHap(0,i));
//        cerr << "\tbloincing0.5\t" << i << endl;
//        cerr << "refhap: 1\tsite: "<< i << "\t" << lp.TestRefHap(1,i) << endl;
//        cerr << "\t"<< lp.TestRefHap(1,i);
//        EXPECT_EQ(0, lp.TestRefHap(1,i));
//        cerr << "\tbloincing1\t" << i << endl;
//        EXPECT_EQ(1,lp.TestRefHap(1,i));
//      cerr << "\tbloincing2\t" << i << endl;
        EXPECT_EQ(0,lp.TestRefHap(2,i));
//        cerr << "\t"<< lp.TestRefHap(2,i);
//        cerr << "\tbloincing3\t" << i << endl;
        EXPECT_EQ(1,lp.TestRefHap(3,i));
//        cerr << "\t"<< lp.TestRefHap(3,i)<<endl;
    }

//    cerr << "BLOINC3\n";
    for(unsigned i = 601; i != 1024; i++){
        EXPECT_EQ(1,lp.TestRefHap(0,i));
        EXPECT_EQ(0,lp.TestRefHap(1,i));
        EXPECT_EQ(0,lp.TestRefHap(2,i));
        EXPECT_EQ(1,lp.TestRefHap(3,i));
    }

}


