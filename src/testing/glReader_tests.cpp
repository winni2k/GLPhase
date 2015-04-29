#include "gtest/gtest.h"
#include "glReader.hpp"

using namespace std;

string sampleDir = "../../samples";
string sampleBin = sampleDir + "/20_011976121_012173018.bin.onlyThree.bin";
string glReaderDir = sampleDir + "/glReader";
string sampleBase = glReaderDir + "/20_011976121_012173018.bin.onlyThree";
string sampleVCF = sampleBase + ".vcf";
string sampleVCFGZ = sampleVCF + ".gz";
string sampleBCF = sampleBase + ".bcf";
string sampleBCFGZ = sampleBCF + ".gz";

namespace test {

// for all tests applicable to files derived from sampleBin
void sampleBinTests(Bio::GLReader &reader) {
  // test names
  auto names = reader.GetNames();
  ASSERT_EQ(names.size(), 3);
  EXPECT_EQ("samp1", names[0]);
  EXPECT_EQ("samp2", names[1]);
  EXPECT_EQ("samp3", names[2]);
  const float abs = 0.001;

  // test GLs
  // see if site exists
  {
    auto gls = reader.GetGLs();
    Bio::snp searchSNP("20", 11977595, "G", "A");
    EXPECT_TRUE(gls.second.exists(searchSNP));
    EXPECT_TRUE(*(gls.second.at(7)) == searchSNP);

    // check that site's GLs
    EXPECT_FLOAT_EQ(0, gls.first.at(3 * 3 * 7));
    EXPECT_FLOAT_EQ(gls.first.at(3 * 3 * 7 + 1), 0);
    EXPECT_FLOAT_EQ(gls.first.at(3 * 3 * 7 + 2), 1);

    EXPECT_FLOAT_EQ(gls.first.at(3 * 3 * 7 + 3), 0);
    EXPECT_FLOAT_EQ(gls.first.at(3 * 3 * 7 + 4), 1);
    EXPECT_FLOAT_EQ(gls.first.at(3 * 3 * 7 + 5), 0);

    EXPECT_FLOAT_EQ(gls.first.at(3 * 3 * 7 + 6), 1);
    EXPECT_FLOAT_EQ(gls.first.at(3 * 3 * 7 + 7), 0);
    EXPECT_FLOAT_EQ(gls.first.at(3 * 3 * 7 + 8), 0);

    // check first site's GLs
    EXPECT_FLOAT_EQ(gls.first.at(0), 0);
    EXPECT_FLOAT_EQ(gls.first.at(1), 0);
    EXPECT_FLOAT_EQ(gls.first.at(2), 1);

    EXPECT_NEAR(0.8, gls.first.at(3), abs);
    EXPECT_NEAR(0.2, gls.first.at(4), abs);
    EXPECT_FLOAT_EQ(gls.first.at(5), 0);

    EXPECT_FLOAT_EQ(gls.first.at(6), 1);
    EXPECT_FLOAT_EQ(gls.first.at(7), 0);
    EXPECT_FLOAT_EQ(gls.first.at(8), 0);
  }
  {
    reader.SetRetGLType(Bio::GLHelper::gl_ret_t::ST_DROP_FIRST);
    auto gls = reader.GetGLs();

    EXPECT_FLOAT_EQ(gls.first.at(0), 0);
    EXPECT_FLOAT_EQ(gls.first.at(1), 1);

    EXPECT_NEAR(0.2, gls.first.at(2), abs);
    EXPECT_FLOAT_EQ(gls.first.at(3), 0);

    EXPECT_FLOAT_EQ(gls.first.at(4), 0);
    EXPECT_FLOAT_EQ(gls.first.at(5), 0);
  }
}
}

TEST(GLReader, loadsSTBin) {

  Bio::GLHelper::init init;
  init.nameFile = sampleBin;
  init.glFile = sampleBin;
  init.glType = Bio::GLHelper::gl_t::STBIN;
  init.glRetType = Bio::GLHelper::gl_ret_t::STANDARD;

  Bio::GLReader reader(init);

  test::sampleBinTests(reader);
}

TEST(GLReader, loadsBCF) {
  vector<string> glFiles{sampleVCF, sampleVCFGZ, sampleBCF, sampleBCFGZ};

  for (auto file : glFiles) {
    Bio::GLHelper::init init;
    init.nameFile = file;
    init.glFile = file;
    init.glType = Bio::GLHelper::gl_t::BCF;
    init.glRetType = Bio::GLHelper::gl_ret_t::STANDARD;

    Bio::GLReader reader(init);

    test::sampleBinTests(reader);
  }
}
