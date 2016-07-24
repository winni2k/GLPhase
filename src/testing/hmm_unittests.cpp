#include "gtest/gtest.h"
#include "hmm.hpp"
#include <gsl/gsl_rng.h>

using namespace std;

gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

class HMMTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    init.rng = rng;
    init.num_sites = 64;
    init.hap_size_in_blocks = 1;
    init.num_possible_emissions = 4 * init.num_sites;

    transitions.resize(3 * init.num_sites, 1);
    z_emissions.resize(init.num_possible_emissions, 0);
    emission_probs.resize(3 * init.num_sites, 0);
    haps.resize(4, 0);

    for (unsigned i = 0; i < init.num_sites; ++i) {
      ImputeHelper::set1(&haps[0], i);
      ImputeHelper::set0(&haps[1], i);
      ImputeHelper::set1(&haps[2], i);
      ImputeHelper::set0(&haps[3], i);

      // let's only allow emissions from hom ref
      z_emissions.at(i * 4) = 1;
      emission_probs.at(i * 3) = 1;
    }

    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        mutation_mat[i][j] = 0;
        if (i == j)
          mutation_mat[i][j] = 1;
      }
    }
  }

  HMMHelper::Init init;
  std::vector<fast> transitions;
  std::vector<fast> z_emissions;
  std::vector<fast> emission_probs;
  std::vector<uint64_t> haps;
  vector<unsigned> prop_hap_numbers = {0, 1, 2, 3};
  fast mutation_mat[4][4]; // mutation matrix
};

TEST_F(HMMTest, alphas) {

  HMM hmm(init, transitions, z_emissions, haps);
  auto alphas = hmm.alphas(0, prop_hap_numbers.data());

  ASSERT_EQ(init.num_sites * 4, alphas.size());
  for (unsigned site = 0; site < init.num_sites; ++site) {
    EXPECT_FLOAT_EQ(0, alphas.at(site * 4));
    EXPECT_FLOAT_EQ(0, alphas.at(site * 4 + 1));
    EXPECT_FLOAT_EQ(0, alphas.at(site * 4 + 2));
    ASSERT_FLOAT_EQ(0.25, alphas.at(site * 4 + 3));
  }
}

TEST_F(HMMTest, sample_z) {

  HMM hmm(init, transitions, z_emissions, haps);
  auto alphas = hmm.alphas(0, prop_hap_numbers.data());

  int previous_z = -1;
  for (int site = init.num_sites - 1; site != -1; --site) {
    size_t sampled_z = hmm.sample_z(alphas.begin(), alphas.begin() + 4,
                                    &transitions.at(site * 3), previous_z);
    previous_z = static_cast<int>(sampled_z);
    ASSERT_FLOAT_EQ(3, sampled_z);
  }
}

TEST_F(HMMTest, sample_y) {

  HMM hmm(init, transitions, z_emissions, haps);
  std::vector<fast> probs(3, 0.25);

  for (int i = 0; i < 100; ++i) {
    ASSERT_FLOAT_EQ(0, hmm.sample_y(0, probs.data(), mutation_mat));
    ASSERT_FLOAT_EQ(1, hmm.sample_y(1, probs.data(), mutation_mat));
    ASSERT_FLOAT_EQ(2, hmm.sample_y(2, probs.data(), mutation_mat));
    ASSERT_FLOAT_EQ(3, hmm.sample_y(3, probs.data(), mutation_mat));
  }
}

TEST_F(HMMTest, sample_y_flat_mut_mat) {

  HMM hmm(init, transitions, z_emissions, haps);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      mutation_mat[i][j] = 1;
    }
  }
  std::vector<fast> probs = {1, 0, 0};

  for (int i = 0; i < 100; ++i) {
    ASSERT_FLOAT_EQ(0, hmm.sample_y(0, probs.data(), mutation_mat));
    ASSERT_FLOAT_EQ(0, hmm.sample_y(1, probs.data(), mutation_mat));
    ASSERT_FLOAT_EQ(0, hmm.sample_y(2, probs.data(), mutation_mat));
    ASSERT_FLOAT_EQ(0, hmm.sample_y(3, probs.data(), mutation_mat));
  }
}

TEST_F(HMMTest, roulette_wheel_sample) {
  std::vector<fast> probs = {0, 0, 0.25, 0};
  for (int i = 0; i < 100; ++i) {
    ASSERT_EQ(2, HMMHelper::roulette_wheel_sample(probs, rng));
  }
}

TEST_F(HMMTest, sample_y_sequence) {

  HMM hmm(init, transitions, z_emissions, haps);
  auto sampled_y_seq = hmm.sample_y_sequence(0, prop_hap_numbers.data(),
                                             mutation_mat, emission_probs);

  ASSERT_EQ(init.num_sites, sampled_y_seq.size());
  for (size_t i = 0; i < init.num_sites; ++i) {
    ASSERT_EQ(0, sampled_y_seq[i]);
  }
}

TEST_F(HMMTest, sample_y_sequence_tricky) {

  vector<size_t> other_sites = {3, 7, 27};
  for (auto site : other_sites) {
    z_emissions[site * 4] = 0;
    z_emissions[site * 4 + 3] = 1;
    emission_probs[site * 3] = 0;
    emission_probs[site * 3 + 2] = 1;
  }

  HMM hmm(init, transitions, z_emissions, haps);
  auto sampled_y_seq = hmm.sample_y_sequence(0, prop_hap_numbers.data(),
                                             mutation_mat, emission_probs);

  ASSERT_EQ(init.num_sites, sampled_y_seq.size());
  for (size_t site = 0; site < init.num_sites; ++site) {
    if (std::find(other_sites.begin(), other_sites.end(), site) !=
        other_sites.end())
      ASSERT_EQ(3, sampled_y_seq[site]);
    else
      ASSERT_EQ(0, sampled_y_seq[site]);
  }
}
