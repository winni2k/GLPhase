/* @(#)hmm.hpp
 */

#ifndef _HMM_HPP
#define _HMM_HPP 1

#include <algorithm>
#include <vector>
#include <gsl/gsl_rng.h>
#include "impute.hpp"

namespace HMMHelper {

struct Init {
  gsl_rng *rng = nullptr;
  unsigned num_sites = 0;
  unsigned hap_size_in_blocks = 0;
  unsigned num_possible_emissions = 0;
};

// calculate number of recombinations necessary for state transition
inline int num_recombinations(unsigned from, unsigned to) {
  unsigned mismatches = from ^ to;
  return (mismatches & 1) + ((mismatches >> 1) & 1);
}

inline size_t roulette_wheel_sample(const std::vector<fast> &values,
                                    gsl_rng *rng) {
  assert(!values.empty());
  const fast sum = std::accumulate(values.begin(), values.end(), 0.0);
  assert(sum > 0);
  fast sampled_sum = gsl_rng_uniform(rng) * sum;
  size_t val_idx = 0;
  do {
    sampled_sum -= values.at(val_idx);
    if (sampled_sum <= 0)
      return val_idx;
    ++val_idx;
  } while (true);
}
}

class HMM {

private:
  HMMHelper::Init m_init;
  const std::vector<fast> &m_transitions;
  const std::vector<fast> &m_emissions;
  const std::vector<uint64_t> &m_haps;

public:
  HMM() = delete;
  HMM(HMMHelper::Init init, const std::vector<fast> &transitions,
      const std::vector<fast> &emissions, const std::vector<uint64_t> &haps);

  std::vector<fast> alphas(unsigned sample_number, const unsigned *proposal_haps);
  size_t sample_z(std::vector<fast>::iterator alpha_begin,
                  std::vector<fast>::iterator alpha_end,
                  const fast *transition_triple, int previous_z);
  size_t sample_y(size_t sampled_z_alleles, const fast *prob_pointer,
                  const fast(&mutation_mat)[4][4]);
  std::vector<size_t>
  sample_y_sequence(unsigned individual, const unsigned *parental_hap_idxs,
                    const fast(&mutation_mat)[4][4],
                    const std::vector<fast> &emission_probs);
};

#endif /* _HMM_HPP */
