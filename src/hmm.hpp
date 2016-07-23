/* @(#)hmm.hpp
 */

#ifndef _HMM_HPP
#define _HMM_HPP 1

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
}

class HMM {

private:
  HMMHelper::Init m_init;
  const std::vector<fast> &transitions;
  const std::vector<fast> &emissions;
  const std::vector<uint64_t> &haps;

public:
  HMM() = delete;
  HMM(HMMHelper::Init init, const std::vector<fast> &transitions_,
      const std::vector<fast> &emissions_, const std::vector<uint64_t> &haps_);

  std::vector<fast> alphas(unsigned sample_number, uint *proposal_haps);
  size_t sample_z(std::vector<fast>::iterator alpha_begin,
                  std::vector<fast>::iterator alpha_end,
                  unsigned first_transission, int previous_z);
  size_t sample_y(size_t sampled_z_alleles, const fast *prob_pointer);
};

#endif /* _HMM_HPP */
