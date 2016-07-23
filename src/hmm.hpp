/* @(#)hmm.hpp
 */

#ifndef _HMM_HPP
#define _HMM_HPP 1

#include <vector>

namespace HMMHelper {

struct Init {
  gsl_rng *rng = nullptr;
  unsigned num_sites = 0;
  unsigned hap_size_in_blocks = 0;
  std::vector<fast> *transitions = nullptr;
  std::vector<fast> *emissions = nullptr;
  std::vector<uint64_t> *haps = nullptr;
}
}

class HMM {

private:
  HMMHelper::Init m_init;

public:
  HMM() = delete;
  HMM(HMMHelper::Init init);
}

#endif /* _HMM_HPP */
