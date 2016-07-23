#include "hmm.hpp"

using namespace std;

namespace HMMHelper {

size_t roulette_wheel_sample(const vector<fast> &values, gsl_rng *rng) {
  assert(!values.empty());
  fast sum = accumulate(values.begin(), values.end(), 0);
  fast sampled_sum = gsl_rng_uniform(rng) * sum;
  size_t val_idx = 0;
  do {
    sampled_sum -= values.at(val_idx);
    if (sampled_sum > 0)
      return val_idx;
    ++val_idx;
  } while (true);
}
}

HMM::HMM(HMMHelper::Init init) : m_init(move(init)) {
  assert(m_init.num_sites > 0);
  assert(m_init.hap_size_in_blocks > 0);
  assert(m_init.transitions != nullptr);
  assert(m_init.emissions != nullptr);
}

vector<fast> HMM::alphas(unsigned sample_number, uint *proposal_haps) {

  // setup the different haplotypes
  uint64_t *f0 = &m_init.haps[proposal_haps[0] * m_init.hap_size_in_blocks],
           *f1 = &m_init.haps[proposal_haps[1] * m_init.hap_size_in_blocks],
           *m0 = &m_init.haps[proposal_haps[2] * m_init.hap_size_in_blocks],
           *m1 = &m_init.haps[proposal_haps[3] * m_init.hap_size_in_blocks];

  vector<fast> alphas;
  alphas.reserve(mn * 4);

  // pull out phase emission and transition probabilities
  fast *emit_pointer = &emissions[sample_number * en],
       *tran_pointer = &transitions[0], sum = 0;

  // l00 = prob of 0|0 phase, etc.
  // all set to 1/4 * emission probability
  fast l00 = 0.25f * emit_pointer[(test(f0, 0) << 1) | test(m0, 0)],
       l01 = 0.25f * emit_pointer[(test(f0, 0) << 1) | test(m1, 0)];
  fast l10 = 0.25f * emit_pointer[(test(f1, 0) << 1) | test(m0, 0)],
       l11 = 0.25f * emit_pointer[(test(f1, 0) << 1) | test(m1, 0)];

  // l00 = prob of 0|0 phase, etc.
  // all set to 1/4 * emission probability
  for (int state = 0; state < 4; ++state) {
    p_haps = parent_hap_pairs[i];
    alphas.push_back(
        0.25f *
        emit_pointer[(test(p_haps.first, 0) << 1) | test(p_haps.second, 0)]);
  }
  alphas.push_back(l00);
  alphas.push_back(l01);
  alphas.push_back(l10);
  alphas.push_back(l11);

  // bxx = backward probabilities of being in phase xx
  fast b00, b01, b10, b11;

  // move to next site for e and t
  emit_pointer += 4;
  tran_pointer += 3;

  assert(alphas.size() == m * 4);
  fast l00 = alphas[m * 4 - 4], l01 = alphas[m * 4 - 3],
       l10 = alphas[m * 4 - 2], l11 = alphas[m * 4 - 1];
  // calculate total probability of model given the four haplotypes
  // passed in, and return as score
  for (unsigned m = 1; m < mn; m++, emit_pointer += 4, tran_pointer += 3) {
    b00 = l00 * t[0] + (l01 + l10) * t[1] + l11 * t[2];
    b01 = l01 * t[0] + (l00 + l11) * t[1] + l10 * t[2];
    b10 = l10 * t[0] + (l00 + l11) * t[1] + l01 * t[2];
    b11 = l11 * t[0] + (l01 + l10) * t[1] + l00 * t[2];
    l00 = b00 * e[(test(f0, m) << 1) | test(m0, m)];
    l01 = b01 * e[(test(f0, m) << 1) | test(m1, m)];
    l10 = b10 * e[(test(f1, m) << 1) | test(m0, m)];
    l11 = b11 * e[(test(f1, m) << 1) | test(m1, m)];

    alphas.push_back(l00);
    alphas.push_back(l01);
    alphas.push_back(l10);
    alphas.push_back(l11);

    // rescale probabilities if they become too small
    if ((sum = l00 + l01 + l10 + l11) < norm) {
      sum = 1.0f / sum;
      score -= logf(sum); // add sum to score
      l00 *= sum;
      l01 *= sum;
      l10 *= sum;
      l11 *= sum;
    }
  }
  assert(alphas.size() == mn * 4);
  return alphas;
}

size_t HMM::sample_z(std::vector<fast>::iterator alpha_begin,
                     std::vector<fast>::iterator alpha_end,
                     unsigned first_transission, int previous_z) {

  const *fast tran_pointer = &m_init.transissions[first_transission];

  // initialization
  assert(distance(alpha_begin, alpha_end) == 4);
  vector<fast> probs;
  vector<fast> transition_probs(4, 1);
  for (unsigned state = 0; state < 4; ++state) {
    if (previous_z != -1) {
      transition_probs[state] =
          tran_pointer[num_recombinations(static_cast(previous_z), state);];
    }
    probs.push_back(transition_probs[state] * *(alpha_begin + state));
  }

  return HMMHelper::roulette_wheel_sample(probs, m_init.rng);
}

size_t HMM:sample_y(size_t sampled_z, 
