#include "hmm.hpp"

using namespace std;
using namespace ImputeHelper;

HMM::HMM(HMMHelper::Init init, const std::vector<fast> &transitions,
         const std::vector<fast> &emissions, const std::vector<uint64_t> &haps)
    : m_init(move(init)), m_transitions(transitions), m_emissions(emissions),
      m_haps(haps) {
  assert(m_init.num_sites > 0);
  assert(m_init.hap_size_in_blocks > 0);
  assert(m_init.rng != nullptr);
}

vector<fast> HMM::alphas(unsigned sample_number, const unsigned *proposal_haps) {

  // setup the different haplotypes
  const uint64_t *f0 = &m_haps[proposal_haps[0] * m_init.hap_size_in_blocks],
                 *f1 = &m_haps[proposal_haps[1] * m_init.hap_size_in_blocks],
                 *m0 = &m_haps[proposal_haps[2] * m_init.hap_size_in_blocks],
                 *m1 = &m_haps[proposal_haps[3] * m_init.hap_size_in_blocks];

  vector<fast> alphas;
  alphas.reserve(m_init.num_sites * 4);

  // pull out phase emission and transition probabilities
  const fast *z_emit =
                 &m_emissions[sample_number * m_init.num_possible_emissions],
             *tran = &m_transitions[0];

  // l00 = prob of 0|0 phase, etc.
  // all set to 1/4 * emission probability
  fast l00 = 0.25f * z_emit[(test(f0, 0) << 1) | test(m0, 0)],
       l01 = 0.25f * z_emit[(test(f0, 0) << 1) | test(m1, 0)];
  fast l10 = 0.25f * z_emit[(test(f1, 0) << 1) | test(m0, 0)],
       l11 = 0.25f * z_emit[(test(f1, 0) << 1) | test(m1, 0)];

  alphas.push_back(l00);
  alphas.push_back(l01);
  alphas.push_back(l10);
  alphas.push_back(l11);

  // move to next site for e and t
  z_emit += 4;
  tran += 3;

  // calculate total probability of model given the four haplotypes
  // passed in, and return as score
  for (unsigned m = 1; m < m_init.num_sites; m++, z_emit += 4, tran += 3) {
    const fast sum00 = l00 * tran[0] + (l01 + l10) * tran[1] + l11 * tran[2];
    const fast sum01 = l01 * tran[0] + (l00 + l11) * tran[1] + l10 * tran[2];
    const fast sum10 = l10 * tran[0] + (l00 + l11) * tran[1] + l01 * tran[2];
    const fast sum11 = l11 * tran[0] + (l01 + l10) * tran[1] + l00 * tran[2];
    l00 = sum00 * z_emit[(test(f0, m) << 1) | test(m0, m)];
    l01 = sum01 * z_emit[(test(f0, m) << 1) | test(m1, m)];
    l10 = sum10 * z_emit[(test(f1, m) << 1) | test(m0, m)];
    l11 = sum11 * z_emit[(test(f1, m) << 1) | test(m1, m)];

    // rescale probabilities if they become too small
    fast sum = 0;
    if ((sum = l00 + l01 + l10 + l11) < norm) {
      sum = 1.0f / sum;
      l00 *= sum;
      l01 *= sum;
      l10 *= sum;
      l11 *= sum;
    }

    alphas.push_back(l00);
    alphas.push_back(l01);
    alphas.push_back(l10);
    alphas.push_back(l11);
  }
  assert(alphas.size() == m_init.num_sites * 4);
  return alphas;
}

size_t HMM::sample_z(std::vector<fast>::iterator alpha_begin,
                     std::vector<fast>::iterator alpha_end,
                     const fast *transition_triple, int previous_z) {
  // initialization
  assert(distance(alpha_begin, alpha_end) == 4);
  vector<fast> probs;
  vector<fast> transition_probs(4, 1);
  for (unsigned state = 0; state < 4; ++state) {
    if (previous_z != -1) {
      assert(previous_z >= 0);
      transition_probs[state] = transition_triple[HMMHelper::num_recombinations(
          static_cast<unsigned>(previous_z), state)];
    }
    probs.push_back(transition_probs[state] * *(alpha_begin + state));
  }

  return HMMHelper::roulette_wheel_sample(probs, m_init.rng);
}

size_t HMM::sample_y(size_t sampled_z_alleles, const fast *prob_pointer,
                     const fast(&mutation_mat)[4][4]) {

  vector<fast> vals;
  for (int y = 0; y < 4; ++y) {
    vals.push_back(prob_pointer[to_gt(y)] * mutation_mat[sampled_z_alleles][y]);
  }
  return HMMHelper::roulette_wheel_sample(vals, m_init.rng);
}

vector<size_t> HMM::sample_y_sequence(unsigned individual,
                                      const unsigned *parental_hap_idxs,
                                      const fast(&mutation_mat)[4][4],
                                      const std::vector<fast> &emission_probs) {
  // setup the different haplotypes
  const uint64_t *father_hap_0 =
                     &m_haps[parental_hap_idxs[0] * m_init.hap_size_in_blocks],
                 *father_hap_1 =
                     &m_haps[parental_hap_idxs[1] * m_init.hap_size_in_blocks],
                 *mother_hap_0 =
                     &m_haps[parental_hap_idxs[2] * m_init.hap_size_in_blocks],
                 *mother_hap_1 =
                     &m_haps[parental_hap_idxs[3] * m_init.hap_size_in_blocks];

  auto alpha_vals = alphas(individual, parental_hap_idxs);
  assert(alpha_vals.size() == m_init.num_sites * 4);
  int previous_z = -1;
  vector<size_t> sampled_y_sequence(m_init.num_sites);
  for (int site = m_init.num_sites - 1; site != -1; --site) {

    size_t sampled_z =
        sample_z(alpha_vals.begin() + site * 4, alpha_vals.begin() + site * 4 + 4,
                 &m_transitions[site * 3], previous_z);
    previous_z = static_cast<int>(sampled_z);

    const vector<uint64_t> parents_to_hap_pair = {
        (test(father_hap_0, site) << 1) | test(mother_hap_0, site),
        (test(father_hap_0, site) << 1) | test(mother_hap_1, site),
        (test(father_hap_1, site) << 1) | test(mother_hap_0, site),
        (test(father_hap_1, site) << 1) | test(mother_hap_1, site)};
    size_t sampled_y = sample_y(
        parents_to_hap_pair.at(sampled_z),
        &emission_probs[individual * 3 * m_init.num_sites] + 3 * site, mutation_mat);
    sampled_y_sequence.at(site) = sampled_y;
  }
  return sampled_y_sequence;
}
