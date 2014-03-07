
#include "MHSampler.hpp"

unsigned MHSampler::SampleHap(unsigned hapNum, double proposal) {

  m_accepted = false;
  assert(proposal <= 0);
  assert(m_firstLike <= 0);
  if (proposal >= m_firstLike) {
    return accept(proposal, hapNum);
  }

  // simple metropolis hastings sampling
  if (m_samplingType == MHType::MH) {
    if (gsl_rng_uniform(m_rng) < exp((proposal - m_firstLike) * m_pen)) {
      return accept(proposal, hapNum);
    } else
      return m_firstHapNum;
  }
  // delayed rejection MH sampler
  else if (m_samplingType == MHType::DRMH) {
    if (m_inDR) {

      // whatever happens next, we are not in DR anymore
      m_inDR = false;
      assert(m_secondLike < m_firstLike);

      // short out if the proposal is clearly going to fail the next test
      if (proposal < m_secondLike) {

        // just a test to make sure things are going as expected
        m_secondLike = 1;
        return m_firstHapNum;
      }

      // calculate complete acceptance probability
      double first2prop = exp((proposal - m_firstLike) * m_pen);
      double accProbProp2sec = 1 - min(1, exp(m_secondLike - proposal));
      double accProbFirst2sec = 1 - min(1, exp(m_secondLike - m_firstLike));

      assert(accProbProp2sec >= 0 && accProbProp2sec <= 1);
      assert(accProbFirst2sec >= 0 && accProbFirst2sec <= 1);
      if (gsl_rng_uniform(m_rng) <
          first2prop * accProbProp2sec / accProbFirst2sec) {
        return accept(proposal, hapNum);
      }

      m_secondLike = 1;
      return m_firstHapNum;
    }
    // we have not delayed rejection before
    else {
      if (gsl_rng_uniform(m_rng) < exp((proposal - m_firstLike) * m_pen)) {
        return accept(proposal, hapNum);
      }
      // delay rejection
      else {
        m_inDR = true;
        m_secondLike = proposal;
        return hapNum;
      }
    }
  }

  // we never want to end up here!
  assert(false);
}
