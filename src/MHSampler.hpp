/* @(#)MHSampler.hpp
 */

#ifndef _MHSampler_H
#define _MHSampler_H 1

#include <gsl/gsl_rng.h>
#include <cassert>
#include <algorithm>
#include <cmath>

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

// MH= Metropolis hastings
// DRMH = Delayed Rejection MH
enum class MHType { MH, DRMH };
class MHSampler {

private:
  MHType m_samplingType;
  gsl_rng *m_rng;
  double m_pen; // scaling constant for likelihoods

  // all likelihoods are log scaled
  double m_firstLike = 1;
  double m_secondLike = 1;
  unsigned m_firstHapNum;
  bool m_inDR = false;
  bool m_accepted = true;

  unsigned accept(double proposal, unsigned hapNum) {
    m_accepted = true;
    m_secondLike = 1;
    assert(proposal <= 0);
    m_firstLike = proposal;
    m_firstHapNum = hapNum;
    return hapNum;
  }

public:
  MHSampler(gsl_rng *rng, double like, unsigned hapNum, double pen = 1)
      : m_samplingType(MHType::MH), m_rng(rng), m_pen(pen), m_firstLike(like),
        m_firstHapNum(hapNum) {};
  MHSampler(gsl_rng *rng, double like, unsigned hapNum, MHType sampler,
            double pen = 1)
      : m_samplingType(sampler), m_rng(rng), m_pen(pen), m_firstLike(like),
        m_firstHapNum(hapNum) {};

  unsigned SampleHap(unsigned hapNum, double proposal);
  bool accepted() {
    return m_accepted;
  };
};

#endif /* _KNN_H */
