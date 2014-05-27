/* @(#)MHSampler.hpp
 */

#ifndef _MHSampler_H
#define _MHSampler_H 1

#include <gsl/gsl_rng.h>
#include <cassert>
#include <algorithm>
#include <cmath>
#include "utils.hpp"

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

// MH= Metropolis hastings
// DRMH = Delayed Rejection MH
enum class MHType { MH, DRMH };

template <class T> class MHSampler {

private:
  MHType m_samplingType;
  gsl_rng *m_rng;
  double m_pen; // scaling constant for likelihoods

  // all likelihoods are log scaled
  double m_firstLike = 1;
  double m_secondLike = 1;
  T m_firstObj;
  T *m_firstObjPtr = NULL;
  bool m_inDR = false;
  bool m_accepted = true;

  void clear() {
    m_secondLike = 1;
    m_firstObjPtr = NULL;
    m_inDR = false;
  }
  bool delayReject(T *object, T prevVal, double proposal) {
    assert(object);
    m_accepted = false;
    m_inDR = true;
    m_firstObjPtr = object;
    m_firstObj = prevVal;
    m_secondLike = proposal;
    return false;
  }
  bool reject(T &object, T prevVal) {
    m_accepted = false;
    object = prevVal;
    if (m_inDR) {
      assert(m_firstObjPtr);
      if (m_firstObjPtr)
        *m_firstObjPtr = m_firstObj;
      else
        exit(1); // logic error
    }
    clear();
    return false;
  }
  bool accept(double proposal) {
    assert(proposal <= 0);
    m_accepted = true;
    m_firstLike = proposal;
    clear();
    return true;
  }

public:
  MHSampler(gsl_rng *rng, double like, double pen = 1)
      : m_samplingType(MHType::MH), m_rng(rng), m_pen(pen), m_firstLike(like),
        m_firstObjPtr(NULL) {};
  MHSampler(gsl_rng *rng, double like, MHType sampler, double pen = 1)
      : m_samplingType(sampler), m_rng(rng), m_pen(pen), m_firstLike(like),
        m_firstObjPtr(NULL) {};

  bool SampleHap(T &object, T prevVal, double proposal);
  bool accepted() {
    return m_accepted;
  };
};

// putting this in header for extra flexibility
template <class T>
bool MHSampler<T>::SampleHap(T &object, T prevVal, double proposal) {

  assert(proposal <= 0);
  assert(m_firstLike <= 0);
  if (proposal >= m_firstLike) {
    return accept(proposal);
  }

  // simple metropolis hastings sampling
  if (m_samplingType == MHType::MH) {
    if (gsl_rng_uniform(m_rng) < std::exp((proposal - m_firstLike) * m_pen)) {
      return accept(proposal);
    } else
      return reject(object, prevVal);
  }
  // delayed rejection MH sampler
  else if (m_samplingType == MHType::DRMH) {
    if (m_inDR) {

      assert(m_secondLike < m_firstLike);

      // short out if the proposal is clearly going to fail the next test
      if (proposal < m_secondLike) {
        return reject(object, prevVal);
      }

      // calculate complete acceptance probability
      double first2prop = std::exp((proposal - m_firstLike) * m_pen);
      double accProbProp2sec =
          1 - std::min(1.0, std::exp(m_secondLike - proposal));
      double accProbFirst2sec =
          1 - std::min(1.0, std::exp(m_secondLike - m_firstLike));

      assert(accProbProp2sec >= 0 && accProbProp2sec <= 1);
      assert(accProbFirst2sec >= 0 && accProbFirst2sec <= 1);
      if (gsl_rng_uniform(m_rng) <
          first2prop * accProbProp2sec / accProbFirst2sec) {
        return accept(proposal);
      }

      // in this case we rejected the second proposal
      return reject(object, prevVal);
    }
    // we have not delayed rejection before
    else {
      if (gsl_rng_uniform(m_rng) < std::exp((proposal - m_firstLike) * m_pen)) {
        return accept(proposal);
      }

      // delay rejection
      else {
        return delayReject(&object, prevVal, proposal);
      }
    }
  }

  // we never want to end up here!
  assert(false);
  exit(1); // logic error
}

#endif /* _KNN_H */
