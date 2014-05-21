/* @(#)glVQ.hpp
 */

#ifndef _GLVQ_HPP
#define _GLVQ_HPP 1

#include <vector>
#include <algorithm>
#include <iterator>

namespace GLVQHelper {

// functor class to compare gls to codebook entries
class EuclidianDist {
public:
  EuclidianDist(float x, float y) : RR(x), Het(y) {}
  bool operator()(pair<float, float> &entryA, pair<float, float> &entryB) {

    return sqrt(pow(abs(entryA.first - RR), 2) +
                pow(abs(entryA.second - Het), 2)) <
           sqrt(pow(abs(entryB.first - RR), 2) +
                pow(abs(entryB.second - Het), 2))
  }

private:
  float RR, Het;
};
}

class GLVQ {
private:
  std::vector<pair<float, float> > m_codeBook;
  gsl_rng &rng;

public:
  char FindGLCode(float RR, float Het);
  std::vector<pair<float, float> > GetCodeBook() const { return m_codeBook; }
};

// create code book a la wikipedia
// https://en.wikipedia.org/wiki/Vector_quantization
GLVQ(const std::vector<float> &inGls, unsigned codeBookSize,
     const float sensInc = 0.1) {

  assert(inGLs.size() % 2 == 0);
  const numSamps = inGLs.size() / 2;
  assert(numSamps > 0);

  assert(codeBookSize > 0);
  m_codeBook.reserve(codeBookSize);

  vector<float> entrySens;
  entrySens.reserve(codeBookSize);

  // initialize code book centroids randomly
  for (unsigned entryNum = 0; entryNum < codeBookSize; ++entryNum) {
    unsigned sampIdx = gsl_rng_uniform_int(&m_rng, numSamps);
    m_codeBook.push_back(make_pair(inGLs[sampIdx * 2], inGLs[sampIdx * 2 + 1]));
    entrySens.push_back(sensInc);
  }

  // update train codebook
  
}

// use code book to find code corresponding to a pair of floats
GLVQ::FindGLCode(float RR, float Het) {

  // the gl code is the distance between the codebook element with the least
  // distance between it and the input gls, and the first element
  // yay for lambdas ...
  EuclidianDist comp(RR, Het);
  return distance(min_element(m_codeBook.begin(), m_codeBook.end(), comp),
                  m_codeBook.begin());
}

#endif /* _GLVQ_H */
