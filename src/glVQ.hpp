/* @(#)glVQ.hpp
 */

#ifndef _GLVQ_HPP
#define _GLVQ_HPP 1

#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
#include <tuple>
#include <gsl/gsl_rng.h>
#include <cassert>
#include <limits>
#include <exception>
#include <stdexcept>
#include <string>
#include <cmath>

namespace GLVQHelper {

// functor class to compare gls to codebook entries
class EuclidianDist {
public:
  EuclidianDist(double x, double y) : RR(x), Het(y) {}
  bool operator()(const std::pair<float, float> &entryA,
                  const std::pair<float, float> &entryB) {
    return sqrt(pow(entryA.first - RR, 2) + pow(entryA.second - Het, 2)) <
           sqrt(pow(entryB.first - RR, 2) + pow(entryB.second - Het, 2));
  };

private:
  float RR, Het;
};

std::pair<float, float> GLTrio2RRHet(std::vector<float>::const_iterator gl);
}

class GLVQ {
private:
  std::vector<std::pair<float, float>> m_codeBook;
  gsl_rng &m_rng;

  void BuildCodeBook(unsigned codeBookSize, const std::vector<float> &inGLs);
  double AssignPoints(
      std::vector<std::tuple<float, float, unsigned char>> &points) const;
  std::vector<std::pair<float, float>> UpdateCodeBook(
      const std::vector<std::tuple<float, float, unsigned char>> &points);

public:
  GLVQ(const std::vector<float> &inGLs, gsl_rng &rng, unsigned codeBookSize);
  unsigned char FindGLCode(float RR, float Het) const;
  std::vector<std::pair<float, float>> GetCodeBook() const {
    return m_codeBook;
  }
};

#endif /* _GLVQ_H */
