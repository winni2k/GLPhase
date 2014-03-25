#ifndef _GLPACK_HPP
#define _GLPACK_HPP

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

#include <vector>
#include <cassert>
#include <cmath>
#include <limits>
#include "utils.hpp"

class GLPack {

  const unsigned m_numSamps;
  const std::vector<float> &m_inGLs;
  const unsigned m_sampleStride;
  const unsigned m_numSites;

  unsigned m_nextSampIdx = 0;
  unsigned m_lastSampIdx = 0;

  char GLTrio2Char(const std::vector<float> &GLs, unsigned idx);
  char GL2HalfChar(float GL);

public:
  GLPack(const std::vector<float> &inGLs, unsigned numSamps,
         unsigned sampleStride);
  std::vector<char> GetPackedGLs();
  unsigned GetSampleStride() { return m_sampleStride; }
  unsigned GetNextSampIdx() { return m_nextSampIdx; }
  unsigned GetLastSampIdx() { return m_lastSampIdx; }
};

#endif
