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

  char GLTrio2Char(unsigned idx) const;
  char GL2HalfChar(float GL) const;

public:
  GLPack(const std::vector<float> &inGLs, unsigned numSamps,
         unsigned sampleStride);
  std::vector<char> GetPackedGLs();
  unsigned GetSampleStride() const { return m_sampleStride; }
  unsigned GetNextSampIdx() const { return m_nextSampIdx; }
  unsigned GetLastSampIdx() const { return m_lastSampIdx; }
  unsigned GetNumSites() const { return m_numSites; }
  unsigned GetNumSamps() const { return m_numSamps; }
};

#endif
