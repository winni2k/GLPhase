#ifndef _GLPACK_HPP
#define _GLPACK_HPP

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

#include <vector>
#include <cassert>
#include <cmath>
#include <limits>
#include "utils.hpp"
#include <vector_types.h>
#include "glVQ.hpp"

namespace GLPackHelper {

struct Init {
  std::vector<float> &inGLs;
  unsigned numSamps = 0;
  unsigned sampleStride = 0;
  bool useVQ = false;
  unsigned numBitsPerGL = 8;
  Init(std::vector<float> &GLs) : inGLs(GLs) {};
};
}

class GLPack {

private:
  const unsigned m_numSamps;
  const std::vector<float> &m_inGLs;
  const unsigned m_sampleStride;
  const unsigned m_numSites;
  const bool m_useVQ;
  unsigned m_numBitsPerGL;
  unsigned m_numGLsPeruchar4;

  unsigned m_nextSampIdx = 0;
  unsigned m_lastSampIdx = 0;
  GLVQ m_VQ;

  char GLTrio2Char(unsigned idx) const;
  char GL2HalfChar(float GL) const;
  void GenerateCodeBook() {};
  char FindGLCode(float RR, float Het);

public:
  // set VQBits > 0 for vector quantization
  GLPack(GLPackHelper::Init &init);

  std::vector<uchar4> GetPackedGLs();
  unsigned GetSampleStride() const { return m_sampleStride; }
  unsigned GetNextSampIdx() const { return m_nextSampIdx; }
  unsigned GetLastSampIdx() const { return m_lastSampIdx; }
  unsigned GetNumSites() const { return m_numSites; }
  unsigned GetNumSamps() const { return m_numSamps; }

  // return the number of VQ bits otherwise 8
  unsigned GetNumBitsPerGL() const { return m_numBitsPerGL; }

  // vector of size bits per sample * 2 (hom ref, alt)
  // gl of hom alt is max(1 - hom ref - alt,1)
  std::vector<pair<float, float> > GetCodeBook() const {
    return m_VQ.GetCodeBook();
  }
};

#endif
