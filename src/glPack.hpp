#ifndef _GLPACK_HPP
#define _GLPACK_HPP

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

#include <vector>
#include <cassert>
#include <cmath>
#include <limits>
#include "utils.hpp"
#include "glVQ.hpp"
#include <limits>
#include <utility>
#include <gsl/gsl_rng.h>
#include "globals.h"

#ifndef NCUDA
#include <vector_types.h>
#endif

// we want our chars to be unsigned
namespace GLPackHelper {

struct Init {
  const std::vector<float> &inGLs;
  gsl_rng &rng;
  unsigned numSamps = 0;
  unsigned sampleStride = 0;
  bool useVQ = false;

  // 32 % numBits needs to be zero
  unsigned numBitsPerGL = 8;

  Init(const std::vector<float> &GLs, gsl_rng &rng) : inGLs(GLs), rng(rng) {};
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
  unsigned m_numGLsPeruint32_t;

  unsigned m_nextSampIdx = 0;
  unsigned m_lastSampIdx = 0;
  GLVQ m_VQ;

  uint32_t GLIdxs2uint32_t(const std::vector<unsigned> &glIdxs) const;
  unsigned char GLs2VQChar(const std::vector<unsigned> &glIdxs, size_t glIdx,
                           size_t numGLsPerChar) const;
  unsigned char GLTrio2Char(unsigned idx) const;
  unsigned char GL2HalfChar(float GL) const;

public:
  // set VQBits > 0 for vector quantization
  GLPack(const GLPackHelper::Init &init);

  std::vector<uint32_t> GetPackedGLs();
  unsigned GetSampleStride() const { return m_sampleStride; }
  unsigned GetNextSampIdx() const { return m_nextSampIdx; }
  unsigned GetLastSampIdx() const { return m_lastSampIdx; }
  unsigned GetNumSites() const { return m_numSites; }
  unsigned GetNumSamps() const { return m_numSamps; }

  // return the number of VQ bits otherwise 8
  unsigned GetNumBitsPerGL() const { return m_numBitsPerGL; }

  // vector of size bits per sample * 2 (hom ref, alt)
  // gl of hom alt is max(1 - hom ref - alt,1)
  std::vector<std::pair<float, float> > GetCodeBook() const {
    return m_VQ.GetCodeBook();
  }
};

#endif
