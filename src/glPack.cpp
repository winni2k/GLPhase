#include "glPack.hpp"

using namespace std;

GLPack::GLPack(const vector<float> &inGLs, unsigned numSamps,
               unsigned sampleStride)
    : m_numSamps{ numSamps }, m_inGLs(inGLs), m_sampleStride{ sampleStride },
      m_numSites(m_inGLs.size() / numSamps / 3) {
  assert(m_inGLs.size() % 3 == 0);
  assert((m_inGLs.size() / 3) % numSamps == 0);
}

vector<char> GLPack::GetPackedGLs() {

  vector<char> packedGLs;

  // Create appropriately strided GLs
  m_lastSampIdx = m_nextSampIdx + m_sampleStride <= m_numSamps
                      ? m_nextSampIdx + m_sampleStride -1
                      : m_numSamps - 1;

  const unsigned numSamps = m_lastSampIdx - m_nextSampIdx + 1;

  // pull out the GLs we are about to work on and repackage them in chars
  packedGLs.resize(numSamps * m_numSites);
  for (unsigned siteIdx = 0; siteIdx < m_numSites; ++siteIdx) {
    unsigned sampIdx = m_nextSampIdx;
    for (unsigned packSampIdx = 0; packSampIdx < numSamps;
         ++sampIdx, ++packSampIdx) {
      packedGLs[packSampIdx + numSamps * siteIdx] =
          GLTrio2Char(m_inGLs, sampIdx * 3 * m_numSites + siteIdx * 3);
    }
  }

  return packedGLs;
}

char GLPack::GLTrio2Char(const vector<float> &GLs, unsigned idx) {

  // make sure I got the right set of GLs!
  if (!(abs(1 - GLs[idx] - GLs[idx + 1] - GLs[idx + 2]) < 0.01))
    throw myException("GLs do not add up to 1\nSample index: " +
                      sutils::uint2str(idx) + "\nGLs: " +
                      sutils::double2str(GLs[idx]) + " " +
                      sutils::double2str(GLs[idx + 1]) + " " +
                      sutils::double2str(GLs[idx + 2]));
  char GLRR = GL2HalfChar(GLs[idx]);
  char GLHet = GL2HalfChar(GLs[idx + 1]);
  return ((GLRR << 4) | GLHet);
}

char GLPack::GL2HalfChar(float GL) {

  // make sure the float is within [0,1)
  assert(GL >= 0);
  GL = max(GL - numeric_limits<float>::min(), 0.0f);
  assert(GL < 1);

  // mash that number into one of numPoss possibilities
  const unsigned numPoss = 16;
  unsigned numerator = floor(GL * numPoss);
  return static_cast<char>(numerator);
}
