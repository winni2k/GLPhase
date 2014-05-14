#include "glPack.hpp"

using namespace std;

GLPack::GLPack(const vector<float> &inGLs, unsigned numSamps,
               unsigned sampleStride)
    : m_numSamps(numSamps), m_inGLs(inGLs), m_sampleStride(sampleStride),
      m_numSites(m_inGLs.size() / m_numSamps / 3) {

  // make sure input GLs have right dimension
  assert(m_inGLs.size() % 3 == 0);
  assert((m_inGLs.size() / 3) % numSamps == 0);

  // make sure sampleStride is valid
  assert(m_sampleStride > 0);
  assert(m_sampleStride <= m_numSamps);
}

vector<char> GLPack::GetPackedGLs() {

  vector<char> packedGLs;

  // Create appropriately strided GLs
  m_lastSampIdx = m_nextSampIdx + m_sampleStride <= m_numSamps
                      ? m_nextSampIdx + m_sampleStride - 1
                      : m_numSamps - 1;

  const unsigned numSamps = m_lastSampIdx - m_nextSampIdx + 1;

  // pull out the GLs we are about to work on and repackage them in chars
  packedGLs.resize(numSamps * m_numSites);
  for (unsigned siteIdx = 0; siteIdx < m_numSites; ++siteIdx) {
    unsigned sampIdx = m_nextSampIdx;
    for (unsigned packSampIdx = 0; packSampIdx < numSamps;
         ++sampIdx, ++packSampIdx) {

      /*
        We are rearranging the order of GLs here.
        inGLs, which GLTrio2Char extracts from, has its GLs ordered by sample:
        samp1site1, samp1site2, samp1site3, samp2site1, etc.
        packegGLs has its GLs ordered by site:
        samp1site1, samp2site1, samp3site1, samp1site2, etc.
      */
      packedGLs[packSampIdx + numSamps * siteIdx] =
          GLTrio2Char(siteIdx * 3 + sampIdx * 3 * m_numSites);
    }
  }

  // move next sampId along
  m_nextSampIdx =
      m_lastSampIdx == m_numSamps - 1 ? 0 : m_nextSampIdx + m_sampleStride;
  return packedGLs;
}

char GLPack::GLTrio2Char(unsigned glIdx) const {

  // normalize GLs to 1
  float sum = m_inGLs[glIdx] + m_inGLs[glIdx + 1] + m_inGLs[glIdx + 2];
  assert(sum > 0);

  char GLRR = GL2HalfChar(m_inGLs[glIdx] / sum);
  char GLHet = GL2HalfChar(m_inGLs[glIdx + 1] / sum);
  return ((GLRR << 4) | GLHet);
}

char GLPack::GL2HalfChar(float GL) const {

  // make sure the float is within [0,1)
  assert(GL >= 0);
  GL = max(GL - numeric_limits<float>::epsilon(), 0.0f);
  assert(GL < 1);

  // mash that number into one of numPoss possibilities
  const unsigned numPoss = 16;
  unsigned numerator = floor(GL * numPoss);
  return static_cast<char>(numerator);
}
