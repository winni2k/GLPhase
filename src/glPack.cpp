#include "glPack.hpp"

using namespace std;

GLPack::GLPack(GLPackHelper::Init &init)
    : m_numSamps(init.numSamps), m_inGLs(init.inGLs),
      m_sampleStride(init.sampleStride),
      m_numSites(m_inGLs.size() / m_numSamps / 3), m_useVQ(init.useVQ),
      m_numBitsPerGL(init.numBitsPerGL), m_VQ(init.inGLs) {

  // make sure input GLs have right dimension
  assert(m_inGLs.size() % 3 == 0);
  assert((m_inGLs.size() / 3) % m_numSamps == 0);

  // make sure sampleStride is valid
  assert(m_sampleStride > 0);
  assert(m_sampleStride <= m_numSamps);

  // because we are using uchar4, we need to make sure that every uchar4 is
  // filled
  assert(m_numSites % 4 == 0);

  if (m_useVQ) {
    assert(m_numBitsPerGL > 0);
    assert(m_numBitsPerGL <= 8);
    assert(8 % m_numBitsPerGL == 0);
  } else
    assert(m_numBitsPerGL == 8);
  m_numGLsPeruchar4 = 4 * 8 / m_numBitsPerGL;
  assert(m_numSites % m_numGLsPeruchar4 == 0);
}

vector<uchar4> GLPack::GetPackedGLs() {

  // Create appropriately strided GLs
  m_lastSampIdx = m_nextSampIdx + m_sampleStride <= m_numSamps
                      ? m_nextSampIdx + m_sampleStride - 1
                      : m_numSamps - 1;

  const unsigned numSamps = m_lastSampIdx - m_nextSampIdx + 1;

  // pull out the GLs we are about to work on and repackage them in chars
  vector<uchar4> packedGLs;
  packedGLs.reserve(numSamps * m_numSites);
  for (unsigned siteIdx = 0; siteIdx < m_numSites; siteIdx += 4) {
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
      vector<unsigned> sampGLIdxs;
      sampGLIdxs.reserve(m_numGLsPeruchar4);
      for (unsigned glNum = 0; glNum < m_numGLsPeruchar4; ++glNum)
        sampGLIdxs.push_back((siteIdx + glNum) * 3 + sampIdx * 3 * m_numSites);
      packedGLs.push_back(GLidxs2uchar4(sampGLIdxs);
    }
  }

  // move next sampId along
  m_nextSampIdx =
      m_lastSampIdx == m_numSamps - 1 ? 0 : m_nextSampIdx + m_sampleStride;
  return packedGLs;
}

uchar4 GLidxs2uchar4(vector<unsigned> &glIxs) {

  uchar4 gls;
  assert(glIdxs.size() == m_numGLsPeruchar4)
  if (!m_useVQ) {
    assert(glIdxs.size() == 4)

    gls.x = GLTrio2Char(glIdxs[0]);
    gls.y = GLTrio2Char(glIdxs[1]);
    gls.z = GLTrio2Char(glIdxs[2]);
    gls.w = GLTrio2Char(glIdxs[3]);
  } else {
    gls.x = GLs2VQChar(glIdxs.data());
    gls.y = GLs2VQChar(glIdxs.data() + m_numGLsPeruchar4 / 4);
    gls.z = GLs2VQChar(glIdxs.data() + 2 * m_numGLsPeruchar4 / 4);
    gls.w = GLs2VQChar(glIdxs.data() + 3 * m_numGLsPeruchar4 / 4);
  }
  return gls;
}

char FindGLCode(float RR, float Het) {

    for(
}
// convert gl trios into VQ code
char GLs2VQChar(vector<unsigned> *glIdx) {
  assert(glIdx);
  char gls;
  for (unsigned glNum = 0; glNum < m_numGLsPeruchar4 / 4; ++glNum) {

    // pulling out the GLs of interest
    float GLRR = m_inGLs[*glIdx];
    float GLHet = m_inGLs[*(++glIdx)];
    ++glIdx;

    // find the code for those two GLs
    char glCode = m_VQ.FindGLCode(GLRR, GLHet);
  }
  return gls;
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
