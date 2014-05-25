#include "glPack.hpp"

using namespace std;

GLPack::GLPack(const GLPackHelper::Init &init)
    : m_numSamps(init.numSamps), m_inGLs(init.inGLs),
      m_sampleStride(init.sampleStride),
      m_numSites(m_inGLs.size() / m_numSamps / 3), m_useVQ(init.useVQ),
      m_numBitsPerGL(init.numBitsPerGL),
      m_numGLsPeruint32_t(UINT32T_SIZE / m_numBitsPerGL),
      m_VQ(m_inGLs, init.rng, pow(2, m_numBitsPerGL)) {

  assert(m_numSamps > 0);
  assert(m_sampleStride > 0);

  // make sure input GLs have right dimension
  assert(m_inGLs.size() % 3 == 0);
  assert((m_inGLs.size() / 3) % m_numSamps == 0);

  // make sure sampleStride is valid
  assert(m_sampleStride > 0);
  assert(m_sampleStride <= m_numSamps);

  if (m_useVQ) {
    assert(m_numBitsPerGL > 0);
    assert(m_numBitsPerGL <= CHAR_BIT);
    assert(8 % m_numBitsPerGL == 0);
  } else
    assert(m_numBitsPerGL == 8);
}

vector<uint32_t> GLPack::GetPackedGLs() {

  // Create appropriately strided GLs
  m_lastSampIdx = m_nextSampIdx + m_sampleStride <= m_numSamps
                      ? m_nextSampIdx + m_sampleStride - 1
                      : m_numSamps - 1;

  const unsigned numSamps = m_lastSampIdx - m_nextSampIdx + 1;

  // pull out the GLs we are about to work on and repackage them in uint32s
  vector<uint32_t> packedGLs;
  packedGLs.reserve(max(numSamps, numSamps * m_numSites / m_numGLsPeruint32_t));
  for (unsigned siteIdx = 0; siteIdx < m_numSites; siteIdx += m_numGLsPeruint32_t) {
    unsigned sampIdx = m_nextSampIdx;
    for (unsigned packSampIdx = 0; packSampIdx < numSamps;
         ++sampIdx, ++packSampIdx) {

      /*
        We are rearranging the order of GLs here.
        inGLs, which GLTrio2Char extracts from, has its GLs ordered by sample:
        samp1site1, samp1site2, samp1site3, samp2site1, etc.
        packegGLs has its GLs ordered by site:
        samp1site1, samp2site1, samp3site1, samp1site2, etc.
        however, we are doing this in chunks of sites at a time as well...
      */
      vector<unsigned> sampGLIdxs;
      sampGLIdxs.reserve(m_numGLsPeruint32_t);
      for (unsigned glNum = 0; glNum < m_numGLsPeruint32_t; ++glNum) {

        // just use the first gl if the GL we want does not exist
        if (siteIdx + glNum >= m_numSites)
          sampGLIdxs.push_back(0);
        else
          sampGLIdxs.push_back((siteIdx + glNum) * 3 +
                               sampIdx * 3 * m_numSites);
      }
      assert(sampGLIdxs.size() == m_numGLsPeruint32_t);
      packedGLs.push_back(GLIdxs2uint32_t(sampGLIdxs));
    }
  }

  // move next sampId along
  m_nextSampIdx =
      m_lastSampIdx == m_numSamps - 1 ? 0 : m_nextSampIdx + m_sampleStride;

  assert(packedGLs.size() ==
         max(numSamps, numSamps * m_numSites / m_numGLsPeruint32_t));
  return packedGLs;
}

uint32_t GLPack::GLIdxs2uint32_t(const vector<unsigned> &glIdxs) const {

  uint32_t gls = 0;
  assert(glIdxs.size() == m_numGLsPeruint32_t);
  if (!m_useVQ) {
    assert(glIdxs.size() == 4);

    gls ^= GLTrio2Char(glIdxs[0]) << (UINT32T_SIZE - CHAR_BIT);
    gls ^= GLTrio2Char(glIdxs[1]) << (UINT32T_SIZE - 2 * CHAR_BIT);
    gls ^= GLTrio2Char(glIdxs[2]) << (UINT32T_SIZE - 3 * CHAR_BIT);
    gls ^= GLTrio2Char(glIdxs[3]);
  } else {
    assert(m_numGLsPeruint32_t % (UINT32T_SIZE / CHAR_BIT) == 0);
    const size_t numGLsPerChar =
        m_numGLsPeruint32_t / (UINT32T_SIZE / CHAR_BIT);
    // pack a set of gl indices equal to the number of GLs that will fit into a
    // char into a char, and then move the char to the right place
    gls ^= GLs2VQChar(glIdxs, 0, numGLsPerChar) << (UINT32T_SIZE - CHAR_BIT);
    gls ^= GLs2VQChar(glIdxs, numGLsPerChar, numGLsPerChar)
           << (UINT32T_SIZE - CHAR_BIT * 2);
    gls ^= GLs2VQChar(glIdxs, numGLsPerChar * 2, numGLsPerChar)
           << (UINT32T_SIZE - CHAR_BIT * 3);
    gls ^= GLs2VQChar(glIdxs, numGLsPerChar * 3, numGLsPerChar);
  }
  return gls;
}

// convert gl trios into VQ code
unsigned char GLPack::GLs2VQChar(const vector<unsigned> &glIdxs, size_t glIdx,
                                 size_t numGLsPerChar) const {
  unsigned char gls = 0;
  //  const unsigned numGLsPerSite = 3;
  assert(glIdxs.size() >= glIdx + numGLsPerChar);
  for (size_t glNum = 0; glNum < numGLsPerChar; ++glNum) {

    // make space for the next GL set
    gls <<= m_numBitsPerGL;

    // pulling out the GLs of interest
    float GLRR = m_inGLs[glIdxs[glIdx + glNum]];
    float GLHet = m_inGLs[glIdxs[glIdx + glNum] + 1];

    // don't forget to normalize to 0... silly me.
    float sum = GLRR + GLHet + m_inGLs[glIdxs[glIdx + glNum] + 2];
    GLRR /= sum;
    GLHet /= sum;
    
    // find the code for those two GLs
    unsigned char glCode = m_VQ.FindGLCode(GLRR, GLHet);
    gls ^= glCode;
  }
  return gls;
}

unsigned char GLPack::GLTrio2Char(unsigned glIdx) const {

  // normalize GLs to 1
  float sum = m_inGLs[glIdx] + m_inGLs[glIdx + 1] + m_inGLs[glIdx + 2];
  assert(sum > 0);

  unsigned char GLRR = GL2HalfChar(m_inGLs[glIdx] / sum);
  unsigned char GLHet = GL2HalfChar(m_inGLs[glIdx + 1] / sum);
  return ((GLRR << 4) | GLHet);
}

unsigned char GLPack::GL2HalfChar(float GL) const {

  // make sure the float is within [0,1)
  assert(GL >= 0);
  GL = max(GL - numeric_limits<float>::epsilon(), 0.0f);
  assert(GL < 1);

  // mash that number into one of numPoss possibilities
  const unsigned numPoss = 16;
  unsigned numerator = floor(GL * numPoss);
  return static_cast<char>(numerator);
}
