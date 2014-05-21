#include "hapPanel.hpp"

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;
using namespace Bio;

void HapPanel::Init(vector<vector<char> > &inHaps, vector<snp> &inSites,
                    vector<string> &inSampleIDs) {

  assert(m_initialized == false);

  // assume sample IDs are OK and work from there
  std::swap(m_sampleIDs, inSampleIDs);

  if (inSites.size() != inHaps.size())
      throw myException("inSites size (" + sutils::uint2str(inSites.size()) +
                      ") is not equal to inHaps size (" +
                        sutils::uint2str(inHaps.size()) + ")");

  // fill the sites of the ref panel
  assert(inSites.size() > 0);
  std::swap(m_sites, inSites);

  // fill the haplotypes of the ref panel
  assert(inHaps[0].size() > 0);
  assert(m_sampleIDs.size() > 0);
  m_numHaps = inHaps[0].size();
  m_numWordsPerHap = static_cast<unsigned>(
      ceil(static_cast<double>(inHaps.size()) / m_wordSize));
  vector<uint64_t> storeHaps =
      Char2BitVec(inHaps, m_numWordsPerHap, m_wordSize);
  std::swap(m_haps, storeHaps);

  assert(m_haps.size() % m_numHaps == 0);
  m_initialized = true;
}

// converts a vector of vectors of chars (0,1) to a bitvector
vector<uint64_t> HapPanel::Char2BitVec(const vector<vector<char> > &inHaps,
                                       unsigned numWords, unsigned wordSize) {

  assert(ceil(static_cast<double>(inHaps.size()) / wordSize) == numWords);

  unsigned numSites = inHaps.size();
  unsigned numHaps = inHaps[0].size();
  vector<uint64_t> storeHaps;
  storeHaps.resize(numWords * numHaps);

  for (unsigned siteNum = 0; siteNum < numSites; siteNum++) {
    for (unsigned hapNum = 0; hapNum < numHaps; hapNum++) {
      if (inHaps[siteNum][hapNum] == 0)
        set0(&storeHaps[hapNum * numWords], siteNum);
      else if (inHaps[siteNum][hapNum] == 1)
        set1(&storeHaps[hapNum * numWords], siteNum);
      else {
        cout << "programming error";
        throw 1;
      }
    }
  }
  return storeHaps;
}
