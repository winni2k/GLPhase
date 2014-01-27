#include "hapPanel.hpp"

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

void HapPanel::Init(vector<vector<char>> & inHaps){

    assert(m_initialized == false);
    m_numSites = inHaps.size();
    assert(m_numSites > 0);
    assert(inHaps[0].size() > 0);
    assert(m_IDs.size() > 0);
    m_numHaps = inHaps[0].size();
    m_numWordsPerHap = static_cast<unsigned>(ceil(static_cast<double>(inHaps.size()) / m_wordSize));
    vector< uint64_t > storeHaps = Char2BitVec(inHaps, m_numWordsPerHap, m_wordSize);
    std::swap(m_haps, storeHaps);

    assert(m_haps.size() % m_numHaps == 0);
    m_initialized = true;

}

// converts a vector of vectors of chars (0,1) to a bitvector
vector< uint64_t > HapPanel::Char2BitVec(const vector<vector<char> > & inHaps,
                                         unsigned numWords, unsigned wordSize)
{

    assert(ceil(static_cast<double>(inHaps.size()) / wordSize) == numWords);

    unsigned numSites = inHaps.size();
    unsigned numHaps = inHaps[0].size();
    vector<uint64_t> storeHaps;
    storeHaps.resize(numWords * numHaps);

    for (unsigned siteNum = 0; siteNum < numSites; siteNum++) {
        for (unsigned hapNum = 0; hapNum < numHaps; hapNum ++) {
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
