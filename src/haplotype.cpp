#include "haplotype.hpp"

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

void Haplotype::Set(unsigned uSite, bool bAllele) {

  assert(uSite < m_uNumAlleles);

  // set uSite to 1
  if (bAllele)
    m_hap[uSite >> WORDSHIFT] |= static_cast<uint64_t>(1)
                                   << (uSite & WORDMOD);
  // set uSite to 0
  else
    m_hap[uSite >> WORDSHIFT] &=
        ~(static_cast<uint64_t>(1) << (uSite & WORDMOD));
}

unsigned Haplotype::HammingDist(const Haplotype &oCompareHap) const {

  assert(oCompareHap.m_uNumAlleles == m_uNumAlleles);

  unsigned uHammingDist = 0;
  for (unsigned uWord = 0; uWord < m_hap.size(); uWord++) {
    bitset<s_wordSize> diff = m_hap[uWord] xor oCompareHap.GetWordBitset(uWord);
    uHammingDist += diff.count();
  }

  return uHammingDist;
}

unsigned Haplotype::HammingDist(const uint64_t *upHap) const {

  unsigned uHammingDist = 0;
  for (unsigned uWord = 0; uWord < m_hap.size(); uWord++, upHap++) {
    static_assert(
        s_wordSize == 64,
        "can't use uint64_t for storage if bitset is not 64 bits in size");
    uint64_t uDiff =
        static_cast<uint64_t>(m_hap[uWord].to_ullong()) xor * upHap;
    uHammingDist += __builtin_popcountll(uDiff);
  }

  return uHammingDist;
}

unsigned Haplotype::HammingDist(const uint64_t *upHap1,
                                const uint64_t *upHap2) const {

  unsigned uHammingDist = 0;
  for (unsigned uWord = 0; uWord < m_hap.size(); uWord++, upHap1++, upHap2++) {
    uint64_t uDiff = *upHap2 xor * upHap1;
    uHammingDist += __builtin_popcountll(uDiff);
  }

  return uHammingDist;
}

unsigned Haplotype::MaxTractLen(const Haplotype &compHap) const {

  assert(compHap.m_uNumAlleles == m_uNumAlleles);

  unsigned maxTracLen = 0;
  unsigned currTracL = 0;
  for (unsigned uWord = 0; uWord < m_hap.size(); uWord++) {

    // number of sites in word may be less than s_wordSize if we are in the last
    // word
    unsigned numSitesInWord = uWord == m_hap.size() - 1
                                  ? m_uNumAlleles - uWord * s_wordSize
                                  : s_wordSize;
    bitset<s_wordSize> diff = m_hap[uWord] xor compHap.GetWordBitset(uWord);

    for (unsigned siteNum = 0; siteNum < numSitesInWord; ++siteNum) {
      if (diff[siteNum] == 1) {
        if (currTracL > maxTracLen)
          maxTracLen = currTracL;
        currTracL = 0;
      } else
        ++currTracL;
    }
  }
  if (currTracL > maxTracLen)
    maxTracLen = currTracL;

  return maxTracLen;
}
