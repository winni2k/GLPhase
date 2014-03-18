#include "haplotype.hpp"

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

void Haplotype::Set(unsigned uSite, bool bAllele) {

  assert(uSite < m_uNumAlleles);

  // set uSite to 1
  if (bAllele)
    m_utHap[uSite >> h_WordShift] |= static_cast<uint64_t>(1)
                                     << (uSite & h_WordMod);
  // set uSite to 0
  else
    m_utHap[uSite >> h_WordShift] &=
        ~(static_cast<uint64_t>(1) << (uSite & h_WordMod));
}

unsigned Haplotype::HammingDist(const Haplotype &oCompareHap) const {

  assert(oCompareHap.m_uNumAlleles == m_uNumAlleles);

  unsigned uHammingDist = 0;
  for (unsigned uWord = 0; uWord < m_utHap.size(); uWord++) {
    uint64_t uDiff = m_utHap[uWord] xor oCompareHap.GetWord(uWord);
    uHammingDist += __builtin_popcountll(uDiff);
  }

  return uHammingDist;
}

unsigned Haplotype::HammingDist(const uint64_t *upHap) const {

  unsigned uHammingDist = 0;
  for (unsigned uWord = 0; uWord < m_utHap.size(); uWord++, upHap++) {
    uint64_t uDiff = m_utHap[uWord] xor * upHap;
    uHammingDist += __builtin_popcountll(uDiff);
  }

  return uHammingDist;
}

unsigned Haplotype::HammingDist(const uint64_t *upHap1,
                                const uint64_t *upHap2) const {

  unsigned uHammingDist = 0;
  for (unsigned uWord = 0; uWord < m_utHap.size();
       uWord++, upHap1++, upHap2++) {
    uint64_t uDiff = *upHap2 xor * upHap1;
    uHammingDist += __builtin_popcountll(uDiff);
  }

  return uHammingDist;
}

unsigned Haplotype::MaxTractLen(const Haplotype &compHap) const {

  assert(compHap.m_uNumAlleles == m_uNumAlleles);

  unsigned maxTracLen = 0;
  unsigned currTracL = 0;
  for (unsigned uWord = 0; uWord < m_utHap.size(); uWord++) {
    uint64_t diff = m_utHap[uWord] xor oCompareHap.GetWord(uWord);
    for (unsigned siteNum = 0; siteNum <= h_WordMod; ++siteNum) {
      if (diff & static_cast<uint64_t>(1)) {
        if (currTracL > maxTracLen)
          maxTracLen = currTracL;
        currTracL = 0;
      } else
        ++currTracL;
      diff >> 1;
    }
  }
  if (currTracL > maxTracLen)
    maxTracLen = currTracL;

  return maxTracLen;
}
