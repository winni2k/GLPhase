/* @(#)haplotype.hpp
 */

#ifndef _HAPLOTYPE_H
#define _HAPLOTYPE_H 1

#include <vector>
#include <stdint.h>
#include <cmath>
#include <assert.h>
#include <bitset>
#include <iostream>

// require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

#define h_WordMod 63  // word size -1
#define h_WordShift 6 // 2^6 = 64

class Haplotype {
public:
  static constexpr unsigned s_wordSize = std::pow(2, h_WordShift);

private:
  std::vector<std::bitset<s_wordSize> > m_hap;

public:
  unsigned m_uNumAlleles;

  // initialize all haps to all 0s
  Haplotype(unsigned uNumAlleles) : m_uNumAlleles(uNumAlleles) {
    unsigned uSize =
        ceil(static_cast<float>(uNumAlleles) / static_cast<float>(s_wordSize));
    m_hap.reserve(uSize);
    for (unsigned uWord = 0; uWord < uSize; uWord++)
      m_hap.push_back(std::bitset<s_wordSize>());
  }

  // site is 0 based site number to be set
  // allele is either 1 or 0 and is the bit to set the site to
  void Set(unsigned uSite, bool bAllele);

  bool TestSite(unsigned uSite) const {
    assert(uSite < m_uNumAlleles);

    // need to rework all of this
    unsigned wordIdx = uSite >> h_WordShift;
    unsigned siteIdx = uSite - s_wordSize * wordIdx;
    return m_hap[wordIdx][siteIdx] == 1;
  };

  // test if bit uSite is 1
  uint64_t TestSite(unsigned uSite, const uint64_t *P) const {
    return (P[uSite >> h_WordShift] >> (uSite & h_WordMod)) &
           static_cast<uint64_t>(1);
  };

  std::bitset<s_wordSize> GetWordBitset(unsigned uWordNum) const {
    return m_hap[uWordNum];
  };
  uint64_t GetWord(unsigned uWordNum) const {
    return static_cast<uint64_t>(m_hap[uWordNum].to_ullong());
  };

  // using a haplotype object
  unsigned HammingDist(const Haplotype &oCompareHap) const;

  // using a vector pointer
  unsigned HammingDist(const uint64_t *upHap) const;

  unsigned HammingDist(const uint64_t *upHap1, const uint64_t *upHap2) const;

  // calculate the max shared tract length like B.Howie
  // using a haplotype object
  unsigned MaxTractLen(const Haplotype &compHap) const;
};

#endif /* _HAPLOTYPE_H */
