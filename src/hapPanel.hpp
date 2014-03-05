/* @(#)hapPanel.hpp
 */

#ifndef _HAPPANEL_H
#define _HAPPANEL_H 1

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

#include <math.h>
#include <vector>
#include <stdint.h>
#include <string>
#include <cassert>
#include <iostream>
#include "snp.hpp"
#include "utils.hpp"

#define WordShift 6
#define WordMod 63

class HapPanel {

private:
  unsigned m_wordSize = 64;
  std::vector<uint64_t> m_haps;
  std::vector<snp> m_sites;
  std::vector<std::string> m_sampleIDs;
  unsigned m_numHaps = 0;
  unsigned m_numWordsPerHap;
  bool m_initialized = false;

  void set1(uint64_t *P, unsigned I) {

    // I >> Uint64_TShift is moving along the array according to which uint64_t
    // I is in
    // e.g. I <= 63 is first uint64_t, and so forth in blocks of 64 bits
    P[I >> WordShift] |= static_cast<uint64_t>(1) << (I & WordMod);
  }

  void set0(uint64_t *P, unsigned I) {
    P[I >> WordShift] &= ~(static_cast<uint64_t>(1) << (I & WordMod));
  }

public:
  void Init(std::vector<std::vector<char> > &inHaps, std::vector<snp> &inSites,
            std::vector<std::string> &inSampleIDs);

  std::string GetID(unsigned idx) {
    assert(idx < m_sampleIDs.size());
    return m_sampleIDs[idx];
  }
  unsigned NumHaps() {
    assert(m_initialized);
    return m_numHaps;
  }
  unsigned NumWordsPerHap() {
    assert(m_initialized);
    return m_numWordsPerHap;
  }
  std::vector<uint64_t> *Haplotypes() {
    assert(m_initialized);
    return &m_haps;
  }
  uint64_t *Hap(unsigned hapNum) {
    assert(hapNum < m_numHaps);
    return &m_haps[hapNum * m_numWordsPerHap];
  }
  unsigned MaxSites() {
    assert(Initialized());
    return NumWordsPerHap() * m_wordSize;
  }
  unsigned NumSites() {
    assert(m_initialized);
    return m_sites.size();
  }
  bool Initialized() { return m_initialized; }
  std::vector<uint64_t>
  Char2BitVec(const std::vector<std::vector<char> > &inHaps, unsigned numWords,
              unsigned wordSize);

  unsigned Position(unsigned idx) { return m_sites[idx].pos; }
};

#endif
