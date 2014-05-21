/* @(#)bio.hpp
 */

#ifndef _BIO_HPP
#define _BIO_HPP 1

#include <string>

// this was written with c++11 in mind
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

namespace Bio {

class Region {
private:
  std::string m_chrom;
  unsigned m_startBP = 0, m_endBP = 0;

public:
  Region() {};
  Region(std::string &chrom, unsigned startBP, unsigned endBP)
      : m_chrom(chrom), m_startBP(startBP), m_endBP(endBP) {}

  // returns region as tabix compatible region specifier of the form
  // chr:start-end
  std::string AsString() const {
    std::string region = m_chrom + ":" + std::to_string(m_startBP) + "-" +
                         std::to_string(m_endBP);
    return region;
  }
};

class snp {
public:
  unsigned pos;
  std::string chr, ref, alt;

  snp(std::string chr, int pos, std::string ref, std::string alt) {
    this->chr = chr;
    this->pos = pos;
    this->ref = ref;
    this->alt = alt;
  }

  bool strand(std::string &r, std::string &a) {
    if (ref == r && alt == a)
      return true;

    return false;
  }
};
}

#endif /* _BIO_HPP */
