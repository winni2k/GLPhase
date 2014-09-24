/* @(#)bio.hpp
 */

#ifndef _BIO_HPP
#define _BIO_HPP 1

#include <string>
#include <boost/functional/hash.hpp>

// this was written with c++11 in mind
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

namespace Bio {

class Region {
private:
  std::string m_chrom;
  unsigned m_startBP = 0, m_endBP = 0;

public:
  Region() {};
  Region(std::string chrom, unsigned startBP, unsigned endBP)
      : m_chrom(std::move(chrom)), m_startBP(startBP), m_endBP(endBP) {}

  // returns region as tabix compatible region specifier of the form
  // chr:start-end
  std::string AsString() const {
    std::string region = m_chrom + ":" + std::to_string(m_startBP) + "-" +
                         std::to_string(m_endBP);
    return region;
  }

  std::string GetChrom() const { return m_chrom; }
  unsigned GetStartBP() const { return m_startBP; }
  unsigned GetEndBP() const { return m_endBP; }
  bool empty() const { return m_chrom.empty(); }
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

  bool operator==(const snp &rhs) const {
    return (chr == rhs.chr && pos == rhs.pos && ref == rhs.ref &&
            alt == rhs.alt);
  }

  bool strand(std::string &r, std::string &a) {
    if (ref == r && alt == a)
      return true;

    return false;
  }
};

struct snpKeyHasher {
  std::size_t operator()(const Bio::snp &k) const {
    using boost::hash_value;
    using boost::hash_combine;

    // Start with a hash value of 0    .
    std::size_t seed = 0;

    // Modify 'seed' by XORing and bit-shifting in
    // one member of 'Key' after the other:
    hash_combine(seed, hash_value(k.chr));
    hash_combine(seed, hash_value(k.pos));
    hash_combine(seed, hash_value(k.ref));
    hash_combine(seed, hash_value(k.alt));

    // Return the result.
    return seed;
  }
};

struct snpPosComp {
  bool operator()(const Bio::snp &a, const Bio::snp &b) {
    return a.pos < b.pos;
  }
};
}

#endif /* _BIO_HPP */
