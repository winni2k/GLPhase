/* @(#)bio.hpp
 */

#ifndef _BIO_HPP
#define _BIO_HPP 1

#include <cassert>
#include <cmath>
#include <limits>
#include <regex>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

// this was written with c++11 in mind
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

namespace Bio {

class Region {
private:
  std::string m_chrom;
  unsigned m_startBP = 0, m_endBP = 0;

public:
  Region(){};
  Region(std::string chrom, unsigned startBP, unsigned endBP)
      : m_chrom(std::move(chrom)), m_startBP(startBP), m_endBP(endBP) {}
  Region(const std::string &region);

  // returns region as tabix compatible region specifier of the form
  // chr:start-end
  // start and end may be omitted
  std::string AsString() const {
    std::string out = m_chrom + ":";
    if (m_startBP)
      out += std::to_string(m_startBP);
    out += "-";
    if (m_endBP)
      out += std::to_string(m_endBP);
    return out;
  }

  std::string Chrom() const { return m_chrom; }

  // return sensible values if the BPs are 0 (i.e. undefined)
  unsigned startBP() const {
    if (m_startBP)
      return m_startBP;
    else
      return 1;
  }
  unsigned endBP() const {
    if (m_endBP)
      return m_endBP;
    else
      return std::numeric_limits<unsigned>::max();
  }
  bool empty() const { return m_chrom.empty(); }
  void clear() {
    m_chrom = "";
    m_startBP = m_endBP = 0;
  }
  bool chrom_eq(const std::string &comp) const { return comp == m_chrom; }
};

class snp {
public:
  unsigned pos;
  std::string chr, ref, alt;

  snp(std::string chr, unsigned pos, std::string ref, std::string alt) {
    this->chr = chr;
    this->pos = pos;
    this->ref = ref;
    this->alt = alt;
  }

  std::string to_string() const {
    return std::string(chr + ":" + std::to_string(pos) + ":" + ref + ":" + alt);
  }

  bool strand(std::string &r, std::string &a) {
    if (ref == r && alt == a)
      return true;

    return false;
  }

  bool operator==(const snp &rhs) const {
    if (pos == rhs.pos && chr == rhs.chr && ref == rhs.ref && alt == rhs.alt)
      return true;
    return false;
  }
};

// this returns a vector filled with the first n_max_tokens
std::vector<std::string> tokenize_partial(const std::string &str,
                                          size_t n_max_tokens,
                                          const std::string &search);

class snp_storage {
protected:
  std::vector<Bio::snp> m_sites;
  // string representation of SNP -> index of SNP in m_sites
  std::unordered_map<std::string, size_t> m_site_set;

public:
  void push_back(Bio::snp input);
  size_t size() const { return m_sites.size(); }
  bool empty() const { return m_sites.empty(); }
  bool exists(const Bio::snp &search) const {
    return m_site_set.count(search.to_string()) > 0;
  }
  size_t index(const Bio::snp &search) const;

  void clear() {
    m_sites.clear();
    m_site_set.clear();
  }

  std::vector<Bio::snp>::const_iterator at(size_t idx) const {
    if (!(idx < m_sites.size()))
      throw std::range_error("Idx [" + std::to_string(idx) + "] out of range");
    return m_sites.cbegin() + idx;
  }
};

class snp_storage_ordered : public Bio::snp_storage {
private:
  std::string m_chrom;

public:
  void clear() {
    m_chrom.clear();
    Bio::snp_storage::clear();
  }
  void push_back(Bio::snp input);
  std::pair<unsigned, unsigned> pos_range() const {
    if (m_sites.empty())
      throw std::range_error("Container is empty");
    return std::make_pair(m_sites.front().pos, m_sites.back().pos);
  }

  std::string chrom() const { return m_chrom; }
  bool eq_chrom(const std::string chrom) { return chrom == m_chrom; }
};

// template functions
template <class F> F gl2prob(F val) {
  static_assert(std::is_floating_point<F>::value,
                "F needs to be a floating point type");
  F tran = std::pow(static_cast<F>(10), val);
  tran = std::max(tran, static_cast<F>(0));
  tran = std::min(tran, static_cast<F>(1));
  return tran;
}

template <class F, class I> F phred2prob(I phred) {
  if (phred > std::numeric_limits<F>::min_exponent10 * -10)
    return static_cast<F>(0);
  F ret = static_cast<F>(pow(10, -phred / 10));
  assert(ret <= 1);
  assert(ret >= 0);
  return ret;
}
}

#endif /* _BIO_HPP */
