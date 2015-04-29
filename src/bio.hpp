/* @(#)bio.hpp
 */

#ifndef _BIO_HPP
#define _BIO_HPP 1

#include <limits>
#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_set>
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

  // returns region as tabix compatible region specifier of the form
  // chr:start-end
  std::string AsString() const {
    std::string region = m_chrom + ":" + std::to_string(m_startBP) + "-" +
                         std::to_string(m_endBP);
    return region;
  }
  std::string Chrom() const { return m_chrom; }
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
  std::unordered_set<std::string> m_site_set;

public:
  void push_back(Bio::snp input) {
    m_site_set.insert(input.to_string());
    m_sites.push_back(std::move(input));
  }
  size_t size() const { return m_sites.size(); }
  bool empty() const { return m_sites.empty(); }
  bool exists(const Bio::snp &search) const {
    return m_site_set.find(search.to_string()) != m_site_set.end();
  }
  void clear() {
    m_sites.clear();
    m_site_set.clear();
  }
  std::vector<Bio::snp>::const_iterator at(size_t idx) const {
    if (m_sites.size() <= idx)
      throw std::range_error("Idx [" + std::to_string(idx) + "]out of range");
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
  void push_back(Bio::snp input) {
    if (!m_sites.empty() && input.pos < m_sites.back().pos)
      throw std::range_error("Tried to insert position [" +
                             std::to_string(input.pos) +
                             "] that was less than last position [" +
                             std::to_string(m_sites.back().pos) + "]");
    if (m_sites.empty())
      m_chrom = input.chr;
    else if (m_chrom != input.chr)
      throw std::range_error("Input snp is wrong chromosome [" + input.chr +
                             "]");
    Bio::snp_storage::push_back(std::move(input));
  }

  std::pair<unsigned, unsigned> pos_range() const {
    if (m_sites.empty())
      throw std::range_error("Container is empty");
    return std::make_pair(m_sites.front().pos, m_sites.back().pos);
  }
  std::string chrom() const { return m_chrom; }
  bool eq_chrom(const std::string chrom) { return chrom == m_chrom; }
};

template <class F> F glToProb(F val) {
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
