/* @(#)bio.cpp
 */

#include "bio.hpp"

using namespace std;
using namespace Bio;

namespace Bio {
typedef std::vector<std::pair<char, char> > compBase_t;
compBase_t ChipSite::s_complementBases;

ChipSite::ChipSite(std::string c, unsigned p, char ina0, char ina1, char strand)
    : m_chr{ c }, m_pos{ p } {

  // initialize static member variable
  if (s_complementBases.empty()) {
    s_complementBases.push_back(make_pair('A', 'T'));
    s_complementBases.push_back(make_pair('C', 'G'));
    s_complementBases.push_back(make_pair('G', 'C'));
    s_complementBases.push_back(make_pair('T', 'A'));
    assert(is_sorted(s_complementBases.begin(), s_complementBases.end()));
  }

  assert(ina0 >= 'A' && ina0 <= 'T');
  assert(ina1 >= 'A' && ina1 <= 'T');
  if (strand == '-') {
    auto got = lower_bound(s_complementBases.begin(), s_complementBases.end(),
                           ina0, complementBasesLowerBoundComp());
    assert(got->first == ina0);
    m_a0 = got->second;
    got = lower_bound(s_complementBases.begin(), s_complementBases.end(), ina1,
                      complementBasesLowerBoundComp());

    assert(got->first == ina1);
    m_a1 = got->second;
  } else if (strand == '+') {
    m_a0 = ina0;
    m_a1 = ina1;
  } else
    throw std::runtime_error("[Bio::ChipSite] Unexpected strand char: " +
                             strand);
  // always store the smaller char first
  if (m_a0 > m_a1)
    std::swap(m_a0, m_a1);
};

bool ChipSite::operator<(const Bio::ChipSite &rhs) const {
  if (m_chr < rhs.chr())
    return true;
  if (m_chr > rhs.chr())
    return false;
  if (m_pos < rhs.pos())
    return true;
  if (m_pos > rhs.pos())
    return false;
  if (m_a0 < rhs.a0())
    return true;
  if (m_a0 > rhs.a0())
    return false;
  if (m_a1 < rhs.a1())
    return true;
  if (m_a1 > rhs.a1())
    return false;
  return false;
}
}
