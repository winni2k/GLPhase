#include "bio.hpp"

using namespace std;
namespace Bio {

Region::Region(const string &region) {
  if (region.empty())
    return;
  regex rgx("(\\S+):(\\d*)-(\\d*)");
  smatch sm;
  regex_match(region, sm, rgx);
  if (sm.empty())
    throw runtime_error("Could not parse region: [" + region + "]");
  assert(sm.size() == 4);
  m_chrom = sm[1].str();

  // allow character sequence to be empty
  if (sm[2].length() != 0)
    m_startBP = stoul(sm[2].str());
  if (m_startBP == 0)
    throw std::range_error("0 is not a valid genomic coordinate");
  if (sm[3].length() != 0)
    m_endBP = stoul(sm[3].str());
}

// this returns a vector filled with the first n_max_tokens
std::vector<std::string> tokenize_partial(const std::string &str,
                                          size_t n_max_tokens,
                                          const std::string &search) {
  std::vector<std::string> tokens;
  tokens.reserve(n_max_tokens);
  std::string::size_type p_last = str.find_first_not_of(search, 0);
  std::string::size_type p_curr = str.find_first_of(search, p_last);

  if (n_max_tokens > 1)
    while ((std::string::npos != p_curr || std::string::npos != p_last) &&
           tokens.size() + 1 != n_max_tokens) {
      tokens.push_back(str.substr(p_last, p_curr - p_last));
      p_last = str.find_first_not_of(search, p_curr);
      p_curr = str.find_first_of(search, p_last);
    }

  tokens.push_back(str.substr(p_last, p_curr - p_last));
  return tokens;
}

void snp_storage::push_back(Bio::snp input) {
  std::string inStr = input.to_string();
  if (m_site_set.count(inStr) > 0)
    throw std::runtime_error("Cannot push SNP, it already exists: [" + inStr +
                             "]");
  m_site_set.insert(make_pair(inStr, m_sites.size()));
  m_sites.push_back(std::move(input));
}

size_t snp_storage::index(const Bio::snp &search) const {
  auto site = m_site_set.find(search.to_string());
  if (site == m_site_set.end())
    throw std::runtime_error("Could not find site [" + search.to_string() +
                             "] in site set");
  return site->second;
}

void snp_storage_ordered::push_back(Bio::snp input) {
  if (!m_sites.empty() && input.pos < m_sites.back().pos)
    throw std::range_error("Tried to insert position [" +
                           std::to_string(input.pos) +
                           "] that was less than last position [" +
                           std::to_string(m_sites.back().pos) + "]");
  if (m_sites.empty())
    m_chrom = input.chr;
  else if (m_chrom != input.chr)
    throw std::range_error("Input snp is wrong chromosome [" + input.chr + "]");
  Bio::snp_storage::push_back(std::move(input));
}
}
