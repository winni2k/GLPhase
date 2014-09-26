/* @(#)hapPanel.hpp
 */

#ifndef _HAPPANEL_H
#define _HAPPANEL_H 1

static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

#include "globals.h"
#include <math.h>
#include <stdint.h>
#include <algorithm>
#include <cassert>
#include <string>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "bio.hpp"
#include "utils.hpp"
#include "tabix.hpp"
#include "vcf_parser.hpp"

namespace HapPanelHelper {

enum class HapFileType { WTCCC, IMPUTE2, VCF, WTCCC_TBX };

std::vector<uint64_t> Char2BitVec(const std::vector<std::vector<char> > &inHaps,
                                  unsigned numWords, unsigned wordSize);
void set1(uint64_t *P, unsigned I);
void set0(uint64_t *P, unsigned I);

// HAPS/SAMPLE format
std::vector<std::string> OpenSample(const std::string &sampleFile);
void OpenHaps(const std::string &hapsFile, const Bio::Region &region,
              std::vector<std::vector<char> > &loadHaps,
              std::vector<Bio::snp> &sites);
void OpenTabHaps(const std::string &hapsFile, const Bio::Region &region,
                 std::vector<std::vector<char> > &loadHaps,
                 std::vector<Bio::snp> &loadSites);

// HAP/LEG/SAMP format
std::vector<std::vector<char> > OpenHap(const std::string &hapFile);
std::vector<Bio::snp> OpenLegend(const std::string &legendFile, string chrom);
std::vector<std::string> OpenSamp(const std::string &sampFile);

// BCF
void OpenVCFGZ(const std::string &vcf, const std::string &region,
               std::vector<std::vector<char> > &loadHaps,
               std::vector<Bio::snp> &loadSites, std::vector<std::string> &ids);

struct Init {

  HapFileType hapFileType = HapFileType::IMPUTE2;

  // find a site in the haplotypes for every site in the
  // keep site list
  bool matchAllSites = false;

  // keep all sites that were found within region
  bool keepAllSitesInRegion = false;

  // filter on region
  Bio::Region keepRegion;
  // or filter on site list
  std::vector<Bio::snp> keepSites;

  // keep only specific samples
  std::vector<std::string> keepSampleIDs;
};
}

class HapPanel {

private:
  // input variables
  std::vector<std::string> m_keepSampIDsList;

  std::unordered_map<std::string, size_t> m_keepSampIDsUnorderedMap;
  std::unordered_map<Bio::snp, size_t, Bio::snpKeyHasher> m_keepSites;
  Bio::Region m_keepRegion;

  HapPanelHelper::Init m_init;
  unsigned m_wordSize = WORDSIZE;

  // internal storage / invariants
  std::vector<std::vector<char> > m_haps;
  std::vector<Bio::snp> m_sites;
  std::vector<std::string> m_sampleIDs;
  std::unordered_map<Bio::snp, size_t, Bio::snpKeyHasher> m_sitesUnorderedMap;

  void OrderSamples(std::vector<std::string> &loadIDs,
                    std::vector<std::vector<char> > &loadHaps);
  void SubsetSamples(std::vector<std::string> &loadIDs,
                     std::vector<std::vector<char> > &loadHaps);
  void FilterSites(std::vector<std::vector<char> > &loadHaps,
                   std::vector<Bio::snp> &loadSites);

  void CheckPanel();

public:
  explicit HapPanel() {};
  explicit HapPanel(const std::vector<std::string> &inFiles,
                    HapPanelHelper::Init init);

  unsigned NumWordsPerHap() const {
    return static_cast<unsigned>(
        ceil(static_cast<double>(m_haps.size()) / m_wordSize));
  }

  void FilterSitesOnAlleleFrequency(double upperBound, double lowerBound,
                                    bool useReferenceAlleleFrequency);

  std::vector<uint64_t> Haplotypes_uint64_t() const {
    return HapPanelHelper::Char2BitVec(m_haps, NumWordsPerHap(), m_wordSize);
  };
  std::vector<uint64_t> Haplotypes_uint64_t(size_t numWordsPerHap) const {
    assert(numWordsPerHap >= NumWordsPerHap());
    return HapPanelHelper::Char2BitVec(m_haps, numWordsPerHap, m_wordSize);
  };

  std::string GetID(size_t idx) const {
    return m_sampleIDs.at(idx);
  };
  size_t NumHaps() const {
    return m_haps[0].size();
  };
  size_t NumSites() const {
    return m_sites.size();
  };

  size_t MaxSites() const { return NumWordsPerHap() * m_wordSize; }

  bool empty() const { return m_sites.empty(); }
  /*
  vector<uint64_t> Hap(unsigned hapNum) const {
      assert(hapNum < m_numHaps);
  return &m_haps[hapNum * m_numWordsPerHap];
}
  */
  bool snp_is_equal(const Bio::snp &lhs, size_t siteIdx) {
    return lhs == m_sites.at(siteIdx);
  }

  int FindSiteIndex(const Bio::snp &test) {
    auto got = m_sitesUnorderedMap.find(test);
    if (got == m_sitesUnorderedMap.end())
      return -1;
    else
      return got->second;
  }

  char GetAllele(size_t siteIdx, size_t haplotypeIdx) {
    return m_haps.at(siteIdx).at(haplotypeIdx);
  }
};

#endif
